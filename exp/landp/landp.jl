using Logging
using LinearAlgebra
using SparseArrays

using ArgParse
using DelimitedFiles
using TimerOutputs

using JuMP

# Solvers
import CPLEX
import Gurobi
import Mosek
import MosekTools

# CPLEX cut parameters
# We set all these to 0 to de-activate all cuts
const CPLEX_CUT_PARAMS = [
    "CPXPARAM_MIP_Cuts_BQP",
    "CPXPARAM_MIP_Cuts_Cliques",
    "CPXPARAM_MIP_Cuts_Covers", 
    "CPXPARAM_MIP_Cuts_Disjunctive",
    "CPXPARAM_MIP_Cuts_FlowCovers",
    "CPXPARAM_MIP_Cuts_PathCut",
    "CPXPARAM_MIP_Cuts_Gomory",
    "CPXPARAM_MIP_Cuts_GUBCovers",
    "CPXPARAM_MIP_Cuts_Implied",
    "CPXPARAM_MIP_Cuts_LocalImplied",
    "CPXPARAM_MIP_Cuts_LiftProj",
    "CPXPARAM_MIP_Cuts_MIRCut",
    # "CPXPARAM_MIP_Cut_MCFCut",  # For some reason this one raises an error
    "CPXPARAM_MIP_Cuts_RLT",
    "CPXPARAM_MIP_Cuts_ZeroHalfCut"
]

using CLaP

import Base.RefValue

function landp_closure(
    finst::String,
    micp_optimizer, cgcp_optimizer,
    time_limit::Float64, nrm::Symbol, max_rounds::Int;
    verbose::Bool=true,
    refine_oa::Bool=true, conic_feas_tol::Float64=1e-1,
    timer=TimerOutput()
)
    # Read model from file
    cbf = MOI.FileFormats.CBF.Model()
    MOI.read_from_file(cbf, finst)

    # Convert to standard form
    @timeit timer "Build SF" sf = build_standard_form(micp_optimizer, cbf, bridge_type=Float64)

    # TODO: add pre-processing option
    @timeit timer "Prepross" extract_implied_integer(sf)
    @info "Pre-processing: $(sum(sf.vartypes)) integer variables ($(sum(sf.vartypes_implied)) implied integers)"

    # TODO: add cut validity checker option
    x_micp = readdlm(joinpath(@__DIR__, "../cplex_micp/res/$(basename(finst)[1:end-4]).sol"))
    @info "MICP objective value: $(dot(sf.c, x_micp))"

    # Set some parameters
    m, n = size(sf.A)
    micp = sf.model
    MOI.set(micp, MOI.TimeLimitSec(), time_limit)
    MOI.set(micp, MOI.Silent(), !verbose)

    # Book-keeping
    ncalls = 0
    ncuts_tot = 0
    nkcuts_tot = 0
    nrounds = 0

    # Set callback
    x = sf.var_indices
    function cblandp(cbdata)
        ncalls += 1

        nrounds > max_rounds && return nothing

        # Get fractional point
        x_ = MOI.get(micp, MOI.CallbackVariablePrimal(cbdata), x)

        # Refine outer approximation before separating cuts
        @timeit timer "Refine OA" if refine_oa
            H = CLaP.extract_conic_infeasibilities(x_, sf.cones) ./ 2
            kflag = false
            for ((kidx, k), η) in zip(sf.cones, H)
                isa(k, CLaP.POLYHEDRAL_CONE) && continue
                @assert isa(k, MOI.SecondOrderCone) "Only SOC are supported (is $k)"

                if η <= -conic_feas_tol / 2
                    # Add K* cut
                    λ = copy(x_[kidx])
                    λ[1] = norm(λ[2:end])
                    λ[2:end] .*= -1
                    λ ./= λ[1]

                    MOI.submit(
                        micp,
                        MOI.UserCut(cbdata),
                        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(λ, x[kidx]), 0.0),
                        MOI.GreaterThan(0.0)
                    )
                    kflag = true
                end
            end

            kflag && return nothing
        end

        nrounds += 1
        # Compute cuts
        Kcuts, Scuts = CLaP.lift_and_project(
            x_, sf, cgcp_optimizer,
            normalization=nrm,
            kcut_pre_check=true,
            kcut_post_check=true,
            strengthen_flag=true,
            timer=timer
        )

        # Submit cuts
        for (js, α, β) in Kcuts
            MOI.submit(
                micp,
                MOI.UserCut(cbdata),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(α, x[js]), 0.0),
                MOI.GreaterThan(β)
            )
        end

        for (js, α, β) in Scuts
            MOI.submit(
                micp,
                MOI.UserCut(cbdata),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(α, x[js]), 0.0),
                MOI.GreaterThan(β)
            )
        end

        # Book-keeping
        nkcuts = length(Kcuts)
        nscuts = length(Scuts)
        @info "Stats for round $nrounds" nkcuts nscuts

        nkcuts_tot += nkcuts
        ncuts_tot += nscuts

        return nothing
    end

    tstart = time()
    MOI.set(micp, MOI.UserCutCallback(), cblandp)

    # Solve
    @timeit timer "MICP" MOI.optimize!(micp)

    # Result log
    @info "User callback was called $(ncalls[]) times" ncuts_tot[] nkcuts_tot[]

    @info MOI.get(micp, MOI.TerminationStatus())
    @info "Final bound: $(MOI.get(micp, MOI.ObjectiveBound()))"

    return micp
end

function parse_commandline(cl_args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--MICPSolver"
            help = "MICP solver"
            arg_type = Symbol
            default = :CPLEX
        "--CGCPSolver"
            help = "CGCP solver"
            arg_type = Symbol
            default = :CPLEX
        "--Normalization"
            help = "Normalization condition"
            arg_type = Symbol
            default = :Conic
        "--Rounds"
            help = "Maximum number of cutting plane rounds"
            arg_type = Int
            default = typemax(Int)
        "--TimeLimit"
            help = "Time limit (in seconds)"
            arg_type = Float64
            default = 1e30
        "finst"
            help = "Instance file (in CBF format)"
            required = true
    end

    return parse_args(cl_args, s)
end

function main(CLARGS)
    cl_args = parse_commandline(CLARGS)
    @info("User options", cl_args)

    # Set solvers
    if cl_args["MICPSolver"] == :CPLEX
        micp_optimizer = MOI.OptimizerWithAttributes(
            CPLEX.Optimizer,
            "CPX_PARAM_PREIND" => 0,
            "CPX_PARAM_THREADS" => 1,
            # disable all cuts
            [cutparam => -1 for cutparam in CPLEX_CUT_PARAMS]...,
            "CPXPARAM_MIP_Limits_CutsFactor" => 1e30,
            "CPXPARAM_MIP_Strategy_HeuristicFreq" => -1,
            "CPXPARAM_MIP_Limits_Nodes" => 0,
            "CPXPARAM_MIP_Strategy_MIQCPStrat" => 2  # OA algorithm
        )
    elseif cl_args["MICPSolver"] == :Gurobi
        GRB_ENV = Gurobi.Env()  # To avoid checking out multiple Gurobi licenses
        micp_optimizer = MOI.OptimizerWithAttributes(
            () -> Gurobi.Optimizer(GRB_ENV),
            "Threads" => 1,
            "Presolve" => 0,
            "Heuristics" => 0,
            "Cuts" => 0,
            "MIQCPMethod" => 1,
            "NodeLimit" => 1,
            "QCPDual" => 0
        )
    else
        error("MICP solver $(cl_args["MICPSolver"]) is not supported.\n Possible values are `CPLEX` and `Gurobi`.")
    end

    if cl_args["CGCPSolver"] == :CPLEX
        cgcp_optimizer = MOI.OptimizerWithAttributes(
            CPLEX.Optimizer,
            "CPX_PARAM_THREADS" => 1,
            "CPX_PARAM_PREIND" => 0,
            "CPX_PARAM_SCRIND" => 0
        )

    elseif cl_args["CGCPSolver"] == :Gurobi
        GRB_ENV = Gurobi.Env()  # To avoid checking out multiple Gurobi licenses
        cgcp_optimizer = MOI.OptimizerWithAttributes(
            () -> Gurobi.Optimizer(GRB_ENV),
            "Threads" => 1,
            "OutputFlag" => 0,
            "Presolve" => 0,
            "QCPDual" => 0
        )

    elseif cl_args["CGCPSolver"] == :Mosek
        cgcp_optimizer = MOI.OptimizerWithAttributes(
            Mosek.Optimizer,
            "NUM_THREADS" => 1,
            "LOG" => 0
        )
    else
        error("CGCP solver $(cl_args["CGCPSolver"]) is not supported.\n Possible values are `CPLEX`, `Gurobi` and `Mosek`.")
    end

    # Warm-up run
    to = TimerOutput()
    @info "Warming up..."
    with_logger(Logging.NullLogger()) do  # this will de-activate logs
        landp_closure(
            cl_args["finst"],
            micp_optimizer, cgcp_optimizer,
            30.0, cl_args["Normalization"], cl_args["Rounds"],
            verbose=false, timer = to
        )
    end

    # Real run
    to = TimerOutput()
    landp_closure(
        cl_args["finst"],
        micp_optimizer, cgcp_optimizer,
        cl_args["TimeLimit"], cl_args["Normalization"], cl_args["Rounds"],
        timer=to
    )

    print_timer(to)
    println()

    return nothing
end

# Run the main method if script is called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end