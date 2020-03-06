include("landp.jl")

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

function main()
    cl_args = parse_commandline(ARGS)
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
        landp(
            cl_args["finst"],
            micp_optimizer, cgcp_optimizer,
            30.0, cl_args["Normalization"], cl_args["Rounds"],
            verbose=false, timer = to
        )
    end

    # Real run
    to = TimerOutput()
    landp(
        cl_args["finst"],
        micp_optimizer, cgcp_optimizer,
        cl_args["TimeLimit"], cl_args["Normalization"], cl_args["Rounds"],
        timer=to
    )

    display(to)
    println()

    return nothing
end

# Run the experiment
main()