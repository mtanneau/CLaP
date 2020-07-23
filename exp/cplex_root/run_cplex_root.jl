using LinearAlgebra
using SparseArrays
using Logging

import CPLEX

using CLaP

const CUT_PARAMS = [
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

function main(finst, maxrounds=200, maxtime=3600.0; verbose::Bool=true)
    # Read model from file
    cbf = MOI.FileFormats.CBF.Model()
    MOI.read_from_file(cbf, finst)

    # Create standard form problem
    micp_optimizer = MOI.OptimizerWithAttributes(CPLEX.Optimizer,
        "CPX_PARAM_PREIND" => 0,                        # Disable presolve (0)
        "CPX_PARAM_THREADS" => 1,                       # Single thread
        "CPXPARAM_MIP_Strategy_MIQCPStrat" => 2,        # 1: NL-B&B, 2:OA
        [cutparam => -1 for cutparam in CUT_PARAMS]..., # Disable all cuts
        "CPXPARAM_MIP_Limits_Nodes" => 0,               # Root node only
        "CPXPARAM_MIP_Strategy_HeuristicFreq" => -1     # Disable heuristics
    )

    sf = build_standard_form(micp_optimizer, cbf, bridge_type=Float64)

    # Solve root node with aggressive L-a-P cuts
    micp = sf.model
    MOI.set(micp, MOI.RawParameter("CPXPARAM_MIP_Limits_CutPasses"), maxrounds)
    MOI.set(micp, MOI.RawParameter("CPXPARAM_MIP_Limits_CutsFactor"), 1e30)
    MOI.set(micp, MOI.RawParameter("CPXPARAM_MIP_Cuts_LiftProj"), 3)

    # Time limit
    MOI.set(micp, MOI.Silent(), !verbose)
    MOI.set(micp, MOI.TimeLimitSec(), maxtime)  # 1-hour time limit

    # Callback for tracking bound progress
    tstart = time()
    ncalls = 0
    function cbinfo(cbdata)
        x_ = MOI.get(micp, MOI.CallbackVariablePrimal(cbdata), sf.var_indices)
        z_ = dot(sf.c, x_)
        ncalls += 1
        t_ = time() - tstart
        @info "Current stats" z_ t_ ncalls
        return nothing
    end

    MOI.set(micp, MOI.UserCutCallback(), cbinfo)

    tstart = time()
    MOI.optimize!(micp)

    verbose && @show MOI.get(micp, MOI.TerminationStatus())    # Termination status
    verbose && @show MOI.get(micp, MOI.ObjectiveBound())       # Lower bound

    return micp
end

if abspath(PROGRAM_FILE) == @__FILE__
    nrounds = parse(Int, ARGS[2])
    maxtime = parse(Float64, ARGS[3])
    
    with_logger(Logging.NullLogger()) do
        main(ARGS[1], 2, 5.0, verbose=false)
    end

    
    main(ARGS[1], nrounds, maxtime)
end