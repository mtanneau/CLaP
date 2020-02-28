using DelimitedFiles
using LinearAlgebra
using SparseArrays

import CPLEX

include(joinpath(@__DIR__, "../../src/standard_form.jl"))

# Read model from file
cbf = MOI.FileFormats.CBF.Model()
fname = ARGS[1]
fdir = dirname(fname)
finst = basename(fname)[1:end-4]

MOI.read_from_file(cbf, fname)

# Create standard form problem
micp_optimizer = MOI.OptimizerWithAttributes(CPLEX.Optimizer,
    "CPX_PARAM_THREADS" => 1,                       # Single thread
    "CPXPARAM_MIP_Strategy_MIQCPStrat" => 2,        # 1: NL-B&B, 2:OA
)

sf = build_standard_form(micp_optimizer, cbf, bridge_type=Float64)

# Solve root node with aggressive L-a-P cuts
micp = sf.model

# Time limit
MOI.set(micp, MOI.TimeLimitSec(), 3600.0)  # 1-hour time limit

MOI.optimize!(micp)

@show MOI.get(micp, MOI.TerminationStatus())    # Termination status
@show MOI.get(micp, MOI.ObjectiveBound())       # Lower bound

# TODO: export solution
x_micp = MOI.get(micp, MOI.VariablePrimal(), [MOI.VariableIndex(j) for j in 1:length(sf.c)])

writedlm(joinpath(@__DIR__, "res/$finst.sol"), x_micp)