using DelimitedFiles
using LinearAlgebra
using SparseArrays

import CPLEX

using CLaP

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
# TODO: provide as input
MOI.set(micp, MOI.TimeLimitSec(), 3600.0)  # 1-hour time limit

MOI.optimize!(micp)

@show MOI.get(micp, MOI.TerminationStatus())    # Termination status
@show MOI.get(micp, MOI.ObjectiveBound())       # Lower bound

# Export solution
# TODO: check for solution status first
x_micp = MOI.get(micp, MOI.VariablePrimal(), sf.var_indices)

writedlm(joinpath(@__DIR__, "res/$finst.sol"), x_micp)