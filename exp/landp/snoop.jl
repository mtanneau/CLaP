include("landp.jl")

# Run the experiment on a small instance
# This one contains Second Order Cones
main([
    "--MICPSolver", "CPLEX",
    "--CGCPSolver", "Gurobi",
    "--Rounds", "200",
    "--TimeLimit", "120.0",
    "--Normalization", "Conic",
    joinpath(@__DIR__, "../../dat/cblib/nvs03.cbf")
])

main([
    "--MICPSolver", "CPLEX",
    "--CGCPSolver", "Gurobi",
    "--Rounds", "200",
    "--TimeLimit", "120.0",
    "--Normalization", "PureConic",
    joinpath(@__DIR__, "../../dat/cblib/sssd-strong-15-8.cbf")
])

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

GRB_ENV = Gurobi.Env()  # To avoid checking out multiple Gurobi licenses
cgcp_optimizer = MOI.OptimizerWithAttributes(
    () -> Gurobi.Optimizer(GRB_ENV),
    "Threads" => 1,
    "OutputFlag" => 0,
    "Presolve" => 0,
    "QCPDual" => 0
)

to = TimerOutput()
landp(
    joinpath(@__DIR__, "../../dat/cblib/nvs03.cbf"),
    micp_optimizer, cgcp_optimizer,
    120.0, :PureConic, 10,
    verbose=true, timer=to
)