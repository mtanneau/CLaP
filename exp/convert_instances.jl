import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIF = MOI.FileFormats

# Use this as cache
MOIU.@model(ClassicConeOptimizer,
    (),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
    (MOI.Reals, MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,),
    true
)

cache = ClassicConeOptimizer{Float64}()
bridge = MOI.Bridges.full_bridge_optimizer(cache, Float64)

# Read initial instance
finst = ARGS[1]
cbf = MOIF.CBF.Model()
MOI.read_from_file(cbf, finst)
MOI.copy_to(bridge, cbf)

# 
cbf2 = MOIF.CBF.Model()
MOI.copy_to(cbf2, cache)
MOI.write_to_file(cbf2, ARGS[2])