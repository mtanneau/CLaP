__precompile__()

module CLaP

import MathOptInterface
const MOI = MathOptInterface

export MOI,
    NONLINEAR_CONE, POLYHEDRAL_CONE, MOI_CONES, SCALAR_SETS, SUPPORTED_SETS,
    StandardProblem, build_standard_form,
    extract_implied_integer,
    build_cgcp,
    strengthen,
    lift_and_project


include("standard_form.jl")
include("prepross.jl")
include("cgcp.jl")
include("strengthening.jl")

# Lift-and-project callback
include("lift_and_project.jl")

end  # module