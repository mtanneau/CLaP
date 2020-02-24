include(joinpath(@__DIR__, "src/standard_form.jl"))
include(joinpath(@__DIR__, "src/cgcp.jl"))
# Read model from file
cbf = MOI.FileFormats.CBF.Model()
MOI.read_from_file(cbf, ARGS[1])

# Create standard form problem
import CPLEX
sf = build_standard_form(MOI.OptimizerWithAttributes(CPLEX.Optimizer), cbf)

# Do stuff
m, n = size(sf.A)
micp = sf.model

# Set some parameters
MOI.set(micp, MOI.RawParameter("CPX_PARAM_PREIND"), 0)
MOI.set(micp, MOI.RawParameter("CPX_PARAM_THREADS"), 1)

MOI.set(micp, MOI.TimeLimitSec(), 30.0)

# TODO: disable all cuts and disable callback
MOI.set(micp, MOI.RawParameter("CPXPARAM_MIP_Strategy_HeuristicFreq"), -1)  # This should disable heuristics
MOI.set(micp, MOI.RawParameter("CPXPARAM_MIP_Limits_EachCutLimit"), 0)  # Should disable all cuts except GMI
MOI.set(micp, MOI.UserCutCallback(),
    cb_data -> begin
        x = [MOI.VariableIndex(j) for j in 1:n]
        x_ = MOI.get(micp, MOI.CallbackVariablePrimal(cb_data), x)
        
        ncuts = 0
        
        for (j, flag) in enumerate(sf.int_flags)
            flag || continue

            # Check if variable is fractional
            f = max(ceil(x_[j]) - x_[j], x_[j] - floor(x_[j]))
            f >= 1e-4 || continue

            # @info "Splitting on x[$j] = $(x_[j])"

            # Look for a split cut
            pi = zeros(n)
            pi[j] = 1.0
            pi0 = floor(x_[j])

            cgcp = build_cgcp(
                MOI.OptimizerWithAttributes(CPLEX.Optimizer),
                m, n, sf.A, sf.b, sf.cones,
                x_, pi, pi0
            )
            MOI.set(cgcp, MOI.Silent(), true)
            MOI.set(cgcp, MOI.NumberOfThreads(), 1)

            # Solve CGCP
            optimize!(cgcp)

            # TODO: check status
            st = termination_status(cgcp)
            if st != MOI.OPTIMAL
                # @warn st
                continue
            end

            # Check if found violated cut
            v = objective_value(cgcp)

            if v < -1e-4
                # Add the cut
                α = value.(cgcp[:α])
                β = value(cgcp[:β])

                # Submit the cut
                MOI.submit(
                    micp,
                    MOI.UserCut(cb_data),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(α, x), 0.0),
                    MOI.GreaterThan(β)
                )

                ncuts += 1
                break
            end
        end

        # @info "User cut callback was called." ncuts 
        # cgcp = build_cgcp(
        #     MOI.OptimizerWithAttributes(CPLEX.Optimizer),
        #     m, n,
        #     sf.A, sf.b, sf.cones,
        #     x_, pi, pi0
        # )
    end
)

MOI.optimize!(micp)

x = [
    MOI.get(micp, MOI.VariablePrimal(), MOI.VariableIndex(j))
    for j in 1:m
]
@info "Final solution" x

# Try separating a point
# x_ = zeros(n)
# x_[1:3] .= [1.0, sqrt(2)/2, sqrt(2)/2]
# pi = zeros(n)
# pi[2] = 1.0
# pi0 = 0.0

# cgcp = build_cgcp(
#     MOI.OptimizerWithAttributes(CPLEX.Optimizer),
#     m, n,
#     sf.A, sf.b, sf.cones,
#     x_, pi, pi0
# )

# # Solve separation problem
# optimize!(cgcp)

# @info "CGCP value" objective_value(cgcp)
# @info value.(cgcp[:α])
