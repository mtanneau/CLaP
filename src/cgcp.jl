using JuMP
using LinearAlgebra

"""

Original problem writes
```
```
and the CGCP writes
```
min  
```
"""
function build_cgcp(
    optimizer::MOI.OptimizerWithAttributes,
    m::Int, n::Int, A, b, cones,
    x_, π, π0
)
    cgcp = Model(optimizer)

    # Create variables
    @variable(cgcp, α[1:n])
    @variable(cgcp, β)

    @variable(cgcp, v[1:m])
    @variable(cgcp, u0 >= 0)
    @variable(cgcp, v0 >= 0)

    # Conic Farkas multipliers
    λ = Vector{JuMP.VariableRef}(undef, n)
    μ = Vector{JuMP.VariableRef}(undef, n)
    
    for (kidx, k) in cones
        d = MOI.dimension(k)
        λ[kidx] .= @variable(cgcp, [1:d] in MOI.dual_set(k))
        μ[kidx] .= @variable(cgcp, [1:d] in MOI.dual_set(k))
    end

    # Constraints
    @constraint(cgcp, f1a, α .==        λ .- (u0 .* π))
    @constraint(cgcp, f1b, α .== A'v .+ μ .+ (v0 .* π))
    @constraint(cgcp, f2a, β  <=           - u0 .* π0)
    @constraint(cgcp, f2b, β  <= b'v + v0 .* (π0 + 1))

    # TODO: normalization
    # @constraint(cgcp, normalization, [1;α] in SecondOrderCone())  # Alpha2 normalization

    # Conic normalization
    # TODO: normalize v as well (?)
    @constraint(cgcp, normalization, u0 + v0 <= 1)
    for (kidx, k) in cones
        kd = MOI.dual_set(k)
        if isa(kd, MOI.Nonnegatives)
            set_normalized_coefficient(normalization, λ[kidx[1]], 1.0)
            set_normalized_coefficient(normalization, μ[kidx[1]], 1.0)
        elseif isa(kd, MOI.Nonpositives)
            set_normalized_coefficient(normalization, λ[kidx[1]], -1.0)
            set_normalized_coefficient(normalization, μ[kidx[1]], -1.0)
        elseif isa(kd, MOI.SecondOrderCone)
            set_normalized_coefficient(normalization, λ[kidx[1]], 1.0)
            set_normalized_coefficient(normalization, μ[kidx[1]], 1.0)
        elseif isa(kd, MOI.Zeros)
            # Nothing to add here
        else
            error("Conic normalization for $(typeof(k)) is not supported")
        end
    end


    # Objective
    @objective(cgcp, Min, dot(x_, α) - β)

    return cgcp
end