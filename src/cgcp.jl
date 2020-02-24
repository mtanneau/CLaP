using JuMP
using LinearAlgebra

"""
    build_cgcp(optimizer, ...)

Build a conic cut-generating problem to separate a split cut.

Original problem writes
```
min     c'x
s.t.    A x = b
          x ∈ K
```
and, for the split `{π'x ⩽ π0} ∧ {π'x ⩾ π0 +1}`, the CGCP writes
```
min     α'x - β
s.t.    α =       λ - u0 * π
        α = A'v + μ + v0 * π
        β ⩽         - u0 * π0
        β ⩽ b'v     + v0 * (π0 + 1)
        λ, μ ∈ K*, u0, v0 ⩾ 0
```

# Arguments
* `optimizer`: Solver for the CGCP. Must support conic constraints
* `m::Int`: Number of equality constraints
* `n::Int`: Number of variables
* `A::AbstractMatrix`: the constraint matrix
* `b::Vector`: The right-hand side of equality constraints
* `cones::Vector{Tuple{Vector{Int}, MOI.AbstractSet}}`: Cones
* `x_::Vector`: the fractional point to separate
* `π::Vector, π0::Float64`: Coefficients of the split
* `nrm::Symbol`: The normalization used in CGCP
    * `:Alpha2`: α-normalization `|α| ⩽ 1`
    * `:Conic`: conic normalization `ρ'λ + ρ'μ + u0 + v0 ⩽ 1` with `ρ ∈ int(K)`
"""
function build_cgcp(
    optimizer::MOI.OptimizerWithAttributes,
    m::Int, n::Int, A::SparseMatrixCSC, b::Vector{Float64}, cones,
    x_::Vector{Float64}, π::Vector{Float64}, π0::Float64;
    nrm::Symbol = :Conic
)
    cgcp = JuMP.direct_model(MOI.instantiate(optimizer))

    # Create variables
    @variable(cgcp, α[1:n])
    @variable(cgcp, β)
    
    @variable(cgcp, -1 <= v[1:m] <= 1)  # Put bounds on v to avoid numerical issues
    @variable(cgcp, u0 >= 0)
    @variable(cgcp, v0 >= 0)

    # Conic Farkas multipliers
    λ = Vector{JuMP.VariableRef}(undef, n)
    μ = Vector{JuMP.VariableRef}(undef, n)
    
    for (kidx, k) in cones
        d = MOI.dimension(k)
        if isa(k, MOI.Nonnegatives)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_lower_bound(λ_, 0.0)
            set_lower_bound(μ_, 0.0)
        elseif isa(k, MOI.Nonpositives)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_upper_bound(λ_, 0.0)
            set_upper_bound(μ_, 0.0)
        elseif isa(k, MOI.Reals)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_lower_bound(λ_, 0.0)
            set_lower_bound(μ_, 0.0)
            set_upper_bound(λ_, 0.0)
            set_upper_bound(μ_, 0.0)
        elseif isa(k, MOI.SecondOrderCone)
            λ_ = @variable(cgcp, [1:d] in SecondOrderCone())
            μ_ = @variable(cgcp, [1:d] in SecondOrderCone())
        else
            error("Cone $(typeof(k)) is not supported in CGCP")
        end
        λ[kidx] .= λ_
        μ[kidx] .= μ_
    end
    cgcp[:λ] = λ
    cgcp[:μ] = μ

    # Constraints
    @constraint(cgcp, f1a, α .==        λ .- (u0 .*  π))
    @constraint(cgcp, f1b, α .== A'v .+ μ .+ (v0 .*  π))
    @constraint(cgcp, f2a, β  <=           -  u0  *  π0)
    @constraint(cgcp, f2b, β  <= b'v       +  v0  * (π0 + 1))

    # Normalization
    if nrm == :Alpha2
        # Alpha normalization
        # |α| ⩽ 1
        @constraint(cgcp, normalization, sum(α .^ 2) <= 1)

    elseif nrm == :Conic
        # Conic normalization
        # ρ'λ + ρ'μ + u0 + v0 ⩽ 1
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

    else
        error("Normalization $nrm is not supported.")
    end

    # Objective
    @objective(cgcp, Min, dot(x_, α) - β)

    return cgcp
end