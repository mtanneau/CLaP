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
and, for the split `{π'x ⩽ π0} ∧ {π'x ⩾ π0 +1}`.

The CGCP in full writes (without the normalization)
```
min     α'x - β
s.t.    α =       λ - u0 * π
        α = A'v + μ + v0 * π
        β ⩽         - u0 * π0
        β ⩽ b'v     + v0 * (π0 + 1)
        λ, μ ∈ K*, u0, v0 ⩾ 0
```
and in compact form
```
min     x'λ - u0 * (x'π - π0) + s
s.t.    A'v + (μ - λ) + (u0 + v0) * π       = 0
        b'v + (s - t) + (u0 + v0) * π0 + v0 = 0
        λ, μ ∈ K*, s, t, u0, v0 ⩾ 0
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
    nrm::Symbol = :Conic,
    compact::Bool=true,
    bridge_type::Union{Nothing, Type}=nothing,
    xref::Vector{Float64}=Float64[]
)
    # Some solvers do not support Rotated Second Order cones directly,
    # so we enable bridges as a (hopefully short-term) work-around
    cgcp = JuMP.direct_model(MOI.instantiate(optimizer, with_bridge_type=bridge_type))

    # Create variables
    if !compact
        @variable(cgcp, α[1:n])
        @variable(cgcp, β)
    end

    @variable(cgcp, v[1:m])

    @variable(cgcp, u0 >= 0)
    @variable(cgcp, v0 >= 0)

    # Conic Farkas multipliers
    λ = Vector{JuMP.VariableRef}(undef, n)
    μ = Vector{JuMP.VariableRef}(undef, n)
    
    for (kidx, k) in cones
        kd = MOI.dual_set(k)
        d = MOI.dimension(k)
        if isa(kd, MOI.Nonnegatives)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_lower_bound(λ_, 0.0)
            set_lower_bound(μ_, 0.0)

        elseif isa(kd, MOI.Nonpositives)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_upper_bound(λ_, 0.0)
            set_upper_bound(μ_, 0.0)

        elseif isa(kd, MOI.Zeros)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)
            set_lower_bound(λ_, 0.0)
            set_lower_bound(μ_, 0.0)
            set_upper_bound(λ_, 0.0)
            set_upper_bound(μ_, 0.0)

        elseif isa(kd, MOI.Reals)
            @assert d == 1
            λ_ = @variable(cgcp)
            μ_ = @variable(cgcp)

        elseif isa(kd, MOI.SecondOrderCone)
            λ_ = @variable(cgcp, [1:d] in SecondOrderCone())
            μ_ = @variable(cgcp, [1:d] in SecondOrderCone())

        elseif isa(kd, MOI.RotatedSecondOrderCone)
            @assert d >= 2
            λ_ = @variable(cgcp, [1:d] in RotatedSecondOrderCone())
            μ_ = @variable(cgcp, [1:d] in RotatedSecondOrderCone())

        else
            error("Primal/dual cones $((typeof(k), typeof(kd))) not supported in CGCP")
        end
        λ[kidx] .= λ_
        μ[kidx] .= μ_
    end
    cgcp[:λ] = λ
    cgcp[:μ] = μ

    # Objective and constraints
    if compact
        @variable(cgcp, η1 >= 0)
        @variable(cgcp, η2 >= 0)

        @constraint(cgcp, A'v .+ ( μ .- λ) .+ (u0 + v0) .* π       .== 0.0)
        @constraint(cgcp, b'v  - (η2 - η1)  + (u0 + v0)  * π0 + v0  == 0.0)

        @objective(cgcp, Min, dot(x_, λ) - u0 * (dot(x_, π) - π0) + η1)
        
    else
        # Full form
        @constraint(cgcp, f1a, α .==        λ .- (u0 .*  π))
        @constraint(cgcp, f1b, α .== A'v .+ μ .+ (v0 .*  π))
        @constraint(cgcp, f2a, β  <=           -  u0  *  π0)
        @constraint(cgcp, f2b, β  <= b'v       +  v0  * (π0 + 1))

        @objective(cgcp, Min, dot(x_, α) - β)
    end
    
    # Normalization
    if nrm == :Alpha2
        add_alpha_normalization(m, n, cgcp, cones, π, p=2)

    elseif nrm == :Interior
        length(xref) == length(x_) || error("xref has length $(length(xref))")

        add_interior_normalization(m, n, cgcp, cones, xref - x_, π)

    elseif nrm == :Trivial
        add_trivial_normalization(m, n, cgcp, cones)

    elseif nrm == :Standard
        add_standard_normalization(m, n, cgcp, cones)

    elseif nrm == :Uniform
        add_uniform_normalization(m, n, cgcp, cones)

    elseif nrm == :Combined
        add_combined_normalization(m, n, cgcp, cones, π)

    else
        error("Normalization $nrm is not supported")
    end
    return cgcp
end

function add_alpha_normalization(m::Int, n::Int, cgcp, cones, π; p::Real=2)

    λ, u0 = cgcp[:λ], cgcp[:u0]
    if p == 2
        @constraint(cgcp, normalization, dot(λ + u0 .* π, λ + u0 .* π) <= 1)
    else
        error("Alpha normalization only supported for p={2} (is $p)")
    end

    return cgcp
end

"""
    Interior normalization


"""
function add_interior_normalization(m::Int, n::Int, cgcp, cones, γ, π)
    λ, u0 = cgcp[:λ], cgcp[:u0]
    @constraint(cgcp, normalization, dot(γ, λ) - u0 * dot(γ, π) <= 1)

    return cgcp
end

"""
    add_trivial_normalization

`u0 + v0 <= 1`
"""
function add_trivial_normalization(m::Int, n::Int, cgcp, cones)
    @constraint(cgcp, normalization, cgcp[:u0] + cgcp[:v0] <= 1)
    return cgcp
end

"""
    add_standard_normalization

`ρ'λ + ρ'μ + u0 + v0 ⩽ 1`
"""
function add_standard_normalization(m::Int, n::Int, cgcp, cones)
    ρ = zeros(n)
    for (kidx, k) in cones
        kd = MOI.dual_set(k)

        if isa(kd, MOI.Nonnegatives)
            ρ[kidx] .= 1
        elseif isa(kd, MOI.Nonpositives)
            ρ[kidx] .= -1
        elseif isa(kd, MOI.Zeros)
            # Nothing to do
        elseif isa(kd, MOI.Reals)
            # Should be treated as equality constraints, i.e.,
            # set one of the multipliers to 0 and let the other one free
            @warn "Some variables are fixed. Should treat them as equality constraints" maxlog = 1
        elseif isa(kd, MOI.SecondOrderCone)
            ρ[kidx[1]] = 1
        elseif isa(kd, MOI.RotatedSecondOrderCone)
            ρ[kidx[1:2]] .= 1
        else
            error("Standard normalization for cone $(typeof(k)) is not supported")
        end
    end

    λ, μ = cgcp[:λ], cgcp[:μ]
    @constraint(cgcp, normalization, dot(ρ, λ) + dot(ρ, μ) + cgcp[:u0] + cgcp[:v0] <= 1)
    
    return cgcp
end

"""
    add_uniform_normalization

`ρ'λ + ρ'μ ⩽ 1`
"""
function add_uniform_normalization(m::Int, n::Int, cgcp, cones)
    ρ = zeros(n)
    for (kidx, k) in cones
        kd = MOI.dual_set(k)

        if isa(kd, MOI.Nonnegatives)
            ρ[kidx] .= 1
        elseif isa(kd, MOI.Nonpositives)
            ρ[kidx] .= -1
        elseif isa(kd, MOI.Zeros)
            # Nothing to do
        elseif isa(kd, MOI.Reals)
            # Should be treated as equality constraints, i.e.,
            # set one of the multipliers to 0 and let the other one free
            
        elseif isa(kd, MOI.SecondOrderCone)
            ρ[kidx[1]] = 1
        elseif isa(kd, MOI.RotatedSecondOrderCone)
            ρ[kidx[1:2]] .= 1
        else
            error("Uniform normalization for cone $(typeof(k)) is not supported")
        end
    end

    λ, μ = cgcp[:λ], cgcp[:μ]
    @constraint(cgcp, normalization, dot(ρ, λ) + dot(ρ, μ) <= 1)

    return cgcp
end

function add_combined_normalization(m::Int, n::Int, cgcp, cones, π)
    ρ = zeros(n)
    for (kidx, k) in cones
        kd = MOI.dual_set(k)

        if isa(kd, MOI.Nonnegatives)
            ρ[kidx] .= 1
        elseif isa(kd, MOI.Nonpositives)
            ρ[kidx] .= -1
        elseif isa(kd, MOI.Zeros)
            # Nothing to do
        elseif isa(kd, MOI.Reals)
            # Should be treated as equality constraints, i.e.,
            # set one of the multipliers to 0 and let the other one free
        elseif isa(kd, MOI.SecondOrderCone)
            ρ[kidx[1]] = 1
        elseif isa(kd, MOI.RotatedSecondOrderCone)
            ρ[kidx[1:2]] .= 1
        else
            error("Uniform normalization for cone $(typeof(k)) is not supported")
        end
    end

    λ, μ = cgcp[:λ], cgcp[:μ]
    u0, v0 = cgcp[:u0], cgcp[:v0]

    if iszero(dot(ρ, π))
        @constraint(cgcp, normalization, dot(ρ, λ) <= 1)
    else
        @constraint(cgcp, normalization, dot(ρ, λ) + u0 <= 1)
    end

    @constraint(cgcp, normalization2, dot(ρ, λ) + dot(ρ, μ) + u0 + v0 <= 100)

    return cgcp
end