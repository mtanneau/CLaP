using SparseArrays
using TimerOutputs

"""
    lift_and_project

Compute lift-and-project cuts.

Returns a list of violated split cuts and K* cuts.
"""
function lift_and_project(
    x_, sf::StandardProblem, cgcp_optimizer;
    nb_rounds_max::Int=typemax(Int),
    normalization::Symbol=:Conic,
    kcut_pre_check::Bool=true,
    kcut_post_check::Bool=true,
    strengthen_flag::Bool=true,
    ϵ_integer::Float64=1e-4,
    ϵ_cut_violation::Float64=1e-4,
    timer=TimerOutput(),
    xref::Vector{Float64}=Float64[]
)
    Kcuts = Tuple{Vector{Int}, Vector{Float64}, Float64}[]  # K*    cuts
    Scuts = Tuple{Vector{Int}, Vector{Float64}, Float64}[]  # Split cuts

    # Problem size
    m, n = size(sf.A)

    # Current lower bound
    z_ = dot(x_, sf.c)

    # Clean ̄x
    @timeit timer "Clean x_" for j in 1:n
        if abs(x_[j]) <= 1e-7
            x_[j] = 0.0
        end
    end

    # Compute fractional parts
    f = min.(ceil.(x_) .- x_, x_ .- floor.(x_))
    f .*= sf.vartypes
    p = sortperm(f, rev=true)
    @debug extrema(f)

    # Check if ̄x is conic feasible
    @timeit timer "Conic infeas." H = extract_conic_infeasibilities(x_, sf.cones) ./ 2
    η = minimum(H)

    # Cut generation loop
    nkcuts_pre = 0
    nkcuts_post = 0
    nscuts = 0
    @timeit timer "Cut-generation" for j in p

        # TODO: check stopping criterion

        # Skip non-fractional variables
        f[j] >= ϵ_integer || continue

        kcut_detected = false
        # Early check for K* cut
        @timeit timer "K*-cut pre-check" if kcut_pre_check && normalization == :Standard
            η_ = min(-f[j], f[j] - 1) / 2

            kcut_detected = (η <= -f[j] / 2) && (η <= (f[j]-1)/2)
            kcut_detected && @debug("K*-cut detected", η)

            if kcut_detected
                for (ηi, (kidx, k)) in zip(H, sf.cones)
                    if ηi <= η
                        # Polyhedral cones are 
                        @debug "K* cut found" x_[kidx] k kidx ηi
                        isa(k, POLYHEDRAL_CONE) && continue
                        
                        @assert isa(k, MOI.SecondOrderCone) "Only SOC are supported (is $k)"

                        # Most violated K* cut according to normalization
                        λ = copy(x_[kidx])
                        λ[1] = norm(λ[2:end])
                        λ[2:end] .*= -1
                        λ ./= λ[1]

                        # Record K* cut
                        push!(Kcuts, (kidx, λ, 0.0))
                    end
                end
                @timeit timer "K*-pre" nkcuts_pre += 1
            end
        end

        kcut_detected && continue

        # Look for a split cut
        pi = zeros(n)
        pi[j] = 1.0
        pi0 = floor(x_[j])

        # Build CGCP
        @timeit timer "CGCP-build" cgcp = build_cgcp(
            cgcp_optimizer,
            m, n, sf.A, sf.b, sf.cones,
            x_, pi, pi0,
            nrm=normalization,
            bridge_type=Float64,
            xref=xref
        )
        cgcp_moi = backend(cgcp)

        # Solve CGCP
        @timeit timer "CGCP-solve" MOI.optimize!(cgcp_moi)

        # TODO: track number of solves && IPM iterations

        # Check termination status
        st = MOI.get(cgcp_moi, MOI.TerminationStatus())
        pst = MOI.get(cgcp_moi, MOI.PrimalStatus())
        if !(
            (pst == MOI.FEASIBLE_POINT)
            || (st == MOI.DUAL_INFEASIBLE && pst == MOI.INFEASIBILITY_CERTIFICATE)
        )
            @warn "CGCP exited with status" st pst
            continue
        end

        # Check if found violated cut
        δ = objective_value(cgcp)
        @debug "Cut violation: $δ"
        if δ >= -ϵ_cut_violation && st != MOI.DUAL_INFEASIBLE
            # Cut is not violated
            @debug "No violated cut" δ st pst
            continue
        end

        # Get Farkas multipliers
        @timeit timer "Get Farkas" begin
            v = value.(cgcp[:v])
            norm(v, Inf) >= 1e2 && @warn "|v| is large: $(extrema(v))"
            u0 = value(cgcp[:u0])
            v0 = value(cgcp[:v0])
            λ = value.(cgcp[:λ])
            μ = value.(cgcp[:μ])
        end

        # Check that u0, v0 are not "too negative", otherwise reject the cut
        (u0 >= -1e-5 && v0 >= -1e-5) || continue

        # Clean u0, v0
        u0 = max(0.0, u0)
        v0 = max(0.0, v0)
        (abs(u0) <= 1e-6) && (u0 = 0.0)
        (abs(v0) <= 1e-6) && (v0 = 0.0)

        # Clean conic multipliers
        @timeit timer "Clean λ" for (kidx, k) in sf.cones
            kd = MOI.dual_set(k)
            if isa(kd, MOI.Nonnegatives)
                λ[kidx[1]] = max(0.0, λ[kidx[1]])
                μ[kidx[1]] = max(0.0, μ[kidx[1]])

            elseif isa(kd, MOI.Nonpositives)
                λ[kidx[1]] = min(0.0, λ[kidx[1]])
                μ[kidx[1]] = min(0.0, μ[kidx[1]])

            elseif isa(kd, MOI.Reals)
                # Nothing to do
                # α[kidx[1]] = 0.0

            elseif isa(kd, MOI.Zeros)
                # TODO
                # abs(a1[kidx[1]] - a2[kidx[1]]) <= 1e-5 || @warn "|a1 - a2| = $(abs(a1[kidx[1]] - a2[kidx[1]])) but should be 0"
                λ[kidx[1]] = 0.0
                μ[kidx[1]] = 0.0

            elseif isa(kd, MOI.SecondOrderCone)
                # Shift dual multipliers if needed
                λ[kidx[1]] = max(λ[kidx[1]], norm(λ[kidx[2:end]]))
                μ[kidx[1]] = max(μ[kidx[1]], norm(μ[kidx[2:end]]))

                # TODO: clean multipliers

            elseif isa(kd, MOI.RotatedSecondOrderCone)
                # start by shifting initial coordinates
                λ[kidx[1]] = max(0.0, λ[kidx[1]])
                λ[kidx[2]] = max(0.0, λ[kidx[2]])
                # Check norm
                η = 2 * λ[kidx[1]] * λ[kidx[2]] - sum(λ[kidx[3:end]] .^ 2)
                if η < 0.0
                    # TODO
                    # Increase first two coordinates
                end
                @assert 2 * λ[kidx[1]] * λ[kidx[2]] - sum(λ[kidx[3:end]] .^ 2) >= -1e-7
                @assert 2 * μ[kidx[1]] * μ[kidx[2]] - sum(μ[kidx[3:end]] .^ 2) >= -1e-7

            else
                @warn "Cone $(typeof(k)) not recognized in feasibility check"
            end
        end

        α = λ .- (u0 .* pi)
        β = min(-u0 * pi0, dot(sf.b, v) + v0 * (pi0 + 1))
        δ = dot(α, x_) - β

        # Check for K* cut
        kcut_flag = (iszero(u0) && β <= 1e-6) || (iszero(v0) && (β - dot(sf.b, v)) <= 1e-6)
        @timeit timer "K*-cut post-check" if kcut_post_check && kcut_flag
            @debug "K*-cut" u0 v0 η kcut_detected δ f[j] j
            # Dis-aggregate the cut
            for (kidx, k) in sf.cones
                if iszero(u0)
                    λ_ = λ[kidx]
                else
                    λ_ = μ[kidx]
                end

                # Skip small components
                norm(λ_) >= 1e-5 || continue

                # Record K* cut
                # TODO: add at most one K* cut per cone
                push!(Kcuts, (kidx, λ_, 0.0))
            end
            @timeit timer "count" nkcuts_post += 1
            continue
        end

        # TODO: strengthening
        @timeit timer "Strengthen" if strengthen_flag && (u0+v0) >= 1e-4
            
            Atv = sf.A'v
            for (kidx, k) in sf.cones

                isa(k, NONLINEAR_CONE) && all(sf.vartypes_implied[kidx] .*  iszero.(x_[kidx])) >= 1 && @info "Could have strengthened a cone"
                
                isa(k, Union{MOI.Nonpositives, MOI.Nonnegatives}) || continue
                j_ = kidx[1]

                # Only strengthen integer variables that are at zero
                (sf.vartypes_implied[j_] && iszero(x_[j_])) || continue

                α[j_] = strengthen(Atv[j_], u0, v0, k)
            end
        end

        # TODO: re-scale the cut

        # TODO: book-keeping

        # TODO: validity check (make it optional)

        # Record the cut
        for j_ in 1:n
            if abs(α[j_]) <= 1e-7 
                α[j_] = 0.0
            end
        end
        abs(β) <= 1e-8 && (β = 0.0)

        # Record split cut
        α_ = sparsevec(α)
        push!(Scuts, (α_.nzind, α_.nzval, β))

        nscuts += 1
    end

    return Kcuts, Scuts
end

function extract_conic_infeasibilities(x, cones)
    # TODO: also return corresponding K* cuts
    ncones = length(cones)
    H = zeros(ncones)

    for (i, (kidx, k)) in enumerate(cones)
        isa(k, MOI.Zeros) && continue  # Skip zero cones
        H[i] = epigraph_violation(x[kidx], k)
    end

    return H
end

function epigraph_violation(x, k::MOI.Nonnegatives)
    length(x) == MOI.dimension(k) || throw(DimensionMismatch(
        "x has length $(length(x)) but cone has dimension $(MOI.dimension(k))"
    ))

    return minimum(x)
end

function epigraph_violation(x, k::MOI.Nonpositives)
    length(x) == MOI.dimension(k) || throw(DimensionMismatch(
        "x has length $(length(x)) but cone has dimension $(MOI.dimension(k))"
    ))

    return -maximum(x)
end

function epigraph_violation(x, k::MOI.Reals)
    length(x) == MOI.dimension(k) || throw(DimensionMismatch(
        "x has length $(length(x)) but cone has dimension $(MOI.dimension(k))"
    ))

    return zero(eltype(x))
end

function epigraph_violation(x, k::MOI.Zeros)
    length(x) == MOI.dimension(k) || throw(DimensionMismatch(
        "x has length $(length(x)) but cone has dimension $(MOI.dimension(k))"
    ))
    error("Cannot compute epigraph violation: zero cone has empty interior")
end

function epigraph_violation(x, k::MOI.SecondOrderCone)
    length(x) == MOI.dimension(k) || throw(DimensionMismatch(
        "x has length $(length(x)) but cone has dimension $(MOI.dimension(k))"
    ))

    return x[1] - norm(x[2:end], 2)
end