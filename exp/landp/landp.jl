using Logging
using LinearAlgebra
using SparseArrays

using ArgParse
using DelimitedFiles
using TimerOutputs

# Solvers
import CPLEX
import Gurobi
import Mosek
import MosekTools

# CPLEX cut parameters
# We set all these to 0 to de-activate all cuts
const CPLEX_CUT_PARAMS = [
    "CPXPARAM_MIP_Cuts_BQP",
    "CPXPARAM_MIP_Cuts_Cliques",
    "CPXPARAM_MIP_Cuts_Covers", 
    "CPXPARAM_MIP_Cuts_Disjunctive",
    "CPXPARAM_MIP_Cuts_FlowCovers",
    "CPXPARAM_MIP_Cuts_PathCut",
    "CPXPARAM_MIP_Cuts_Gomory",
    "CPXPARAM_MIP_Cuts_GUBCovers",
    "CPXPARAM_MIP_Cuts_Implied",
    "CPXPARAM_MIP_Cuts_LocalImplied",
    "CPXPARAM_MIP_Cuts_LiftProj",
    "CPXPARAM_MIP_Cuts_MIRCut",
    # "CPXPARAM_MIP_Cut_MCFCut",  # For some reason this one raises an error
    "CPXPARAM_MIP_Cuts_RLT",
    "CPXPARAM_MIP_Cuts_ZeroHalfCut"
]

# Code
include(joinpath(@__DIR__, "../../src/standard_form.jl"))
include(joinpath(@__DIR__, "../../src/prepross.jl"))
include(joinpath(@__DIR__, "../../src/strengthening.jl"))
include(joinpath(@__DIR__, "../../src/cgcp.jl"))

import Base.RefValue

function landp_callback(
    cbdata, micp, x_micp, sf::StandardProblem, cgcp_optimizer,
    ncalls::RefValue{Int}, ncuts_tot::RefValue{Int}, ncgcp_solve::RefValue{Int},
    max_rounds, nrm,
    timer, tstart, time_limit
)

    m, n = size(sf.A)

    x = sf.var_indices
    x_ = MOI.get(micp, MOI.CallbackVariablePrimal(cbdata), x)

    # Compute objective value and log
    tnow = time()
    @info "Objective value after $(ncalls[]) rounds and $(tnow - tstart) seconds: $(dot(sf.c, x_))"

    # Update number of rounds
    ncalls[] += 1
    ncalls[] <= max_rounds || return nothing

    # Clean x_
    @timeit timer "Cleaning" for j in 1:n
        if abs(x_[j]) <= 1e-7
            x_[j] = 0.0
        end
    end

    # We try to separate most fractional coordinates first
    f = min.(ceil.(x_) .- x_, x_ .- floor.(x_))
    f .*= sf.vartypes
    p = sortperm(f, rev=true)

    @timeit timer "Cut-generation" for j in p

        # Time limit
        time() - tstart <= time_limit || break

        # Skip non-fractional variables
        f[j] >= 1e-4 || continue
        
        # Look for a split cut
        pi = zeros(n)
        pi[j] = 1.0
        pi0 = floor(x_[j])

        # Build CGCP
        # TODO: make this faster
        @timeit timer "CGCP-build" cgcp = build_cgcp(
            cgcp_optimizer,
            m, n, sf.A, sf.b, sf.cones,
            x_, pi, pi0,
            nrm=nrm,
            bridge_type=Float64
        )
        cgcp_moi = backend(cgcp)
        ncgcp_solve[] += 1

        # Solve CGCP
        @timeit timer "CGCP-solve" MOI.optimize!(cgcp_moi)

        # Check termination status
        st = MOI.get(cgcp_moi, MOI.TerminationStatus())
        if st != MOI.OPTIMAL
            @warn "CGCP exited with status" st
            continue
        end

        # Check if found violated cut
        δ = objective_value(cgcp)

        if δ >= -1e-4
            # Cut is not violated
            continue
        end

        v = value.(cgcp[:v])
        norm(v, Inf) >= 1e2 && @warn "|v| is large: $(extrema(v))"
        u0 = value(cgcp[:u0])
        v0 = value(cgcp[:v0])
        λ = value.(cgcp[:λ])
        μ = value.(cgcp[:μ])

        # Check that u0, v0 are ⩾ 0
        # If not, we reject the cut
        u0 >= 0.0 && v0 >= 0.0 || continue

        # @info "Norm of multipliers" u0 + v0 norm(v, 2) norm(λ, 2) norm(μ, 2) dot(λ, x_)

        a1 = -u0 .* pi
        a2 = sf.A'v .+ v0 .* pi

        # Check validity of conic multipliers
        @timeit timer "CGCP-clean λ" for (kidx, k) in sf.cones
            kd = MOI.dual_set(k)

            if isa(kd, MOI.Nonnegatives)
                # α[kidx[1]] = max(a1[kidx[1]], a2[kidx[1]])
                λ[kidx[1]] = max(0.0, λ[kidx[1]])
                μ[kidx[1]] = max(0.0, μ[kidx[1]])

            elseif isa(kd, MOI.Nonpositives)
                # α[kidx[1]] = min(a1[kidx[1]], a2[kidx[1]])
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
                    # Increase first two coordinates
                    
                end
                @assert 2 * λ[kidx[1]] * λ[kidx[2]] - sum(λ[kidx[3:end]] .^ 2) >= -1e-7
                @assert 2 * μ[kidx[1]] * μ[kidx[2]] - sum(μ[kidx[3:end]] .^ 2) >= -1e-7

            else
                @warn "Cone $(typeof(k)) not recognized in feasibility check"
            end
        end

        β = min(-u0 * pi0, dot(sf.b, v) + v0 * (pi0 + 1))
        α = λ .- (u0 .* pi)
        σ = α - a2

        # Check that \mu ∈ K*
        @timeit timer "CGCP-clean μ" for (kidx, k) in sf.cones
            kd = MOI.dual_set(k)

            if isa(kd, MOI.Nonnegatives)
                σ[kidx[1]] >= -1e-5 || @warn "σ = $(σ[kidx[1]]) but should be ⩾ 0"
                
            elseif isa(kd, MOI.Nonpositives)
                σ[kidx[1]] <=  1e-5 || @warn "σ = $(σ[kidx[1]]) but should be ⩽ 0"

            elseif isa(kd, MOI.Reals)
                # Nothing to do

            elseif isa(kd, MOI.Zeros)
                # TODO
                # abs(σ[kidx[1]]) <= 1e-5 || @warn "σ = $(σ[kidx[1]]) but should be 0"

            elseif isa(kd, MOI.SecondOrderCone)
                # Shift dual multipliers if needed
                σ[kidx[1]] >= norm(σ[kidx[2:end]], 2) - 1e-6 || @warn "σ1 - |σ| = $(σ[kidx[1]] - norm(σ[kidx[2:end]], 2)) but should be ⩾ 0"

                # α[kidx[1]] += max(0.0, norm(σ[kidx[2:end]]) - σ[kidx[1]])

            elseif isa(kd, MOI.RotatedSecondOrderCone)
                # TODO
                (
                    (σ[kidx[1]] >= - 1e-8) && (σ[kidx[2]] >= - 1e-8)
                    && (2 * σ[kidx[1]] * σ[kidx[2]] >= (sum(σ[kidx[3:end]] .^ 2) - 1e-6))
                ) || @warn "RSOC multiplier invalid.\nσ1 = $(σ[kidx[1]]) (should be ⩾ 0)\nσ1 = $(σ[kidx[2]]) (should be ⩾ 0)\n2*σ1*σ1 - |σ|^2 = $(2 * σ[kidx[1]] * σ[kidx[2]] - sum(σ[kidx[3:end]] .^ 2)) (should be ⩾ 0)"

            else
                @warn "Cone $(typeof(k)) not recognized in feasibility check"
            end
        end

        δ = dot(α, x_) - β

        # Strengthening
        # TODO: add option to strengthen or not
        @timeit timer "Strengthen" if u0 + v0 > 1e-4
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

        # TODO: Proper cleaning and re-scaling for numerics
        @timeit timer "CGCP-clean α" for j_ in 1:n
            if abs(α[j_]) <= 1e-7 
                α[j_] = 0.0
            end
        end
        abs(β) <= 1e-8 && (β = 0.0)
        r = norm(α, 2)

        # TODO: book-keeping

        # TODO: check that cut does not cut optimal solution
        z = (dot(α, x_micp) - β) / r
        z >= -1e-6 || (@error("Optimal solution is cut by $z", norm(x_, 2), norm(α, 2), β, u0, v0, j, x_[j], f[j]); continue)

        # Submit the cut
        @timeit timer "Submit" MOI.submit(
            micp,
            MOI.UserCut(cbdata),
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(α ./ r, x), 0.0),
            MOI.GreaterThan(β / r)
        )

        ncuts_tot[] += 1  # book-keeping

    end  # cut-generation loop
    return nothing
end  # callback function


function landp(
    finst::String,
    micp_optimizer, cgcp_optimizer,
    time_limit::Float64, nrm::Symbol, max_rounds::Int;
    verbose::Bool=true, timer=TimerOutput()
)
    # Read model from file
    cbf = MOI.FileFormats.CBF.Model()
    MOI.read_from_file(cbf, finst)

    # Convert to standard form
    @timeit timer "Build SF" sf = build_standard_form(micp_optimizer, cbf, bridge_type=Float64)

    # TODO: add pre-processing option
    @timeit timer "Prepross" extract_implied_integer(sf)
    @info "Pre-processing: $(sum(sf.vartypes)) integer variables ($(sum(sf.vartypes_implied)) implied integers)"

    # TODO: add cut validity checker option
    x_micp = readdlm(joinpath(@__DIR__, "../cplex_micp/res/$(basename(finst)[1:end-4]).sol"))
    @info "MICP objective value: $(dot(sf.c, x_micp))"

    # Set some parameters
    m, n = size(sf.A)
    micp = sf.model
    MOI.set(micp, MOI.TimeLimitSec(), time_limit)
    MOI.set(micp, MOI.Silent(), !verbose)

    # Book-keeping
    ncalls = Ref(0)
    ncuts_tot = Ref(0)
    ncgcp_solve = Ref(0)

    # Set callback
    tstart = time()
    MOI.set(micp, MOI.UserCutCallback(),
        cbdata -> landp_callback(
            cbdata, micp, x_micp, sf, cgcp_optimizer,
            ncalls, ncuts_tot, ncgcp_solve,
            max_rounds, nrm, timer, tstart, time_limit
        )
    )

    # Solve
    @timeit timer "MICP" MOI.optimize!(micp)

    # Result log
    @info "User callback was called $(ncalls[]) times" ncuts_tot[] ncgcp_solve[]

    @info MOI.get(micp, MOI.TerminationStatus())
    @info MOI.get(micp, MOI.ObjectiveBound())
end