using SparseArrays

import MathOptInterface
const MOI = MathOptInterface

# TODO: use a trait instead
const NONLINEAR_CONE = Union{
    MOI.SecondOrderCone
}
const POLYHEDRAL_CONE = Union{
    MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives
}
const MOI_CONES = Union{POLYHEDRAL_CONE, NONLINEAR_CONE}
const SCALAR_SETS = Union{
    MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}
}
const SUPPORTED_SETS = Union{
    SCALAR_SETS, MOI_CONES, MOI.Integer
}

mutable struct StandardProblem
    model::MOI.ModelLike

    A::SparseMatrixCSC{Float64, Int}  # constraint matrix
    b::Vector{Float64}  # Right-hand sides
    c::Vector{Float64}  # Objective
    c0::Float64  # Objective constant

    cones::Vector{Tuple{Vector{Int}, MOI.AbstractSet}}  # cones
    int_flags::Vector{Bool}  # Integer flags
end

"""
    convert_to_standard_form(m::MOI.ModelLike)

Convert instance to standard form.
"""
function build_standard_form(
    optimizer::MOI.OptimizerWithAttributes,
    src::MOI.ModelLike
)
    dst = MOI.instantiate(optimizer)

    b = Float64[]
    arows = Int[]
    acols = Int[]
    avals = Float64[]
    
    cones = Tuple{Vector{Int}, MOI.AbstractSet}[]

    # First, copy all variables
    x_old = MOI.get(src, MOI.ListOfVariableIndices())
    nvar = length(x_old)  # Original number of variables
    x = MOI.add_variables(dst, nvar)  # variables in new model
    int_flags = zeros(Bool, nvar)  # integer flags
    @info "Original problem has $nvar variables"
    
    varindices = Dict(x_old .=> 1:nvar)  # Map from old indices to new ones

    # Objective coefficients
    c = zeros(Float64, nvar)
    fobj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c0 = fobj.constant
    for t in fobj.terms
        j = varindices[t.variable_index]
        c[j] = t.coefficient
    end

    MOI.set(dst,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), c0)
    )

    # Objective sense
    MOI.set(dst, MOI.ObjectiveSense(), MOI.get(src, MOI.ObjectiveSense()))

    @info "Objective imported successfully"

    # Now, constraints
    ncon = 0
    contypes = MOI.get(src, MOI.ListOfConstraints())
    b = Float64[]  # This will be the right-hand side

    for (F, S) in contypes
        S <: SUPPORTED_SETS || error("Set $S is not supported.")

        # Get list of corresponding constraint indices
        conindices = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())

        for cidx in conindices
            fun = MOI.get(src, MOI.ConstraintFunction(), cidx)
            set = MOI.get(src, MOI.ConstraintSet(), cidx)

            # Update standard form
            add_to_standard_form!(dst,
                fun, set,
                varindices,
                c, b, arows, acols, avals, cones, int_flags
            )
        end
    end

    nvar = length(c)
    ncon = length(b)

    @info "Problem stats:\n\tVariables  : $nvar\n\tConstraints: $ncon\n\tCones      : $(length(cones))"
    
    return StandardProblem(dst, sparse(arows, acols, avals, ncon, nvar), b, c, c0, cones, int_flags)
end

"""

Convert integer restriction
"""
function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.SingleVariable,
    ::MOI.Integer,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)

    j = varindices[fun.variable]
    int_flags[j] = true

    MOI.add_constraint(dst, MOI.SingleVariable(MOI.VariableIndex(j)), MOI.Integer())

    return nothing
end

"""

Convert vector polyhedral cone into series of smaller cones
"""
function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.VectorOfVariables,
    set::POLYHEDRAL_CONE,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)
    if isa(set, MOI.Zeros)
        set_ = MOI.EqualTo(0.0)
        k = MOI.Zeros(1)
    elseif isa(set, MOI.Nonnegatives)
        set_ = MOI.GreaterThan(0.0)
        k = MOI.Nonnegatives(1)
    elseif isa(set, MOI.Nonpositives)
        set_ = MOI.LessThan(0.0)
        k = MOI.Nonpositives(1)
    else
        error("")
    end

    for j_old in fun.variables
        j = varindices[j_old]

        MOI.add_constraint(dst,
            MOI.SingleVariable(MOI.VariableIndex(j)),
            set_
        )

        push!(cones, ([j], k))
    end

    return nothing
end

"""

Convert vector polyhedral constraint into series of smaller cones.

This adds `m` variables, where `m` is the dimension of the conic constraint
"""
function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.VectorAffineFunction,
    set::POLYHEDRAL_CONE,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)
    nvar = length(c)
    ncon = length(b)

    m = MOI.dimension(set)

    if isa(set, MOI.Nonnegatives)
        set_ = MOI.GreaterThan(0.0)
        k = MOI.Nonnegatives(1)
    elseif isa(set, MOI.Nonpositives)
        set_ = MOI.LessThan(0.0)
        k = MOI.Nonpositives(1)
    else
        error("")
    end

    # Add slack variables
    s = MOI.add_variables(dst, m)
    # Set bounds on slack variables
    for i in 1:m
        MOI.add_constraint(dst,
            MOI.SingleVariable(s[i]),
            set_
        )
        push!(cones, ([nvar + i], k))
    end
    append!(c, zeros(m))
    append!(int_flags, zeros(Bool, m))

    # Update matrix coefficients
    append!(arows, ncon .+ collect(1:m))
    append!(acols, nvar .+ collect(1:m))
    append!(avals, -ones(m))

    # Extract terms of VectorAffineFunction
    funs_ = [
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(-1.0, s[i])], 0.0)
        for i in 1:m
    ]
    for t in fun.terms
        i = t.output_index
        j = varindices[t.scalar_term.variable_index]
        v = t.scalar_term.coefficient
        push!(funs_[i].terms, MOI.ScalarAffineTerm(v, MOI.VariableIndex(j)))
        
        # Update coefficients of A
        push!(arows, ncon + i)
        push!(acols, j)
        push!(avals, v)
    end

    # Add each constraint separately
    for i in 1:m
        # We add each equality constraint separately
        rhs = -fun.constants[i]
        push!(b, rhs)

        # Add equality constraint
        MOI.add_constraint(dst, funs_[i], MOI.EqualTo(rhs))
    end

    return nothing
end

function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.VectorAffineFunction,
    set::MOI.Zeros,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)
    nvar = length(c)
    ncon = length(b)

    m = MOI.dimension(set)

    set_ = MOI.EqualTo(0.0)
    k = MOI.Zeros(1)

    # No need to add slack variables

    # Extract terms of VectorAffineFunction
    funs_ = [
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
        for i in 1:m
    ]
    for t in fun.terms
        i = t.output_index
        j = varindices[t.scalar_term.variable_index]
        v = t.scalar_term.coefficient
        push!(funs_[i].terms, MOI.ScalarAffineTerm(v, MOI.VariableIndex(j)))

        # Update coefficients of A
        push!(arows, ncon + i)
        push!(acols, j)
        push!(avals, v)
    end

    # Add each constraint separately
    for i in 1:m
        # We add each equality constraint separately
        rhs = -fun.constants[i]
        push!(b, rhs)

        # Add equality constraint
        MOI.add_constraint(dst, funs_[i], MOI.EqualTo(rhs))
    end

    return nothing
end


# TODO: Add option to use extended formulation for second order cone
"""

Convert vector non-linear conic constraint into standard form

This adds `m` variables, where `m` is the dimension of the conic constraint.
"""
function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.VectorAffineFunction,
    set::NONLINEAR_CONE,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)
    nvar = length(c)
    ncon = length(b)

    m = MOI.dimension(set)

    # Add slack variables
    s = MOI.add_variables(dst, m)
    append!(c, zeros(m))
    append!(int_flags, zeros(Bool, m))

    # Update matrix coefficients
    append!(arows, ncon .+ collect(1:m))
    append!(acols, nvar .+ collect(1:m))
    append!(avals, -ones(m))

    # Set conic domain on slacks
    MOI.add_constraint(dst, MOI.VectorOfVariables(s), set)
    push!(cones, (nvar .+ collect(1:m), set))

    # Extract terms of VectorAffineFunction
    funs_ = [
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(-1.0, s[i])], 0.0)
        for i in 1:m
    ]
    for t in fun.terms
        i = t.output_index
        j = varindices[t.scalar_term.variable_index]
        v = t.scalar_term.coefficient
        push!(funs_[i].terms, MOI.ScalarAffineTerm(v, MOI.VariableIndex(j)))

        # Update coefficients of A
        push!(arows, ncon + i)
        push!(acols, j)
        push!(avals, v)
    end

    # Add each constraint separately
    for i in 1:m
        # We add each equality constraint separately
        rhs = -fun.constants[i]
        push!(b, rhs)

        # Add equality constraint
        MOI.add_constraint(dst, funs_[i], MOI.EqualTo(rhs))
    end

    return nothing
end

function add_to_standard_form!(
    dst::MOI.ModelLike,
    fun::MOI.VectorOfVariables,
    set::NONLINEAR_CONE,
    varindices,
    c, b, arows, acols, avals, cones, int_flags
)
    m = MOI.dimension(set)

    x_ = [
        MOI.VariableIndex(varindices[v])
        for v in fun.variables
    ]

    # Add conic constraint
    MOI.add_constraint(dst, MOI.VectorOfVariables(x_), set)

    # Record cone
    push!(cones, ([varindices[v] for v in fun.variables], set))

    return nothing
end