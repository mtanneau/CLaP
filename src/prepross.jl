using SparseArrays

is_int(x::Float64) = x == floor(x)

"""
    extract_implied_integer(sf)

Detect implied integer variables.

Given a problem
```
min     c'x
s.t.    A * x = b
        x ∈ K
```
detects all components of `x` that must be integer.

A variable `xj` is implied integer if it appears in a constraint
```
± x_j0 = b_i - Σ a_i,j x_j
```
with ``x_j``, ``j ≠ j0`` (implied) integer and ``b_i, a_i,j`` integers.
"""
function extract_implied_integer(sf::StandardProblem)
    m, n = size(sf.A)
    A = sf.A
    b = sf.b

    At = sparse(sf.A')

    rhs_integer_flag = is_int.(sf.b)
    
    nz_per_row = zeros(Int, m)
    n_intvar_per_row = zeros(Int, m)
    integer_row_flags = copy(rhs_integer_flag)

    arows, acols, avals = findnz(sf.A)

    # Initialize number of integer 
    for (i, j, v) in zip(arows, acols, avals)
        nz_per_row[i] += 1
        n_intvar_per_row[i] += sf.vartypes_implied[j]
        integer_row_flags[i] &= is_int(v)
    end

    has_changed = true
    while has_changed
        has_changed = false

        # Go through rows
        for i in 1:m
            integer_row_flags[i] || continue  # Skip non-integer rows

            # Go through coefficients
            row = At[:, i]

            j0 = 0
            has_implied_integer = true
            for (j, v) in zip(row.nzind, row.nzval)
                # Check if row writes ± x0 = b - Σ aj xj
                # where b, aj, are integers, and xj, j ≠ j0 are integer variables

                if sf.vartypes_implied[j]
                    # Check if coefficient is integer
                    if !is_int(v)
                        # Not good
                        has_implied_integer = false
                        break
                    end
                else
                    # This variable is fractional
                    if j0 == 0 && (abs(v) == 1.0)
                        j0 = j
                    else
                        # Two fractional variables, or wrong coefficient
                        has_implied_integer = false
                        break
                    end
                end

            end

            has_implied_integer && (j0 != 0) || continue

            # Flag variable as implied integer
            sf.vartypes_implied[j0] = true
            # @info "Variable $(j0) flagged as integer"
            has_changed = true
        end
    end

    @info "Model pre-processing" m n sum(integer_row_flags) sum(sf.vartypes_implied)

    return sf
end