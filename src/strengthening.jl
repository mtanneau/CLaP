import MathOptInterface
const MOI = MathOptInterface

"""
    strengthen(a, u0, v0, ::MOI.Nonnegatives)

Compute strenghtened coefficient for non-negative variable.

The strengthened coefficient is given by
```
α = min(-u0 * ⌊π⌋, a + v0 * ⌊π⌋) 
```
where ``π = -a / (u0 + v0)``
"""
function strengthen(a::Float64, u0::Float64, v0::Float64, k::MOI.Nonnegatives)
    
    # Sanity checks
    @assert MOI.dimension(k) == 1 "Trying to strengthen a $(typeof(k)) of dimension $(MOI.dimension(k)) (should be 1)"
    @assert u0 + v0 > 1e-4  "u0 + v0 = $(u0 + v0) (should be > 0)"  # Otherwise cut is essentially a K* cut

    π = -a / (u0 + v0)
    return min(
        -u0 * floor(π),
        a + v0 * ceil(π)
    )
end

"""
    strengthen(a, u0, v0, ::MOI.Nonpositive)

Compute strenghtened coefficient for non-positive variable.

The strengthened coefficient is given by
```
α = min(-u0 * ⌊π⌋, a + v0 * ⌊π⌋) 
```
where ``π = -a / (u0 + v0)``
"""
function strengthen(a::Float64, u0::Float64, v0::Float64, k::MOI.Nonpositives)
    
    # Sanity checks
    @assert MOI.dimension(k) == 1 "Trying to strengthen a $(typeof(k)) of dimension $(MOI.dimension(k)) (should be 1)"
    @assert u0 + v0 > 1e-4  # Otherwise cut is essentially a K* cut
    
    π = -a / (u0 + v0)
    return max(
        -u0 * ceil(π),
        a + v0 * floor(π)
    )
end

# TODO: strengthening for conic components