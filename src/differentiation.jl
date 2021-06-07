export differentiate

"""
    differentiate(p::AbstractPolynomialLike, v::AbstractVariable, deg::Union{Int, Val}=1)

Differentiate `deg` times the polynomial `p` by the variable `v`.


    differentiate(p::AbstractPolynomialLike, vs, deg::Union{Int, Val}=1)

Differentiate `deg` times the polynomial `p` by the variables of the vector or
tuple of variable `vs` and return an array of dimension `deg`. It is recommended
to pass `deg` as a `Val` instance when the degree is known at compile time, e.g.
`differentiate(p, v, Val{2}())` instead of `differentiate(p, x, 2)`, as this
will help the compiler infer the return type.

    differentiate(p::AbstractArray{<:AbstractPolynomialLike, N}, vs, deg::Union{Int, Val}=1) where N

Differentiate the polynomials in `p` by the variables of the vector or tuple of variable `vs` and return an array of dimension `N+deg`.
If `p` is an `AbstractVector` this returns the Jacobian of `p` where the i-th row containts the partial
derivaties of `p[i]`.

### Examples

```julia
p = 3x^2*y + x + 2y + 1
differentiate(p, x) # should return 6xy + 1
differentiate(p, x, Val{1}()) # equivalent to the above
differentiate(p, (x, y)) # should return [6xy+1, 3x^2+1]
differentiate( [x^2+y, z^2+4x], [x, y, z]) # should return [2x 1 0; 4 0 2z]
```
"""
function differentiate end

# Fallback for everything else
differentiate(α::T, v::AbstractVariable) where T = zero(T)
differentiate(v1::AbstractVariable, v2::AbstractVariable) = v1 == v2 ? 1 : 0
differentiate(t::AbstractTermLike, v::AbstractVariable) = coefficient(t) * differentiate(monomial(t), v)
# The polynomial function will take care of removing the zeros
function differentiate(p::APL, v::AbstractVariable)
    if iszero(p)
        # As `terms(p)` is empty, `differentiate.(terms(p), v)` gives `Any[]`.
        T = typeof(differentiate(one(termtype(p)), v))
        polynomial!(T[], SortedUniqState())
    else
        polynomial!(differentiate.(terms(p), v), SortedState())
    end
end
differentiate(p::RationalPoly, v::AbstractVariable) = (differentiate(p.num, v) * p.den - p.num * differentiate(p.den, v)) / p.den^2

const ARPL = Union{APL, RationalPoly}

function differentiate(ps::AbstractArray{PT}, xs::AbstractArray) where {PT <: ARPL}
    differentiate.(reshape(ps, (size(ps)..., 1)), reshape(xs, 1, :))
end

function differentiate(ps::AbstractArray{PT}, xs::Tuple) where {PT <: ARPL}
    differentiate(ps, collect(xs))
end


# TODO: this signature is probably too wide and creates the potential
# for stack overflows
differentiate(p::ARPL, xs) = [differentiate(p, x) for x in xs]

# differentiate(p, [x, y]) with TypedPolynomials promote x to a Monomial
differentiate(p::ARPL, m::AbstractMonomial) = differentiate(p, variable(m))

# The `R` argument indicates a desired result type. We use this in order
# to attempt to preserve type-stability even though the value of `deg` cannot
# be known at compile time. For scalar `p` and `x`, we set R to be the type
# of differentiate(p, x) to give a stable result type regardless of `deg`. For
# vectors p and/or x this is impossible (since differentiate may return an array),
# so we just set `R` to `Any`
function (_differentiate_recursive(p, x, deg::Int, ::Type{R})::R) where {R}
    if deg < 0
        throw(DomainError(deg, "degree is negative"))
    elseif deg == 0
        return p
    else
        return differentiate(differentiate(p, x), x, deg-1)
    end
end

differentiate(p, x, deg::Int) = _differentiate_recursive(p, x, deg, Base.promote_op(differentiate, typeof(p), typeof(x)))
differentiate(p::AbstractArray, x,                              deg::Int) = _differentiate_recursive(p, x, deg, Any)
differentiate(p,                x::Union{AbstractArray, Tuple}, deg::Int) = _differentiate_recursive(p, x, deg, Any)
differentiate(p::AbstractArray, x::Union{AbstractArray, Tuple}, deg::Int) = _differentiate_recursive(p, x, deg, Any)


# This is alternative, Val-based interface for nested differentiation.
# It has the advantage of not requiring an conversion or calls to
# Base.promote_op, while maintaining type stability for any argument
# type.
differentiate(p, x, ::Val{0}) = p
differentiate(p, x, ::Val{1}) = differentiate(p, x)

@static if VERSION < v"v0.7.0-"
    # Marking this @pure helps julia v0.6 figure this out
    Base.@pure _reduce_degree(::Val{N}) where {N} = Val{N - 1}()
    function differentiate(p, x, deg::Val{N}) where N
        if N < 0
            throw(DomainError(deg))
        else
            differentiate(differentiate(p, x), x, _reduce_degree(deg))
        end
    end
else
    # In Julia v0.7 and above, we can remove the _reduce_degree trick
    function differentiate(p, x, deg::Val{N}) where N
        if N < 0
            throw(DomainError(deg))
        else
            differentiate(differentiate(p, x), x, Val{N - 1}())
        end
    end
end
