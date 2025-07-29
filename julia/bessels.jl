module Bessels

export sph_besselj, sph_besselj_prime

using SpecialFunctions

"""
    sph_besselj(l, x)

Compute the spherical Bessel function of the first kind of order `l` at `x` using
`SpecialFunctions.sph_besselj` for accuracy.
"""
sph_besselj(l::Integer, x::Real) = SpecialFunctions.sph_besselj(l, x)

"""
    sph_besselj_prime(l, x)

Derivative of the spherical Bessel function of the first kind with respect to
`x`. Uses the recurrence relation
`j_l'(x) = j_{l-1}(x) - (l + 1) / x * j_l(x)` for efficient evaluation.
"""
function sph_besselj_prime(l::Integer, x::Real)
    if l == 0
        return SpecialFunctions.sph_besselj(1, x) - 1 / x * SpecialFunctions.sph_besselj(0, x)
    else
        return SpecialFunctions.sph_besselj(l - 1, x) - (l + 1) / x * SpecialFunctions.sph_besselj(l, x)
    end
end

end # module
