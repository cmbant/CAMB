module Reionization

export TanhReionization, electron_fraction, optical_depth

using ..Constants
using ..: Cosmology, E
using QuadGK

"""
    TanhReionization(; z_re=8.0, delta_z=0.5, fraction=1.08)

Simple hyperbolic tangent parameterization of the reionization history.
`z_re` sets the mid-point redshift and `delta_z` controls the width of the
transition. `fraction` is the asymptotic ionized fraction after reionization.
"""
struct TanhReionization
    z_re::Float64
    delta_z::Float64
    fraction::Float64
end

TanhReionization(; z_re=8.0, delta_z=0.5, fraction=1.08) =
    TanhReionization(z_re, delta_z, fraction)

"""
    electron_fraction(model, z)

Ionization fraction at redshift `z`.
"""
function electron_fraction(m::TanhReionization, z::Float64)
    x_e = m.fraction / 2 * (1 + tanh((m.z_re - z) / m.delta_z))
    return max(0.0, min(m.fraction, x_e))
end

"""
    optical_depth(c, model; zmax=40.0)

Approximate Thomson scattering optical depth due to reionization
for cosmology `c` and reionization model `model`.
"""
function optical_depth(c::Cosmology, m::TanhReionization; zmax::Float64=40.0)
    integrand(z) = electron_fraction(m, z) / ((1 + z) * E(c, 1 / (1 + z)))
    tau, _ = quadgk(integrand, 0.0, zmax)
    return tau
end

end # module
