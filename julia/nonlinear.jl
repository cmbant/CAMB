module NonLinear

export nonlinear_power_spectrum

using ..Constants: const_pi
using ..: Cosmology, matter_power_spectrum, sigma_R
using ..MathUtils: brentq

"""
    nonlinear_power_spectrum(c, k; z=0.0)

Approximate nonlinear matter power spectrum using a simplified Halofit-like
fitting formula. This implementation estimates the scale where
`\sigma(R=1/k_sigma,z) = 1` and boosts the linear spectrum accordingly.
"""
function nonlinear_power_spectrum(c::Cosmology, k::Float64; z::Float64=0.0)
    f(x) = sigma_R(c, 1/x; z=z) - 1
    k_sigma = brentq(f, 1e-3, 1e1)
    delta_lin = matter_power_spectrum(c, k; z=z) * k^3 / (2 * const_pi^2)
    a = 0.482
    alpha = 3.0
    delta_nl = delta_lin * (1 + (k / k_sigma)^alpha)^a
    return delta_nl / k^3 * 2 * const_pi^2
end

end # module
