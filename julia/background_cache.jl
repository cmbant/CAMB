module BackgroundCache

export BackgroundCache, precompute_background

using Interpolations
using ..Constants
using ..: Cosmology, E, hubble_parameter, comoving_radial_distance,
            angular_diameter_distance, luminosity_distance,
            distance_modulus, lookback_time, age_at_z, age_universe

struct BackgroundCache{I1,I2,I3}
    cosmo::Cosmology
    H_itp::I1
    chi_itp::I2
    t_itp::I3
end

function precompute_background(c::Cosmology; a_min=1e-4, n=800)
    a_vals = collect(range(a_min, 1.0; length=n))
    H_vals = similar(a_vals)
    chi_vals = similar(a_vals)
    t_vals = similar(a_vals)

    integrand_chi(x) = 1 / (x^2 * E(c, x))
    integrand_t(x) = 1 / (x * E(c, x))

    chi_vals[1] = quadgk(integrand_chi, 0.0, a_vals[1])[1]
    t_vals[1] = quadgk(integrand_t, 0.0, a_vals[1])[1]
    H_vals[1] = hubble_parameter(c, 1 / a_vals[1] - 1)

    for i in 2:n
        a0 = a_vals[i-1]
        a1 = a_vals[i]
        H_vals[i] = hubble_parameter(c, 1 / a1 - 1)
        chi_vals[i] = chi_vals[i-1] + quadgk(integrand_chi, a0, a1)[1]
        t_vals[i] = t_vals[i-1] + quadgk(integrand_t, a0, a1)[1]
    end

    H_itp = CubicSplineInterpolation(a_vals, H_vals; extrapolation_bc=Line())
    chi_itp = CubicSplineInterpolation(a_vals, chi_vals; extrapolation_bc=Line())
    t_itp = CubicSplineInterpolation(a_vals, t_vals; extrapolation_bc=Line())

    return BackgroundCache{typeof(H_itp), typeof(chi_itp), typeof(t_itp)}(c, H_itp, chi_itp, t_itp)
end

function hubble_parameter(bc::BackgroundCache, z::Float64)
    a = 1 / (1 + z)
    return bc.H_itp(a)
end

function comoving_radial_distance(bc::BackgroundCache, z::Float64)
    a = 1 / (1 + z)
    return (Constants.c / 1000 / bc.cosmo.H0) * bc.chi_itp(a)
end

function angular_diameter_distance(bc::BackgroundCache, z::Float64)
    chi = comoving_radial_distance(bc, z)
    Ok = bc.cosmo.Omega_k
    if abs(Ok) < 1e-8
        dm = chi
    else
        rcurv = (Constants.c / 1000 / bc.cosmo.H0) / sqrt(abs(Ok))
        if Ok > 0
            dm = rcurv * sinh(chi / rcurv)
        else
            dm = rcurv * sin(chi / rcurv)
        end
    end
    return dm / (1 + z)
end

luminosity_distance(bc::BackgroundCache, z::Float64) = angular_diameter_distance(bc, z) * (1 + z)^2

distance_modulus(bc::BackgroundCache, z::Float64) = distance_modulus(bc.cosmo, z)

function lookback_time(bc::BackgroundCache, z::Float64)
    a = 1 / (1 + z)
    t_H = 9.778 / (bc.cosmo.H0 / 100.0)
    return t_H * (bc.t_itp(1.0) - bc.t_itp(a))
end

function age_at_z(bc::BackgroundCache, z::Float64)
    a = 1 / (1 + z)
    t_H = 9.778 / (bc.cosmo.H0 / 100.0)
    return t_H * bc.t_itp(a)
end

age_universe(bc::BackgroundCache) = age_at_z(bc, 0.0)

end # module
