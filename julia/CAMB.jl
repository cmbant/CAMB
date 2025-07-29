module CAMB

export Cosmology,
       distance_modulus,
       lookback_time,
       age_at_z,
       age_universe,
       transfer_nw,
       transfer_eh,
       growth_factor,
       growth_rate,
       matter_power_spectrum,
       nonlinear_power_spectrum,
       comoving_radial_distance,
       angular_diameter_distance,
       luminosity_distance,
       hubble_parameter,
       recombination_history,
       angular_power_spectrum,
       sigma_R,
       sigma8,
       critical_density,
       hubble_distance,
       sound_horizon,
       drag_redshift,
       sound_horizon_drag,
       Omega_nu,
       Omega_de,
       Omega_m_z,
       Omega_de_z,
       Constants,
       Config,
       tensor_param_indeptilt,
       tensor_param_rpivot,
       tensor_param_AT,
       InitialPowerLaw,
       scalar_power,
       tensor_power,
       effective_ns,
       TransferParams,
       AccuracyParams,
       CAMBparams,
       default_params,
       BackgroundOutputs,
       CAMBdata,
       compute_background,
       update_background!,
       BackgroundCache,
       precompute_background,
       TanhReionization,
       electron_fraction,
       optical_depth,
       sph_besselj,
       sph_besselj_prime,
       integrate_romberg,
       brentq,
       gauss_legendre

include("constants.jl")
using .Constants
include("config.jl")
using .Config
include("initial_power.jl")
using .InitialPower
include("params.jl")
using .Params
include("results.jl")
using .Results
include("reionization.jl")
using .Reionization
include("bessels.jl")
using .Bessels
include("math_utils.jl")
using .MathUtils
include("transfer.jl")
using .Transfer
include("nonlinear.jl")
using .NonLinear
include("background_cache.jl")
using .BackgroundCache
using QuadGK
using LinearAlgebra
using ForwardDiff

const k_pivot = 0.05

"""
    Cosmology(; H0, Omega_b, Omega_c, Omega_k, ns, A_s, Tcmb,
              w0=-1.0, wa=0.0, Neff=3.046, mnu=0.0)

Container for cosmological parameters used by the Julia implementation.
`w0` and `wa` define the CPL dark energy equation-of-state
``w(a) = w0 + wa(1-a)``, and `Neff` specifies the effective number of
massless neutrino species used when computing the radiation density.
`mnu` gives the sum of neutrino masses in eV.
"""
struct Cosmology
    H0::Float64      # Hubble parameter [km/s/Mpc]
    Omega_b::Float64 # baryon density fraction today
    Omega_c::Float64 # cold dark matter density fraction today
    Omega_k::Float64 # curvature fraction today
    ns::Float64      # scalar spectral index
    A_s::Float64     # scalar amplitude at k_pivot
    Tcmb::Float64    # CMB temperature [K]
    w0::Float64      # dark energy equation of state parameter today
    wa::Float64      # dark energy evolution parameter
    Neff::Float64    # effective number of massless neutrino species
    mnu::Float64     # total mass of massive neutrinos [eV]
end

w(c::Cosmology, a::Float64) = c.w0 + c.wa * (1 - a)

"""
    Omega_nu(c)

Density fraction of massive neutrinos today assuming a total mass `mnu`
(eV) distributed equally among species.
"""
function Omega_nu(c::Cosmology)
    return Constants.neutrino_mass_fac * c.mnu / h(c)^2
end

Omega_m(c::Cosmology) = c.Omega_b + c.Omega_c + Omega_nu(c)
h(c::Cosmology) = c.H0 / 100.0

"""
    Omega_de(c)

Dark-energy density fraction today assuming a flat or curved universe.
"""
Omega_de(c::Cosmology) = 1 - Omega_m(c) - c.Omega_k - Omega_r(c)

"""
    Omega_m_z(c, z)

Matter density fraction at redshift `z`.
"""
function Omega_m_z(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    return Omega_m(c) / (a^3 * E(c, a)^2)
end

"""
    Omega_de_z(c, z)

Dark-energy density fraction at redshift `z`.
"""
function Omega_de_z(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    Ol0 = Omega_de(c)
    wa_term = exp(-3 * c.wa * (1 - a))
    return Ol0 * a^(-3 * (1 + c.w0 + c.wa)) * wa_term / E(c, a)^2
end

"""
    growth_rate(c, z)

Logarithmic derivative of the growth factor with respect to scale factor.
"""
function growth_rate(c::Cosmology, z::Float64)
    D = z -> growth_factor(c, z)
    dDdz = ForwardDiff.derivative(D, z)
    return -(1 + z) * dDdz / D(z)
end

const Tcmb_fid = 2.7255

"""
    Omega_gamma(c)

Physical photon density fraction today.
"""
Omega_gamma(c::Cosmology) = 2.472e-5 * (c.Tcmb / Tcmb_fid)^4 / h(c)^2

"""
    Omega_r(c)

Total radiation density fraction including photons and massless neutrinos.
"""
Omega_r(c::Cosmology) = Omega_gamma(c) * (1 + 0.2271 * c.Neff)

"""
    E(c, a)

Dimensionless Hubble parameter at scale factor `a` for a \u03bbCDM cosmology.
"""
function E(c::Cosmology, a::Float64)
    Or = Omega_r(c)
    Om = Omega_m(c)
    Ok = c.Omega_k
    Ol = 1 - Or - Om - Ok
    wa_term = exp(-3 * c.wa * (1 - a))
    return sqrt(Or / a^4 + Om / a^3 + Ok / a^2 +
                Ol * a^(-3 * (1 + c.w0 + c.wa)) * wa_term)
end

"""
    growth_factor(c, z)

Linear growth factor normalized to unity today.
"""
function growth_factor(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    integrand(a) = 1 / (a^3 * E(c, a)^3)
    integral, _ = quadgk(integrand, 0.0, a)
    d = 5 * Omega_m(c) * E(c, a) * integral
    norm, _ = quadgk(integrand, 0.0, 1.0)
    d_today = 5 * Omega_m(c) * E(c, 1.0) * norm
    return d / d_today
end


const ckms = Constants.c / 1000

"""
    hubble_parameter(c, z)

Hubble parameter in km/s/Mpc at redshift `z`.
"""
hubble_parameter(c::Cosmology, z::Float64) = c.H0 * E(c, 1 / (1 + z))

"""
    comoving_radial_distance(c, z)

Comoving radial distance from us to redshift `z` in Mpc.
"""
function comoving_radial_distance(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    integrand(a) = 1 / (a^2 * E(c, a))
    chi, _ = quadgk(integrand, a, 1.0)
    return (ckms / c.H0) * chi
end

"""
    angular_diameter_distance(c, z)

Angular diameter distance to redshift `z` in Mpc.
"""
function angular_diameter_distance(c::Cosmology, z::Float64)
    chi = comoving_radial_distance(c, z)
    Ok = c.Omega_k
    if abs(Ok) < 1e-8
        dm = chi
    else
        rcurv = (ckms / c.H0) / sqrt(abs(Ok))
        if Ok > 0
            dm = rcurv * sinh(chi / rcurv)
        else
            dm = rcurv * sin(chi / rcurv)
        end
    end
    return dm / (1 + z)
end

"""
    luminosity_distance(c, z)

Luminosity distance to redshift `z` in Mpc.
"""

luminosity_distance(c::Cosmology, z::Float64) = angular_diameter_distance(c, z) * (1 + z)^2

"""
    distance_modulus(c, z)

Distance modulus at redshift `z`.
"""
function distance_modulus(c::Cosmology, z::Float64)
    dl = luminosity_distance(c, z) * 1e6
    return 5 * (log10(dl) - 1)
end

"""
    lookback_time(c, z)

Lookback time in gigayears from redshift `z` to today.
"""
function lookback_time(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    integrand(a) = 1 / (a * E(c, a))
    t, _ = quadgk(integrand, a, 1.0)
    t_H = 9.778 / h(c)
    return t_H * t
end

"""
    age_at_z(c, z)

Age of the universe in gigayears at redshift `z`.
"""
function age_at_z(c::Cosmology, z::Float64)
    a = 1 / (1 + z)
    integrand(a) = 1 / (a * E(c, a))
    t, _ = quadgk(integrand, 0.0, a)
    t_H = 9.778 / h(c)
    return t_H * t
end

"""
    age_universe(c)

Age of the universe today in gigayears.
"""
age_universe(c::Cosmology) = age_at_z(c, 0.0)

"""
    sound_horizon(c; zdec=1059.0)

Comoving sound horizon at photon decoupling. `zdec` defaults to the Planck
2018 value for the redshift of decoupling.
"""
function sound_horizon(c::Cosmology; zdec::Float64=1059.0)
    a_dec = 1 / (1 + zdec)
    R0 = 3 * c.Omega_b / (4 * Omega_gamma(c))
    integrand(a) = 1 / sqrt(3 * (1 + R0 * a)) / (a^2 * E(c, a))
    r, _ = quadgk(integrand, 0.0, a_dec)
    return (ckms / c.H0) * r
end

"""
    drag_redshift(c)

Approximate redshift of the drag epoch using the Eisenstein & Hu (1998)
formula.
"""
function drag_redshift(c::Cosmology)
    wm = Omega_m(c) * h(c)^2
    wb = c.Omega_b * h(c)^2
    b1 = 0.313 * wm^(-0.419) * (1 + 0.607 * wm^0.674)
    b2 = 0.238 * wm^0.223
    return 1291 * wm^0.251 / (1 + 0.659 * wm^0.828) * (1 + b1 * wb^b2)
end

"""
    sound_horizon_drag(c)

Comoving sound horizon at the drag epoch.
"""
sound_horizon_drag(c::Cosmology) = sound_horizon(c; zdec=drag_redshift(c))

"""
    recombination_history(c)

Compute a simplified recombination history. Returns an array of redshifts
and the corresponding ionization fraction. This is a crude approximation
and serves as a placeholder for a full Boltzmann solver.
"""
function recombination_history(c::Cosmology)
    z = collect(0.0:10.0:2000.0)
    x_e = @. 1 / (1 + (z / 1100)^1.5)
    return z, x_e
end

"""
    angular_power_spectrum(c, \ellmax)

Compute an approximate CMB temperature power spectrum up to multipole
`\ellmax`. The calculation uses the matter power spectrum as a proxy and
does not reproduce the accuracy of the full CAMB code.
"""
function angular_power_spectrum(c::Cosmology, \ellmax::Integer)
    ells = collect(2:\ellmax)
    cl = similar(ells, Float64)
    for (i, \ell) in enumerate(ells)
        k = \ell / 14000  # heuristic projection
        cl[i] = matter_power_spectrum(c, k) / \ell^(\ell + 1)
    end
    return ells, cl
end

"""
    sigma_R(c, R; z=0, kmin=1e-4, kmax=1e2)

RMS linear density fluctuation in spheres of radius `R` (Mpc) at redshift `z`.
The integration limits `kmin` and `kmax` are in units of \(h\,\mathrm{Mpc}^{-1}\).
"""
function sigma_R(c::Cosmology, R::Float64; z::Float64=0.0, kmin::Float64=1e-4,
                  kmax::Float64=1e2)
    W(x) = 3 * (sin(x) - x * cos(x)) / x^3
    integrand(k) = k^2 * matter_power_spectrum(c, k; z=z) * W(k * R)^2
    val, _ = quadgk(integrand, kmin, kmax)
    return sqrt(val / (2 * const_pi^2))
end

"""
    sigma8(c; z=0)

Convenience wrapper for `sigma_R` with `R = 8 / h(c)` in Mpc.
"""
sigma8(c::Cosmology; z::Float64=0.0) = sigma_R(c, 8 / h(c); z=z)

"""
    critical_density(c; z=0.0)

Critical density in kg/m^3 at redshift `z` (defaults to today).
"""
function critical_density(c::Cosmology; z::Float64=0.0)
    H = hubble_parameter(c, z) * 1000 / Constants.Mpc
    return 3 * H^2 / Constants.kappa
end

"""
    hubble_distance(c; z=0.0)

Hubble distance in Mpc at redshift `z`.
"""
function hubble_distance(c::Cosmology; z::Float64=0.0)
    return ckms / hubble_parameter(c, z)
end

end # module
