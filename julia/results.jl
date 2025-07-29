module Results

export BackgroundOutputs, CAMBdata, compute_background, update_background!,
       lookback_time, age_at_z, age_universe,
       drag_redshift, sound_horizon_drag,
       growth_rate, Omega_de, Omega_m_z, Omega_de_z,
       electron_fraction, optical_depth, critical_density,
       hubble_distance,
       nonlinear_power_spectrum

using ..CAMB: Cosmology, comoving_radial_distance, angular_diameter_distance,
                luminosity_distance, hubble_parameter, distance_modulus,
                lookback_time, age_at_z, age_universe,
                sigma_R, sigma8, sound_horizon,
                drag_redshift, sound_horizon_drag,
                growth_rate, Omega_de, Omega_m_z, Omega_de_z,
                electron_fraction, optical_depth, critical_density, hubble_distance,
                nonlinear_power_spectrum

struct BackgroundOutputs
    H::Vector{Float64}
    DA::Vector{Float64}
    rs_by_D_v::Vector{Float64}
end

struct CAMBdata
    cosmology::Cosmology
    transfer_times::Vector{Float64}
    sigma_8::Vector{Float64}
    background::BackgroundOutputs
end

function CAMBdata(c::Cosmology)
    CAMBdata(c, Float64[], Float64[], BackgroundOutputs(Float64[], Float64[], Float64[]))
end

function compute_background(c::Cosmology, a_vals::Vector{Float64})
    Hvals = similar(a_vals)
    DAvals = similar(a_vals)
    rsDv = similar(a_vals)
    @inbounds @simd for i in eachindex(a_vals)
        a = a_vals[i]
        z = 1 / a - 1
        Hvals[i] = hubble_parameter(c, z)
        DAvals[i] = angular_diameter_distance(c, z)
        rsDv[i] = comoving_radial_distance(c, z) / sqrt(1 + z)
    end
    return BackgroundOutputs(Hvals, DAvals, rsDv)
end

function update_background!(data::CAMBdata, a_vals::Vector{Float64})
    data.background = compute_background(data.cosmology, a_vals)
    return data.background
end

hubble_parameter(data::CAMBdata, z::Float64) = hubble_parameter(data.cosmology, z)
comoving_radial_distance(data::CAMBdata, z::Float64) = comoving_radial_distance(data.cosmology, z)
angular_diameter_distance(data::CAMBdata, z::Float64) = angular_diameter_distance(data.cosmology, z)
luminosity_distance(data::CAMBdata, z::Float64) = luminosity_distance(data.cosmology, z)
distance_modulus(data::CAMBdata, z::Float64) = distance_modulus(data.cosmology, z)
lookback_time(data::CAMBdata, z::Float64) = lookback_time(data.cosmology, z)
age_at_z(data::CAMBdata, z::Float64) = age_at_z(data.cosmology, z)
age_universe(data::CAMBdata) = age_universe(data.cosmology)
sigma_R(data::CAMBdata, R; z=0.0, kmin=1e-4, kmax=1e2) = sigma_R(data.cosmology, R; z=z, kmin=kmin, kmax=kmax)
sigma8(data::CAMBdata; z=0.0) = sigma8(data.cosmology; z=z)
sound_horizon(data::CAMBdata; zdec=1059.0) = sound_horizon(data.cosmology; zdec=zdec)
drag_redshift(data::CAMBdata) = drag_redshift(data.cosmology)
sound_horizon_drag(data::CAMBdata) = sound_horizon_drag(data.cosmology)
growth_rate(data::CAMBdata, z::Float64) = growth_rate(data.cosmology, z)
Omega_de(data::CAMBdata) = Omega_de(data.cosmology)
Omega_m_z(data::CAMBdata, z::Float64) = Omega_m_z(data.cosmology, z)
Omega_de_z(data::CAMBdata, z::Float64) = Omega_de_z(data.cosmology, z)
electron_fraction(data::CAMBdata, z::Float64, model) =
    electron_fraction(model, z)
optical_depth(data::CAMBdata, model; zmax=40.0) =
    optical_depth(data.cosmology, model; zmax=zmax)
hubble_distance(data::CAMBdata; z=0.0) = hubble_distance(data.cosmology; z=z)
critical_density(data::CAMBdata; z=0.0) = critical_density(data.cosmology; z=z)
nonlinear_power_spectrum(data::CAMBdata, k; z=0.0) =
    NonLinear.nonlinear_power_spectrum(data.cosmology, k; z=z)

end # module
