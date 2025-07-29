using Test


push!(LOAD_PATH, "../src")

using CAMB

c = Cosmology(
    H0=67.5,
    Omega_b=0.049,
    Omega_c=0.264,
    Omega_k=0.0,
    ns=0.965,
    A_s=2.1e-9,
    Tcmb=2.7255,
    w0=-1.0,
    wa=0.0,
    Neff=3.046,
    mnu=0.06,
)

@testset "Basic cosmology" begin
    @test hubble_parameter(c, 0.0) ≈ c.H0
    @test Omega_m_z(c, 0.0) ≈ CAMB.Omega_m(c)
    @test growth_factor(c, 0.0) ≈ 1.0 atol=1e-6
    @test age_universe(c) > 10
end

@testset "Distances" begin
    @test comoving_radial_distance(c, 0.5) > 0
    @test angular_diameter_distance(c, 0.5) > 0
    @test luminosity_distance(c, 0.5) > 0
    @test distance_modulus(c, 1.0) > 0
end

@testset "Power spectrum" begin
    k = 0.1
    @test matter_power_spectrum(c, k) > 0
    @test nonlinear_power_spectrum(c, k) >= matter_power_spectrum(c, k)
    @test sigma8(c) > 0
end

@testset "Misc" begin
    @test sound_horizon(c) > 0
    @test drag_redshift(c) > 0
    @test critical_density(c) > 0
    @test hubble_distance(c) > 0
end
