using CAMB
using Plots

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

z = range(0, 3; length=200)
chi = [comoving_radial_distance(c, zi) for zi in z]
da = [angular_diameter_distance(c, zi) for zi in z]
dl = [luminosity_distance(c, zi) for zi in z]
mu = [distance_modulus(c, zi) for zi in z]
lb = [lookback_time(c, zi) for zi in z]
Hvals = [hubble_parameter(c, zi) for zi in z]
D = [growth_factor(c, zi) for zi in z]

plot(z, Hvals, xlabel="z", ylabel="H(z) [km/s/Mpc]", label="H(z)")
savefig("H_z.png")

plot(z, chi, xlabel="z", ylabel="Distance [Mpc]", label="Comoving")
plot!(z, da, label="Angular diameter")
plot!(z, dl, label="Luminosity")
savefig("distances.png")

plot(z, mu, xlabel="z", ylabel="μ(z)", label="distance modulus")
savefig("distance_modulus.png")

plot(z, lb, xlabel="z", ylabel="Lookback time [Gyr]", label="t_L")
savefig("lookback_time.png")

plot(z, D, xlabel="z", ylabel="Growth factor", label="D(z)")
savefig("growth_factor.png")

k = range(1e-3, 10; length=200)
Pk = [matter_power_spectrum(c, ki) for ki in k]
Pk_nl = [nonlinear_power_spectrum(c, ki) for ki in k]
plot(k, Pk, xlabel="k [h/Mpc]", ylabel="P(k) [(Mpc/h)^3]", xscale=:log10, yscale=:log10, label="linear")
plot!(k, Pk_nl, label="nonlinear")
savefig("power_spectrum.png")

zrec, xe = recombination_history(c)
plot(zrec, xe, xlabel="z", ylabel="x_e(z)", label="e fraction")
savefig("recombination.png")

ells, cl = angular_power_spectrum(c, 1000)
plot(ells, cl, xlabel="ℓ", ylabel="C_ℓ", xscale=:log10, yscale=:log10, label="TT approx")
savefig("angular_power.png")

println("Plots saved to current directory")
