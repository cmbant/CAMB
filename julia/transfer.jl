module Transfer

export transfer_nw, transfer_eh, matter_power_spectrum

using ..: Cosmology, Omega_m, h
using ..Constants: const_pi

"""
    transfer_nw(c, k)

Eisenstein & Hu (1998) no-wiggle transfer function. `k` is in units of
\(h\,\mathrm{Mpc}^{-1}\).
"""
function transfer_nw(c::Cosmology, k::Float64)
    wm = Omega_m(c) * h(c)^2
    wb = c.Omega_b * h(c)^2
    rb = c.Omega_b / Omega_m(c)
    s = 44.5 * log(9.83 / wm) / sqrt(1 + 10 * wb^0.75)
    alpha = 1 - 0.328 * log(431 * wm) * rb + 0.38 * log(22.3 * wm) * rb^2
    gamma = Omega_m(c) * h(c) * (alpha + (1 - alpha) / (1 + (0.43 * k * s * h(c))^4))
    q = k * (c.Tcmb / 2.7)^2 / gamma
    L = log(2 * exp(1) + 1.8 * q)
    C = 14.2 + 731 / (1 + 62.5 * q)
    return L / (L + C * q^2)
end

"""
    transfer_eh(c, k)

Full Eisenstein & Hu (1998) transfer function including baryon oscillations.
This is still an approximation and does not capture all effects of massive
neutrinos or detailed recombination physics.
"""
function transfer_eh(c::Cosmology, k::Float64)
    wm = Omega_m(c) * h(c)^2
    wb = c.Omega_b * h(c)^2
    f_b = c.Omega_b / Omega_m(c)
    theta = c.Tcmb / 2.7

    zeq = 2.50e4 * wm / theta^4
    Req = 31.5 * wb / theta^4 / zeq
    b1 = 0.313 * wm^(-0.419) * (1 + 0.607 * wm^0.674)
    b2 = 0.238 * wm^0.223
    zd = 1291 * wm^0.251 / (1 + 0.659 * wm^0.828) * (1 + b1 * wb^b2)
    Rd = 31.5 * wb / theta^4 / zd
    keq = 7.46e-2 * wm / theta^2
    s = 2 / (3 * keq) * sqrt(6 / Req) * log((sqrt(1 + Rd) + sqrt(Rd + Req)) / (1 + sqrt(Req)))

    alpha_c = (46.9 * wm)^(0.670) * (1 + (32.1 * wm)^(-0.532))
    alpha_c = alpha_c^(-f_b)
    beta_c = 1 / (1 + 0.944 * f_b) + (1 - 1 / (1 + 0.944 * f_b)) / (1 + (459 * wm)^(-1))

    k_silk = 1.6 * wb^0.52 * wm^0.73 * (1 + (10.4 * wm)^(-0.95))

    alpha_b = 2.07 * keq * s * (1 + Rd)^(-0.75)
    beta_b = 0.5 + f_b + (3 - 2f_b) * sqrt((17.2 * wm)^2 + 1)

    q = k / (13.41 * keq)
    T_c_ln = log(exp(1) + 1.8 * beta_c * q)
    T_c_ln_beta = T_c_ln / (T_c_ln + (14.2 / alpha_c + 386 / (1 + 69.9 * q^1.08)) * q^2)
    T_c = T_c_ln_beta / (1 + (k * s / 5.4)^4)

    s_tilde = s / (1 + (beta_b / (k * s))^3)^(1/3)
    T_b = (T_c_ln_beta / (1 + (k * s / 5.2)^2) + alpha_b / (1 + (beta_b / (k * s))^3) * exp(-(k / k_silk)^1.4)) * sin(k * s_tilde) / (k * s_tilde)

    return f_b * T_b + (1 - f_b) * T_c
end

"""
    matter_power_spectrum(c, k; z=0, wiggles=true)

Linear matter power spectrum using the Eisenstein & Hu transfer function.
Set `wiggles=false` to use the no-wiggle approximation.
"""
function matter_power_spectrum(c::Cosmology, k::Float64; z::Float64=0.0, wiggles::Bool=true)
    g = CAMB.growth_factor(c, z)
    T = wiggles ? transfer_eh(c, k) : transfer_nw(c, k)
    return c.A_s * (k / CAMB.k_pivot)^(c.ns - 1) * T^2 * g^2
end

end # module
