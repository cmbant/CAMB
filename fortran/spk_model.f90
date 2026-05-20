module SPkModel
use precision
implicit none
private

real(dl), parameter, public :: SPk_calibrated_z_min = 0.0_dl
real(dl), parameter, public :: SPk_calibrated_z_max = 3.0_dl
real(dl), parameter, public :: SPk_calibrated_k_min = 1e-12_dl
real(dl), parameter, public :: SPk_calibrated_k_max = 12.0_dl
real(dl), parameter, public :: SPk_min_suppression = 1e-6_dl

integer, parameter, public :: SPk_rel_power_law = 1
integer, parameter, public :: SPk_rel_cosmo_power_law = 2
integer, parameter, public :: SPk_rel_double_power_law = 3

public :: SPk_Suppression

contains

pure function SPk_Poly2(x, c0, c1, c2) result(y)
real(dl), intent(in) :: x, c0, c1, c2
real(dl) :: y

y = c2 * x * x + c1 * x + c0

end function SPk_Poly2

pure subroutine SPk_GetParams(SO, z, spk_a_out, spk_b_out, spk_g_out, lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c)
integer, intent(in) :: SO
real(dl), intent(in) :: z
real(dl), intent(out) :: spk_a_out, spk_b_out, spk_g_out, lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c
real(dl) :: x

x = 1.0_dl + z

if (SO == 500) then
    spk_a_out = SPk_Poly2(x, 14.783423122120318_dl, -0.999062404857228_dl, 0.12062854541689262_dl)
    spk_b_out = SPk_Poly2(x, 14.620528368613265_dl, -0.9136466201011957_dl, 0.10835389086945699_dl)
    spk_g_out = SPk_Poly2(x, 0.9671320682693298_dl, -0.03185388045484575_dl, 0.02650236152450093_dl)
    lambda_a = SPk_Poly2(x, 0.019349810078190303_dl, -0.007410668383424459_dl, 0.0008334762393555539_dl)
    lambda_b = SPk_Poly2(x, 2.9566773924238143_dl, 0.6205340408676114_dl, -0.001928273640110775_dl)
    mu_a = SPk_Poly2(x, 0.715853343781141_dl, -0.19276613600825665_dl, 0.04948240117059147_dl)
    mu_b = SPk_Poly2(x, 3.385355123440431_dl, 0.9658906605139421_dl, -0.06825861100375574_dl)
    mu_c = SPk_Poly2(x, 4.457257708010122_dl, -2.191853871334233_dl, 0.45457701107254733_dl)
    nu_a = SPk_Poly2(x, 478.86477329610375_dl, 429.88795783439946_dl, 249.25655627821902_dl)
    nu_b = SPk_Poly2(x, -11.227459319819815_dl, -0.5581080204509223_dl, 0.4489962047114509_dl)
    nu_c = SPk_Poly2(x, 3.499449440557995_dl, -0.08488559389068073_dl, -0.0923847866118189_dl)
else
    spk_a_out = SPk_Poly2(x, 15.24311120000861_dl, -1.2436699435560352_dl, 0.14837558774401766_dl)
    spk_b_out = SPk_Poly2(x, 14.969187892657688_dl, -1.0993025612653198_dl, 0.12905587245129102_dl)
    spk_g_out = SPk_Poly2(x, 0.8000441576980428_dl, -0.01715621131893159_dl, 0.06131887249968379_dl)
    lambda_a = SPk_Poly2(x, 0.02178116280689233_dl, -0.0077564325654746955_dl, 0.0007915576054589781_dl)
    lambda_b = SPk_Poly2(x, 3.0878286643613437_dl, 0.4529677646796634_dl, 0.001552571083240605_dl)
    mu_a = SPk_Poly2(x, 0.6930259177449359_dl, -0.16913553700233935_dl, 0.04263185199898842_dl)
    mu_b = SPk_Poly2(x, 3.161914061444856_dl, 0.8616834297321924_dl, 0.011346427353554053_dl)
    mu_c = SPk_Poly2(x, 5.532188503256583_dl, -3.0864672185252537_dl, 0.5083422518560442_dl)
    nu_a = SPk_Poly2(x, 413.00988701513904_dl, 311.63957063032285_dl, 37.89105940901369_dl)
    nu_b = SPk_Poly2(x, -11.243859405779181_dl, -0.34421412616421965_dl, 0.3343548325485801_dl)
    nu_c = SPk_Poly2(x, 3.476463891168505_dl, -0.018333059687988575_dl, -0.08276237963970698_dl)
end if

end subroutine SPk_GetParams

pure function SPk_Suppression(SO, kh, z, relation_kind, fb_a, fb_pow, fb_pivot, rel_alpha, rel_beta, rel_gamma, &
    rel_epsilon, rel_m_pivot, e_ratio) result(sup)
integer, intent(in) :: SO, relation_kind
real(dl), intent(in) :: kh, z, fb_a, fb_pow, fb_pivot
real(dl), intent(in) :: rel_alpha, rel_beta, rel_gamma, rel_epsilon, rel_m_pivot
real(dl), intent(in) :: e_ratio
real(dl) :: sup
real(dl) :: spk_a, spk_b, spk_g
real(dl) :: lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c
real(dl) :: best_mass, m_opt, fb, x, x0, x1, x2
real(dl) :: e_ratio_eff

call SPk_GetParams(SO, z, spk_a, spk_b, spk_g, lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c)

best_mass = spk_a - (spk_a - spk_b) * (kh ** spk_g)
m_opt = 10.0_dl ** best_mass
e_ratio_eff = max(e_ratio, SPk_calibrated_k_min)

select case (relation_kind)
case (SPk_rel_cosmo_power_law)
    fb = (exp(rel_alpha) / 100.0_dl) * ((m_opt / 1.0e14_dl) ** (rel_beta - 1.0_dl)) * (e_ratio_eff ** rel_gamma)
case (SPk_rel_double_power_law)
    fb = 0.5_dl * rel_epsilon * (((m_opt / rel_m_pivot) ** rel_alpha) + ((m_opt / rel_m_pivot) ** rel_beta)) * &
        (e_ratio_eff ** rel_gamma)
case default
    fb = fb_a * ((m_opt / fb_pivot) ** fb_pow)
end select

x = log10(kh)
x0 = 1.0_dl + lambda_a * exp(lambda_b * x)
x1 = mu_a + ((1.0_dl - mu_a) / (1.0_dl + exp(mu_b * x + mu_c)))
x2 = nu_a * exp(-0.5_dl * ((x - nu_b) / nu_c) ** 2)

sup = x0 - (x0 - x1) * exp(-x2 * fb)
sup = max(sup, SPk_min_suppression)

end function SPk_Suppression

end module SPkModel
