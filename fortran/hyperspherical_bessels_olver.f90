    module HypersphericalBesselOlver
    ! Approximate calculation of ultraspherical Bessel functions for non-flat universe
    ! Uses the Olver approximation to relate to normal spherical bessels,
    ! with fallback where not reliable to next order Olver/Airy approx or recursion.
    ! Precision target 1e-4 of peak, with max error < 2e-4.
    use precision
    use constants, only: const_pi
    use FlatBessels, only: BJL
    use SpherBessels, only: phi_recurs
    use HypersphericalBesselUtils, only: normalize_chi, turning_point, curved_radius, qintegral_exact, CACHE_EPS
    use HypersphericalBesselAiry, only: airy_u_normalized, airy_ok
    use HypersphericalBesselSmallNu, only: open_smallnu_u
    implicit none
    private

    ! Pointwise raw-Olver gate calibrated to keep peak-normalized errors below
    ! about 1e-4 on the open/closed validation grid.
    real(dl), parameter :: OLVER_OPEN_ALPHA_FLOOR = 0.095_dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_JOIN = 0.12_dl
    real(dl), parameter :: OLVER_OPEN_L_JOIN = 500._dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_LOW_EXP = 0.70_dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_HIGH_EXP = 0.14_dl
    real(dl), parameter :: OLVER_GATE_OPEN_EPS = 3.0e-3_dl
    real(dl), parameter :: OLVER_GATE_CLOSED_EPS = 7.0e-3_dl
    real(dl), parameter :: SMALLCHI_GATE_ALPHA = 2.5_dl
    integer, parameter :: SMALLCHI_GATE_LOW_ALPHA_L_MIN = 50
    real(dl), parameter :: SMALLCHI_GATE_METRIC = 5.0e-2_dl
    integer, parameter :: AIRY_FALLBACK_OPEN_L_MIN = 100
    integer, parameter :: AIRY_FALLBACK_CLOSED_L_MIN = 150

    public :: phi_olver, u_olver

    contains

    function phi_olver(l, K, nu, chi) result(phi)
    ! phi_l^nu(chi).  Assumes l >= 0, chi >= 0, and for K=1 an integer nu > l
    ! (so nu > sqrt(l(l+1)) and the turning point asin(sqrt(l(l+1))/nu) exists);
    ! the Python wrapper validates both.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    real(dl) :: phi

    phi = olver_value(l, K, nu, chi, reduced=.false.)

    end function phi_olver

    function u_olver(l, K, nu, chi) result(u)
    ! Reduced function u = S_K(chi) phi_l^nu(chi).  Same input assumptions as phi_olver.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    real(dl) :: u

    u = olver_value(l, K, nu, chi, reduced=.true.)

    end function u_olver

    function olver_value(l, K, nu, chi, reduced) result(val)
    ! Shared driver for phi_olver/u_olver: fold chi into the fundamental domain,
    ! dispatch the special cases, and convert between phi and u using the single
    ! S_K(achi) evaluation that the Olver map needs anyway.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    logical, intent(in) :: reduced
    real(dl) :: val
    real(dl) :: achi, symm, sin_k, j_l

    call normalize_chi(l, K, nu, chi, achi, symm)

    if (K == 0) then
        call BJL(l, nu*achi, j_l)
        if (reduced) then
            val = achi*j_l
        else
            val = j_l
        end if
        return
    end if

    if (l <= 2) then
        val = symm*phi_recurs(l, K, nu, achi)
        if (reduced) val = val*curved_radius(K, achi)
        return
    end if

    if (achi <= CACHE_EPS) then
        val = 0._dl
        return
    end if

    sin_k = curved_radius(K, achi)
    val = olver_reduced(l, K, nu, achi, sin_k, symm)
    if (.not. reduced) val = val/sin_k

    end function olver_value

    function olver_reduced(l, K, nu, achi, sin_k, symm) result(u)
    ! Reduced u for K = +-1 and l >= 3: pick the near-flat small-chi map, the
    ! recurrence-class fallback, or the leading Olver map, then map the flat
    ! solution z j_l(nu z) back with the Liouville-Green amplitude.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, achi, sin_k, symm
    real(dl) :: u
    real(dl) :: alpha_gate, z, amp, j_l

    alpha_gate = nu/real(l, dl)

    if (use_smallchi_map(l, nu, achi, alpha_gate)) then
        call compute_olver_z_amp_smallchi(l, K, nu, achi, z, amp)
        if (amp <= 0._dl) then
            u = 0._dl
            return
        end if
    else
        ! Empirical recurrence fallback for corners where the leading Olver map
        ! is not pointwise reliable at the 1e-4 peak-normalized level.
        ! The raw open-space approximation is accepted over the full sampled
        ! interval when alpha ~= nu/l clears a smooth, grid-calibrated cutoff
        ! in l. Below that cutoff, and for closed-space endpoint tails, the
        ! accepted region is limited by the small endpoint parameter
        !
        !   open:   eps_peak = O(chi / (2 nu))
        !   closed: eps_peak = O(chi / (2 (nu-l))).
        !
        ! Constants are grid-calibrated against phi_recurs, including low-l
        ! cases, high-l cases to l=10000, dense turning samples, very small chi,
        ! and open-space tails about 80 oscillations past the turning region.
        !
        ! Both metrics are cheap, so they are tested before open_alpha_cut,
        ! which costs a real-exponent power.
        if (K == 1) then
            ! nu is an integer mode > l for K=1, so the denominator is >= 2.
            if (achi/(2._dl*(nu - real(l, dl))) > OLVER_GATE_CLOSED_EPS) then
                u = fallback_reduced(l, K, nu, achi, sin_k, symm)
                return
            end if
        else if (K == -1) then
            if (achi/(2._dl*max(nu, tiny(1._dl))) > OLVER_GATE_OPEN_EPS &
                .and. alpha_gate < open_alpha_cut(l)) then
                u = fallback_reduced(l, K, nu, achi, sin_k, symm)
                return
            end if
        end if

        call compute_olver_z_amp(l, K, nu, achi, sin_k, z, amp)
    end if

    call BJL(l, nu*z, j_l)
    u = symm*amp*z*j_l

    end function olver_reduced

    function fallback_reduced(l, K, nu, achi, sin_k, symm) result(u)
    ! The Olver gates above are unchanged; within their recursive fallback
    ! region, use validated faster approximations where available before
    ! falling back to phi_recurs.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, achi, sin_k, symm
    real(dl) :: u
    logical :: ok

    if (use_airy_fallback(l, K, nu, achi)) then
        u = airy_u_normalized(l, K, nu, achi, ok)
        if (ok) then
            u = symm*u
            return
        end if
    end if

    if (K == -1) then
        u = open_smallnu_u(l, nu, achi, ok)
        if (ok) then
            u = symm*u
            return
        end if
    end if

    u = symm*phi_recurs(l, K, nu, achi)*sin_k

    end function fallback_reduced

    logical function use_airy_fallback(l, K, nu, achi) result(use_airy)
    ! Second-order Olver/Airy patch: validated only well away from low l.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, achi

    use_airy = .false.
    if (.not. airy_ok(l, K, nu, achi)) return

    select case (K)
    case (-1)
        use_airy = l >= AIRY_FALLBACK_OPEN_L_MIN
    case (1)
        use_airy = l >= AIRY_FALLBACK_CLOSED_L_MIN
    end select

    end function use_airy_fallback

    elemental real(dl) function open_alpha_cut(l) result(alpha_cut)
    ! Smooth grid-calibrated cutoff in alpha = nu/l below which the raw open-space
    ! Olver map is not trusted pointwise; a broken power law in l with a floor.
    integer, intent(in) :: l
    real(dl) :: ell

    ell = real(l, dl)
    if (ell < OLVER_OPEN_L_JOIN) then
        alpha_cut = OLVER_OPEN_ALPHA_JOIN*(OLVER_OPEN_L_JOIN/ell)**OLVER_OPEN_ALPHA_LOW_EXP
    else
        alpha_cut = OLVER_OPEN_ALPHA_JOIN*(OLVER_OPEN_L_JOIN/ell)**OLVER_OPEN_ALPHA_HIGH_EXP
    end if
    alpha_cut = max(OLVER_OPEN_ALPHA_FLOOR, alpha_cut)

    end function open_alpha_cut

    elemental logical function use_smallchi_map(l, nu, achi, alpha_gate) result(use_smallchi)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, achi, alpha_gate

    ! Pointwise version of the near-flat small-chi integration gate. The
    ! full integration uses 0.3 with chi_max; using the current achi as the
    ! local endpoint needs a smaller threshold to preserve the pointwise
    ! phi_olver envelope. The gate uses nu/l; the map below still uses the
    ! more accurate sqrt(l(l+1)) curvature scale.  Callers guarantee l >= 3.
    use_smallchi = (alpha_gate > SMALLCHI_GATE_ALPHA .or. &
        (l >= SMALLCHI_GATE_LOW_ALPHA_L_MIN .and. alpha_gate > 1._dl)) .and. &
        real(l, dl)**2*achi**7/nu < SMALLCHI_GATE_METRIC

    end function use_smallchi_map

    pure subroutine compute_olver_z_amp(l, K, nu, achi, sin_k, z, amp)
    ! Leading Olver map chi -> z: the flat coordinate whose Liouville-Green
    ! action from the turning point matches the curved one, with the amplitude
    ! amp = (dz/dchi)^(-1/2) so that u(chi) = amp * z j_l(nu z).
    ! sin_k = S_K(achi) is supplied by the caller, which needs it anyway.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, achi, sin_k
    real(dl), intent(out) :: z
    real(dl), intent(out), optional :: amp

    real(dl) :: ell, alpha, turn_z, turn_chi, action, eps, c1, c2, turn_scale

    if (K == 0) then
        z = achi
        if (present(amp)) amp = 1._dl
        return
    end if

    ell = sqrt(real(l, dl)*real(l + 1, dl))
    alpha = nu/ell
    turn_z = ell/nu
    turn_chi = turning_point(ell, nu, K)

    if (achi <= min(1.0e-6_dl, 1.0e-4_dl*max(turn_chi, 1._dl))) then
        ! z -> chi only up to a factor 1 + O(K/alpha^2); wherever that factor is
        ! not negligible alpha <~ 1, and then nu*achi <= 1e-6 nu << l makes phi
        ! utterly negligible compared with its peak.
        z = achi
    else if (abs(achi - turn_chi) <= 1.0e-4_dl*max(turn_chi, turn_z)) then
        ! Near the turning point analytic_amplitude below would be evaluated as a
        ! ratio of two quantities that both vanish there, so expand the map
        ! directly in eps = chi - chi_t instead.  Matching
        ! (dz/dchi)^2 (alpha^2 - 1/z^2) = alpha^2 - 1/S_K^2 order by order gives
        ! z = z_t + eps (c1 + c2 eps) with c1^3 = cos_K(chi_t); carrying c2 costs
        ! no transcendentals and cuts the amplitude error here by ~1e4.
        turn_scale = turning_scale(K, turn_z)
        c1 = turn_scale**(1._dl/3._dl)
        c2 = (3._dl*alpha**2*(c1**4 - turn_scale**2) - real(K, dl))/(10._dl*alpha*c1**2)
        eps = achi - turn_chi
        z = turn_z + eps*(c1 + c2*eps)
        if (present(amp)) amp = 1._dl/sqrt(c1 + 2._dl*c2*eps)
        return
    else
        action = qintegral_exact(sin_k, alpha, K)
        z = invert_flat_action(action, turn_z, achi < turn_chi)
    end if

    if (present(amp)) then
        amp = analytic_amplitude(achi, sin_k, z, K, alpha, turn_chi)
    end if

    end subroutine compute_olver_z_amp

    elemental real(dl) function turning_scale(K, turn_z) result(turn_scale)
    ! cos_K(chi_t) = sqrt(1 - K S_K(chi_t)^2) at the turning point S_K(chi_t) = turn_z;
    ! cos(asin(x)) and cosh(asinh(x)) both reduce to sqrt(1 - K x^2).
    integer, intent(in) :: K
    real(dl), intent(in) :: turn_z

    turn_scale = sqrt(1._dl - real(K, dl)*turn_z*turn_z)

    end function turning_scale

    elemental subroutine compute_olver_z_amp_smallchi(l, K, nu, chi, z, amp)
    ! Small-chi curvature expansion for the Olver action map.
    !
    ! This solves the differentiated Liouville-Green map perturbatively for
    ! small curvature across the interval,
    !
    !   (dz/dchi)^2 * (alpha^2 - 1/z^2)
    !       = alpha^2 - 1/S_K(chi)^2,
    !
    ! through O((K/alpha^2)^3), writing z = chi * F(alpha*chi,K/alpha^2).
    ! For this local curvature expansion the exact centrifugal coefficient
    ! sqrt(l(l+1)) is more accurate than the Langer parameter at low ell, and
    ! is indistinguishable from it at high ell. The full cached Olver map uses
    ! the same coefficient in compute_olver_z_amp.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    real(dl), intent(out) :: z, amp

    real(dl) :: ell2, nu2, chi2
    real(dl) :: rk, h, h2, a, a2, ha
    real(dl) :: F, D

    real(dl), parameter :: c6 = 1._dl/6._dl
    real(dl), parameter :: c360 = 1._dl/360._dl
    real(dl), parameter :: c45360 = 1._dl/45360._dl

    if (K == 0 .or. abs(chi) <= CACHE_EPS) then
        z = chi
        amp = 1._dl
        return
    end if

    rk = real(K, dl)

    ! ell^2 = l(l+1), avoiding sqrt(l(l+1))
    ell2 = real(l, dl)*real(l + 1, dl)

    nu2 = nu*nu
    chi2 = chi*chi

    ! h = K / alpha^2 = K * ell^2 / nu^2.
    h = rk*ell2/nu2

    ! a = h * t^2 = K * chi^2.
    a = rk*chi2
    a2 = a*a

    h2 = h*h
    ha = h*a

    F = 1._dl - h*( &
        c6 &
        + c360*(4._dl*a + 13._dl*h) &
        + c45360*(48._dl*a2 + 148._dl*ha + 737._dl*h2))

    D = 1._dl - h*( &
        c6 &
        + c360*(12._dl*a + 13._dl*h) &
        + c45360*(240._dl*a2 + 444._dl*ha + 737._dl*h2))

    z = chi*F

    if (D > 0._dl) then
        amp = 1._dl/sqrt(D)
    else
        amp = 0._dl
    end if

    end subroutine compute_olver_z_amp_smallchi

    elemental real(dl) function analytic_amplitude(chi, sin_k, z, K, alpha, turn_chi)
    ! Liouville-Green amplitude (dz/dchi)^(-1/2) = |(alpha^2 - 1/z^2)/(alpha^2 - 1/S_K^2)|^(1/4).
    ! Both terms vanish at the turning point, so the ratio is only used away from it;
    ! compute_olver_z_amp expands the map directly inside its near-turning window and
    ! the guards here are a backstop for the residual cancellation.
    real(dl), intent(in) :: chi, sin_k, z, alpha, turn_chi
    integer, intent(in) :: K

    real(dl) :: alpha2, flat_term, curved_term

    if (K == 0) then
        analytic_amplitude = 1._dl
        return
    end if

    if (chi <= 1.0e-8_dl .or. z <= 1.0e-8_dl) then
        analytic_amplitude = 1._dl
        return
    end if

    alpha2 = alpha*alpha
    if (abs(chi - turn_chi) <= 1.0e-10_dl*max(1._dl, turn_chi)) then
        analytic_amplitude = turning_scale(K, 1._dl/alpha)**(-1._dl/6._dl)
        return
    end if

    flat_term = alpha2 - 1._dl/z**2
    curved_term = alpha2 - 1._dl/sin_k**2

    if (abs(flat_term) + abs(curved_term) <= 100._dl*CACHE_EPS*max(1._dl, alpha2)) then
        analytic_amplitude = turning_scale(K, 1._dl/alpha)**(-1._dl/6._dl)
        return
    end if

    analytic_amplitude = sqrt(sqrt(abs(flat_term/curved_term)))

    end function analytic_amplitude

    elemental real(dl) function invert_flat_action(action, z_turn, below_turn)
    ! Invert the flat action q for u = z/z_turn: q = t - tanh(t), u = sech(t) below
    ! the turning point, q = tan(theta) - theta, u = sec(theta) above it.
    real(dl), intent(in) :: action, z_turn
    logical, intent(in) :: below_turn

    real(dl) :: q, s, u
    real(dl) :: A, e
    real(dl) :: x, qplus, invqplus, invqplus2

    ! Branch cuts in q equivalent to p = (3q)^(1/3) < 1.72 and < 1.97, so the
    ! cube root is only taken on the near-turning polynomial branches.  Both are
    ! set where the polynomial fit and the asymptotic form have equal error
    ! against the exact inverse (worst 1.5e-6 evanescent, 8.5e-9 oscillatory).
    real(dl), parameter :: Q_CUT_EVAN = 1.72_dl**3/3._dl
    real(dl), parameter :: Q_CUT_OSC = 1.97_dl**3/3._dl

    q = max(action, 0._dl)

    ! The near-turning coefficients below are numerical polynomial fits for u(p^2),
    ! not Taylor coefficients. This avoids a root solve in the hot Olver map path.
    if (below_turn) then

        ! Evanescent branch:
        !   q = t - tanh(t),  u = sech(t)
        !
        ! Near the turning point use a polynomial fit in p^2.
        ! Farther out use the large-q asymptotic inversion for sech(t).
        if (q < Q_CUT_EVAN) then
            s = (3._dl*q)**(2._dl/3._dl)

            u = (((((( &
                2.2687533176976695e-05_dl*s &
                - 1.7210470362266376e-04_dl)*s &
                - 3.1231653154203167e-04_dl)*s &
                + 2.1586618823186233e-04_dl)*s &
                + 7.5055874727872618e-02_dl)*s &
                - 5.0000720190143755e-01_dl)*s &
                + 1.0000000000000000_dl)

        else
            A = q + 1._dl
            x = exp(-A)
            e = x*x

            u = 2._dl*x*(1._dl + e*(1._dl + 3._dl*e))
        end if

    else

        ! oscillatory inverse:
        !
        !   q = tan(theta) - theta
        !   u = sec(theta)
        !
        ! Uses a polynomial fit in p^2 near the turning point and a high-order
        ! asymptotic expansion in q + pi/2 farther out.

        if (q < Q_CUT_OSC) then
            s = (3._dl*q)**(2._dl/3._dl)

            u = ((((((( &
                3.4828644185221886e-07_dl*s &
                - 8.3097833638549810e-06_dl)*s &
                + 8.8500265224195657e-05_dl)*s &
                - 4.8620948599887724e-04_dl)*s &
                - 3.5056578112088423e-04_dl)*s &
                + 7.4998103457285289e-02_dl)*s &
                + 5.0000018483290443e-01_dl)*s &
                + 1.0000000000000000_dl)

        else
            qplus = q + const_pi/2._dl
            invqplus = 1._dl/qplus
            invqplus2 = invqplus*invqplus

            u = qplus - invqplus*( &
                0.5_dl + invqplus2*( &
                7._dl/24._dl + invqplus2*( &
                83._dl/240._dl + invqplus2*( &
                6949._dl/13440._dl + invqplus2*( &
                23399._dl/26880._dl + invqplus2* &
                266317._dl/168960._dl)))))
        end if

    end if

    invert_flat_action = z_turn*u

    end function invert_flat_action

    end module HypersphericalBesselOlver
