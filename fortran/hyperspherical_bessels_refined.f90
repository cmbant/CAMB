    module HypersphericalBesselRefined
    ! Short exact series and explicit Airy corrections for folded curved modes.
    ! All public evaluations return the reduced function u = S_K(chi)*Phi.
    ! A false ok means that the caller must use its existing evaluator.
    use precision, only: dl
    use constants, only: const_pi
    use MathUtils, only: airy_ai_fast, airy_ai_reference
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: refined_u
    real(dl), parameter :: tol = 2.e-12_dl
    integer, parameter :: max_terms = 40
    contains

    elemental real(dl) function log1p(x) result(v)
    real(dl), intent(in) :: x
    real(dl) :: y
    y = 1._dl+x
    if (y == 1._dl) then
        v = x
    else
        v = log(y)*(x/(y-1._dl))
    end if
    end function log1p

    pure complex(dl) function stirling(z) result(v)
    ! R(z) in log Gamma(z+1) = (z+1/2) log z-z+log(2*pi)/2+R(z).
    ! All callers have |z| >= 8 and Re z >= 0.
    complex(dl), intent(in) :: z
    complex(dl) :: w, w2
    w = 1._dl/z
    w2 = w*w
    v = w*(1._dl/12+w2*(-1._dl/360+w2*(1._dl/1260+w2*(-1._dl/1680+w2*(1._dl/1188+ &
        w2*(-691._dl/360360+w2*(1._dl/156+w2*(-3617._dl/122400))))))))
    end function stirling

    pure real(dl) function half_logratio(z) result(v)
    ! log Gamma(z+1/2)/Gamma(z), avoiding subtraction of large log-gammas.
    real(dl), intent(in) :: z
    if (z < 10._dl) then
        v = log_gamma(z+0.5_dl)-log_gamma(z)
    else
        v = 0.5_dl*log(z)+z*log1p(0.5_dl/z)-0.5_dl+ &
            real(stirling(cmplx(z+0.5_dl,0._dl,dl))-stirling(cmplx(z,0._dl,dl)),dl)
    end if
    end function half_logratio

    pure subroutine closed_series(l, nu, s, c, qmax, u, ok)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, s, c, qmax
    real(dl), intent(out) :: u
    logical, intent(out) :: ok
    integer :: p, m, parity, j
    real(dl) :: a, y, t, total, asum, q, qnext, tail, lp, logs
    ok = .false.
    u = 0._dl
    p = nint(nu)-l-1
    if (p < 0) return
    m = p/2
    parity = mod(p,2)
    a = real(l+1,dl)+m+parity
    y = c*c
    ! The absolute term ratios decrease monotonically. Reject expensive cases early.
    if (m > 8 .and. y*m*a/(parity+0.5_dl) > qmax) return
    t = 1._dl
    total = t
    asum = abs(t)
    do j = 0, min(m,max_terms)
        q = 0._dl
        tail = 0._dl
        if (j < m) then
            q = y*((m-j)/real(j+1,dl))*((a+j)/(j+parity+0.5_dl))
            qnext = y*((m-j-1)/real(j+2,dl))*((a+j+1)/(j+parity+1.5_dl))
            if (qnext >= 1._dl) then
                tail = huge(1._dl)
            else
                tail = abs(t*q)/(1._dl-qnext)
            end if
        end if
        ! Tail bound is exact in real arithmetic; the roundoff allowance is heuristic.
        if (tail+32*epsilon(1._dl)*(j+1)*asum <= tol*abs(total)) then
            if (abs(c) < 0.5_dl) then
                logs = 0.5_dl*log1p(-y)
            else
                logs = log(s)
            end if
            lp = real(l+1,dl)*logs-0.5_dl*(log(2*nu)+ &
                half_logratio(real(l+1+m,dl))+half_logratio(m+0.5_dl))
            if (parity == 1) total = total*c*sqrt(2._dl*(real(l+1,dl)+m)*(2._dl*m+1))
            if (mod(m,2) == 1) total = -total
            if (total == 0._dl) then
                u = 0._dl
            else if (abs(lp) < 650._dl) then
                u = exp(lp)*total
            else
                u = sign(exp(lp+log(abs(total))),total)
            end if
            ok = ieee_is_finite(u)
            return
        end if
        if (j == m .or. j == max_terms) return
        t = -t*q
        total = total+t
        asum = asum+abs(t)
        if (asum > 1.e100_dl) return
    end do
    end subroutine closed_series

    pure real(dl) function phase_over_nu(l, nu) result(v)
    ! sum atan(nu/j)/nu via a stable difference of log-gamma phases.
    integer, intent(in) :: l
    real(dl), intent(in) :: nu
    real(dl) :: x, w, w2, t, n2, pow, wp, tail, shift
    integer :: j, p
    real(dl), parameter :: zeta_odd(6) = [1.2020569031595942854_dl,1.0369277551433699263_dl, &
        1.0083492773819228268_dl,1.0020083928260822144_dl,1.0004941886041194646_dl,1.0001227133475784891_dl]
    x = real(l,dl)
    if (l <= 12) then
        v = 0._dl
        do j = 1,l
            t = nu/j
            if (abs(t) < 1.e-4_dl) then
                v = v+(1-t*t/3+t**4/5)/j
            else
                v = v+atan(t)/nu
            end if
        end do
    else if (nu <= 0.125_dl) then
        w = 1/x
        w2 = w*w
        v = log(x)+0.5772156649015328606_dl+w/2-w2*(1._dl/12-w2*(1._dl/120-w2* &
            (1._dl/252-w2*(1._dl/240-w2/132))))
        n2 = nu*nu
        pow = -n2
        wp = w2
        do j = 1,6
            p = 2*j+1
            tail = wp*(1._dl/(p-1)-w/2+p*w2*(1._dl/12-(p+1)*(p+2)*w2* &
                (1._dl/720-(p+3)*(p+4)*w2/30240)))
            v = v+pow*(zeta_odd(j)-tail)/p
            pow = -pow*n2
            wp = wp*w2
        end do
    else if (nu >= 10._dl) then
        t = x/nu
        v = 0.5_dl*log1p(t*t)+((x+0.5_dl)*atan2(nu,x)-const_pi/4)/nu+ &
            aimag(stirling(cmplx(x,nu,dl))-stirling(cmplx(0._dl,nu,dl)))/nu
    else
        shift = 0._dl
        do j = 1,10
            shift = shift+atan(nu/j)
        end do
        v = 0.5_dl*log((x*x+nu*nu)/(100+nu*nu))+ &
            ((x+0.5_dl)*atan(nu/x)-10.5_dl*atan(nu/10)+shift+ &
            aimag(stirling(cmplx(x,nu,dl))-stirling(cmplx(10._dl,nu,dl))))/nu
    end if
    end function phase_over_nu

    pure subroutine open_jost(l, nu, chi, z, qmax, u, ok)
    ! Exact outgoing wave: exp(i nu chi) (1-z)^(-l) 2F1(-l,-l-i nu;1-i nu;z).
    ! Store the imaginary parts divided by nu, including the limit nu=0.
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, chi, z, qmax
    real(dl), intent(out) :: u
    logical, intent(out) :: ok
    integer :: j
    real(dl) :: ell, nu2, d, theta, us, cs, tr, ti, sr, si, aa, bb, fac, fr, fi, newr
    real(dl) :: mag, q, h, hn, mn, s0, s1, ar, ai, rs, ri, scale, err, rnd, amp
    u = 0._dl
    ok = .false.
    ell = real(l,dl)
    nu2 = nu*nu
    q = z*ell*sqrt((ell*ell+nu2)/(1+nu2))
    ! Keep the amplification (1-z)^(-l) bounded as well as the term ratios.
    if (z > 0.25_dl .or. q > qmax) return
    amp = exp(-ell*log1p(-z))
    d = chi-phase_over_nu(l,nu)
    theta = nu*d
    if (abs(theta) < 1.e-4_dl) then
        us = d*(1-theta*theta/6+theta**4/120)
    else
        us = sin(theta)/nu
    end if
    cs = cos(theta)
    tr = 1._dl
    ti = 0._dl
    sr = tr
    si = ti
    mag = 1._dl
    h = 0._dl
    rs = mag
    ri = 0._dl
    do j = 0, min(l,max_terms)
        u = us*sr+cs*si
        ! Apply the absolute error floor to the final, amplified reduced function.
        scale = max(abs(u),1/(amp*max(1._dl,nu)))
        err = 0._dl
        aa = ell-j
        bb = real(j+1,dl)
        mn = 0._dl
        hn = h
        if (j < l) then
            q = z*(aa/bb)*sqrt((aa*aa+nu2)/(bb*bb+nu2))
            mn = mag*q
            hn = h+1/aa+1/bb
            if (q >= 1._dl) then
                err = huge(1._dl)
            else
                s0 = 1/(1-q)
                s1 = q*s0*s0
                ar = mn*s0
                ! Each coefficient phase increment divided by nu is <=1/a+1/b<=2.
                ai = mn*(hn*s0+2*s1)
                if (nu >= 1._dl) ai = min(ai,ar/nu)
                err = abs(us)*ar+abs(cs)*ai
            end if
        end if
        rnd = 8*epsilon(1._dl)*(j+1)*(abs(us)*rs+abs(cs)*ri)
        if (err+rnd <= tol*scale) then
            u = amp*u
            ok = ieee_is_finite(u)
            return
        end if
        if (j == l .or. j == max_terms) return
        fac = z*(aa/bb)/(bb*bb+nu2)
        fr = fac*(aa*bb-nu2)
        fi = fac*(ell+1)
        newr = tr*fr-nu2*ti*fi
        ti = tr*fi+ti*fr
        tr = newr
        sr = sr+tr
        si = si+ti
        mag = mn
        h = hn
        rs = rs+mag
        ri = ri+mag*h
    end do
    end subroutine open_jost

    pure real(dl) function gauss24_x(i) result(v)
    integer, intent(in) :: i
    real(dl), parameter :: tab(24) = [ &
        2.40639000148934468e-03_dl, 1.26357220143452631e-02_dl, 3.08627239986336011e-02_dl, &
        5.67922364977994643e-02_dl, 8.99990070130485265e-02_dl, 1.29937904210722821e-01_dl, &
        1.75953174031512227e-01_dl, 2.27289264305580219e-01_dl, 2.83103246186977464e-01_dl, &
        3.42478660151918302e-01_dl, 4.04440566263191859e-01_dl, 4.67971553568697185e-01_dl, &
        5.32028446431302759e-01_dl, 5.95559433736808197e-01_dl, 6.57521339848081698e-01_dl, &
        7.16896753813022536e-01_dl, 7.72710735694419837e-01_dl, 8.24046825968487773e-01_dl, &
        8.70062095789277179e-01_dl, 9.10000992986951474e-01_dl, 9.43207763502200480e-01_dl, &
        9.69137276001366343e-01_dl, 9.87364277985654737e-01_dl, 9.97593609998510655e-01_dl ]
    v = tab(i)
    end function gauss24_x

    pure real(dl) function gauss24_w(i) result(v)
    integer, intent(in) :: i
    real(dl), parameter :: tab(24) = [ &
        6.17061489999434516e-03_dl, 1.42656943144667796e-02_dl, 2.21387194087097061e-02_dl, &
        2.96492924577181800e-02_dl, 3.66732407055400805e-02_dl, 4.30950807659766025e-02_dl, &
        4.88093260520569630e-02_dl, 5.37221350579827825e-02_dl, 5.77528340268627580e-02_dl, &
        6.08352364639016471e-02_dl, 6.29187281734141235e-02_dl, 6.39690976733760103e-02_dl, &
        6.39690976733760103e-02_dl, 6.29187281734141235e-02_dl, 6.08352364639016471e-02_dl, &
        5.77528340268627580e-02_dl, 5.37221350579827825e-02_dl, 4.88093260520569630e-02_dl, &
        4.30950807659766025e-02_dl, 3.66732407055400805e-02_dl, 2.96492924577181800e-02_dl, &
        2.21387194087097061e-02_dl, 1.42656943144667796e-02_dl, 6.17061489999434516e-03_dl ]
    v = tab(i)
    end function gauss24_w

    pure real(dl) function gauss32_x(i) result(v)
    integer, intent(in) :: i
    real(dl), parameter :: tab(32) = [ &
        1.36806907525921506e-03_dl, 7.19424422736580915e-03_dl, 1.76188722062468051e-02_dl, &
        3.25469620311301666e-02_dl, 5.18394221169739544e-02_dl, 7.53161931337150148e-02_dl, &
        1.02758102016028807e-01_dl, 1.33908940629855144e-01_dl, 1.68477866534892384e-01_dl, &
        2.06142121379618848e-01_dl, 2.46550045533885320e-01_dl, 2.89324361934682361e-01_dl, &
        3.34065698858936166e-01_dl, 3.80356318873931454e-01_dl, 4.27764019208601742e-01_dl, &
        4.75846167156130817e-01_dl, 5.24153832843869183e-01_dl, 5.72235980791398258e-01_dl, &
        6.19643681126068491e-01_dl, 6.65934301141063889e-01_dl, 7.10675638065317639e-01_dl, &
        7.53449954466114624e-01_dl, 7.93857878620381152e-01_dl, 8.31522133465107616e-01_dl, &
        8.66091059370144856e-01_dl, 8.97241897983971137e-01_dl, 9.24683806866285041e-01_dl, &
        9.48160577883025990e-01_dl, 9.67453037968869833e-01_dl, 9.82381127793753195e-01_dl, &
        9.92805755772634191e-01_dl, 9.98631930924740785e-01_dl ]
    v = tab(i)
    end function gauss32_x

    pure real(dl) function gauss32_w(i) result(v)
    integer, intent(in) :: i
    real(dl), parameter :: tab(32) = [ &
        3.50930500473525288e-03_dl, 8.13719736545287162e-03_dl, 1.26960326546310120e-02_dl, &
        1.71369314565108823e-02_dl, 2.14179490111134178e-02_dl, 2.54990296311880457e-02_dl, &
        2.93420467392677825e-02_dl, 3.29111113881808415e-02_dl, 3.61728970544241690e-02_dl, &
        3.90969478935351139e-02_dl, 4.16559621134733535e-02_dl, 4.38260465022018916e-02_dl, &
        4.55869393478818899e-02_dl, 4.69221995404022554e-02_dl, 4.78193600396373542e-02_dl, &
        4.82700442573638297e-02_dl, 4.82700442573638297e-02_dl, 4.78193600396373542e-02_dl, &
        4.69221995404022554e-02_dl, 4.55869393478818899e-02_dl, 4.38260465022018916e-02_dl, &
        4.16559621134733535e-02_dl, 3.90969478935351139e-02_dl, 3.61728970544241690e-02_dl, &
        3.29111113881808415e-02_dl, 2.93420467392677825e-02_dl, 2.54990296311880457e-02_dl, &
        2.14179490111134178e-02_dl, 1.71369314565108823e-02_dl, 1.26960326546310120e-02_dl, &
        8.13719736545287162e-03_dl, 3.50930500473525288e-03_dl ]
    v = tab(i)
    end function gauss32_w

    pure real(dl) function open_log_norm(l, nu) result(v)
    ! log(prod sqrt(1+(nu/j)^2)), using the infinite product for sinh(pi*nu).
    integer, intent(in) :: l
    real(dl), intent(in) :: nu
    real(dl) :: x, t, t2, p, r
    x = real(l,dl)
    t = nu/x
    t2 = t*t
    p = const_pi*nu
    if (p < 0.1_dl) then
        r = p*p*(1._dl/6+p*p*(-1._dl/180+p*p*(1._dl/2835-p*p/37800)))
    else
        r = p-log(2*p)+log1p(-exp(-2*p))
    end if
    v = 0.5_dl*(r+(x+0.5_dl)*log1p(t2))-nu*atan(t)+ &
        real(stirling(cmplx(x,nu,dl))-stirling(cmplx(x,0._dl,dl)),dl)
    end function open_log_norm

    pure real(dl) function compact_integrand(y, ell, b, nu) result(v)
    real(dl), intent(in) :: y, ell, b, nu
    real(dl) :: z
    z = sqrt(b)*y
    v = exp(ell*log1p(-y*y/ell))*cos(2*nu*asinh(z))/sqrt(1+z*z)
    end function compact_integrand

    pure subroutine open_compact(l, nu, chi, s, c, u, ok)
    ! Real integral after y=sqrt(l)*sinh(t/2)/sinh(chi/2).
    ! In the compact box the integrand is bounded by exp(-y^2).
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, chi, s, c
    real(dl), intent(out) :: u
    logical, intent(out) :: ok
    real(dl) :: ell, b, ymax, j24, j32, lp, z
    integer :: i
    ok = .false.
    u = 0._dl
    if (l < 8) return
    ell = real(l,dl)
    b = (s/(1+c))*s/(2*ell)
    if (b > 0.25_dl .or. nu*nu*b > 5.25_dl .or. b <= 0._dl .or. nu > 4.6_dl) return
    ! The omitted tail above y=7 is <= exp(-49)/14. The two rules
    ! are an empirical convergence check, not a certified quadrature error.
    ymax = min(sqrt(ell),7._dl)
    j24 = 0._dl
    j32 = 0._dl
    do i = 1,24
        j24 = j24+gauss24_w(i)*compact_integrand(ymax*gauss24_x(i),ell,b,nu)
    end do
    do i = 1,32
        j32 = j32+gauss32_w(i)*compact_integrand(ymax*gauss32_x(i),ell,b,nu)
    end do
    if (abs(j24-j32)*ymax > 2.e-9_dl) return
    z = exp(-chi)
    if (z < 0.25_dl) then
        lp = ell*(log1p(-z)-log1p(z))
    else
        lp = ell*log(tanh(chi/2))
    end if
    u = 2*sqrt(b)*exp(open_log_norm(l,nu)+lp)*ymax*j32
    ok = ieee_is_finite(u)
    end subroutine open_compact

    pure real(dl) function effective_mu(l, k, nu) result(mu)
    integer, intent(in) :: l,k
    real(dl), intent(in) :: nu
    real(dl) :: la,x,y
    la = l+0.5_dl
    if (k == -1) then
        x = min(la,nu)
        y = max(la,nu)
        mu = x/sqrt(1+(x/y)**2)
    else
        x = la/nu
        y = 1+x*x
        mu = la*((nu-la)/nu)*(1+x)/(y*sqrt(y))
    end if
    end function effective_mu

    pure real(dl) function map_norm(la,h,oneh,k,nu,second) result(v)
    ! Stirling expansion after cancellation of the origin-action normalization.
    real(dl), intent(in) :: la,h,oneh,nu
    integer, intent(in) :: k
    logical, intent(in) :: second
    real(dl) :: a,x,y,x2,y2,a2,logr
    a = 1/la
    a2 = a*a
    y = a/oneh
    x = h*y
    x2 = x*x
    y2 = y*y
    logr = -a*a2*a2/1260+(31._dl/40320)*x*x2*(5*y2+10*x*y+x2) &
        +a*a2**3/1680-(127._dl/215040)*x2*x2*(7*y*y2+35*x*y2+21*x2*y+x*x2) &
        -a*a2**4/1188+(511._dl/608256)*x*x2*x2*(9*y2*y2+84*x*y*y2+126*x2*y2+36*x*x2*y+x*x2)
    if (.not.second) logr = logr+a*a2/360-(7._dl/2880)*x2*(3*y+x)
    v = sqrt(const_pi/(nu*la))*(1+logr*(1+0.5_dl*logr))
    if (k == -1 .and. nu < 6) v = v*sqrt(1-exp(-2*const_pi*nu))
    end function map_norm

    pure real(dl) function mapped_ai(x, action) result(ai)
    ! Keep the accurate action in the asymptotic tails. CAMB's fast Ai is
    ! peak-accurate, but its positive-tail truncation is not relative-accurate.
    real(dl), intent(in) :: x, action
    real(dl) :: z, z2, a, b, term
    integer :: j
    real(dl), parameter :: coef(0:6) = [1._dl,5._dl/72,385._dl/10368,85085._dl/2239488, &
        37182145._dl/644972544,5391411025._dl/46438023168._dl,5849680962125._dl/20061226008576._dl]
    if (x >= -12._dl .and. x < 2.09_dl) then
        ai = airy_ai_fast(x)
    else if (x >= 2.09_dl .and. x < 10._dl) then
        ai = airy_ai_reference(x)
    else if (x > 0._dl) then
        if (action > 745._dl) then
            ai = 0._dl
            return
        end if
        z = 1/action
        a = 1._dl
        term = 1._dl
        do j = 1,6
            term = -term*z
            a = a+coef(j)*term
        end do
        ai = exp(-action)*a/(2*sqrt(const_pi*sqrt(x)))
    else
        z = 1/action
        z2 = z*z
        a = 1+z2*(-coef(2)+z2*(coef(4)-z2*coef(6)))
        b = z*(coef(1)+z2*(-coef(3)+z2*coef(5)))
        ai = (sin(action+const_pi/4)*a-cos(action+const_pi/4)*b)/sqrt(const_pi*sqrt(-x))
    end if
    end function mapped_ai

    pure subroutine corrected_airy(l,k,nu,s,c,u,ok,second)
    ! Uniform Airy map g=tau+b+c1, including the optional second correction.
    ! tau=lambda^(2/3)*zeta; all shift derivatives below are wrt tau.
    integer, intent(in) :: l,k
    real(dl), intent(in) :: nu,s,c
    real(dl), intent(out) :: u
    logical, intent(out) :: ok
    logical, intent(in) :: second
    real(dl) :: la,b,b2,h,oneh,s2,d2,d,x,t,act,zeta,rt,g,r,rp,dd,ddp,coef(0:16),hn,hs, &
        v,b0,bp,psi,amp,tau,sc,c0,cp,pp(0:4),sg,zz,cv,cvp,tt,jac,iv,hh, &
        qpow(0:16),ppow(0:8),ff(5:16),fac,cor,xx
    integer :: j,m,nt
    ! Generated by tests/generate_refined_airy.py using exact rational polynomial arithmetic.
    real(dl),parameter::F_local(0:8,5:16)=reshape([ &
    -6.45907692307692294e+00_dl, 5.16923076923076934e+01_dl, -7.07692307692307736e+01_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -1.67378062300511274e+01_dl, 1.41516228895820745e+02_dl, -2.67012561135010117e+02_dl, &
    1.06673734428836468e+02_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -2.81114863687732850e+01_dl, 2.54991877430292789e+02_dl, -6.06928935330375907e+02_dl, &
    4.06508152111993638e+02_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -3.90175901703309975e+01_dl, 3.81867506746916263e+02_dl, -1.09452125187505999e+03_dl, &
    1.04346702633770656e+03_dl, -1.72686464009576980e+02_dl, 0.00000000000000000e+00_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -4.87497848683082182e+01_dl, 5.15251855644299894e+02_dl, -1.72592068029721668e+03_dl, &
    2.16021804154328174e+03_dl, -7.71079692019535855e+02_dl, 0.00000000000000000e+00_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -5.70826378685569509e+01_dl, 6.50578859659298018e+02_dl, -2.49327431458207093e+03_dl, &
    3.88800192837051873e+03_dl, -2.18186059467484847e+03_dl, 2.06396354904949931e+02_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -6.40301364338917267e+01_dl, 7.84866366367451405e+02_dl, -3.38697387861232392e+03_dl, &
    6.34419147108669222e+03_dl, -4.91578390128887986e+03_dl, 1.10152289154608252e+03_dl, &
    0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -6.97120922551179660e+01_dl, 9.16210328955103591e+02_dl, -4.39690209732313815e+03_dl, &
    9.63235292192156703e+03_dl, -9.59081711256942981e+03_dl, 3.56190796479534038e+03_dl, &
    -2.21344893001411236e+02_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -7.42852066852422865e+01_dl, 1.04343994524934874e+03_dl, -5.51311734008427993e+03_dl, &
    1.38432728341934890e+04_dl, -1.69157839857810504e+04_dl, 8.96377502517564790e+03_dl, &
    -1.38715692538893359e+03_dl, 0.00000000000000000e+00_dl, 0.00000000000000000e+00_dl, &
    -7.79094971770334581e+01_dl, 1.16588168557171025e+03_dl, -6.72621383985410466e+03_dl, &
    1.90563210173502666e+04_dl, -2.76760503877355150e+04_dl, 1.92703128090917635e+04_dl, &
    -5.10656218755381087e+03_dl, 2.25890281356459241e+02_dl, 0.00000000000000000e+00_dl, &
    -8.07332587675702626e+01_dl, 1.28319740171656781e+03_dl, -8.02749266522812923e+03_dl, &
    2.53408669611733421e+04_dl, -4.27214337723826320e+04_dl, 3.71027476470379552e+04_dl, &
    -1.43490488920329735e+04_dl, 1.63138126384244424e+03_dl, 0.00000000000000000e+00_dl, &
    -8.28874426383395786e+01_dl, 1.39527308076575628e+03_dl, -9.40902257678715614e+03_dl, &
    3.27576218214679539e+04_dl, -6.29561495806573366e+04_dl, 6.57981760312559491e+04_dl, &
    -3.39909667117589051e+04_dl, 6.76749513268668215e+03_dl, -2.24724494132639592e+02_dl ], [9,12])
    ok = .false.
    u = 0._dl
    la = l+0.5_dl
    b = nu/la
    b2 = b*b
    h = k/b2
    oneh = 1-h
    if (k == 1) oneh = (nu-la)*(nu+la)/(nu*nu)
    s2 = s*s
    d2 = (1-b*s)*(1+b*s)
    d = sqrt(abs(d2))
    if (c == 0._dl) then
        x = -huge(1._dl)
    else
        x = d2/(c*c)
    end if
    c0 = 0._dl
    cp = 0._dl
    sc = max(1._dl,abs(h))
    if (abs(x)*sc <= 0.2_dl) then
        ! R(x)=sum 3(1+h+...+h^j)x^j/(2j+3), regular at the turn.
        t = x*sc
        hn = 1._dl
        hs = 1._dl
        coef(0) = 1._dl
        iv = 1/sc
        hh = h*iv
        nt = 16
        if (abs(t)<0.1_dl) nt = 13
        if (abs(t)<0.03_dl) nt = 10
        if (abs(t)<0.008_dl) nt = 8
        if (.not.second) then
            if (abs(t)<0.03_dl) nt = 8
            if (abs(t)<0.008_dl) nt = 6
        end if
        do j = 1,nt
            hn = hn*hh
            hs = hs*iv+hn
            coef(j) = 3*hs/(2*j+3)
        end do
        r = coef(nt)
        rp = 0._dl
        do j = nt-1,0,-1
            rp = rp*t+r
            r = r*t+coef(j)
        end do
        dd = coef(nt)-coef(1)*coef(nt-1)
        ddp = 0._dl
        do j = nt-1,2,-1
            ddp = ddp*t+dd
            dd = dd*t+coef(j)-coef(1)*coef(j-1)
        end do
        g = (la*oneh*r/2)**(1._dl/3)
        zeta = x*g*g
        b0 = 5._dl/48*dd*(sc/g**2)**2
        bp = 5._dl/48*(sc/g**2)**3*(ddp-4*dd*rp/(3*r))/(1+2*t*rp/(3*r))
        if (second) then
            ! Divided numerator coefficients, in symmetric homogeneous form.
            qpow(0) = 1._dl
            ppow(0) = 1._dl
            do j = 1,nt
                qpow(j) = qpow(j-1)*(iv+hh)
            end do
            do j = 1,nt/2
                ppow(j) = ppow(j-1)*(iv*hh)
            end do
            do j = 5,nt
                ff(j) = 0._dl
                do m = 0,j/2
                    ff(j) = ff(j)+F_local(m,j)*qpow(j-2*m)*ppow(m)
                end do
            end do
            cv = ff(nt)
            cvp = 0._dl
            do j = nt-1,5,-1
                cvp = cvp*t+cv
                cv = cv*t+ff(j)
            end do
            c0 = cv*(sc/g**2)**5/46080
            cp = (sc/g**2)**6/46080*(cvp-10*cv*rp/(3*r))/(1+2*t*rp/(3*r))
        end if
        act = la*abs(x)*sqrt(abs(x))*oneh*r/3
        amp = sqrt(s*g/c)
    else
        ! Stable elementary expressions for the curved Liouville-Green action.
        if (d2 > 0) then
            if (k == -1) then
                act = asinh(d/(s*sqrt(1+b2)))-b*atan2(d,b*c)
            else
                act = asinh(d/(s*(b*sqrt(oneh))))-b*asinh(d/(b*sqrt(oneh)))
            end if
        else
            if (k == -1) then
                act = b*asinh(d/sqrt(1+b2))-atan2(d,c)
            else
                act = ((nu-la)/la)*atan2(d,b*c)-atan2(((nu-la)/la)*d*c,b*c*c+d*d)
            end if
        end if
        if (act <= 0._dl) return
        act = la*act
        zeta = sign((1.5_dl*act)**(2._dl/3),d2)
        rt = sqrt(abs(zeta))
        v = c/d
        sg = sign(1._dl,d2)
        b0 = ((5._dl/3)*v*v*sg-(1+h))*v*sg/(8*la*oneh*rt)-5/(48*zeta*zeta)
        psi = -zeta*s2*(b2*(b2-4*k)*s2+4*b2-k)/(4*d2**3*la*la)+5/(16*zeta*zeta)
        bp = (psi-b0)/(2*zeta)
        if (second) then
            zz = sg*v*v
            pp = [5525._dl,-9945*(1+h),4779+16983*h+4779*h*h, &
                -375*(1+h**3)-7395*h*(1+h),405*h*(1+h*h)+2610*h*h]
            cv = pp(0)
            cvp = 9*pp(0)
            do j = 1,4
                cv = cv*zz+pp(j)
                cvp = cvp*zz+(9-2*j)*pp(j)
            end do
            tt = sg*v*cv/(5760*(la*oneh)**3*rt)
            jac = s2*b2*oneh/(d2*d2*la)
            c0 = tt-1105._dl/(9216*zeta**5)-b0*b0/(4*zeta)+5*b0/(32*zeta**3)
            cp = -cvp*jac/(5760*(la*oneh)**3)-tt/(2*zeta)+5525._dl/(9216*zeta**6) &
                -b0*bp/(2*zeta)+b0*b0/(4*zeta*zeta)+5*bp/(32*zeta**3)-15*b0/(32*zeta**4)
        end if
        amp = sqrt(s*rt/d)
    end if
    cor = b0+c0
    tau = zeta+cor
    fac = 1+bp+cp
    if (fac <= 0._dl .or. .not.ieee_is_finite(fac) .or. .not.ieee_is_finite(tau)) return
    if (zeta < -12._dl .and. abs(cor) < 0.1_dl*abs(zeta)) then
        xx = sqrt(-zeta)
        rt = sqrt(-tau)
        act = act-(2._dl/3)*cor*(xx*xx+xx*rt+rt*rt)/(xx+rt)
    else if (zeta > 10._dl .and. abs(cor) < 0.1_dl*abs(zeta)) then
        xx = sqrt(zeta)
        rt = sqrt(tau)
        act = act+(2._dl/3)*cor*(xx*xx+xx*rt+rt*rt)/(xx+rt)
    else
        act = (2._dl/3)*abs(tau)*sqrt(abs(tau))
    end if
    u = map_norm(la,h,oneh,k,nu,second)*amp*mapped_ai(tau,act)/sqrt(fac)
    ok = ieee_is_finite(u)
    end subroutine corrected_airy

    pure subroutine refined_u(l,k,nu,chi,s,u,ok,method)
    ! Folded chi and S_K(chi) supplied by the caller. Intended for the
    ! expensive fallback region; leave the fast leading Olver path unchanged.
    integer, intent(in) :: l,k
    real(dl), intent(in) :: nu,chi,s
    real(dl), intent(out) :: u
    logical, intent(out) :: ok
    integer, optional, intent(out) :: method
    real(dl) :: c, mu
    ok = .false.
    u = 0._dl
    if (present(method)) method = 0
    if (l < 3 .or. chi <= 0._dl .or. s <= 0._dl .or. nu < 0._dl) return
    if (k /= 1 .and. k /= -1) return
    ! Avoid very large absolute phase arguments outside the validated range.
    if (.not.ieee_is_finite(nu) .or. .not.ieee_is_finite(chi) .or. .not.ieee_is_finite(s)) return
    if (nu > 1.e10_dl/max(1._dl,chi)) return
    if (k == 1) then
        if (nu <= l .or. nu > real(huge(l),dl)) return
    end if
    mu = effective_mu(l,k,nu)
    if (k == 1) then
        c = cos(chi)
        call closed_series(l,nu,s,c,merge(0.25_dl,8._dl,mu >= 4._dl),u,ok)
        if (ok) then
            if (present(method)) method = 1
            return
        end if
    else if (k == -1) then
        call open_jost(l,nu,chi,exp(-2*chi),merge(0.5_dl,4._dl,mu >= 4._dl),u,ok)
        if (ok) then
            if (present(method)) method = 2
            return
        end if
        if (chi > 100._dl) return
        c = cosh(chi)
    else
        return
    end if
    if (mu >= 4._dl .and. nu > 0._dl) then
        call corrected_airy(l,k,nu,s,c,u,ok,mu < 16._dl)
        if (ok .and. present(method)) method = 3
    else if (k == -1) then
        call open_compact(l,nu,chi,s,c,u,ok)
        if (ok .and. present(method)) method = 4
    end if
    end subroutine refined_u
    end module HypersphericalBesselRefined
