# Olver Hyperspherical Bessel Approximation

This note records the current implementation and accuracy checks for the
Olver-style hyperspherical Bessel approximation. The target use case is
near-flat open/closed models, with small radial distances, high multipoles,
and comoving wavenumbers up to the usual CAMB transfer range.

## Equation And Variables

CAMB writes the hyperspherical radial function as

$$
\phi_\ell(\chi) = \frac{u_\ell(\chi)}{S_K(\chi)},
\qquad
S_K(\chi)=
\begin{cases}
\sin\chi & K=+1,\\
\chi & K=0,\\
\sinh\chi & K=-1 .
\end{cases}
$$

The reduced radial function satisfies

$$
u_\ell''(\chi)
+ \left[\beta^2 - \frac{\ell(\ell+1)}{S_K^2(\chi)}\right]u_\ell(\chi)=0 .
$$

The current full Olver map uses the exact centrifugal scale

$$
\lambda=\sqrt{\ell(\ell+1)},\qquad \alpha = \frac{\beta}{\lambda}.
$$

The curved and flat turning points are

$$
S_K(\chi_t)=\frac{1}{\alpha}, \qquad z_t=\frac{1}{\alpha}.
$$

For example, in an open model $\chi_t=\operatorname{asinh}(1/\alpha)$.

## Full Olver Map

The approximation maps the curved radial equation to a flat spherical Bessel
problem by equating Liouville-Green actions. In the oscillatory region,

$$
\int_{\chi_t}^{\chi}
  \sqrt{\alpha^2-\frac{1}{S_K^2(x)}}\,dx
=
\int_{z_t}^{z}
  \sqrt{\alpha^2-\frac{1}{y^2}}\,dy .
$$

The evanescent side uses the corresponding positive action from the point to
the turning point. After solving for $z(\chi)$, the approximation is

$$
\phi_\ell(\chi) \simeq
\frac{z}{S_K(\chi)}
\left|
\frac{\alpha^2-z^{-2}}{\alpha^2-S_K(\chi)^{-2}}
\right|^{1/4}
j_\ell(\beta z).
$$

At the turning point the ratio is evaluated by its limiting value,
$(dS_K/d\chi|_{\chi_t})^{-1/6}$.

## Fast Implementation Details

The production implementation is in `fortran/olver_hyperspherical_bessel.f90`.

- `qintegral_exact` evaluates the curved action analytically for open and
  closed models. The inverse trigonometric branches use `atan2` forms to avoid
  quadrant/sign mistakes.
- `invert_flat_action` inverts the universal flat action. The action is written
  as $q$, and the cache coordinate is
  $p=(3q)^{1/3}$ so the turning point is regular.
- The evanescent cache stores
  $u = \alpha z = \operatorname{sech} t$, where $q=t-\tanh t$.
- The oscillatory cache stores
  $\theta$, where $q=\tan\theta-\theta$, and evaluates
  $u=\alpha z=\sec\theta$ after spline lookup.
- Very close to the turning point, the code uses the analytic local limit of
  the action map,
  $z-z_t\simeq (dS_K/d\chi|_{\chi_t})^{1/3}(\chi-\chi_t)$, with the matching
  amplitude limit. This avoids cancellation in the closed-form action before
  the universal inverse-action spline is used.

Splining $\theta(p)$, rather than $u(p)=\sec\theta(p)$, is important near the
turning point. The earlier broad spline of $u$ was too coarse just above the
turning point and produced errors such as $\beta\,\Delta z\simeq 5\times10^{-2}$
in high-$L$ cases. With the current $\theta$ spline, the same checks give
$|\beta\,\Delta z|\sim 5\times10^{-8}$.

The cache is universal: it depends only on whether the point is below or above
the turning point, not on $K$, $\ell$, or $\beta$. Runtime evaluation is one
analytic action plus one direct-indexed spline lookup.

## Small-Chi Approximation

For small curvature across the relevant interval, a series approximation to the
Olver map is useful for tests and possible fast paths. Define

$$
\lambda = \sqrt{\ell(\ell+1)},\qquad
\alpha=\frac{\beta}{\lambda},
\qquad
t=\alpha\chi,
\qquad
h=\frac{K}{\alpha^2}.
$$

The tested third-order map is

$$
z=\chi F(t,h), \qquad A=(dz/d\chi)^{-1/2},
$$

with

$$
F =
1-\frac{h}{6}
-\frac{h^2(4t^2+13)}{360}
-\frac{h^3(48t^4+148t^2+737)}{45360},
$$

and

$$
\frac{dz}{d\chi} =
1-\frac{h}{6}
-\frac{h^2(12t^2+13)}{360}
-\frac{h^3(240t^4+444t^2+737)}{45360}.
$$

The approximation is then

$$
\phi_\ell(\chi) \simeq
\frac{A z}{S_K(\chi)}\,j_\ell(\beta z).
$$

This is a series approximation to the full Olver action map. When the final
flat spherical Bessel is evaluated accurately, the full Olver map is normally
more accurate than the small-chi series. Apparent small-chi wins in the current
Fortran comparison are mostly cancellations against the present `bjl` transition
error, not evidence that the series map is better than the full map.

## Shifted-q Approximation

The simpler shifted approximation used for comparison is

$$
\phi_\ell(\chi)\approx
\sqrt{\frac{q}{\beta}}\frac{\chi}{S_K(\chi)}j_\ell(q\chi),
\qquad
q^2=\beta^2-\frac{K\lambda^2}{3},\qquad
\lambda=\sqrt{\ell(\ell+1)}.
$$

It can tie the small-chi expansion at extremely small $\chi_t$ when it is
gated tightly enough, where all tested approximations sit near the same
reference/Bessel error floor. At larger $\chi_t$, especially high $L$, it
rapidly becomes much worse than either Olver map unless a phase-error gate is
applied.

The current diagnostic gates use

$$
q_\lambda^2=\beta^2-\frac{K\lambda^2}{3},
$$

and the endpoint estimates

$$
\alpha=\frac{\beta}{\lambda},\qquad
t_{\max}=\alpha\chi_{\max},\qquad
\epsilon_\phi=\frac{\lambda^2\chi_{\max}^3}{90\beta},
$$

with the shifted-correction estimates

$$
\epsilon_x=
\frac{\lambda t_{\max}(t_{\max}^2+2)}{90\alpha^4},\qquad
\epsilon_A=\frac{|t_{\max}^2-2|}{180\alpha^4},\qquad
\epsilon_\mathrm{shift}=\frac{1}{2}\epsilon_x+\epsilon_A.
$$

## Accuracy Checks

The main harnesses are:

- `fortran/tests/olver_phi_compare.f90`
- `accuracy_plots/olver_hyperspherical/compare_smallchi_olver.py`
- `accuracy_plots/olver_hyperspherical/scan_olver_chi_bins.py`
- `accuracy_plots/olver_hyperspherical/phi_qapprox_peaknorm_cases.f90`

### Full Olver Cached Tests

The full cached action map is tested against `phi_recurs_stable`; it does not
call `compute_olver_z_amp_smallchi`. The focused open/flat/closed grid uses
$\chi_t<0.1$, representative $1\le\ell\le6000$, and samples $\alpha$ from just
above the $\chi_t<0.1$ boundary up to $2\times10^4$. The output is recorded in
`accuracy_plots/olver_hyperspherical/qapprox_fullolver_l1_gates_summary.csv`,
with the later relaxed-gate check in
`accuracy_plots/olver_hyperspherical/qapprox_err1e3_gate_l1_summary.csv`.

For open, near-flat cases with $\chi<0.2$, $\chi_t<0.2$, and
$k\le 2\,\mathrm{Mpc}^{-1}$, the current full Olver implementation removes the
previous turning-point spline error. In the Fortran-end-to-end comparison the
remaining visible floor is now about $4.4\times10^{-5}$ of the peak, consistent
with direct checks of the current Fortran `bjl` transition floor against SciPy
`spherical_jn`.

On the focused grid, the full cached result stays at that same floor over all
tested endpoints:

| $\chi_{\max}$ | worst full-Olver cached error | worst cached-smallchi delta |
| ---: | ---: | ---: |
| 0.1 | $4.4\times10^{-5}$ | $8.8\times10^{-6}$ |
| 0.2 | $4.4\times10^{-5}$ | $8.6\times10^{-6}$ |
| 0.3 | $4.4\times10^{-5}$ | $8.1\times10^{-6}$ |
| 0.5 | $4.4\times10^{-5}$ | $2.9\times10^{-5}$ |
| 0.75 | $4.4\times10^{-5}$ | $3.8\times10^{-4}$ |
| 1.0 | $4.4\times10^{-5}$ | $2.4\times10^{-3}$ |
| 1.5 | $4.4\times10^{-5}$ | $3.5\times10^{-2}$ |
| 2.0 | $4.4\times10^{-5}$ | $4.6\times10^{-2}$ |

The direct cached-smallchi delta is included to show where the small-chi series
is still a useful independent test of the full action map. In this grid the
delta is below $10^{-5}$ through $\chi_{\max}=0.3$ and below
$3\times10^{-5}$ at $\chi_{\max}=0.5$; at larger endpoints the series, not the
full map, is leaving its validity range.

### Small-Chi Series Tests

The production small-chi gate was also scanned without any separate
$\chi_t$ or $\chi_{\max}$ cut, using open/flat/closed cases with
$1\le\ell\le6000$, endpoints up to $\chi_{\max}=3$, and
$\alpha$ from just above 2 to 5000. The tested variable was

$$
\eta_\chi = \frac{\ell^2\chi_{\max}^7}{\beta}.
$$

The initially tested gate $\alpha>2$ and $\eta_\chi<0.3$ is not sufficient:
the worst accepted error is $6.7\times10^{-4}$, from open $\ell=50$,
$\alpha=2.01$, $\chi_{\max}=0.533$. Failures also occur for low closed
multipoles where $\eta_\chi$ is tiny but $\alpha$ is too close to the
curvature scale. Raising the alpha floor fixes these cases. The accepted
fractions below are over the broad scan grid:

| gate | accepted fraction | worst accepted small-chi error |
| --- | ---: | ---: |
| $\alpha>2,\ \eta_\chi<0.3$ | 0.547 | $6.7\times10^{-4}$ |
| $\alpha>2.2,\ \eta_\chi<0.2$ | 0.453 | $2.6\times10^{-4}$ |
| $\alpha>3,\ \eta_\chi<0.3$ | 0.432 | $2.6\times10^{-4}$ |
| $\alpha>3,\ \eta_\chi<0.35$ | 0.455 | $3.3\times10^{-4}$ |

The adopted gate is therefore

$$
\alpha>3,\qquad \eta_\chi < 0.3,
$$

with the threshold divided by the non-flat accuracy boost in production. This
keeps the worst accepted case below $3\times10^{-4}$ in the broad scan, while
the nearby relaxation $\eta_\chi<0.35$ is already slightly too loose.

### Shifted-q Approximation Tests

For the shifted-$q$ approximation above, the same focused open/flat/closed grid
was compared against `phi_recurs_stable`. Without the phase/amplitude gate,
$\chi_t<0.1$ alone is not a useful accuracy condition:

| $\chi_{\max}$ | worst shifted-$q$ error |
| ---: | ---: |
| 0.1 | $2.5\times10^{-3}$ |
| 0.2 | $1.1\times10^{-2}$ |
| 0.3 | $2.0\times10^{-2}$ |
| 0.5 | $5.1\times10^{-2}$ |
| 1.0 | $8.4\times10^{-2}$ |
| 2.0 | $9.7\times10^{-2}$ |

The proposed shifted-correction gate,
$q_\lambda^2>0$ and $\epsilon_\mathrm{shift}<3\times10^{-4}$, works without an
explicit $\ell>20$ cut in this focused grid. The accepted curved-model
$\alpha$ range is shown for each endpoint:

| $\chi_{\max}$ | accepted curved $\alpha$ range | accepted curved cases | worst curved shifted-$q$ error |
| ---: | ---: | ---: | ---: |
| 0.1 | $10$--$2\times10^4$ | 399 | $1.4\times10^{-4}$ |
| 0.2 | $10$--$2\times10^4$ | 357 | $1.3\times10^{-4}$ |
| 0.3 | $10$--$2\times10^4$ | 305 | $1.2\times10^{-4}$ |
| 0.5 | $10.1$--$2\times10^4$ | 237 | $9.6\times10^{-5}$ |
| 0.75 | $14.8$--$2\times10^4$ | 176 | $5.9\times10^{-5}$ |
| 1.0 | $40$--$2\times10^4$ | 142 | $1.9\times10^{-5}$ |
| 1.5 | $160$--$2\times10^4$ | 96 | $5.9\times10^{-6}$ |
| 2.0 | $320$--$2\times10^4$ | 70 | $3.4\times10^{-6}$ |

The worst accepted case is closed $\ell=200$, $\alpha=12$,
$\chi_{\max}=0.1$, with error $1.4\times10^{-4}$ of peak. For
$\ell\le20$, the worst accepted curved shifted-$q$ case is closed $\ell=20$,
$\alpha=12$, $\chi_{\max}=0.3$, with error $1.2\times10^{-4}$.
This is the gate used by the shifted-$q$ near-flat integration checks.

Relaxing the same uncapped proxy to $\epsilon_\mathrm{shift}<10^{-3}$ includes
more curved cases, but is marginally above a $3\times10^{-4}$ target in this
grid:

| $\chi_{\max}$ | curved cases accepted | worst accepted shifted-$q$ error |
| ---: | ---: | ---: |
| 0.1 | 0.942 | $3.0\times10^{-4}$ |
| 0.2 | 0.861 | $3.4\times10^{-4}$ |
| 0.3 | 0.803 | $3.8\times10^{-4}$ |
| 0.5 | 0.651 | $2.8\times10^{-4}$ |
| 0.75 | 0.530 | $2.3\times10^{-4}$ |
| 1.0 | 0.436 | $2.3\times10^{-4}$ |
| 1.5 | 0.318 | $6.2\times10^{-5}$ |
| 2.0 | 0.246 | $5.6\times10^{-5}$ |

The alternative gate $q_\lambda^2>0$, $\epsilon_\phi<10^{-3}$, and
$\epsilon_A<5\times10^{-5}$ is looser at high multipole and is not a
$3\times10^{-4}$ bound:

| $\chi_{\max}$ | worst curved shifted-$q$ error |
| ---: | ---: |
| 0.1 | $5.1\times10^{-4}$ |
| 0.2 | $3.4\times10^{-4}$ |
| 0.3 | $2.6\times10^{-4}$ |
| 0.5 | $1.6\times10^{-4}$ |
| 0.75 | $1.7\times10^{-4}$ |
| 1.0 | $7.6\times10^{-5}$ |
| 1.5 | $1.7\times10^{-5}$ |
| 2.0 | $1.4\times10^{-5}$ |

## Low Multipoles

For $\ell<20$, the high-$\ell$ Langer choice $\lambda=\ell+1/2$ is not the best
small-chi curvature parameter. A scan over $2\le\ell<20$, $\chi\le0.3$, and
$\chi_t\le0.3$ found that using

$$
\lambda=\sqrt{\ell(\ell+1)}
$$

inside the small-chi series improved every scanned case. The worst errors were:

| range | $\lambda=\ell+1/2$ | $\lambda=\sqrt{\ell(\ell+1)}$ |
| --- | ---: | ---: |
| $\chi_t \le 0.10$ | $8.7\times10^{-5}$ | $3.4\times10^{-7}$ |
| $\chi_t \le 0.20$ | $3.4\times10^{-4}$ | $5.2\times10^{-6}$ |
| $\chi_t \le 0.24$ | $4.9\times10^{-4}$ | $1.1\times10^{-5}$ |
| $\chi_t \le 0.30$ | $8.7\times10^{-4}$ | $2.9\times10^{-5}$ |

The worst low-$\ell$ case is $\ell=2$ at the largest tested $\chi_t$. The
public cached path therefore uses the same $\sqrt{\ell(\ell+1)}$ curvature
scale. It remains the full cached Olver action map for all $\ell\ge1$; the
small-chi series is kept as an independent test/diagnostic path.
