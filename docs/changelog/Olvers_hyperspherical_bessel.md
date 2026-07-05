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

The production implementation is in `fortran/hyperspherical_bessels_olver.f90`.

- `qintegral_exact` evaluates the curved action analytically for open and
  closed models. The inverse trigonometric branches use `atan2` forms to avoid
  quadrant/sign mistakes.
- `invert_flat_action` inverts the universal flat action using hard-coded
  polynomial/asymptotic approximations. The action is written as $q$, with
  $p=(3q)^{1/3}$ so the turning point is regular.
- The evanescent branch uses $q=t-\tanh t$ and $u=\alpha z=\operatorname{sech}t$.
- The oscillatory branch uses $q=\tan\theta-\theta$ and
  $u=\alpha z=\sec\theta$.
- Very close to the turning point, the code uses the analytic local limit of
  the action map,
  $z-z_t\simeq (dS_K/d\chi|_{\chi_t})^{1/3}(\chi-\chi_t)$, with the matching
  amplitude limit. This avoids cancellation in the closed-form action before
  the universal inverse-action approximation is used.
- `phi_olver`, `phi_olver_raw`, and `u_olver` share the same internal reduced
  evaluation path. `u_olver` returns $u_\ell=S_K\phi_\ell$ directly for callers
  that do not need the final trigonometric division.
- The production `phi_olver` path may use the small-$\chi$ map as a fast local
  approximation when its pointwise gate is satisfied. The diagnostic
  `phi_olver_raw` path always uses the full action map, except for the exact
  flat and low-$\ell$ stable-reference paths.

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

## Shifted-ν Approximation

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

### Full Olver Tests

The current broad validation grid compares `phi_olver`, diagnostic
`phi_olver_raw`, `phi_olver_smallchi`, and `phi_langer` against
`phi_recurs`. It uses open/flat/closed cases, every integer
$1\le\ell\le20$, representative larger multipoles up to $\ell=6000$, alpha
targets from near the curvature scale to $2\times10^4$, and
$\chi_{\max}=0.01,0.1,0.3,0.5,1.0,1.57,2.0$. Closed modes are rounded to
allowed integer $\beta$ and skipped when $\beta\le\ell$.

The production `phi_olver` gate uses the simple gate alpha

$$
\alpha_g=\frac{\beta}{\ell}
$$

rather than $\beta/\sqrt{\ell(\ell+1)}$. The action map itself still uses
$\lambda=\sqrt{\ell(\ell+1)}$. With the $\alpha_g$ convention the calibrated
pointwise fallback is

$$
\alpha_g\ge4
\quad\hbox{or}\quad
\frac{\chi}{2\beta}\le2.6\times10^{-2}\quad (K=-1),
\quad\hbox{or}\quad
\frac{\chi}{2(\beta-\ell)}\le6.2\times10^{-3}\quad (K=+1).
$$

The small-$\chi$ fast path inside `phi_olver` uses a stricter pointwise version
of the near-flat integration gate,

$$
\alpha_g>4,\qquad \frac{\ell^2\chi^7}{\beta}<5\times10^{-2}.
$$

This keeps the pointwise `phi_olver` envelope at the original target while still
using the cheaper small-$\chi$ map where it is clearly safe. The broader
integration-level gate remains looser because it is applied to a controlled
near-flat integration interval rather than to arbitrary pointwise calls.

The latest broad-grid comparison gives:

| method | total CPU time | time per evaluation | worst peak-normalized error |
| --- | ---: | ---: | ---: |
| `phi_recurs` | 50.5 s | 9.73 us | reference |
| `phi_olver` | 0.573 s | 0.110 us | $9.9\times10^{-5}$ |
| `phi_olver_smallchi` gated | 0.152 s | 0.072 us | $6.2\times10^{-5}$ |
| `phi_langer` | 0.450 s | 0.087 us | $5.0\times10^{-2}$ |

The raw full action map without any low-alpha fallback is not bounded at this
level on the broad grid: accepting all raw points gives a worst error
$3.0\times10^{-2}$, from a closed low-alpha case around $\ell=30$,
$\beta=31$, and $\chi\simeq1.57$.

Earlier focused full-Olver tests gave smaller errors, about
$4.4\times10^{-5}$, because they were a much narrower near-flat/high-alpha
scan, with $\chi_t<0.1$ and endpoints chosen around the intended near-flat
small-$\chi$ use case. In that restricted region the full action-map error is
mostly hidden below the current Fortran `bjl` transition floor. The broad grid
now deliberately includes low-alpha and large-endpoint cases, where the raw
action map needs either the stable fallback or the stricter pointwise gate.

### Langer And Timing Comparison

The broad-grid timing comparison above also scores the existing `phi_langer`
WKB/Langer approximation from `fortran/bessels.f90`. `phi_langer` is slightly
faster than `phi_olver` on this grid, but is not competitive as a general
replacement: the worst errors are low-alpha curved cases, reaching about
$5.0\times10^{-2}$ of peak at $\ell=2$. Restricting to higher multipoles helps
only gradually:

| multipole cut | worst `phi_langer` error |
| ---: | ---: |
| $\ell\ge20$ | $2.5\times10^{-2}$ |
| $\ell\ge50$ | $2.5\times10^{-2}$ |
| $\ell\ge100$ | $9.7\times10^{-3}$ |
| $\ell\ge500$ | $2.6\times10^{-3}$ |
| $\ell\ge1000$ | $1.3\times10^{-3}$ |

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

### Shifted-ν Approximation Tests

For the shifted-$q$ approximation above, the same focused open/flat/closed grid
was compared against `phi_recurs`. Without the phase/amplitude gate,
$\chi_t<0.1$ alone is not a useful accuracy condition:

| $\chi_{\max}$ | worst shifted-$q$ error |
| ---: | ---: |
| 0.1 | $2.5\times10^{-3}$ |
| 0.2 | $1.1\times10^{-2}$ |
| 0.3 | $2.0\times10^{-2}$ |
| 0.5 | $5.1\times10^{-2}$ |
| 0.75 | $8.4\times10^{-2}$ |
| 1.0 | $8.4\times10^{-2}$ |
| 1.5 | $8.3\times10^{-2}$ |
| 2.0 | $1.1\times10^{-1}$ |

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
public Olver map therefore uses the same $\sqrt{\ell(\ell+1)}$ curvature scale
inside the action and small-$\chi$ maps. The runtime gates use the simpler
$\beta/\ell$ alpha convention, with a retuned alpha cutoff, because that is
sufficient for deciding when to use the small-$\chi$ fast path or stable
fallback.
