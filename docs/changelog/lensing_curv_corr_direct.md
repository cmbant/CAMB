# Direct Curved-Sky Lensing Method 4

## Summary

This change adds a new Fortran lensing option,
`lensing_method_curv_corr_direct = 4`, which evaluates the curved-sky correlation
integrals directly with Gauss-Legendre quadrature in the same overall style as
`camb.correlations.lensed_cls`.

It also adds `lensing_method_optimized = 5`, which uses the old curved-sky
method when `AccurateBB = False` and the direct Gauss-Legendre implementation
when `AccurateBB = True`.

The main implementation work was:

- add method-4 dispatch in `fortran/lensing.f90` for the normal lensing path and
  the custom-spectrum path,
- expose the method constant in Python via `camb._config` and `camb.__init__`,
- allow per-call overrides from Python with
  `CAMBdata.get_lensed_cls_with_spectrum(..., lensing_method=...)` and
  `CAMBdata.get_partially_lensed_cls(..., lensing_method=...)`,
- add a focused regression in `camb/tests/camb_test.py`,
- cache Gauss-Legendre nodes and weights across repeated calls,
- parallelize the direct correlation accumulation over angle points with OpenMP.

## What "custom spectrum" means here

The custom-spectrum path is
`CAMBdata.get_lensed_cls_with_spectrum(clpp, ...)`, where Python passes in a
lensing potential spectrum
`[L(L+1)]^2 C_L^{\phi\phi} / 2\pi` instead of using the one already stored in the
current `CAMBdata` object. This is the cleanest way to time or compare the
lensing step by itself, because the expensive transfer and source calculations
have already been done.

## Method 1 vs Method 4

Method 1, `lensing_method_curv_corr`, is the existing curved-sky correlation
function implementation. It uses the non-perturbative isotropic term together
with a second-order expansion in `C_gl,2`.

Method 4, `lensing_method_curv_corr_direct`, evaluates the correlation function
more directly with Gauss-Legendre sampling over angle and explicit accumulation
of the curved-sky Wigner-d combinations. The Fortran implementation was written
to follow the structure of `camb.correlations.lensed_cls`, but there are still a
few important differences:

- The Fortran code uses CAMB's existing high-`l` template extension when the
  internal lensing convolution needs `l > Params.max_l`.
- The Python reference code used for comparison does not have that Fortran-side
  template machinery unless the input arrays are built with the same internal
  CAMB support.
- Method 4 in Fortran now caches Gauss-Legendre nodes and weights and uses
  OpenMP thread-local accumulation. The Python routine is purely a reference
  calculation here, not a performance target.
- For `AccurateBB = False`, both methods intentionally restrict the angular
  integration range. That keeps TT, EE and TE fast, but low-`l` BB is not meant
  to be high-accuracy in that mode.

## Correct way to compare accuracy

The first round of comparisons used a higher-support reference with requested
`lmax + 1000`, `lens_potential_accuracy = 2`, and then compared ordinary runs at
the original `lmax` against that higher-support result.

That comparison is useful for one question:

- how much total error is present from the full calculation, including the
  physical lensing support cut and convolution margin.

It is not a pure method-accuracy comparison, because the reference contains more
small-scale unlensed and lensing-potential information than the lower-support
run actually had.

For a like-with-like method comparison, the reference has to use the same input
spectra and the same physical scale cut as the Fortran methods. The matched
comparison used here is:

1. build one `CAMBdata` object at support `lmax_out + 1000`,
2. read that object's actual internal support `Params.max_l`,
3. extract the unlensed spectra and `clpp` to exactly that internal `Params.max_l`,
4. run both Fortran methods through `get_lensed_cls_with_spectrum` using that
   same `clpp`,
5. compare them to a high-density Python direct reference built from those same
   unlensed and `clpp` arrays.

This removes the extra-physics mismatch and leaves only the difference between
the numerical lensing algorithms on the same inputs.

The exact benchmark used for the numbers below is left in the working tree as
`scripts/benchmark_lensing_custom.py`. The final tables below were produced with
its `--strict-reference` mode.

## Accuracy Results: Strict Matched Support, AccurateBB = True

Reference setup:

- cosmology: `H0=67.5`, `ombh2=0.0224`, `omch2=0.12`, `tau=0.054`, `mnu=0.06`,
  `As=2.1e-9`, `ns=0.965`,
- support request: `lmax_out + 1000`,
- benchmark CAMB run: `lens_potential_accuracy = 2`, `NonLinear_lens`, `AccurateBB = True`,
- reference CAMB run: separate strict run with
  `AccuracyBoost=lSampleBoost=lAccuracyBoost=3`,
  `min_l_logl_sampling=10000`, and the default validated
  `use_cl_spline_template=True`,
- reference lensing call:
  `correlations.lensed_cls(..., sampling_factor=4.0, theta_max=None, apodize_point_width=14)`.

The BB bands of most interest here are:

- tensor contamination band: `2 <= l <= 399`,
- intermediate BB band: `400 <= l <= 1200`.

The table below shows maximum relative error against the matched-support
reference in those bands, plus the last-500-`l` tail for context.

| `lmax_out` | internal support | method | `BB(l<400)` max | `BB(400..1200)` max | `BB tail` max |
| --- | ---: | --- | ---: | ---: | ---: |
| `2000` | `3150` | method 1 | `1.42e-3` | `6.27e-4` | `1.78e-2` |
| `2000` | `3150` | method 4 | `2.98e-4` | `1.15e-3` | `2.08e-2` |
| `4000` | `5150` | method 1 | `1.48e-3` | `6.65e-4` | `5.12e-2` |
| `4000` | `5150` | method 4 | `4.55e-4` | `6.46e-4` | `5.47e-2` |
| `6000` | `7150` | method 1 | `1.48e-3` | `6.51e-4` | `7.62e-2` |
| `6000` | `7150` | method 4 | `4.59e-4` | `6.42e-4` | `8.24e-2` |

Under this stricter reference, the main robust `AccurateBB=True` result is that
method 4 is clearly better in the tensor-sensitive `l < 400` BB band, while the
`400 <= l <= 1200` BB band is only a near-tie at `lmax=4000` and `6000` and is
slightly worse at `lmax=2000`. For context beyond BB, the full-range `TT`, `EE`
and scaled `TE` maxima remain dominated by the same cutoff-local effect seen in
the BB-tail column above, with method 1 still slightly smaller at the top end
for `lmax=4000` and `6000`.

### What is going wrong near the top cutoff

The large residual in the last few hundred multipoles is mainly a support-margin
problem, not a failure of the core lensing integral.

The reason the extrapolation template does not remove this sensitivity is that
it is only a scaled fiducial tail. CAMB loads a fixed high-`l` template from
`HighLExtrapTemplate_lenspotentialCls.dat`, and the lensing code then matches
its amplitude at `CP%Max_l` before using that shape for `TT`, `EE`, `TE` and
`PP` above the internal cutoff. That is good enough for small numerical
corrections, but it is not the same as having the actual cosmology's true
high-`l` unlensed and lensing-potential spectra available over a wide range.

So when the requested output `l` is too close to the region fed by the scaled
template, the lensed spectrum near the top end inherits the template mismatch.
Increasing the margin matters because it moves the output band farther away from
that approximate tail.

One important correction to the first quick template check is that a default
high-`l` CAMB output is not automatically a better reference than the shipped
template file. At `l > 5000`, the unlensed `TT` and especially `EE` spectra are
already deep in the damping tail, so the dense output is sensitive to CAMB's own
`l` interpolation choices. With the default settings, the high-`l` output can be
off by tens of percent relative to a stricter run, even though the absolute
signal there is tiny.

A more reliable reference was made by increasing
`AccuracyBoost=lSampleBoost=lAccuracyBoost=3`, setting
`min_l_logl_sampling=10000`, and disabling the spline-template interpolation.
That reference is itself well converged: increasing the boosts further to `4`
changes `TT`, `EE`, scaled `TE` and `PP` by only about `1e-4` over
`4000 <= l <= 6000`.

Using that stricter reference, the conclusion is:

- the earlier default-style run was not a safe truth standard for judging the
  template tail;
- with the same high-accuracy settings, turning CAMB's spline-template
  interpolation back on changes the unlensed `TT`, `EE`, scaled `TE` and `PP`
  by only about `1e-5` to `5e-5` over `4000 <= l <= 6000`;
- the shipped `HighLExtrapTemplate_lenspotentialCls.dat` is therefore not the
  limiting error source in that range once the underlying run is itself
  converged.

So the large discrepancies seen in the first quick tail probe were mostly a
comparison against an under-converged default high-`l` output, not evidence that
the shipped template file was poor. The support-margin conclusion above still
holds for the lensed top-end output, but the template itself should be compared
against a stricter reference than the default CAMB run.

In these matched-support tests, the data object was built with an internal
support only about `1150` above the requested output `lmax_out`. That is enough
for the `l < 400` and `400..1200` BB bands, but it is not enough to make the
last-500-`l` tail a clean method comparison.

When the support is increased while keeping the method fixed, the BB-tail error
collapses rapidly:

| `lmax_out` | effective internal margin | method 1 tail max | method 4 tail max |
| --- | ---: | ---: | ---: |
| `4000` | `1150` | `5.31e-2` | `5.66e-2` |
| `4000` | `1750` | `1.28e-2` | `1.50e-2` |
| `4000` | `2350` | `1.72e-3` | `3.19e-3` |
| `6000` | `1150` | `7.93e-2` | `8.55e-2` |
| `6000` | `1750` | `8.95e-3` | `1.30e-2` |
| `6000` | `2350` | `2.98e-3` | `6.58e-6` |

This shows two things.

- The dominant problem near the top cutoff is insufficient support margin.
- Once the support is pushed high enough, method 4 improves very strongly,
  especially when it no longer has to rely on the extrapolated high-`l` tail.

So the right way to improve the top end is first to increase the lensing margin
or internal support. In normal CAMB use, the relevant user-facing control is the
`lens_output_margin` argument to `set_for_lmax`.

## AccurateBB = False

The same matched-support custom-spectrum benchmark was also run for
`AccurateBB = False`, using CPU time, one warm cached run, and four measured
runs.

### Accuracy

In this mode, neither method is acceptable for the tensor-sensitive low-`l` BB
band. The `l < 400` maximum relative error stays around `3e-3` to `4e-3` for
both methods.

For the `400 <= l <= 1200` BB band, the strict-reference runs show only a small
method difference: method 4 is slightly better at `lmax=4000` and `6000`, but
it is worse at `lmax=2000`:

| `lmax_out` | method | `BB(l<400)` max | `BB(400..1200)` max |
| --- | --- | ---: | ---: |
| `2000` | method 1 | `3.89e-3` | `6.23e-4` |
| `2000` | method 4 | `2.22e-3` | `1.15e-3` |
| `4000` | method 1 | `3.98e-3` | `6.77e-4` |
| `4000` | method 4 | `3.35e-3` | `6.47e-4` |
| `6000` | method 1 | `3.98e-3` | `6.64e-4` |
| `6000` | method 4 | `3.71e-3` | `6.43e-4` |

So the current short-range mode is still not suitable for low-`l` BB science,
which is exactly why `AccurateBB=True` exists. The non-accurate mode is mainly a
fast TT/EE/TE option, with only approximate BB. Under the stricter reference,
the `400..1200` BB differences between methods are much smaller than the first
exploratory numbers suggested.

### Timing

One-thread custom-spectrum CPU-time averages are:

| `lmax_out` | method 1 custom | method 4 custom |
| --- | ---: | ---: |
| `2000` | `0.010-0.017 s` | `0.027 s` |
| `4000` | `0.030-0.034 s` | `0.064 s` |
| `6000` | `0.058-0.064 s` | `0.117 s` |

So in `AccurateBB=False`, method 4 is only marginally better in the
`400..1200` band at high `lmax`, but it is still about two to three times slower
than method 1 and does not fix the low-`l` BB problem.

## Short-Range Cut And Apodization Scan

The direct method's `AccurateBB=False` path was checked against the full-range
reference by varying:

- the angular cut `theta_max = pi / range_fac`, with `range_fac` tested at
  `24`, `32` and `40`,
- the apodization width in quadrature points.

Main outcome:

- `range_fac = 24` improves the `400..1200` BB band, but costs noticeably more
  CPU time and worsens the low-`l` tensor band.
- `range_fac = 40` is faster and somewhat better in the tensor band, but it
  degrades the `400..1200` band too much.
- `range_fac = 32` remains the best compromise when both BB bands matter.

For the taper, replacing the point-index half-Gaussian with a C2 smoothstep in
the physical angle gave a cleaner improvement, but only once the transition was
made broader than the nominal Gaussian width. The retained Fortran short-range
direct method uses

```fortran
apodize_theta_width = min(theta_max, 48*pi/lmax)
```

and applies the smoothstep over `theta_max - apodize_theta_width < theta <
theta_max`. Narrower smoothstep widths, including `12*pi/lmax`, did not improve
the old taper: they left the `400..1200` band essentially unchanged and made the
low-`l` tensor band slightly worse. The wider C2 window improves the
tensor-sensitive `l < 400` BB band while leaving `400..1200` BB at the previous
level.


Main points:

- The earlier very large full-range errors were indeed mixing method error with
  physical support and margin error.
- With matched inputs, low and mid `l` agreement is much better than those first
  numbers suggested.
- Method 4 is much better in the tensor-sensitive `l < 400` BB band.
- In the `400 <= l <= 1200` BB band, the strict-reference runs show only a very
  small method difference: method 4 is slightly better at `lmax=4000` and
  `6000`, but not by a large factor, and it is worse at `lmax=2000`.
- Method 1 is slightly better right at the top of the output range for these
  `lmax=4000` and `6000` matched-support tests.
- By `lmax=4000` and `6000`, the dominant residual method disagreement is very
  concentrated near the output cutoff rather than across the full spectrum.

This is the most useful interpretation of the current code state: method 4 gives
clearly better low-`l` BB behavior at the same physical support, has only a
small edge or deficit in the `400..1200` BB band depending on `lmax`, while
method 1 still has a modest edge very close to the high-`l` cutoff.

## Timing Results: Custom-Spectrum Lensing Only

The timings below are only for the isolated lensing calculation,
`get_lensed_cls_with_spectrum`, using the same matched-support setup as above.
No `get_results` timings are included here.

Reported values below are CPU times from `time.process_time()`, after one warmup
call to ensure the Gauss-Legendre cache is populated, averaged over four
measured custom-spectrum runs.

### One thread

| `lmax_out` | method 1 custom | method 4 custom |
| --- | ---: | ---: |
| `2000` | `0.733 s` | `0.841 s` |
| `4000` | `2.410 s` | `1.935 s` |
| `6000` | `4.632 s` | `4.146 s` |

### Four threads

| `lmax_out` | method 1 custom | method 4 custom |
| --- | ---: | ---: |
| `2000` | `0.982 s` | `1.113 s` |
| `4000` | `3.224 s` | `2.652 s` |
| `6000` | `5.799 s` | `4.656 s` |

Main timing points:

- Method 4 is not a speed win at `lmax=2000`.
- At `lmax=4000` and `6000`, method 4 is clearly faster for the isolated
  `AccurateBB` lensing step.
- CPU time is the right measure for work here, but it does not show wall-clock
  parallel speedup. With OpenMP, the four-thread CPU times remain similar to, or
  somewhat larger than, the one-thread CPU times because they sum work across
  threads.
- Gauss-Legendre caching removes the need to regenerate nodes and weights for
  repeated calls at the same quadrature size.

## Parameter Exploration

Several internal tuning changes were explored while bringing method 4 into the
Fortran path.

### 1. Gauss-Legendre caching

Method 4 originally regenerated nodes and weights for every call. A module-level
cache keyed by the quadrature size now reuses them across repeated runs. This is
especially useful for repeated custom-spectrum calls and helps the isolated
lensing timings above.

### 2. OpenMP thread-local accumulation

The direct method now accumulates into thread-local work arrays over angle
points and reduces at the end. This avoids synchronization inside the inner loop
and is the main reason method 4 becomes competitive at high `lmax`.

### 3. Method-1 scratch-array refactor

Method 1 originally wrote each sampled-`l` contribution into a scratch array and
then immediately summed those arrays back into `corr(1:4)`. The disabled spline
interpolation path associated with those scratch arrays was also removed.

The current method-1 code now accumulates the sampled contributions directly in
the `j` loop, with the same exact low-`l` and interpolated high-`l` weights as
before. This is not a large algorithmic change, but it removes a clean piece of
method-1 overhead with no observed regression in the focused test.

### 4. Direct-method sampling factor

The current direct method uses

```fortran
sampling_factor = 1.4_dl*LensAccuracyBoost
```

An additional exploratory bump,

```fortran
if (lmax > 3500) sampling_factor = sampling_factor*1.1_dl
```

was tested and then removed. Under the corrected matched-support benchmark it
did not give any material accuracy improvement, but it did measurably slow the
direct method. For the same A/B, the matched-support accuracy differences were
negligible at the quoted precision.

So the extra high-`l` bump was not retained.

### 5. Method-1 theta-floor retune after the later C2 taper change

This note mainly covers method 4, but the later `AccurateBB=False` work on method 1 changed the
short-range taper enough that the original follow-up `ThetaSampleBoost` floors were worth revisiting.

After switching method 1 to the C2 taper with

```fortran
apodize_width = 0.012_dl
apodize_point_width = nint(apodize_width/dtheta)
```

the earlier high-l `ThetaSampleBoost = 2.6` floor from the follow-up note stopped being the minimal
passing value. With the low-l floor fixed at `1.6`, the current strict-reference
`params_lmax6000.ini` scan gives:

| high-l `ThetaSampleBoost` floor | pass? | score | standard wall time |
| --- | --- | ---: | ---: |
| `2.0` | pass | `0.8687158` | `0.444 s` |
| `2.2` | pass | `0.8687197` | `0.489 s` |
| `2.4` | pass | `0.8687231` | `0.479 s` |
| `2.6` | pass | `0.8687261` | `0.482 s` |

while the no-uplift case (`ThetaSampleBoost = max(., 1.6)` with no `Max_l > 3500` override) fails
TT in the top tail at `L = 5894` with `0.00421` against the `0.003` tolerance. So the current
method-1 short-range path still needs a dedicated high-l uplift, but `2.6` is no longer carrying a
visible accuracy advantage over `2.0` on this check. The retained method-1 setting uses `2.2` as a
small cushion above the minimal tested passing value.

The associated plots are in `accuracy_plots/theta_highl_scan/`, with TT/EE fractional-difference
curves against the strict reference and against the current `2.6` baseline.

## Current Takeaway

Method 4 is now a practical additional option rather than just a direct port of
the Python calculation.

- It is wired through both the normal and custom-spectrum Fortran paths.
- It is selectable from Python, including per-call overrides.
- It is validated by the focused regression test.
- It is faster than method 1 for the high-`lmax`, `AccurateBB=True` custom
  lensing cases that are the main target here.
- In strict matched-support accuracy tests, it clearly improves BB in the
  `l < 400` tensor band, is only marginally different in the `400..1200` band,
  and method 1 still remains slightly better very near the top output cutoff.
- In `AccurateBB=False`, it is still substantially slower than method 1 and does
  not deliver a comparably strong accuracy advantage.

So the present tradeoff is not "method 4 is uniformly more accurate everywhere".
The more accurate statement is:

- method 4 is the faster choice for high-`lmax` `AccurateBB` runs,
- method 4 gives a clear BB gain in the tensor band and only a small change in
  the `400..1200` band away from the extreme top-end cutoff,
- method 1 still has a small accuracy advantage right at the highest output `l`.
