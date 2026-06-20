# How `l_max` is set in CAMB

CAMB has several distinct multipole cut-offs that interact: the user-facing
`l_max` for scalars and tensors, a `lens_output_margin` that controls how
far below `Max_l` the lensed output is guaranteed to reach, an `l_max`
used while extrapolating the unlensed spectrum into the lensing
convolution, and a separate `l_max` for the Boltzmann hierarchies. This
document spells out where each is set and how the values relate.

## 1. The user-facing cut-offs on `CAMBparams`

The Fortran type `CAMBparams` (mirrored in Python as
[`camb.model.CAMBparams`](../camb/model.py)) carries these fields:

- `min_l` — minimum L for the scalar `C_L` (1 or 2; `L=1` dipoles are
  Newtonian gauge).
- `max_l` — `l_max` for the **scalar** `C_L`.
- `max_l_tensor` — `l_max` for the **tensor** `C_L`.
- `lens_output_margin` (default `200`) — number of L below `max_l` down
  to which lensed `C_L` are guaranteed to be output. The unlensed
  spectrum used in the lensing convolution is treated as reliable to
  `max_l − (lens_output_margin − 50)`, so this field also fixes the
  lensed convolution margin (always 50 less than `lens_output_margin`;
  see §3).
- `max_eta_k` — maximum `k η₀` used in the scalar transfer integration;
  this is what actually controls `k_max`.
- `max_eta_k_tensor` — the same for tensors.

`max_l` is the value used internally everywhere as "the scalar `l_max`".
Lensed CMB spectra are *output* below this (see §3), and unlensed spectra
are output up to and including `max_l`.

These fields are exposed both from `.ini` files and from the Python
interface.

## 2. Setting `l_max` from a `.ini` file

The `.ini` reader is in [`fortran/camb.f90`](../fortran/camb.f90) (see
`CAMBparams_ReadParams`, lines ~387–420):

| `.ini` key            | Fortran field             | Default if not set                                   |
|-----------------------|---------------------------|------------------------------------------------------|
| `l_min`               | `P%Min_l`                 | `2`                                                  |
| `l_max_scalar`        | `P%Max_l`                 | *required when `get_scalar_cls = T` or `get_vector_cls = T`* |
| `k_eta_max_scalar`    | `P%Max_eta_k`             | `2.5 * l_max_scalar`                                 |
| `lens_output_margin`  | `P%lens_output_margin`    | `200`                                                |
| `l_max_tensor`        | `P%Max_l_tensor`          | *required when `get_tensor_cls = T`*                 |
| `k_eta_max_tensor`    | `P%Max_eta_k_tensor`      | `max(500, 2 * l_max_tensor)`                         |

The sample [`inifiles/params.ini`](../inifiles/params.ini) sets:

```
l_max_scalar       = 2200
# k_eta_max_scalar = 4000    (commented out, so default 2.5*l_max_scalar used)
# lens_output_margin = 200   (commented out, so default 200 used)
l_max_tensor       = 1500
k_eta_max_tensor   = 3000
```

The comments in that file capture the operational rule of thumb:

- `C_L` within ~50 of `l_max_scalar` are inaccurate (5%-level), so go
  higher than you need.
- Lensed power spectra are guaranteed to be output up to
  `l_max_scalar − lens_output_margin`, and at most up to
  `l_max_scalar − (lens_output_margin − 50)` (see §3).
- For accurate lensed BB you need `l_max_scalar > 2000` and
  `k_eta_max_scalar > 10000`.
- For an accurate lensing potential you also need
  `k_eta_max_scalar > 10000`.
- Otherwise the default `k_eta_max_scalar = 2.5 * l_max_scalar` is enough
  (matching the Python `set_for_lmax` default `k_eta_fac = 2.5`).

When `do_lensing = T`, the `.ini` user is expected to set
`l_max_scalar` high enough that `l_max_scalar − lens_output_margin`
covers the L range they actually care about.

## 3. The lensed output range

The lensed output upper limit lives in `CLout%lmax_lensed` and is
computed identically by all lensing methods. From
[`fortran/lensing.f90`](../fortran/lensing.f90) (lines ~251–256 for
method 1, ~739–743 for method 4):

```fortran
max_lensed_ix = lSamp%nl-1
do while(lSamp%l(max_lensed_ix) > CP%Max_l - (CP%lens_output_margin - lens_convolution_gap))
    max_lensed_ix = max_lensed_ix - 1
end do
CLout%lmax_lensed = max(lSamp%l(max_lensed_ix), CP%Max_l - CP%lens_output_margin)
```

where `lens_convolution_gap = 50` is a parameter in
[`fortran/config.f90`](../fortran/config.f90) — the fixed gap between
the lensed output margin and the unlensed-spectrum convolution margin.
So if we write `M = CP%lens_output_margin` and `G = 50`:

- **Upper bound:** `lmax_lensed ≤ Max_l − (M − G)`, because the
  while-loop walks the sparse L-sampling grid (§6) back to a point that
  is at most `Max_l − (M − G)`.
- **Lower bound (floor):** `lmax_lensed ≥ Max_l − M`.

So `lmax_lensed ∈ [Max_l − M, Max_l − (M − G)]` — the exact value
depends on where the sparse L-sampling falls in that band.

With the default `M = 200` and `G = 50`, that becomes
`[Max_l − 200, Max_l − 150]`.

### Why `set_for_lmax` adds the full margin

`set_for_lmax(..., lens_output_margin=M)` sets both `pars.max_l =
user_lmax + M` and `pars.lens_output_margin = M`. With `Max_l − M =
user_lmax`, the floor in `lmax_lensed` guarantees the lensed output
reaches the user's target `L` even in the worst sparse-sampling case.
If the floor were e.g. `Max_l − (M − G)`, the lensed output would
fall up to `G = 50` L short of the request in the worst case.

The Python and Fortran sides used to use hardcoded margin constants here;
they are now derived from `CP%lens_output_margin` so
that changing the margin from one side (Python or `.ini`) propagates
through to the lensing convolution and the L-sampling logic.

### Extrapolation beyond `Max_l`

When convolving the unlensed `C_L` with the lensing potential, the
unlensed `C_L` is extrapolated past `Max_l` using a fiducial high-L
template up to
`min(lmax_extrap_highl, output_lmax + lens_convolution_gap + extrap_margin)`,
where `output_lmax = Max_l - lens_output_margin`,
`extrap_margin = max(750, ceiling((0.45 * output_lmax + 400) * LensAccuracyBoost))`
when `AccuracyTarget > 0`, and `lmax_extrap_highl = 8000` is a
`parameter` constant in [`fortran/config.f90`](../fortran/config.f90).

## 4. Setting `l_max` from Python

There are two convenience setters; both end up writing to `max_l`,
`lens_output_margin`, `max_l_tensor`, `max_eta_k`, `max_eta_k_tensor`.

### 4a. `CAMBparams.set_for_lmax`

[`camb/model.py`](../camb/model.py), method `set_for_lmax`. Signature:

```python
pars.set_for_lmax(
    lmax,
    max_eta_k=None,
    lens_potential_accuracy=None,
    lens_output_margin=200,
    k_eta_fac=2.5,
    lens_k_eta_reference=18000.0,
    nonlinear=None,
)
```

Behaviour:

- `lmax` is the `l_max` the user wants the **output** lensed scalar
  spectrum to be accurate to.
- `pars.lens_output_margin` is set to `lens_output_margin` so the
  Fortran lensing code uses the matching floor (§3).
- If `DoLensing` is true, the internal `max_l` is set to `lmax +
  lens_output_margin` (default `+200`). If `DoLensing` is false,
  `max_l = lmax`.
- `max_eta_k` is set to `max_eta_k` if given, otherwise
  `k_eta_fac * max_l` (default `k_eta_fac = 2.5`).
- If `lens_potential_accuracy is None`, CAMB uses the automatic high-accuracy
  default `max(4, (lmax - 1500) / 500)`. Set `lens_potential_accuracy=0` to
  reproduce the previous low-accuracy default.
- The automatic rule is calibrated for target accuracy in the lensed CMB
  spectra. The lensing potential spectrum may not be equally converged at every
  multipole, but is expected to be easily good enough for CMB lensing
  reconstruction applications.
- If `lens_potential_accuracy > 0`, `max_eta_k` is further bumped up to at
  least `lens_k_eta_reference * lens_potential_accuracy` (default
  `18000 * lens_potential_accuracy`). This is what controls `k_max` for
  accurate lensing.
- At high `L`, the remaining numerical error in the lensed CMB spectra can be
  much smaller than the theoretical uncertainty from non-linear matter-power
  modelling.
- `max_l_tensor` is **not** touched here — set it directly on `pars` if
  you want tensors.

The docstring warns that this does *not* fix the output `L` range; the
spectra may be returned above `lmax`. Use the `lmax` argument on the
`get_*_cls` accessors (§5) to truncate output.

### 4b. `camb.set_params(..., lmax=..., lens_potential_accuracy=..., max_eta_k=...)`

[`camb/camb.py`](../camb/camb.py), function `set_params`. This is a
convenience wrapper that forwards relevant kwargs to the individual
setters in a fixed order; `lmax`, `lens_potential_accuracy`,
`lens_output_margin`, `k_eta_fac`, `lens_k_eta_reference`, `nonlinear`,
and `max_eta_k` are routed to `CAMBparams.set_for_lmax`. Anything else
with a matching attribute name (`max_l_tensor`, `max_eta_k_tensor`,
`min_l`, `min_l_logl_sampling`, …) is set directly on the matching
field at the end.

### 4c. Direct assignment

For full control, write the fields directly:

```python
pars.max_l = 2500
pars.max_l_tensor = 1500
pars.max_eta_k = 10000.0
pars.max_eta_k_tensor = 3000.0
pars.lens_output_margin = 200
pars.min_l = 2
```

`set_for_lmax` exists to do the bookkeeping (margin, `max_eta_k`
defaults, non-linear lensing flag) on your behalf; direct assignment
skips it.

## 5. `lmax` on the Python output accessors

[`camb/results.py`](../camb/results.py) defines `_lmax_setting`
(lines ~447–458) and the per-spectrum accessors. The rule is:

- For lensed outputs (`get_total_cls`, `get_lensed_scalar_cls`,
  `get_cmb_power_spectra` default, `get_cmb_correlation_functions`,
  `save_cmb_power_spectra`), the "calculated" `l_max` is
  `f_get_lmax_lensed()`, i.e. `CLout%lmax_lensed` from the Fortran
  lensing module (§3).
- For unlensed outputs (`get_unlensed_scalar_cls`, `get_unlensed_total_cls`,
  `get_lens_potential_cls`, `get_tensor_cls`,
  `get_unlensed_scalar_array_cls`), the "calculated" `l_max` is just
  `Params.max_l` (or `max_l_tensor` for tensors — note `get_tensor_cls`
  defaults its `lmax` to `Params.max_l_tensor`).
- If the user passes `lmax` larger than the calculated value, a warning
  is logged and the spectrum is returned padded/zeroed above the
  calculated range.
- `get_cmb_power_spectra` uses the *lensed* calculated `lmax` for **all**
  spectra it returns; if you want the unlensed and lensing-potential
  spectra to the higher unlensed `lmax`, call the individual accessors.

## 6. Internal `l_max` used inside the solver (not user-tunable)

These are not the same as `max_l`; they control the Boltzmann hierarchy
length per `k`:

- [`fortran/equations.f90`](../fortran/equations.f90):
  `max_l_evolve = 256` is the compile-time cap for any one mode's
  hierarchy length. Per-mode `EV%lmaxg`, `EV%lmaxgpol`, `EV%lmaxnr`,
  `EV%lmaxnu` are chosen as functions of `q` and
  `Accuracy%lAccuracyBoost` (lines ~870–940). The aggregate
  `EV%MaxlNeeded` is required to stay below `max_l_evolve`.
- The L-sampling on which `C_L` is interpolated is built in
  `lSamples_init` in [`fortran/results.f90`](../fortran/results.f90)
  (line ~1336). It runs from `CP%Min_l` up to `Max_l`, with denser
  sampling controlled by `Accuracy%lSampleBoost` and the
  `CP%min_l_logl_sampling` parameter (sparser log spacing above that L).
  When `DoLensing` is on, the sampler also tries to insert a point at
  `Max_l − (lens_output_margin − 50)` so that the lensed convolution
  has a sample to interpolate from there.
- Inside lensing (`lensing.f90`), the per-method internal extrapolation
  upper bound is
  `min(lmax_extrap_highl, output_lmax + lens_convolution_gap + extrap_margin)`
  (see §3).

## 7. Quick checklist

If you want the **lensed** scalar `C_L` accurate to some target `L*`:

1. Either call `pars.set_for_lmax(L*, lens_potential_accuracy=k)` —
   which sets `pars.max_l = L* + 200`, `pars.lens_output_margin = 200`,
   and bumps `max_eta_k` — or set `pars.max_l = L* + 200`,
   `pars.lens_output_margin = 200`, and `pars.max_eta_k`
   (≥ `2.5 * pars.max_l`, and ≥ `18000 * k` for accurate lensing) by
   hand. The two margins must agree.
2. Read the lensed spectrum back with
   `results.get_lensed_scalar_cls(lmax=L*)`; values above `L*` are
   zero-padded.

If you want the **unlensed** scalar `C_L` accurate to `L*`:

1. Set `pars.max_l = L* + 50` or so (the last ~50 are inaccurate; see
   the `.ini` comment in §2). Default `max_eta_k = 2 * max_l` suffices.
2. Read back with `results.get_unlensed_scalar_cls(lmax=L*)`.

If you are reading parameters from a `.ini` file, `l_max_scalar` plays
the role of `pars.max_l` directly. There is no automatic `+200` margin
applied for you on the `.ini` path, so set `l_max_scalar` high enough
that `l_max_scalar − lens_output_margin` covers the L range you
actually care about for the lensed output. The default
`lens_output_margin = 200` can be overridden with the same-named key.
