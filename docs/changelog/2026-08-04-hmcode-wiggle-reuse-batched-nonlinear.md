# HMcode BAO wiggle re-use and batched non-linear matter power

## Summary

Three related speed changes:

- `init_wiggle` (BAO wiggle isolation) is now only called for the halofit versions that use the
  de-wiggled spectrum at all, i.e. HMcode-2020 (`mead2020`, `mead2020_feedback`, internally
  `imead==3`). `mead2015`, `mead2016` and `halomodel` never call it.
- For those 2020 versions, the wiggle is extracted once per `HMcode` call, from the available
  spectrum nearest `z = 0`, when the cosmology has nearly scale-independent growth. Otherwise it
  is still extracted separately at every redshift, exactly as before. When it applies is
  controlled by the new `Halofit` parameter `HMcode_wiggle_max_fnu` (default `0.01`), and by a
  new `TDarkEnergyModel` method, `assume_scale_indep_lowz_growth()` (see Gates below).
- Non-linear ratios for the matter power spectrum are now obtained for all redshifts in one
  `GetNonLinRatios` call rather than one call per redshift, so redshift-independent set-up
  (growth table, wiggle extraction, cosmology assignment) is shared.

## Why the wiggle can usually be re-used

`fill_plintab` stores the CAMB linear spectrum divided by the HMcode growth factor, i.e. a
`z = 0`-normalised table, and `p_dewiggle` scales the stored wiggle back up with
`cached_grow(z)**2`. So re-extracting the wiggle at each redshift only makes a difference through
the parts of the growth that are *not* captured by that single scale-independent factor:

- scale-dependent growth from massive neutrinos, and
- growth that does not follow HMcode's own internal growth table.

## Gates

Set in `assign_HM_cosmology`:

```fortran
use_wiggle = halofit_version==halofit_mead2020 .or. halofit_version==halofit_mead2020_feedback
cosm%wiggle_once   = use_wiggle .and. cosm%f_nu < this%HMcode_wiggle_max_fnu .and. &
    CP%DarkEnergy%assume_scale_indep_lowz_growth()
cosm%wiggle_each_z = use_wiggle .and. .not. cosm%wiggle_once
```

The neutrino condition is a new `Halofit` class parameter, with the same `.ini` name, so the
trade-off can be changed without editing the code:

| parameter | default | effect |
| --- | --- | --- |
| `HMcode_wiggle_max_fnu` | `0.01` | largest `f_nu` for which the `z = 0` wiggle is re-used; set to 0 to always re-extract |

The standard `mnu = 0.06 eV` case has `f_nu = om_nu/om_m ~ 0.0045` and qualifies; `mnu = 0.3 eV`
gives `f_nu ~ 0.022` and does not. For example `pars.NonLinearModel.HMCode_wiggle_max_fnu = 0.05`
extends the re-use to `mnu = 0.3 eV`, where it costs `9.6e-5` in the non-linear power (12
redshifts, `0 <= z <= 3`) and saves about 15% of `get_matter_power_spectrum`. Setting
`HMCode_wiggle_max_fnu = 0` reproduces the old per-redshift results to `1e-14`.

The dark energy condition is a new `TDarkEnergyModel` method, `assume_scale_indep_lowz_growth()`,
rather than a `Halofit` parameter, since whether growth is scale-independent is a property of the
dark energy model, not of the non-linear correction being applied to it. It only has to hold for the
low redshifts relevant for late-time structure, not at all `z`. `TDarkEnergyModel`'s base implementation
is conservative and returns `.false.`; overrides in the shipped subclasses are:

- `TDarkEnergyEqnOfState` (covering `TDarkEnergyFluid`, the standard `w`/`wa`/tabulated-`w` fluid
  model, and `TDarkEnergyPPF`) returns `this%cs2_lam > 0.99`. With unit rest-frame sound speed
  (the default, `cs2_lam = 1`), this dark energy is smooth on the sub-horizon scales relevant for
  structure formation regardless of the equation of state, so growth is scale-independent to good
  approximation whether or not `w = -1`. Lowering `cs2_lam` lets the fluid cluster, so the gate
  turns itself off.
- `TQuintessence` (general canonical scalar-field quintessence) returns `.true.`: the field
  perturbation equation has a unit-sound-speed `k^2` term, so for a slowly-rolling potential
  (effective mass `a^2 V''` well below `k^2` at the scales of interest) it is smooth for the same
  reason. `TEarlyQuintessence` and `TAxionEffectiveFluid` (the two early-dark-energy
  implementations, the latter extending `TDarkEnergyModel` directly) also return `.true.`: their
  perturbations are genuinely scale-dependent while oscillating, but by construction that is
  confined to `a` close to their transition scale factor `a_c`, with the early component's density
  fraction negligible by `z = 0`; since the method only needs to hold at the low redshifts
  relevant for late-time structure, well below `zc` for the intended use of these classes, that is
  enough.
- Custom dark energy classes default to `.false.` unless they override the method.

This is a wider default than before: any dark energy model with unit-or-near-unit sound speed and
`f_nu` below the threshold now gets the wiggle re-use automatically, not just an exact
cosmological constant. For `w0 = -0.9, wa = 0.3` (`mead2020`, `mnu = 0.06`) this now costs `2.7e-5`
in the non-linear power and `2.5e-6` in `C_L^phiphi` by default, where it previously cost 0 (exact
match to per-redshift extraction); an `EarlyQuintessence`/`AxionEffectiveFluid` example (`zc =
3000`, `fde_zc = 0.05`) costs `2.9e-5` in the non-linear power. (Python attribute name uses the
`HMCode_` spelling of the existing baryon parameters; the Fortran and `.ini` name is
`HMcode_wiggle_max_fnu`.)

`fill_plintab` also now records which redshift index its table holds (`plin_iz`) and returns
immediately if asked to refill the same one, so the precompute does not duplicate work when
`HMcode` is called with a single redshift.

## Batching the matter power path

`Transfer_GetNonLinRatio_index` previously computed a single-redshift `MatterPowerData` and called
`GetNonLinRatios` on it, so `get_matter_power_spectrum` with N redshifts made N `HMcode` calls,
each of which redid the growth table and (with the gate above) the wiggle extraction. The new
`Transfer_CacheNonLinRatios` fills `State%CAMB_PK` for all transfer redshifts in one call and
caches it, which is what `Transfer_GetUnsplinedNonlinearPower` (used by
`get_nonlinear_matter_power_spectrum` and `get_matter_power_interpolator`) already did; that
routine now shares the same helper. `Transfer_SaveMatterPower` with
`transfer_interp_matterpower = F` goes through `Transfer_GetNonLinRatio_index` for the same reason.

The helper only fills the cache when there is more than one power spectrum redshift, when the
power spectrum redshifts are the full transfer redshift list, and when not `OnlyTransfer`. The
second condition is what excludes non-linear lensing, where the transfer redshifts are the merged
list of power spectrum and non-linear lensing interpolation redshifts. Nothing is lost there:
`MakeNonlinearSources` has already computed the ratios for the whole merged list in a single
`GetNonLinRatios` call before any matter power is requested, and `Transfer_GetNonLinRatio_index`
reads the power spectrum redshifts straight out of that cache. With `NonLinear_both`,
`get_matter_power_spectrum` for three redshifts takes 0.07 ms (a cache lookup), against 7.1 ms
when HMcode has to be run for the same three redshifts with lensing off. So the combined case
already makes exactly one HMcode call, and the only case that still runs HMcode per redshift is
non-linear lensing redshifts being present without `cmbmain` having run (`WantCls = False` with
`NonLinear = lens` or `both`), where batching all the transfer redshifts would compute ~30
interpolation redshifts nobody asked for.

## Accuracy

Non-linear lensing (`NonLinear_both`, `lmax = 2500`, `lens_potential_accuracy = 4`, 30 non-linear
lensing redshifts), maximum relative difference against the previous code for `L >= 2`:

| model | `C_L^phiphi` | lensed TT | EE | BB | TE / sqrt(TT.EE) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `mead2020`, `mnu = 0.06` | `2.5e-6` | `1.0e-7` | `6.8e-8` | `1.5e-6` | `3.2e-8` |
| `mead2020_feedback`, `mnu = 0.06` | `2.6e-6` | `9.9e-8` | `6.7e-8` | `1.6e-6` | `3.2e-8` |
| `mead2020`, `w0 = -0.9`, `wa = 0.3` | `2.5e-6` | `8.6e-8` | `6.0e-8` | `1.5e-6` | `2.8e-8` |
| `mead2020`, `mnu = 0.3` | ~0 | ~0 | ~0 | ~0 | ~0 |
| `mead2015`, `mead2016`, `halomodel` | ~0 | ~0 | ~0 | ~0 | ~0 |

`~0` means float-rounding noise only (`~1e-16` to `~1e-15`, from the underlying Cl calculation
itself): the gate does not apply (`f_nu` above threshold, or a version that never uses the
wiggle), so these fall back to the original per-redshift code. `w0 = -0.9, wa = 0.3` now applies
the gate by default (see Gates above), so it is no longer in the zero-difference group. (Both
builds were compiled from scratch for this comparison.)

Non-linear matter power (`get_matter_power_spectrum`, 12 redshifts `0 <= z <= 3`,
`1e-3 <= k/h <= 10`):

| model | max relative difference |
| --- | ---: |
| `mead2020` / `mead2020_feedback`, `mnu = 0.06` | `2.9e-5` |
| `mead2020`, `w0 = -0.9`, `wa = 0.3` | `2.7e-5` |
| `mead2020`, `EarlyQuintessence` (`zc = 3000`, `fde_zc = 0.05`) | `2.9e-5` |
| `mead2020`, `AxionEffectiveFluid` (`zc = 3000`, `fde_zc = 0.05`) | `2.9e-5` |
| `mead2020` `mnu = 0.3`, `mead2016`, `mead2015`, `takahashi` | `1.4e-14` |
| any version, single redshift | 0 |

The `EarlyQuintessence`/`AxionEffectiveFluid` rows are against the `HMcode_wiggle_max_fnu = 0`
per-redshift result for the same model, rather than against pre-feature code (neither class'
`assume_scale_indep_lowz_growth()` override existed before this whole feature, so there is no
prior behavior to reproduce here). The `~2.9e-5`/`2.7e-5` differences are the growth-scaled
wiggle. They sit on the BAO scales, peaking at `k/h ~ 0.2`, are exactly zero at `z = 0` (where the
wiggle is extracted), and grow to about `3e-5` by `z ~ 1.4`; this is four orders of magnitude
below HMcode's own accuracy. The
`1.4e-14` entries are last-bit rounding: the linear power table for the batched call is built for
all redshifts at once, which the compiler vectorises differently. The `.ini` output path
(`do_nonlinear = 3`, `halofit_version = 9`, redshifts 2, 1, 0) agrees with this, with `3.6e-5`,
`3.3e-5` and `0` for the three redshifts, for both `transfer_interp_matterpower = T` and `F`.

## Timing

`get_matter_power_spectrum` for 12 redshifts, `kmax = 20`, 4 threads, minimum of 15 repeats and
best of three separate runs:

| version | before | after | gain |
| --- | ---: | ---: | ---: |
| `mead2016`, `mnu = 0.06` | 27.8 ms | 24.7 ms | 11% |
| `mead2020`, `mnu = 0.06` | 27.7 ms | 24.6 ms | 11% |
| `mead2020_feedback`, `mnu = 0.06` | 42.4 ms | 39.1 ms | 8% |
| `mead2020`, `mnu = 0.3` | 26.1 ms | 25.9 ms | 1% |
| `takahashi` | 3.6 ms | 3.5 ms | 0 |

`mead2016` gains from the version gate alone (the wiggle it never used); `mead2020` with
`mnu = 0.06` gains from the wiggle being extracted once for all 12 redshifts instead of 12 times.
The `mnu = 0.3` and `takahashi` rows are the control cases where neither gate applies.

For non-linear lensing (30 redshifts in a single `HMcode` call) the wiggle saving is internal to
that one call: HMcode CPU time drops from about 0.32-0.37 s to 0.30-0.32 s, which is under 1% of
a full `get_results` and not separable from run-to-run noise in the total.

A reordering of the `Transfer_GetMatterPowerData` loops (redshift outermost, for contiguous
transfer reads) was tried and reverted: it perturbed every spectrum at the `1e-14` level through
vectorised `log`, and the timing gain it appeared to give turned out to be measurement noise.

## Tests

```text
python -m unittest camb.tests.camb_test      22 tests, OK
pre-commit run --files fortran/halofit.f90 fortran/results.f90 fortran/DarkEnergyInterface.f90 \
    fortran/DarkEnergyFluid.f90 fortran/DarkEnergyQuintessence.f90 camb/nonlinear.py
```

`HMcode_wiggle_max_fnu` also round-trips through `params.write_ini` / `camb.read_ini`.
