# Parallelize halofit/HMcode over redshift

## Summary

`GetNonLinRatios` receives all the matter power redshifts in one call, and the redshifts are
independent, so the loop over redshift is now the outermost OpenMP region for both the standard
halofit branch and HMcode. Previously the redshift loop was serial and only the loops *inside*
each redshift were threaded, which scaled poorly: most of the per-redshift cost is set-up
(`fill_plintab`, `fill_sigtab`, `halomod_init`) made up of small parallel regions separated by
serial stretches.

Results are bit-identical to before, for every halofit version and any number of threads.

## What was blocking it

Per-redshift state lived on the shared `THalofit` object and in the shared `HM_cosmology`:

- `this%imead` was reassigned inside the redshift loop (and, for `mead2020_feedback`, again for
  each of the three response passes), and read by `Delta_v`, `delta_c`, `eta`, `kstar`, `As`,
  `fdamp`, `alpha`, `neff`, `conc_bull`, `halomod`, `p_1h` and `p_2h`.
- `this%om_m`, `om_v`, `fnu`, `omm0`, `acur`, `w_hf`, `wa_hf` were set per redshift for the
  standard halofit fitting formulae.
- `cosm` was rewritten at every redshift by `fill_plintab`, `fill_sigtab` and `init_wiggle`.

## Changes

- `imead` moved from `THalofit` to `HM_tables%imead`, set by `halomod_init(this,imead,z,lut,cosm)`.
  `Delta_v` gains a `lut` argument. The routines that now read only `lut` are `nopass` bindings,
  so their call sites are unchanged.
- The standard halofit background quantities moved into a local `THalofit_zparams`, passed to
  `halofit()`. `THalofit` now holds no per-redshift state and is read-only while a calculation runs.
- Each thread takes its own `HM_cosmology` by intrinsic assignment (`cosm = cosi`), which
  deep-copies the allocatable interpolation tables (~25 kB per thread, against ~4 ms of work per
  redshift), and its own `HM_tables`. The body of one redshift is factored into
  `HMcode_redshift` so it can be called from either the parallel or the serial branch.
- The inner regions (the loops over k and mass, the `PARALLEL SECTIONS` in `halomod_init`,
  `fill_plintab`, `fill_sigtab`) get `IF(HM_par_inner())`, which is `.not. omp_in_parallel()`.
  They therefore still run in parallel for a single redshift, and are serialized when the
  redshift loop is doing the threading.
- For a single redshift, HMcode branches around the outer region in Fortran rather than using an
  `IF` clause on `!$OMP PARALLEL`. Even an inactive region puts the inner regions on libgomp's
  nested-team path, which measured ~2 ms of extra overhead per call.
- `SCHEDULE(DYNAMIC)` on the redshift loops and on `fill_sigtab` (the cost of `sigma_integral`
  varies a lot with R, and of a redshift with how quickly the adaptive integrals converge).

## Fixed: `MatterPowerData_k` index cache

The standard halofit branch called `MatterPowerData_k(CAMB_PK, rk, itf)` with no `index_cache`,
falling back on the `save`d `i_last` in `results.f90`. That is a data race once the redshift loop
is threaded, so it now passes a local `index_cache`, as `wint_pk_table` already did.

## HMcode-2020 feedback response

The response model runs `halomod_init` three times per redshift, for `imead` 3 (HMcode-2020),
4 (denominator) and 5 (numerator). Those three share everything expensive: `sigv`, `sigv100`,
`sig8z`, `sig8z_cold`, the mass/radius/sigma tables, `delta_c` and `Delta_v` (both `_Mead` for
3, 4 and 5), `r_nl`, `neff` (all use `itype=1`), the Bullock collapse redshifts, and the Dolag
correction (`pow=1` for all three). Only the halo-concentration amplitude `A = As(lut,cosm)`
differs, along with the cheap scalars `eta`, `kstar`, `fdamp`, `alpha` and the baryon mass
fraction.

So `halomod_init` is split: `halomod_tables` does the shared work, and is skipped on the second
and third passes (`reuse=`); `fill_conc` fills `lut%c` from the collapse redshifts and the cached
Dolag factors (`lut%dolag_inf`, `lut%dolag_z`), and is cheap enough to redo. Caching the Dolag
factors rather than the scaled concentrations keeps the arithmetic in exactly the original order,
so `lut%c` is bit-identical.

## Timings

Non-linear step only, `mead2020`, `kmax=30`, gfortran, devcontainer, ms per redshift:

| threads | nz=1 before | nz=1 after | nz=96 before | nz=96 after |
|--------:|------------:|-----------:|-------------:|------------:|
| 1       | 3.59        | 3.48       | 4.09         | 4.08        |
| 4       | 2.83        | 1.78       | 2.35         | 1.12        |
| 8       | 3.27        | 2.52       | 3.16         | 0.83        |

At nz=96 that is 2.1x on 4 threads and 3.8x on 8, where before 8 threads were *slower* than 4.
Serial performance is unchanged.

`mead2020_feedback` with `mnu=0.06` at nz=96: 3.60 -> 1.96 ms per redshift on 4 threads (the
set-up re-use accounts for about 9% of that on its own, measured single-threaded).

Standard halofit was entirely serial before; at nz=96 on 4 threads, `takahashi` goes 0.29 ->
0.09 ms per redshift and `casarini` (which does a root solve per redshift in `PKequal`) 0.61 ->
0.19 ms.

## Follow-up

Hoisting `fill_plintab` + `fill_sigtab` out of the redshift loop when growth is scale-independent
removes about half of what is left. It is not bit-identical, and was tried as an opt-in parameter and rejected as only accurate for very low neutrino masses.
