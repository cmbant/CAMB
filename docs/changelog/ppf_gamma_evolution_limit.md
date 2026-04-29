# PPF Gamma Experiment Summary

## Motivation

The PPF dark-energy Gamma equation is numerically fragile for high-redshift, high-`k/H` modes when dark energy is dynamically negligible. With the extra dark-energy source tolerance boost in `cmbmain.f90` commented out:

```fortran
!if (.not. CP%DarkEnergy%is_cosmological_constant) tol1=tol1/5
```

the default PPF branch shows large normalized dark-energy perturbation excursions and increased spectrum sensitivity to `IntTolBoost`. The goal of these experiments is to find a PPF Gamma update that is stable at `IntTolBoost=1`, close to the existing `ckH^2 > 1000` behavior where that behavior is reliable, and avoids the physical/pathological shifts caused by a hard `ckH^2 > 30` cutoff.

## Modes Compared

- `mode 0`: current branch, hard `Gamma=0` only for `ckH^2 > 1000`.
- `mode 1`: for `ckH^2 > 30`, use the algebraic quasi-static fixed point `Gamma_qs = S_Gamma / (1 + ckH^2)^2`, with weak state relaxation controlled by `ppf_test_param`.
- `mode 2`: same fixed point as mode 1, but evolve toward it with a smooth capped relaxation rate controlled by `Gamma_relax_cap`.
- `mode 3`: reference bad case, hard `Gamma=0` for `ckH^2 > 30`.
- `mode 4`: for `ckH^2 > 30`, keep using the evolved `Gamma` state but relax it toward the same quasi-static fixed point with `ppf_test_param`. This avoids the direct algebraic replacement `Gamma = Gamma_qs`.

## Main Test Cases

- Near-LCDM non-crossing model: `w0=-0.95`, `wa=0.15`.
- Near-LCDM crossing model: `w0=-0.998`, `wa=-0.01`.
- Odd stress-test model from the earlier note: `h=0.67`, `mnu=0.2`, `Omegam=0.3`, `Omegab=0.046`, `w0=0.1`, `wa=-0.11`.
- Fluid dark energy is used as a control only for non-crossing cases.

All fractional-difference plots now label their reference explicitly, for example:

- relative to the same PPF mode at `IntTolBoost=4`,
- relative to `fluid(IntTolBoost=4)`,
- or relative to `mode 0(IntTolBoost=4)`.

Plots that include mode 3 hard 30 are clipped where needed so its large failures do not hide the other curves. The clipping is annotated on the affected axes.

## Key Results

For the near-LCDM `k=0.1` mode evolution:

| case | max \|Delta_de\| | max \|v_de\| |
|---|---:|---:|
| mode 0, `IntTolBoost=1` | `3.15e3` | `3.17e5` |
| mode 1, `IntTolBoost=1` | `2.82e-2` | `2.68` |
| mode 2 cap 30, `IntTolBoost=1` | `2.82e-2` | `2.68` |
| mode 3 hard 30, `IntTolBoost=1` | `2.87e-2` | `2.68` |
| mode 4 relax 31, `IntTolBoost=1` | `2.82e-2` | `2.68` |

For near-LCDM spectra sensitivity, comparing `IntTolBoost=1` to `IntTolBoost=4` within the same mode:

| case | lensed TT max frac | lens PP max frac |
|---|---:|---:|
| mode 0 | `6.81e-4` | `8.77e-4` |
| mode 1 | `6.37e-5` | `1.08e-5` |
| mode 2 cap 30 | `6.37e-5` | `1.08e-5` |
| mode 3 hard 30 | `6.37e-5` | `1.08e-5` |
| mode 4 relax 31 | `6.37e-5` | `1.08e-5` |

For the near-LCDM crossing `w0=-0.998`, `wa=-0.01` check, the same stabilization holds:

| case | max \|Delta_de\| | max \|v_de\| |
|---|---:|---:|
| mode 0, `IntTolBoost=1` | `7.56e4` | `1.67e8` |
| mode 0, `IntTolBoost=4` | `2.23e3` | `5.62e6` |
| mode 1, `IntTolBoost=1` | `1.13e-3` | `2.68` |
| mode 2 cap 30, `IntTolBoost=1` | `1.13e-3` | `2.68` |
| mode 3 hard 30, `IntTolBoost=1` | `1.15e-3` | `2.68` |
| mode 4 relax 31, `IntTolBoost=1` | `1.13e-3` | `2.68` |

For that crossing model's spectra sensitivity, comparing `IntTolBoost=1` to `IntTolBoost=4` within the same mode:

| case | lensed TT max frac | lens PP max frac |
|---|---:|---:|
| mode 0 | `4.46e-4` | `9.71e-4` |
| mode 1 | `6.76e-5` | `1.19e-5` |
| mode 2 cap 30 | `6.76e-5` | `1.19e-5` |
| mode 3 hard 30 | `6.76e-5` | `1.19e-5` |
| mode 4 relax 31 | `6.76e-5` | `1.19e-5` |

For the near-LCDM non-crossing fluid control:

- `Delta_cdm` differs from `fluid(IntTolBoost=4)` by only `~5e-6` for mode 1 at `IntTolBoost=1`.
- CMB spectra remain close to fluid away from the lowest multipoles. For mode 1 at `IntTolBoost=1`, max fractional differences for `ell >= 30` are roughly:
  - TT: `3.8e-5`
  - EE: `6.2e-5`
  - lens PP: `4.8e-5`

For the odd stress-test matter power relative to `fluid(IntTolBoost=4)`:

| case | z | k/h | PPF/fluid - 1 |
|---|---:|---:|---:|
| mode 0, `IntTolBoost=4` | 0 | 0.098 | `-5.95e-2` |
| mode 1, `IntTolBoost=1` | 0 | 0.098 | `-5.86e-2` |
| mode 2 cap 30, `IntTolBoost=1` | 0 | 0.098 | `-5.94e-2` |
| mode 3 hard 30, `IntTolBoost=1` | 0 | 0.098 | `-3.27e-2` |
| mode 4 relax 31, `IntTolBoost=1` | 0 | 0.098 | `-5.95e-2` |
| mode 0, `IntTolBoost=4` | 0 | 10 | `-5.07e-3` |
| mode 1, `IntTolBoost=1` | 0 | 10 | `-2.16e-3` |
| mode 2 cap 30, `IntTolBoost=1` | 0 | 10 | `-1.56e-3` |
| mode 3 hard 30, `IntTolBoost=1` | 0 | 10 | `+1.08e1` |
| mode 4 relax 31, `IntTolBoost=1` | 0 | 10 | `-1.55e-3` |
| mode 1, `IntTolBoost=1` | 100 | 10 | `-2.19e-3` |
| mode 2 cap 30, `IntTolBoost=1` | 100 | 10 | `-1.55e-3` |
| mode 3 hard 30, `IntTolBoost=1` | 100 | 10 | `+2.63e-1` |
| mode 4 relax 31, `IntTolBoost=1` | 100 | 10 | `-1.55e-3` |

The odd-model `Delta_cdm` comparison gives mode 1 and mode 2 essentially the same agreement with fluid as mode 0, while mode 3 has a large `k=1` excursion:

| case | k | max \|Delta_cdm/fluid - 1\| |
|---|---:|---:|
| mode 0, `IntTolBoost=4` | 0.1 | `3.29e-2` |
| mode 1, `IntTolBoost=1` | 0.1 | `3.26e-2` |
| mode 2 cap 30, `IntTolBoost=1` | 0.1 | `3.28e-2` |
| mode 3 hard 30, `IntTolBoost=1` | 0.1 | `3.07e-2` |
| mode 4 relax 31, `IntTolBoost=1` | 0.1 | `3.28e-2` |
| mode 0, `IntTolBoost=4` | 1 | `5.03e-3` |
| mode 1, `IntTolBoost=1` | 1 | `5.07e-3` |
| mode 2 cap 30, `IntTolBoost=1` | 1 | `5.03e-3` |
| mode 3 hard 30, `IntTolBoost=1` | 1 | `4.11e-1` |
| mode 4 relax 31, `IntTolBoost=1` | 1 | `5.03e-3` |

For the odd stress-test C_l comparison, the mode 3 hard 30 C_l arrays themselves are all `NaN` for `ell >= 30`. The plot can otherwise look innocuous because matplotlib skips those points, so the legend now marks this case as `all NaN/Inf`. This is the clearest evidence that the hard `ckH^2 > 30` cutoff can hide the near-LCDM instability while causing unacceptable behavior elsewhere.

## Current Takeaway

Mode 4 with `ppf_test_param=31` is now the cleanest replacement candidate. This value matches the original local relaxation rate `1 + ckH^2` exactly at the `ckH^2=30` transition. It captures the same stabilization benefit as mode 1, avoids the hard cutoff's large odd-model matter-power distortions and `NaN` C_l behavior, and does not introduce an algebraic jump by replacing `Gamma` directly with `Gamma_qs`.

Mode 1 used `ppf_test_param=1` in these tests. It is stable, but because it sets `Gamma = Gamma_qs` in the high-`ckH` stress-energy calculation, it is less continuous than mode 4. Mode 2 with `Gamma_relax_cap` behaves very similarly to mode 4 relax 31, but mode 4 has the appealing property that it leaves the original PPF equation unchanged below `ckH^2=30` and matches the original local relaxation rate exactly at the threshold.
