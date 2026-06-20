# Matter-Power Transfer Accuracy Review

Branch scope: the template-residual interpolation changes described here were
implemented on the separate `mpk-interpolation` branch/prototype
(`550f93f4 mpk interpolation`). They are not active in the current `devel_acc`
checkout unless that branch is merged. The hierarchy-tuning discussion below is
still useful context for v2 matter-power accuracy.

This note summarizes the follow-up review of the matter-power interpolation change
and the targeted transfer hierarchy and massive-neutrino accuracy retuning.

## Executive summary

- The templated Python matter-power interpolator is the right default for public
  arbitrary-`k` calls through `CAMBdata.get_matter_power_interpolator`.
- The public regularly-spaced Python matter-power path now uses the same
  templated interpolator by default. The legacy Fortran interpolation path is
  still available with `use_spline_template=False`.
- The earlier high-precision transfer matter-power failures are controlled by
  the photon hierarchy floors plus the massless-neutrino hierarchy depth in the
  affected transfer range. The older `lmaxnr` onset at `q = 0.10` was too late
  for the HMcode low-`k` matter-power case.
- The accepted `lmaxnr` change is a monotone transfer/21cm ramp:
  `8*lAccuracyBoost` at `q = 0.010`, `14*lAccuracyBoost` at `q = 0.05`, and
  `45*lAccuracyBoost` at `q = 0.12`, with the existing `lmaxnr >= lmaxnu`
  consistency floor applied afterward.
- The HMcode low-`k` matter-power mismatch improves from `7.704e-3` to
  `1.171e-3` with timing indistinguishable from noise in the focused transfer
  reproducer.
- `params_tranfer_delta10.ini` exposed a checker-reference issue for fixed
  native transfer grids: the standard `transfer_k_per_logint = 10` output passes
  when compared to a dense boosted native-grid reference, but failed when the
  boosted reference copied the same sparse output grid.
- In the generated `transfer_high_precision=T`, `get_transfer=T` sweep, there
  are no remaining over-tolerance matter-power cases after using a dense boosted
  reference grid for fixed-output transfer checks. Remaining failures in that
  subset are CMB/lensing rows or one interpolation-range error in a curved
  model.
- Tensor-only CMB checks now compare only through `max_l_tensor`, use
  tensor-only TT/TE/EE tolerance bands, and no longer produce the old TE
  `0/0` false failure.

## Python and Fortran interpolation paths

The Python APIs now split into three relevant paths:

- `get_linear_matter_power_spectrum` and
  `get_nonlinear_matter_power_spectrum` return the native unsplined Fortran
  transfer grid.
- `get_matter_power_interpolator` builds a Python `RectBivariateSpline` from the
  native grid. With `use_spline_template=True`, it splines the residual after
  subtracting the bundled high-accuracy template.
- `get_matter_power_spectrum` builds its regular log-spaced output from the
  same public interpolator by default, so arbitrary-`k` Python outputs share the
  same template treatment.
- `get_matter_power_spectrum(..., use_spline_template=False)` preserves the
  legacy Fortran `CAMBdata_GetMatterPower` wrapper, which calls
  `Transfer_GetMatterPowerD`.

The Fortran `Transfer_GetMatterPowerD` path constructs

```text
log(T_i T_j k 2 pi^2 h^3)
```

on the transfer nodes, cubic-splines that quantity in `log(k/h)`, and multiplies
the primordial scalar power after interpolation. It does not subtract the new
template shape. This means Python arbitrary-`k` interpolator calls and Fortran
regular-grid matter-power output are intentionally not the same interpolation
algorithm only when the legacy `use_spline_template=False` path is requested.

Other Fortran-side consumers also remain on the native Fortran matter-power
data structures:

- `Transfer_GetSigmaRArray` integrates using `MatterPowerData_k`.
- nonlinear-ratio setup uses `Transfer_GetMatterPowerData`.
- saved INI matter-power output uses either native unsplined nodes or
  `Transfer_GetMatterPowerS/D`.
- CMB/lensing solver paths use solver-side transfer data, not the Python
  template.

Suggested consistency policy:

- Keep `camb check_accuracy` testing the public Python interpolator for
  automatic-grid arbitrary-`k` comparisons.
- Keep fixed `transfer_k_per_logint != 0` checks on their requested output nodes;
  otherwise the diagnostic mixes requested fixed output with arbitrary dense
  interpolation quality.
- Keep `get_matter_power_spectrum(..., use_spline_template=False)` as the
  compatibility path for non-smooth primordial spectra or consumers that need
  the old Fortran interpolation exactly. A Fortran template implementation is
  only needed if INI-output or direct Fortran regular-grid output must also
  match the Python templated interpolation.

## Template derivative fix

The templated interpolator changes the value returned by
`PK(z, log(k))`, so derivative calls must include the derivative of the added
template. The fix adds the piecewise-linear template slope for `dy=1` calls and
leaves derivatives in `z` unchanged because the template is independent of
redshift. A finite-difference regression test now covers this.

## Hierarchy tuning

The original staged transfer-high-precision hierarchy additions were:

- `lmaxnu >= 28` for `0.03 < q < 0.2`;
- transfer-only `lmaxg >= 7` for `0.03 < q < 0.04`;
- mixed CMB+transfer `lmaxg >= 12` for `0.03 < q < 0.2`;
- `lmaxnr >= 50` for `0.03 < q < 0.2`.

Initial A/B checks, before the HMcode-specific low-`k` follow-up, showed:

| Variant | Key result |
| --- | --- |
| remove transfer-only photon floor | `params_tranfer_highprec` mid-`k` matter power fails at `1.666e-3` |
| remove mixed high-precision `lmaxg` floor | `params_all_nonlinhigh` matter power fails at `1.080e-3` |
| remove `lmaxnu >= 28` | tested target matter rows still pass |
| remove `lmaxnr >= 50` | tested target matter rows still pass and slightly improve in mixed cases |
| use `35*lAccuracyBoost` above the low-`q` cut | mixed high-precision matter rows fail at `1.059e-3` |
| keep low-`q` reduction for `q < 0.10`, use `45*lAccuracyBoost` above | tested target matter rows pass |
| keep low-`q` reduction for `q < 0.12`, use `45*lAccuracyBoost` above | tested target matter rows pass |
| keep low-`q` reduction for `q < 0.13`, use `45*lAccuracyBoost` above | mixed high-precision matter rows fail at `1.061e-3` near `k/h = 0.1855` |
| keep low-`q` reduction for `q < 0.15` | mixed high-precision matter rows fail at `1.096e-3` |
| ramp `lmaxnr` as `min(45, max(14, slope*q))*lAccuracyBoost`, slope `300..420` | tested target rows pass, but this does more low-`q` work than the threshold rule |
| ramp `lmaxnr` from the existing low-`q` value to `45*lAccuracyBoost` over `0.10 <= q < 0.12` | tested target rows pass, avoids the `q = 0.12` discontinuity, and keeps transfer-only timing essentially unchanged |
| broaden transfer-only `lmaxg >= 7` to `0.02 < q < 0.08` | `params_tranfer_highprec` fails at `1.062e-3` near `k/h = 0.0292` |
| remove the pre-`0.04` photon floor, compared to `lAccuracyBoost=3`, `DoLateRadTruncation=F` only | `params_tranfer_highprec` fails at `1.755e-3` near `k/h = 0.0562` |
| ramp `lmaxg` from `7` to `10` over `0.03 < q < 0.04`, same isolated reference | passes, but worsens the maximum to `9.435e-4` |
| ramp `lmaxg` from `7` to `10` over `0.03 < q < 0.05`, same isolated reference | fails at `1.290e-3` |
| use constant `lmaxg >= 8` over `0.03 < q < 0.04`, same isolated reference | fails at `1.272e-3` |
| apply common boost-scaled high-precision transfer floors after `TransferOnly` caps | `params_tranfer_highprec` passes at `9.018e-4` against the isolated reference and `9.395e-4` against the default reference |

The accepted hierarchy shape after the HMcode low-`k` follow-up is therefore:

- apply the same boost-scaled high-precision transfer floors whether or not the
  run is transfer-only: `lmaxg >= 12*lAccuracyBoost` and
  `lmaxgpol >= 5*lAccuracyBoost` for `0.03 < q < 0.2`;
- drop the unneeded `lmaxnu >= 28`;
- drop the unneeded `lmaxnr >= 50`;
- for `transfer_high_precision` or 21cm, make the massless-neutrino hierarchy
  monotone with `q`: ramp from `8*lAccuracyBoost` at `q = 0.010` to
  `14*lAccuracyBoost` at `q = 0.05`, then to `45*lAccuracyBoost` at
  `q = 0.12`;
- keep `45*lAccuracyBoost` for `q >= 0.12`;
- keep the existing final consistency floor `lmaxnr >= lmaxnu`, so
  high-precision massive-neutrino transfer runs retain a sufficiently deep
  massless hierarchy for the relativistic-correction subtraction.

This targets the expensive massless-neutrino hierarchy smoothly with physical
scale: large scales still use a shallower hierarchy, equality/BAO scales get
enough depth before the old high-`q` target, and massive-neutrino runs remain
self-consistent because the `lmaxnu` floor is applied afterward.

Final HMcode low-`k` check:

| Variant | Focused HMcode low-`k` max error | HMcode cpu mean | `params_tranfer_highprec` cpu mean |
| --- | ---: | ---: | ---: |
| baseline | `7.704e-3` | `6.85s` | `0.766s` |
| ramp from `q = 0.020` | `4.076e-3` | `6.92s` | `0.761s` |
| ramp from `q = 0.015` | `2.314e-3` | `6.70s` | `0.781s` |
| final ramp from `q = 0.010` | `1.171e-3` | `6.62s` | `0.772s` |

Final and earlier variant artifacts:

- [HMcode_cos26_hm66_lmaxnr_cases_a0p4_a1p0.png](../../accuracy_plots/lmaxnr_hmcode_isolation/HMcode_cos26_hm66_lmaxnr_cases_a0p4_a1p0.png)
- [HMcode_cos26_hm66_lmaxnr_heatmap_baseline_vs_q010.png](../../accuracy_plots/lmaxnr_hmcode_isolation/HMcode_cos26_hm66_lmaxnr_heatmap_baseline_vs_q010.png)
- [check_accuracy_q010_highprec/matter_power_errors.png](../../accuracy_plots/lmaxnr_hmcode_isolation/check_accuracy_q010_highprec/matter_power_errors.png)
- [check_accuracy_q010_all_nonlinhigh/matter_power_errors.png](../../accuracy_plots/lmaxnr_hmcode_isolation/check_accuracy_q010_all_nonlinhigh/matter_power_errors.png)

Earlier sweep artifacts:

- [lmaxnr_q_profiles.png](../../accuracy_plots/lmaxnr_q_scaling/lmaxnr_q_profiles.png)
- [params_tranfer_highprec_matter_variant_overlay.png](../../accuracy_plots/lmaxnr_q_scaling/params_tranfer_highprec_matter_variant_overlay.png)
- [params_all_nonlinhigh_matter_variant_overlay.png](../../accuracy_plots/lmaxnr_q_scaling/params_all_nonlinhigh_matter_variant_overlay.png)
- [variant_timing_and_worst_rows.csv](../../accuracy_plots/lmaxnr_q_scaling/variant_timing_and_worst_rows.csv)
- [final_onset_block/params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmaxg_isolated_reference/final_onset_block/params_tranfer_highprec/matter_power_errors.png)
- [no_pre04/params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmaxg_isolated_reference/no_pre04/params_tranfer_highprec/matter_power_errors.png)
- [ramped_onset/params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmaxg_isolated_reference/ramped_onset/params_tranfer_highprec/matter_power_errors.png)
- [ramp_to_005/params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmaxg_isolated_reference/ramp_to_005/params_tranfer_highprec/matter_power_errors.png)
- [const8_pre04/params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmaxg_isolated_reference/const8_pre04/params_tranfer_highprec/matter_power_errors.png)
- [isolated_params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmax_flow_common_floor/isolated_params_tranfer_highprec/matter_power_errors.png)
- [default_params_tranfer_highprec/matter_power_errors.png](../../accuracy_plots/lmax_flow_common_floor/default_params_tranfer_highprec/matter_power_errors.png)

The low-`L` EE row near `L = 25` is not a transfer fix. It appears when the
reference changes only `lAccuracyBoost` and remains with `IntTolBoost` or
`lSampleBoost`; raising `AccuracyBoost` to `2` brings that low-`L` EE row below
tolerance, so it belongs to the CMB/polarization accuracy family.

## Fixed-grid `params_tranfer_delta10.ini`

The `params_tranfer_delta10.ini` matter-power failure is not fixed by the
`lmaxnr` ramp and is not a massive-neutrino hierarchy-depth issue. It uses
`transfer_k_per_logint = 10`, so `check_accuracy` compares the standard run on
its requested native transfer nodes against the boosted run interpolated onto
those nodes.

The standard and boosted fixed native grids have a small phase shift. For
`transfer_k_per_logint = 10`, the median boosted/native same-index shift is
`6.876e-3` in `k`, or `0.0685` of a log-grid step. Comparing on the boosted
nodes or on a dense common interpolation grid gives the same qualitative
failure if the boosted reference is forced to use the same sparse native output
grid. However, keeping the standard run at `transfer_k_per_logint = 10` and
raising only the boosted reference grid to `100` samples per log interval gives
a passing comparison. The standard output nodes are therefore stable; the
checker was using an under-resolved reference interpolant.

Same-grid sampling sweep, where the standard and boosted reference both use the
listed `transfer_k_per_logint`:

| `transfer_k_per_logint` | Status | Max at checker nodes | Dense common-grid max |
| ---: | --- | ---: | ---: |
| `0` | PASS | `9.397e-4` | n/a |
| `5` | FAIL | `4.050e-2` | `4.264e-2` |
| `10` | FAIL | `3.754e-3` | `3.215e-3` |
| `15` | FAIL | `6.660e-3` | `5.779e-3` |
| `20` | FAIL | `1.799e-3` | `1.798e-3` |
| `30` | FAIL | `1.088e-3` | `1.214e-3` |
| `50` | PASS | `9.477e-4` | `9.420e-4` |
| `100` | PASS | `9.650e-4` | `9.425e-4` |

Keeping the standard run fixed at `transfer_k_per_logint = 10` and varying only
the boosted reference grid gives:

| Boosted reference `transfer_k_per_logint` | Max at standard nodes |
| ---: | ---: |
| `10` | `3.754e-3` |
| `20` | `1.681e-3` |
| `30` | `1.061e-3` |
| `50` | `8.502e-4` |
| `100` | `8.257e-4` |
| `200` | `8.206e-4` |
| `500` | `8.188e-4` |

The checker fix is to keep fixed-grid comparisons on the standard requested
output nodes, but make the boosted reference native transfer grid at least
`100` samples per log interval when `transfer_k_per_logint != 0`.

Artifacts:

- [delta10_kperlog_sweep.png](../../accuracy_plots/lmaxnr_hmcode_isolation/delta10_kperlog_sweep/delta10_kperlog_sweep.png)
- [grid_phase_summary.csv](../../accuracy_plots/lmaxnr_hmcode_isolation/delta10_kperlog_sweep/grid_phase_summary.csv)
- [check_accuracy_q010_delta10_dense_ref/matter_power_errors.png](../../accuracy_plots/lmaxnr_hmcode_isolation/check_accuracy_q010_delta10_dense_ref/matter_power_errors.png)

## Full referenced-INI validation

Generated INIs:

```bash
python fortran/tests/CAMB_test_files.py accuracy_plots/current_check_accuracy/inis --make_ini --no_run_test
```

Checker command:

```bash
python -c 'import sys; from camb._command_line import _run_check_accuracy; raise SystemExit(_run_check_accuracy(sys.argv[1:]))' INI
```

The earlier sweep also classified lensed-CMB failures within 100 multipoles of
the common lensed endpoint as ignored for now. In this rerun only
`params_nolens.ini` has ignored endpoint rows, and it has no remaining
actionable rows after endpoint ignoring. The previous tensor-only endpoint/TE
failure is resolved by capping the tensor-only lensed comparison at
`max_l_tensor`.

Artifacts:

- [summary CSV](../../accuracy_plots/current_check_accuracy/referenced_ini_summary_ignore_endpoint.csv)
- [failed-row CSV](../../accuracy_plots/current_check_accuracy/referenced_ini_failed_rows.csv)
- [per-ini plots](../../accuracy_plots/current_check_accuracy/plots/)

Earlier full referenced results after ignoring those lensed endpoint rows:

| Ini | Status ignoring endpoint | Actionable worst row | Matter-power worst row |
| --- | --- | --- | --- |
| `params_nolens.ini` | PASS | none; raw TT/EE endpoint rows near `L = 2176` ignored | none |
| `params_lmax6000.ini` | PASS | none | none |
| `params_all_nonlinhigh.ini` | PASS | none | OK, mid-`k` `9.859e-4` |
| `params_all_nonlin1.ini` | PASS | none | OK, broad row `1.527e-3` vs relaxed `3e-3` |
| `params_delta_xe.ini` | PASS | none | OK, mid-`k` `9.860e-4` |
| `params_all_nonlin2.ini` | FAIL | nonlinear high-`k` matter power, `5.393e-3` at `k/h = 7.14286` | FAIL |
| `params_HMcode.ini` | FAIL | nonlinear high-`k` matter power, `3.722e-3` at `k/h = 75.8158` | FAIL |
| `params_halomodel.ini` | FAIL | nonlinear high-`k` matter power, `3.777e-3` at `k/h = 100` | FAIL |
| `params_tensor_tranfer.ini` | PASS | none | OK, mid-`k` `9.394e-4` |
| `params_tranfer_highprec.ini` | PASS | none | OK, mid-`k` `9.394e-4` |
| `params_tranfer_delta10.ini` | PASS with dense boosted fixed-grid reference | none | OK, mid-`k` `8.257e-4` |
| `params_tranfer_only.ini` | FAIL | transfer endpoint matter power, `3.335e-3` at `k/h = 2` | FAIL |
| `params_tranfer_redshifts.ini` | FAIL | transfer endpoint matter power, `3.335e-3` at `z = 1`, `k/h = 2` | FAIL |
| `params_tranfer_redshifts2.ini` | FAIL | transfer endpoint matter power, `3.335e-3` at `z = 0.7`, `k/h = 2` | FAIL |
| `params_tranfer_nonu.ini` | FAIL | transfer endpoint matter power, `3.335e-3` at `k/h = 2` | FAIL |

Those difference plots for remaining failures are in the same artifact
tree. The most useful entry points are:

- [params_all_nonlin2/matter_power_errors.png](../../accuracy_plots/current_check_accuracy/plots/params_all_nonlin2/matter_power_errors.png)
- [params_HMcode/matter_power_errors.png](../../accuracy_plots/current_check_accuracy/plots/params_HMcode/matter_power_errors.png)
- [params_halomodel/matter_power_errors.png](../../accuracy_plots/current_check_accuracy/plots/params_halomodel/matter_power_errors.png)
- [params_tranfer_only/matter_power_errors.png](../../accuracy_plots/current_check_accuracy/plots/params_tranfer_only/matter_power_errors.png)

The Planck 2018 high-`L` lensed accuracy check passes with the corrected
peak-scale definition:

```bash
python -c 'import sys; from camb._command_line import _run_check_accuracy; raise SystemExit(_run_check_accuracy(sys.argv[1:]))' \
  inifiles/planck_2018.ini --set-for-lmax 4000 --lens-potential-accuracy 4
```

Current key rows:

- lensed EE `0 <= L < 600`: `2.294e-3` at `L = 20`, passing `3e-3`;
- lensed EE `600 <= L < 3500`: `4.340e-4` at `L = 1081`, passing `1e-3`;
- matter power `0.005 <= k/h < 2`: `8.815e-4`, passing `1e-3`;
- result: `PASS`.

## High-precision transfer WantTransfer sweep

After the final `q = 0.010` `lmaxnr` ramp, all generated inis under
`accuracy_plots/current_check_accuracy/inis` with both
`transfer_high_precision = T` and `get_transfer = T` were checked with the
default boosted `check_accuracy` reference.

Summary:

- checked: `38` generated inis;
- pass: `24`;
- fail: `13`;
- error: `1`, from `params_omk_-0.030.ini` hitting a matter-power interpolation
  range mismatch;
- over-tolerance matter-power rows: none.

The remaining failures in this subset are CMB/lensing rows, not matter-power
hierarchy failures:

- lensing potential: `params_21cm_base2.ini`, `params_omk_-0.001.ini`;
- high-`L` EE/TE: `params_hubble_62.000.ini`, `params_hubble_67.000.ini`,
  `params_hubble_78.000.ini`, `params_ombh2_0.025.ini`,
  `params_omch2_0.080.ini`, `params_omch2_0.100.ini`,
  `params_omch2_0.150.ini`, `params_omk_0.040.ini`,
  `params_w_-0.750.ini`, `params_w_-1.200.ini`;
- low-`L` EE: `params_re_optical_depth_0.110.ini`;
- interpolation-range error: `params_omk_-0.030.ini`.

Artifact:

- [summary_after_dense_ref.csv](../../accuracy_plots/lmaxnr_hmcode_isolation/high_precision_transfer_sweep/summary_after_dense_ref.csv)

## Massive-neutrino transfer follow-up

`params_mu_masssplit.ini` is CMB-only in the generated standard test set, so its
matter-power stability was isolated with a transfer-only variant using
`transfer_high_precision = T`, `get_transfer = T`, `do_lensing = F`, and
`transfer_kmax = 5`.

The failure was not fixed by transfer k-grid density, time stepping, or
integration tolerance. It is also not intrinsically a mass-split issue:
three equal-mass eigenstates show the same instability over a wide `omnuh2`
range. The old hard behavior around the large-mass cuts was not physical: the
three-degenerate sweep failed below and above `m_nu c^2/k_B T_nu0 = 2500`, and
also changed behavior around the pre-existing `5000` "megadamping" multiplier.
The controlling pieces are:

- the massive-neutrino momentum quadrature once an eigenstate is heavy enough
  to be nonrelativistic through the transfer range: the 3-point optimized rule
  leaves high-precision transfer spikes, the 4-point rule removes the first
  set, and a 5-point rule is needed for the heaviest transfer cases to agree
  with a 5-point reference;
- the active massive-neutrino hierarchy depth and associated
  nonrelativistic/massive-neutrino approximation timing, implemented as local
  multiplicative boosts so non-default accuracy settings continue to scale
  smoothly.

The accepted code path is mass- and scale-gated rather than split-gated:

- the 4-point q rule is used only for high-precision transfer runs with
  `max(m_nu c^2 / k_B T_nu0) > 600`;
- the 5-point q rule is used only for high-precision transfer runs with
  `max(m_nu c^2 / k_B T_nu0) > 4300`;
- `TNuPerturbations_init` now uses an exclusive accuracy-threshold ladder, so
  the selected q rule is not assigned and overwritten by several independent
  `if` statements;
- the hierarchy/timing multiplier is applied only to individual eigenstates above
  `m_nu c^2/k_B T_nu0 > 600`; the timing multiplier is `2.0` and the
  hierarchy multiplier is `2.15`, rather than a hard `max(..., 2)` floor;
- the hierarchy multiplier is additionally limited to
  `q * tau_nonrelativistic_physical < 500`, using an unboosted physical
  nonrelativistic time. The delayed switch time is kept separate for the actual
  approximation handoff, so the gate targets modes whose evolution is still
  controlled by the nonrelativistic free-streaming turnover rather than
  boosting all high-`k` massive-neutrino modes;
- all of these targeted transfer changes are gated by `AccuracyTarget > 0`;
- `EV%lmaxnu` is raised consistently with the per-eigenstate `lmaxnu_tau`
  target, and the existing final `lmaxnr >= lmaxnu` rule keeps the massless
  hierarchy deep enough for the relativistic-correction subtraction;
- global `GetNumEqns` `l_accuracy_boost` remains unchanged, so photon,
  massless-neutrino, and high-k massive-neutrino hierarchies are not broadly
  raised.

Saved test INIs:

The focused reproducer INIs used for this isolation are saved in
[matter_power_transfer_accuracy_inis](matter_power_transfer_accuracy_inis/).
They are transfer-only high-precision matter-power cases with
`transfer_k_per_logint = 0` and `do_lensing = F`. The three-degenerate files
use INI `transfer_kmax = 5` (`k/h = 5`). The `params_mu_*` variants reproduce
the Python `set_matter_power(kmax=5)` checks and therefore serialize as
`transfer_kmax = 7.142857...` for `h = 0.7`.

| File | Purpose |
| --- | --- |
| [params_mu_masssplit_transfer_highprec_k5.ini](matter_power_transfer_accuracy_inis/params_mu_masssplit_transfer_highprec_k5.ini) | Original split-mass MPK reproducer derived from `fortran/testfiles/params_mu_masssplit.ini` |
| [params_mu_masssplit_transfer_highprec_k5_accmassnu.ini](matter_power_transfer_accuracy_inis/params_mu_masssplit_transfer_highprec_k5_accmassnu.ini) | Same split-mass reproducer with `accurate_massive_neutrino_transfers = T` |
| [params_mu_mass0p1_transfer_highprec_k5.ini](matter_power_transfer_accuracy_inis/params_mu_mass0p1_transfer_highprec_k5.ini) | Light single-eigenstate control derived from `params_mu_mass0.1.ini` |
| [three_deg_omnuh2_0p00320.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p00320.ini) | Three equal eigenstates just below the `m/Tnu0 = 600` gate |
| [three_deg_omnuh2_0p00325.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p00325.ini) | Three equal eigenstates just above the `m/Tnu0 = 600` gate |
| [three_deg_omnuh2_0p01340.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p01340.ini) | Three equal eigenstates near the old `m/Tnu0 = 2500` behavior change |
| [three_deg_omnuh2_0p02000.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p02000.ini) | Representative 4-point q case |
| [three_deg_omnuh2_0p02300.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p02300.ini) | Three equal eigenstates at the 5-point q onset |
| [three_deg_omnuh2_0p02400.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p02400.ini) | Representative 5-point q threshold case |
| [three_deg_omnuh2_0p02500.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p02500.ini) | Heavy three-eigenstate q5/q100 comparison case |
| [three_deg_omnuh2_0p03000.ini](matter_power_transfer_accuracy_inis/three_deg_omnuh2_0p03000.ini) | Heavier three-eigenstate q5/q100 comparison case |

Validation:

| Ini | Status | Mid-`k` matter-power max | Standard cpu |
| --- | --- | ---: | ---: |
| `params_mu_masssplit_transfer_highprec_k5.ini` before | FAIL | `4.721e-3` | `0.718s` |
| `params_mu_masssplit_transfer_highprec_k5.ini` after | PASS | `7.295e-4` | `0.993s` |
| `params_mu_masssplit_transfer_highprec_k5_accmassnu.ini` after | PASS | `7.295e-4` | `1.24s` |
| `params_mu_mass0.1_transfer_highprec_k5.ini` after | PASS | `9.587e-4` | unchanged within timing noise |
| `three equal eigenstates, omnuh2=0.00320`, `m/Tnu0 = 598.5` | PASS | `9.04e-4` | no local boost |
| `three equal eigenstates, omnuh2=0.00325`, `m/Tnu0 = 607.9` | PASS | `9.57e-4` | hierarchy/switch boost active |
| `three equal eigenstates, omnuh2=0.01340`, `m/Tnu0 = 2506.3` | PASS | `9.44e-4` | 4-point q active |
| `three equal eigenstates, omnuh2=0.02000`, `m/Tnu0 = 3740.7` | PASS | `7.86e-4` | 4-point q active |
| `three equal eigenstates, omnuh2=0.02400`, `m/Tnu0 = 4488.9` | PASS | `8.76e-4` | 5-point q active |
| `three equal eigenstates, omnuh2=0.03000`, `m/Tnu0 = 5611.1` | PASS | `8.88e-4` | 5-point q active |

Threshold checks:

- `m/Tnu0 = 561.1` and `598.5` pass without the local massive-neutrino boost.
- `m/Tnu0 = 673.3` fails at `1.018e-3` if the mass gate is moved up to
  `700`, so the final `600` gate is the smallest tested clean onset.
- Reducing the hierarchy window to `q * tau_nonrelativistic_physical <= 400` or `450`
  leaves the `omnuh2 = 0.02` three-eigenstate case failing at `1.16e-3` near
  `k/h = 1.88`; the final `500` window is the smallest tested clean
  `q * tau_nonrelativistic_physical` cut.
- With multiplicative scaling, hierarchy multipliers `1.5`, `1.75`, `1.9`,
  `2.0`, and `2.1` were not clean across the sweep; `2.15` is the smallest
  tested value that passed the `omnuh2 = 0.02` and `0.025` tight cases.
- Comparing the 4-point q rule to an explicit 5-point q reference is clean up
  to `m/Tnu0 = 4301.9`, but `m/Tnu0 = 4488.9` fails at `1.04e-3`; this sets
  the final 5-point q onset at `4300`.

Slow full-q reference check:

The slow branch of `TNuPerturbations_init` was checked by using the standard
reference boosts plus `neutrino_q_boost = 5`, so the reference call has
`AccuracyBoost*neutrino_q_boost = 10` and evolves `100` q nodes. Results:

| Ini | Standard q case | Status vs q100 reference | Mid-`k` matter-power max | Note |
| --- | --- | --- | ---: | --- |
| `three equal eigenstates, omnuh2=0.00320`, `m/Tnu0 = 598.5` | 3q | PASS | `9.126e-4` | below the local q boost gate |
| `three equal eigenstates, omnuh2=0.02000`, `m/Tnu0 = 3740.7` | 4q | PASS | `7.831e-4` | 4q transfer regime |
| `params_mu_masssplit_transfer_highprec_k5.ini`, `m/Tnu0,max = 13263.3` | 5q | PASS | `7.099e-4` | original split-mass case |
| `three equal eigenstates, omnuh2=0.02300`, `m/Tnu0 = 4301.9` | 5q | FAIL | `1.010e-3` | just above the q5 transition |
| `three equal eigenstates, omnuh2=0.02400`, `m/Tnu0 = 4488.9` | 5q | FAIL | `1.042e-3` | worst high-`k` row `1.345e-3`, still below high-`k` tolerance |
| `three equal eigenstates, omnuh2=0.02500`, `m/Tnu0 = 4675.9` | 5q | FAIL | `1.070e-3` | heavy equal-mass transition regime |
| `three equal eigenstates, omnuh2=0.03000`, `m/Tnu0 = 5611.1` | 5q | FAIL | `1.180e-3` | heavy equal-mass transition regime |

The q100 reference itself is stable: for the `omnuh2 = 0.02400` case,
q100 compared to q200 gives mid-`k` `2.38e-6` and worst-row `2.64e-6`.
Forcing the standard side onto the slow uniform branch with `Accuracy = 3.001`
(`30` q nodes) makes the heavy equal-mass cases pass the q100 reference
(`9.27e-4` to `9.38e-4` mid-`k`), but raises representative standard CPU time
from about `1.3s` to `3.7--4.9s`. That fallback is therefore not included in
the current minimal patch; it is a possible follow-up if agreement with the
uniform slow-q reference is required for high-precision transfer cases with
several eigenstates in the heavy transition band.

Holding the rest of the boosted reference fixed and changing only the q table
shows that q5 itself is much closer to q100 than the earlier standard-vs-q100
comparison. A temporary diagnostic build disabled only the automatic q-count
promotion, then used `AccuracyBoost = lAccuracyBoost = lSampleBoost =
IntTolBoost = 2` throughout. With `neutrino_q_boost = 0.5`, `1.0`, `1.001`,
and `5.0`, the perturbation initializer used the 3-point, 4-point, 5-point,
and q100 rules respectively. Mid-`k` matter-power changes relative to q100 were:

| Ini | q3 | q4 | q5 |
| --- | ---: | ---: | ---: |
| `three equal eigenstates, omnuh2=0.00320`, `m/Tnu0 = 598.5` | `1.62e-4` | `3.45e-5` | `2.03e-5` |
| `three equal eigenstates, omnuh2=0.02000`, `m/Tnu0 = 3740.7` | `1.90e-3` | `4.35e-4` | `2.25e-4` |
| `three equal eigenstates, omnuh2=0.02400`, `m/Tnu0 = 4488.9` | `2.50e-3` | `5.80e-4` | `2.98e-4` |
| `three equal eigenstates, omnuh2=0.03000`, `m/Tnu0 = 5611.1` | `3.51e-3` | `8.32e-4` | `4.17e-4` |
| `params_mu_masssplit_transfer_highprec_k5.ini`, `m/Tnu0,max = 13263.3` | `4.52e-3` | `1.02e-3` | `5.47e-4` |

So the q100-vs-standard residual in the heavy equal-mass cases is not mainly a
5-point quadrature-table error when all other reference boosts are matched. The
q-table-only error is already below the mid-`k` tolerance for q5 in these
representative cases, while q3 is clearly insufficient above the nonrelativistic
transfer threshold and q4 is marginal for the split-mass case.

After separating the physical nonrelativistic time from the delayed switch time,
the representative high-precision transfer checks still pass:
`omnuh2 = 0.00325`, `0.02000`, `0.02400`, and the split-mass `kmax = 5` case
give mid-`k` matter-power maxima of `9.54e-4`, `8.88e-4`, `8.76e-4`, and
`7.30e-4` respectively.

## Massive-neutrino CMB/lensing cross-check

The same neutrino settings were also checked for lensed CMB spectra and the
lensing potential, with transfer output disabled and high-precision transfer
off. The generated CMB/lensing INIs are saved in
[matter_power_transfer_accuracy_inis/cmb_lensing_lmax4000_lpa6](matter_power_transfer_accuracy_inis/cmb_lensing_lmax4000_lpa6/).
They use `get_transfer = F`, `do_lensing = T`, `l_max_scalar = 4150`, and
`k_eta_max_scalar = 108000`, matching
`set_for_lmax(4000, lens_potential_accuracy=6)` with the standard lens-output
margin. The checker comparison used `lmax = 4000`, `set_for_lmax = 4000`,
`lens_potential_accuracy = 6`, and the default boosted reference accuracy
settings.

At this stage the massive-neutrino transfer changes were still gated by
`AccuracyTarget > 0`, `WantTransfer`, and `Transfer.high_precision`, so these
CMB/lensing runs did not activate the new q-rule promotion, the
massive-neutrino switch delay, or the local hierarchy boost. The failures below
therefore identified a separate CMB accuracy sensitivity in heavy-neutrino
models, rather than showing that the full MPK-specific transfer fix should be
used unchanged for CMB spectra.

| Ini | Status | Worst lensed failure | Worst lensing-potential row |
| --- | --- | --- | --- |
| `params_mu_mass0p1_transfer_highprec_k5` | PASS | none; low-`L` EE `2.41e-3` vs `3e-3` | high-`L` PP `2.39e-3` vs `5e-3` |
| `params_mu_masssplit_transfer_highprec_k5` | FAIL | EE `2.81e-3` vs `1e-3` at `L = 721` | PP `5.53e-3` vs `5e-3` at `L = 1899` |
| `params_mu_masssplit_transfer_highprec_k5_accmassnu` | FAIL | EE `2.81e-3` vs `1e-3` at `L = 721` | PP `5.53e-3` vs `5e-3` at `L = 1899` |
| `three_deg_omnuh2_0p00320` | PASS | none | high-`L` PP `1.85e-3` vs `5e-3` |
| `three_deg_omnuh2_0p00325` | PASS | none | high-`L` PP `1.85e-3` vs `5e-3` |
| `three_deg_omnuh2_0p01340` | PASS | none | mid-`L` PP `2.75e-3` vs `5e-3` |
| `three_deg_omnuh2_0p02000` | FAIL | EE `4.67e-3` vs `3e-3` at `L = 3985` | mid-`L` PP `4.06e-3` vs `5e-3` |
| `three_deg_omnuh2_0p02300` | FAIL | EE `4.75e-3` vs `3e-3` at `L = 3985` | mid-`L` PP `4.64e-3` vs `5e-3` |
| `three_deg_omnuh2_0p02400` | FAIL | EE `4.68e-3` vs `3e-3` at `L = 3985` | mid-`L` PP `4.80e-3` vs `5e-3` |
| `three_deg_omnuh2_0p02500` | FAIL | EE `1.53e-3` vs `1e-3` at `L = 837` | PP `5.01e-3` vs `5e-3` at `L = 845` |
| `three_deg_omnuh2_0p03000` | FAIL | EE `2.99e-3` vs `1e-3` at `L = 600` | PP `5.96e-3` vs `5e-3` at `L = 1016` |

The light and low-mass three-degenerate cases pass cleanly. The failures start
for the heavier three-degenerate models and for the original split-mass model,
mainly in lensed EE around the high-`L` endpoint or around the lensing-kernel
transition, with lensing-potential rows marginally failing only in the split and
heaviest equal-mass cases. The `accurate_massive_neutrino_transfers` variant is
identical to the standard split case here because transfer output is disabled.

A temporary ungated test removed the `WantTransfer` and
`Transfer.high_precision` gates from the new massive-neutrino transfer q-rule
promotion, switch delay, and local hierarchy boost, while keeping the mass and
`q*tau_nonrelativistic_physical` cuts. This is not kept in the code. It shows
that the same changes do improve the lensing-potential rows and remove the
dominant `L ~ 1000` PP/EE failures, but they do not make the CMB/lensing checks
pass: the remaining worst rows move to high-`L` EE near the endpoint.

| Ini | Gated worst failure | Ungated worst failure | Ungated PP row |
| --- | --- | --- | --- |
| `params_mu_masssplit_transfer_highprec_k5` | EE `2.81e-3` vs `1e-3` at `L = 721`; PP `5.53e-3` vs `5e-3` at `L = 1899` | EE `4.38e-3` vs `3e-3` at `L = 3984` | passes, worst PP `1.02e-3` |
| `three_deg_omnuh2_0p02000` | EE `4.67e-3` vs `3e-3` at `L = 3985`; PP `4.06e-3` vs `5e-3` | EE `4.24e-3` vs `3e-3` at `L = 3985` | passes, worst PP `9.98e-4` |
| `three_deg_omnuh2_0p02500` | EE `1.53e-3` vs `1e-3` at `L = 837`; PP `5.01e-3` vs `5e-3` at `L = 845` | EE `4.47e-3` vs `3e-3` at `L = 3985` | passes, worst PP `1.21e-3` |
| `three_deg_omnuh2_0p03000` | EE `2.99e-3` vs `1e-3` at `L = 600`; PP `5.96e-3` vs `5e-3` at `L = 1016` | EE `3.91e-3` vs `3e-3` at `L = 4000` | passes, worst PP `1.42e-3` |

Warmed standard-run timings, with one warm-up calculation discarded and then two
runs averaged, show the ungated version is not cheap for CMB/lensing:

| Ini | Gated CPU | Ungated CPU | Change |
| --- | ---: | ---: | ---: |
| `params_mu_masssplit_transfer_highprec_k5` | `4.13s` | `5.15s` | `+24.6%` |
| `three_deg_omnuh2_0p02000` | `5.00s` | `5.78s` | `+15.6%` |
| `three_deg_omnuh2_0p02500` | `4.59s` | `5.36s` | `+16.6%` |
| `three_deg_omnuh2_0p03000` | `5.33s` | `6.61s` | `+24.2%` |

So the transfer-specific hierarchy boost should remain gated to high-precision
transfer. The PP improvement needs a separate CMB-targeted setting that avoids
paying the full transfer hierarchy cost and leaves the independent high-`L` EE
endpoint residual to a separate fix.

The follow-up tuned CMB setting keeps the high-precision transfer hierarchy
boost at `2.15`, but uses a smaller `1.5` local hierarchy multiplier for CMB
spectra. The q-rule promotion and nonrelativistic switch delay are enabled for
both high-precision transfer and CMB massive-neutrino runs. This fixes the
original lensed middle-`L` EE/TT/TE and lensing-potential rows while still
leaving the high-`L` EE endpoint residual as a separate issue.

| Ini | Middle EE | Middle TT | Middle PP | Remaining failure |
| --- | ---: | ---: | ---: | --- |
| `params_mu_masssplit_transfer_highprec_k5` | `8.08e-4` | `6.32e-4` | `1.45e-3` | high-`L` EE `4.55e-3` vs `3e-3` |
| `three_deg_omnuh2_0p02000` | `9.53e-4` | `7.31e-4` | `1.33e-3` | high-`L` EE `4.41e-3` vs `3e-3` |
| `three_deg_omnuh2_0p02500` | `9.53e-4` | `6.99e-4` | `1.62e-3` | high-`L` EE `4.45e-3` vs `3e-3` |
| `three_deg_omnuh2_0p03000` | `8.61e-4` | `5.97e-4` | `1.90e-3` | high-`L` EE `3.90e-3` vs `3e-3` |

The light control remains fully passing. A q-only change did not pass the
middle EE rows, a switch-only change left PP failing, and hierarchy-only fixed
PP but not middle EE. The combination above is the smallest tested clean CMB
setting for the original lensed middle-`L` failures. Raising the CMB hierarchy
multiplier above this is not monotonic: `2.5` worsened the middle EE rows, so
the smaller `1.5` multiplier is preferred.

Unlensed CMB spectra were checked with the same neutrino settings, `do_lensing
= F`, `Want_CMB_lensing = F`, and `lmax = 4000`. With the CMB boost still
restricted to lensed/lensing-potential runs, the split and heavy equal-mass
unlensed cases have large middle-`L` residuals, e.g. split-mass EE
`2.96e-3`, `omnuh2 = 0.03000` EE `3.00e-3`, and `omnuh2 = 0.02000` TT
`1.47e-3`, all against the `1e-3` tolerance. Enabling the same tuned q-rule,
switch delay, and smaller CMB hierarchy boost for all `WantCls` runs removes
that large heavy-neutrino component:

| Ini | Unlensed middle EE before | Unlensed middle EE after | Note |
| --- | ---: | ---: | --- |
| `params_mu_masssplit_transfer_highprec_k5` | `2.96e-3` | `1.05e-3` | remaining row is just above tolerance |
| `three_deg_omnuh2_0p02000` | `1.47e-3` TT | `1.28e-3` EE | TT now passes; EE near generic residual |
| `three_deg_omnuh2_0p02500` | `1.54e-3` | `1.20e-3` | low-`L` EE also returns below tolerance |
| `three_deg_omnuh2_0p03000` | `3.00e-3` | `1.15e-3` | large heavy-neutrino residual removed |

The unlensed light control, whose mass is below the local massive-neutrino
gates, still has EE `1.12e-3` in the same `600 <= L < 3000` band. That means
the remaining small unlensed EE excess is not caused by the massive-neutrino
boost gates and should be treated as a separate baseline unlensed-CMB accuracy
issue. The massive-neutrino-specific part is therefore not only a lensed-`Cl`
problem, so the final code gates the smaller CMB multiplier on `WantCls`, not on
`DoLensing` or `Want_CMB_lensing`.

Warmed standard-run timings compared to the original transfer-gated code:

| Ini | Original gated CPU | Tuned CPU | Change |
| --- | ---: | ---: | ---: |
| `params_mu_masssplit_transfer_highprec_k5` | `4.13s` | `4.81s` | `+16.4%` |
| `three_deg_omnuh2_0p02000` | `5.00s` | `4.99s` | `-0.1%` |
| `three_deg_omnuh2_0p02500` | `4.59s` | `5.51s` | `+19.9%` |
| `three_deg_omnuh2_0p03000` | `5.33s` | `6.43s` | `+20.7%` |

The high-precision transfer MPK path was rechecked after the split transfer/CMB
settings. It is unchanged within the focused diagnostics: the split-mass
transfer reproducer passes with middle-`k` matter power `7.30e-4`, and the
three-degenerate `omnuh2 = 0.02400` case passes with `8.76e-4`.

Artifacts:

- [massive_neutrino_before_after_mpk_diff.png](../../accuracy_plots/masssplit_isolation/before_after_final/massive_neutrino_before_after_mpk_diff.png)
- [three_degenerate_omnu_sweep_final.png](../../accuracy_plots/masssplit_isolation/three_degenerate_omnu_sweep_final/three_degenerate_omnu_sweep_final.png)
- [5-point q reference masssplit plot](../../accuracy_plots/masssplit_isolation/final_multiplicative_q5_checks/q5/masssplit/matter_power_errors.png)
- [5-point q reference three-degenerate plot](../../accuracy_plots/masssplit_isolation/final_multiplicative_q5_checks/q5/three_deg_0p02400/matter_power_errors.png)
- [q100 reference summary CSV](../../accuracy_plots/masssplit_isolation/q10_reference_checks/summary.csv)
- [q100 reference combined plot](../../accuracy_plots/masssplit_isolation/q10_reference_checks/q10_reference_mpk_errors.png)
- [q100 vs q200 stability plot](../../accuracy_plots/masssplit_isolation/q10_reference_checks/q10_vs_q20/matter_power_errors.png)
- [q30 standard-side probe summary CSV](../../accuracy_plots/masssplit_isolation/q10_reference_checks/q30_standard_probe/summary_extra.csv)
- [q table vs q100 summary CSV](../../accuracy_plots/masssplit_isolation/q_table_vs_q100_reference/summary.csv)
- [q table vs q100 plot](../../accuracy_plots/masssplit_isolation/q_table_vs_q100_reference/q_table_vs_q100_mpk_errors.png)
- [physical tau gate recheck summary CSV](../../accuracy_plots/masssplit_isolation/physical_tau_gate_recheck/summary.csv)
- [CMB/lensing lmax4000 lens-potential-accuracy-6 summary CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/summary.csv)
- [split-mass lensed Cl errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/params_mu_masssplit_transfer_highprec_k5/lensed_cls_errors.png)
- [split-mass lensing-potential errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/params_mu_masssplit_transfer_highprec_k5/lens_potential_errors.png)
- [heavy three-degenerate lensed Cl errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/three_deg_omnuh2_0p03000/lensed_cls_errors.png)
- [heavy three-degenerate lensing-potential errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/three_deg_omnuh2_0p03000/lens_potential_errors.png)
- [ungated transfer-boost CMB/lensing summary CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_ungated_transfer_boosts/summary.csv)
- [ungated timing CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_ungated_transfer_boosts/timing_ungated.csv)
- [gated timing CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6/timing_gated.csv)
- [final CMB-tuned summary CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_final_cmb_tuned/summary.csv)
- [final CMB-tuned timing CSV](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_final_cmb_tuned/timing_final.csv)
- [final split-mass lensed Cl errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_final_cmb_tuned/params_mu_masssplit_transfer_highprec_k5/lensed_cls_errors.png)
- [final split-mass lensing-potential errors](../../accuracy_plots/masssplit_isolation/cmb_lensing_lmax4000_lpa6_final_cmb_tuned/params_mu_masssplit_transfer_highprec_k5/lens_potential_errors.png)
- [unlensed CMB before WantCls boost summary CSV](../../accuracy_plots/masssplit_isolation/cmb_unlensed_lmax4000/summary.csv)
- [unlensed CMB after WantCls boost summary CSV](../../accuracy_plots/masssplit_isolation/cmb_unlensed_lmax4000_wantcls_boost/summary.csv)
- [final transfer recheck summary CSV](../../accuracy_plots/masssplit_isolation/final_cmb_tuned_transfer_recheck/summary.csv)
- [before/after summary CSV](../../accuracy_plots/masssplit_isolation/before_after_final/summary.csv)
- [three-degenerate sweep CSV](../../accuracy_plots/masssplit_isolation/three_degenerate_omnu_sweep_final/summary.csv)
- [default reference summary CSV](../../accuracy_plots/masssplit_isolation/final_multiplicative_q5_checks/default_summary.csv)
- [5-point q reference summary CSV](../../accuracy_plots/masssplit_isolation/final_multiplicative_q5_checks/q5_summary.csv)
- [timing summary CSV](../../accuracy_plots/masssplit_isolation/timing_multiplicative_q5/timing_summary.csv)

## Recommended next fixes

1. Treat the high-precision transfer matter-power hierarchy fix as done: keep
   the photon hierarchy floors and the final monotone `q = 0.010` `lmaxnr`
   ramp.
2. Keep the CMB massive-neutrino adjustment separate from the
   high-precision transfer hierarchy fix: the transfer path uses the stronger
   `2.15` local hierarchy multiplier, while CMB spectra use the tuned `1.5`
   multiplier plus the same mass/q gates. The remaining high-`L` EE endpoint
   residual and the small unlensed EE residual should be treated as separate
   CMB accuracy issues.
3. For fixed `transfer_k_per_logint` cases, keep the fixed-output-node contract
   but use a denser boosted reference grid. The standard
   `params_tranfer_delta10.ini` output with `10` samples per log interval passes
   once the boosted reference is evaluated with at least `100` samples per log
   interval.
4. Keep the new `get_matter_power_spectrum` default routed through the public
   Python templated interpolator, with `use_spline_template=False` as the
   legacy Fortran interpolation escape hatch.
5. Keep high-`k` nonlinear matter-power residuals separate from the spline
   template work. They survive the interpolation fix and are likely tied to
   nonlinear ratio/tail behavior rather than BAO interpolation.

## Endpoint follow-up

`params_tranfer_only.ini` still shows a matter-power endpoint/tail residual:

```text
matter P(k), all z, 0.0001 <= k/h <= 2: 3.335e-3 at k/h = 2
```

Clipping the top of the comparison grid shows that the error rises smoothly into
the upper boundary:

| Top samples clipped | Worst remaining error | Worst `k/h` |
| ---: | ---: | ---: |
| 0 | `3.335e-3` | `2.000` |
| 1 | `3.303e-3` | `1.961` |
| 5 | `3.173e-3` | `1.811` |
| 10 | `3.009e-3` | `1.640` |
| 20 | `2.745e-3` | `0.0596` |

Forcing the public interpolator to build guard support beyond the requested
endpoint with `extrap_kmax = 1.05 * kmax` reduces the exact endpoint error from
`3.335e-3` to `2.847e-3`; larger guards keep the endpoint below tolerance but
leave a marginal interior high-k tail around `k/h ~= 1.43`.

A simple one-point residual ghost padding of the templated spline was tested
and rejected: it left the checker max at the endpoint. The cleaner fix is to
give the interpolator real guard support above requested `kmax`, then evaluate
and advertise only the requested range. Concretely, for automatic-grid public
Python interpolation:

- internally request or construct matter-power support to at least
  `kmax * 1.05` to `kmax * 1.2`;
- keep `PKInterpolator.kmax` and returned public `k` nodes at the user-requested
  maximum, so API behavior does not silently expand;
- keep fixed `transfer_k_per_logint != 0` checks on requested native nodes,
  because those are an output-grid contract rather than an arbitrary-`k`
  interpolation contract;
- if internal recomputation above `Transfer.kmax` is too invasive, use the
  existing power-law extrapolation support as a temporary guard, but treat it as
  a boundary-condition fix, not as physical high-k accuracy.

The endpoint plot for the current tree is:

- [params_tranfer_only_current/matter_power_errors.png](../../accuracy_plots/mpk_end_effects/params_tranfer_only_current/matter_power_errors.png)
