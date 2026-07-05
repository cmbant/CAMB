# Non-flat source integration: direct Olver evaluation for high-substep ranges

## Change

In `DoRangeInt` (`fortran/cmbmain.f90`), when a single source time step would need
`nSubSteps >= 10` Numerov substeps (`direct_olver_min_substeps`), phi is now filled
by direct `phi_olver` evaluation at each time-step point (new helper
`FillDirectOlverPhiVals`) instead of Numerov-stepping the hyperspherical Bessel ODE
on the finer `dchimax ~ 0.3/nu` chi grid.

The former `ujl`/`ujl_vals` variables (and derived names: `minujl`, `MINUJl1`,
`FillShiftedNuUjlVals`, `FillSmallChiUjlVals`, `HypersphericalUjlOnsetChi`) were
renamed to `phi`/`phi_vals` etc. throughout `cmbmain.f90` for consistency with the
hyperspherical Bessel `phi` naming.

`DoRangeIntTensor` is deliberately left unchanged: instrumented runs (default tensor
settings and pushed `max_l_tensor=1500`, `max_eta_k_tensor=4500`, omk=+/-0.05) show
tensor `nSubSteps` never reaches 10 (averages 1.3-3.2). The tensor pass uses the
denser `dtau0 = max(taurst/40, tau0/2000)` time grid (the tensor source oscillates at
frequency ~k all along the line of sight, with late-time evolution kept only for
k*tau < boost*max_eta_k_tensor/10), and the fast-oscillation range skip bounds q*dtau
in the oscillatory region, so the gate would be dead code there.

## Motivation / measurements

Instrumented profiling (devcontainer, 8 threads, lmax=2500, lens_potential_accuracy=1)
showed the non-flat integration cost is dominated by Numerov substepping, not by the
source multiplies (~2% of runtime): one `u_olver` call costs about the same as ~13
Numerov substeps (~93 ns vs ~7.3 ns/substep), while closed models spend most substeps
in ranges with 8-31 substeps per time step.

Timings (min over repeats, before -> after):

| omk   | before | after  | speedup |
|-------|--------|--------|---------|
| -0.10 | 2.07 s | ~1.5 s | ~38%    |
| -0.05 | 2.27 s | 1.52 s | ~33%    |
| -0.02 | 1.05 s | 0.81 s | ~30%    |
| +0.02 | 0.50 s | 0.50 s | ~0%     |
| +0.05 | 0.65 s | 0.63 s | ~3%     |
| +0.10 | 0.90 s | 0.83 s | ~9%     |

Open models gain less because more of their work uses the near-flat fast paths and
low-substep ranges.

## Accuracy

- Default-accuracy scalar Cl shifts vs previous code: max ~2.3e-5 (TT, omk=-0.1),
  EE <= 4e-6; tensor spectra bit-identical (tensor code path unchanged).
- Errors vs an AccuracyBoost=2 reference are unchanged to displayed precision for
  omk=+/-0.05; consistency also verified running both codes at AccuracyBoost=2
  (diffs <= 6e-7).
- Safety argument: ranges with large `nSubSteps` were already re-seeded from `u_olver`
  every 25-200 Numerov steps by the rebootstrap logic, so Olver absolute accuracy was
  already load-bearing there; direct evaluation removes accumulated phase drift rather
  than adding error.
- Threshold sweep: 8 is slightly faster for closed but moves Cls by ~5e-5; 12 gives
  ~1e-6 shifts; 10 adopted as the agreed speed/accuracy point (few 1e-5 acceptable).
- `python -m camb.check_accuracy` PASSes for params.ini with omk=-0.05 and omk=+0.05.
- `python -m unittest camb.tests.camb_test`: 18 tests OK.

## Notes

- The exploration also showed that coarsening the source time sampling itself (taking
  time steps twice as long for low q) cannot help: the source-multiply pass is ~2% of
  runtime and only ~11% of steps are in `nSubSteps==1` ranges where a longer time step
  would also halve Bessel ODE work.
- Exploration worktree with instrumentation and threshold-sweep scripts:
  `/tmp/camb-nf-speed` (disposable).
