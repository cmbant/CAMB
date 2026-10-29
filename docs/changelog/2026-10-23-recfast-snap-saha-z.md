# RECFAST Recombination_saha_z recorded at an output node

- Date: 2026-10-23
- Scope: optional RECFAST setting `snap_saha_z` (ini key `RECFAST_snap_saha_z`), off by default. `lag_corrections`
  turns it on (see `2026-10-23-recfast-lag-corrections.md`).

## Symptom

`Recombination_saha_z` is the redshift below which CAMB evolves the delta x_e perturbation when `Evolve_delta_xe` is
on. With the setting off, RECFAST records it in two places, at the output nodes of the rows the Rosenbrock solver
integrates, and in the right-hand side `ION` wherever the integrator steps inside 0.985 < x_H < 0.99 (x_H is
hydrogen's ionized fraction). The second write wins whenever it happens.

At CAMB's defaults (`rosenbrock_handoff_xH = 0.976`) the right-hand side does not record it, and the value is the node
at z ≈ 1554. With `rosenbrock_handoff_xH = 0.996` or `use_rosenbrock = F`, the value lands between nodes, at z ≈ 1555
in the cosmology of the lag_corrections entry, and it moves with the solver's tolerance. Changing only `rk45_tol`,
from 3e-6 to 1e-8, moves it by up to about 0.01 in z, and not monotonically. `Recombination_saha_z` is not exposed to
Python, so these readings come from a debug print at the end of `TRecfast_init`.

## Change

With the setting on, RECFAST records `Recombination_saha_z` at the last output node with x_H > 0.985, whichever solver
produced that node.

## Effect

- At CAMB's defaults both recorders give the node at z ≈ 1554, so the spectra do not change, even with
  `Evolve_delta_xe` on.
- With `rosenbrock_handoff_xH = 0.996` or `use_rosenbrock = F` the setting records the node at z ≈ 1559. With
  `Evolve_delta_xe` on, lensed TT changes by about 5e-10 and 3e-10 respectively for 2 <= l <= 2500.
- The recorded value no longer depends on `rk45_tol`.

## Validation

- With the new settings at their defaults, outputs match CAMB 2.0.5 bit for bit in the same 392,476 numbers and 541
  files as in `2026-10-23-recfast-heiii-saha-continued.md`.
- `testRecfastSnapSahaZ` checks that the spectra do not change at CAMB's defaults, that they move by less than 1e-8
  with `use_rosenbrock = F`, and that the setting survives an ini round trip.
- `python scripts/recfast_lag_corrections.py --all` derives the node from x_e in its section 9 (1559.14 at
  `rosenbrock_handoff_xH = 0.996`, 1554.25 at the default).
