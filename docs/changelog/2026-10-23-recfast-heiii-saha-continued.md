# RECFAST HeIII Saha root continued to z = 3500

- Date: 2026-10-23
- Scope: optional RECFAST setting `heIII_saha_continued` (ini key `RECFAST_heIII_saha_continued`), off by default.
  `lag_corrections` turns it on (see `2026-10-23-recfast-lag-corrections.md`).

## Symptom

Above z_scale = 5000 RECFAST takes x_e from the He++ <-> He+ Saha root, and between z_scale = 5000 and 3500 it sets
x_e = 1 + fHe. Here z_scale = (1 + z) T_CMB/(2.7255 K) - 1, the variable RECFAST's other fixed redshift cuts use. At
the z = 5000 node of the default output grid (Nz = 2046) the root still has about 1.1e-3 of the helium doubly ionized,
so x_e drops by about 9.3e-5 across that cut.

Read through CAMB's output spline, the drop shows up as the largest x_e error over 400 < z < 8000, 9.486e-5 at
z = 4995 against the reference CAMB of the lag_corrections entry. CAMB's largest error below z = 3400 is 5.790e-5.

## Change

With the setting on, the rows down to z_scale = 3500 keep the same Saha root, as `get_init` already does. There the
root lies within about 4e-13 of 1 + fHe, so x_e reaches 1 + fHe without a step.

## Effect

- RECFAST's x_e changes at the rows with 3500 < z_scale <= 5000. Read through CAMB's output spline, the change
  reaches z ≈ 3480 to 5120, while x_e below z = 3400 and theta* stay the same bit for bit.
- The largest x_e error over 400 < z < 8000 falls from 9.486e-5 to 5.790e-5, CAMB's own error below z = 3400.
- The step itself, taken from the Saha root at z = 5000 with CAMB's constants, is 9.2813e-5. That number does not
  depend on how the output is sampled, while the 9.486e-5 above does.
- Run time stays within timing noise, at 0.997 to 1.015 of the time with the setting off in five interleaved runs.

## Validation

- With the new settings at their defaults, outputs match CAMB 2.0.5 bit for bit in 392,476 numbers from background,
  derived-parameter, CMB-spectrum and matter-power runs in 18 configurations, and in the 541 numerical files
  `fortran/tests/run_tests.py` writes.
- `testRecfastHeIIIContinued` checks that the step is gone and that the setting survives an ini round trip.
- `python scripts/recfast_lag_corrections.py --all` prints the step and both x_e readings in its section 5.
