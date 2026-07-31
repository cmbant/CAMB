# Non-linear lensing redshift grid

## Current implementation

`CAMBdata_GetComputedPKRedshifts` now uses 30 non-linear-lensing (NLL)
interpolation nodes at the default accuracy,

```text
N = max(2, nint(30 * AccuracyBoost * NonlinSourceBoost))
z_i = z_max * ((N - i) / (N - 1))**(3/2),  i = 1, ..., N.
```

The nodes are ordered from high to low redshift. `z_max` remains 10 normally
and 15 when `AccuracyBoost * NonlinSourceBoost >= 2.5`. The endpoint-inclusive
formula puts the first node exactly at `z_max` and the final node at zero. The
non-linear `kmax` rule is unchanged.

The previous origin/master implementation used 50 uniformly spaced nodes but
excluded the nominal high-redshift endpoint:

```text
N = nint(50 * AccuracyBoost * NonlinSourceBoost)
z_i = (N - i) / (N / z_max),  i = 1, ..., N.
```

For the default `AccuracyBoost = NonlinSourceBoost = 1` and `z_max = 10`, the
sampled redshifts were:

```text
9.8, 9.6, 9.4, 9.2, 9.0, 8.8, 8.6, 8.4, 8.2, 8.0,
7.8, 7.6, 7.4, 7.2, 7.0, 6.8, 6.6, 6.4, 6.2, 6.0,
5.8, 5.6, 5.4, 5.2, 5.0, 4.8, 4.6, 4.4, 4.2, 4.0,
3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0,
1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0
```

The current 30-node `z^(3/2)` grid is:

```text
10.000000, 9.487244, 8.983564, 8.489127, 8.004109, 7.528698,
7.063090, 6.607498, 6.162145, 5.727274, 5.303144, 4.890033,
4.488243, 4.098104, 3.719974, 3.354247, 3.001360, 2.661797,
2.336103, 2.024897, 1.728888, 1.448899, 1.185905, 0.941087,
0.715909, 0.512263, 0.332725, 0.181112, 0.064033, 0.000000
```

Thus the new grid spends fewer nodes at high redshift and samples the rapidly
varying low-redshift non-linear correction more densely, while making the
linear/non-linear interpolation boundary explicit at `z_max`.

## Why an absolute merge tolerance was insufficient

The requested power-spectrum redshifts (`PK_redshifts`) and NLL nodes are
merged into one transfer-redshift array. That same array is used by
`MakeNonlinearSources` to spline the non-linear-to-linear power ratio in
conformal time.

The original merge test used a fixed absolute tolerance of `1e-5` in redshift.
That is too small relative to the interpolation spacing: for the 60-node
`z^(3/2)` grid used by the boosted PPF check, an NLL node occurs at
`z = 0.499294708`, only `7.0529e-4` from the requested `z = 0.5`. The fixed
tolerance retained both nodes, creating a very short spline interval. Natural
cubic-spline second derivatives are global, so this near-duplicate perturbed
the source rescaling and caused a large standard-versus-boosted discrepancy.
The underlying PPF non-linear ratio is smooth; the problem was the knot
layout.

The current code instead uses a local tolerance,

```text
merge_tol = max(1e-5, 0.1 * local_NLL_spacing),
```

and keeps the exact requested PK redshift when a nearby NLL node is redundant.
For example, near `z = 0.5` in the 60-node grid the local merge tolerance is
about `9.65e-3`, so the `0.499294708` NLL node is replaced by the requested
`0.5` node without applying a dangerously large absolute tolerance everywhere.

## Accuracy checks

The current compiled code was checked with:

```text
python -m camb.check_accuracy \
  fortran/testfiles/params_cosmomc_theta_ppf_w-0p82_wa0p35_lmax4000.ini
python -m camb.check_accuracy \
  fortran/testfiles/params_cosmomc_theta_lmax4000_lpa_auto.ini
```

Both report `PASS`. The standard runs use the 30-node grid and the boosted
runs use the corresponding higher-accuracy grid.

| model | lens PP max, `5 <= L < 2500` | tolerance | result |
| --- | ---: | ---: | --- |
| PPF (`w=-0.82`, `wa=0.35`) | `5.881e-4` at `L=2093` | `5e-3` | pass |
| LCDM | `4.126e-4` at `L=2499` | `5e-3` | pass |

The other lensed-CMB and lensing-potential ranges also pass. In particular,
the PPF high-`L` lens-PP range has maximum `6.927e-3` against its `1e-2`
tolerance. The previous fixed-tolerance implementation reached `2.843e-2`
for the PPF low-`L` lens-PP comparison; the temporary global `0.01` test
reduced that to `5.337e-3` but still failed, whereas the local-spacing merge
passes with comfortable margin.

The production merge floor remains `1e-5`; no global `0.01` tolerance is used.
