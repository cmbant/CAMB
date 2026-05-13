# Variable Redefinitions

## Summary

Redefined the background-state variable `chi0` as `DMt0`, the comoving angular
diameter distance to the big bang. The new definition is
`curvature_radius * rofChi(tau0 / curvature_radius)`, so it is a physical
distance and reduces to `tau0` in the flat limit.

## Main Conclusions

- The old `chi0` state variable mixed a dimensionless curved-space quantity with
  call sites that expected a physical distance.
- Renaming the state field to `DMt0` makes the intended geometry explicit and
  removes repeated implicit factors of `curvature_radius`.
- The CMB source-grid and integration formulas that depend on the distance to
  last scattering must use `DMt0` directly to stay consistent with the flat
  limit.
- In particular, the `q_cmb` expression should scale as `l / DMt0`; using the
  old dimensionless `chi0` was inconsistent away from flat space.
