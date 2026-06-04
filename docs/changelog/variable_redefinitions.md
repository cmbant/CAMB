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

## `scale`

`State%scale` is not the cosmological scale factor. It stores the CMB
peak-position and multipole-sampling ratio

\[
s_{\rm peak}=\frac{\theta_{\rm drag}^{\rm P18}}{\theta_{\rm drag}},\qquad
\theta_{\rm drag}=\frac{r_s(z_{\rm drag})}{D_M(\tau_0)}.
\]

Here \(\theta_{\rm drag}^{\rm P18}\) is the same expression evaluated for the
fiducial Planck 2018 parameters.

It is initialized to `1` while the background distances are being built, then
updated after `z_drag` and the sound horizon are known.

This ratio is used to decide whether the existing flat-model `l` sampling and
flat spherical-Bessel table machinery can be reused. If
`abs(State%scale - 1) <= near_flat_scale_tol` with
`near_flat_scale_tol = 0.03`, equivalently \(|s_{\rm peak}-1|\le0.03\), CAMB
treats the peak shift as small enough to keep
the near-flat table strategy, enlarging the tabulated Bessel argument range when
needed for shifted-`q` calls. If the ratio is outside this tolerance, the
acoustic scale has moved enough that the `l` sampling is scaled by
`State%scale` instead of assuming the fiducial flat sampling.
