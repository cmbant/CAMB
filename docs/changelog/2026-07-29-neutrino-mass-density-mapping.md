# Neutrino mass to density mapping

- Date: 2026-07-29
- Files: `fortran/massive_neutrinos.f90`, `fortran/model.f90`, `camb/model.py`

## Change

`set_cosmology(mnu=...)` now maps the neutrino mass sum to `Omega_nu h^2` through
the actual thermal density, rather than the non-relativistic proxy
`omnuh2 = mnu/94.07 * (standard_neutrino_neff/3)**0.75` previously computed in
Python. `mnu` is now the physical sum of eigenstate masses in eV, and the density
it implies includes the finite-temperature correction (it still neglects spectral
distortions to the neutrino distribution).

The pieces are:

1. `ThermalNuBackground%omnuh2_from_mass(mass, degeneracy, number, TCMB)` returns
   the present-day density of an eigenstate of `number` neutrinos of physical mass
   `mass`, sharing relativistic degeneracy `degeneracy`:

   ```
   omnuh2 = degeneracy * T_nu * (TCMB/COBE_CMBTemp)**3 * rho(am)
            / (neutrino_mass_fac * nu_fit_rho_scale)
   ```

   with `T_nu` the standard neutrino temperature in eV and `am = mass/(T_nu *
   (degeneracy/number)**0.25)`. Since `rho -> nu_fit_rho_scale*am` when
   non-relativistic, this reduces to the old proxy for heavy eigenstates, and is
   the exact inverse of the `omnuh2 -> find_nu_mass_for_rho -> nu_masses` mapping
   in `results.f90`.

2. `CAMBparams%SetNeutrinoMasses(mnu, omnuh2_sterile, nnu, neutrino_hierarchy,
   num_massive_neutrinos, standard_neutrino_neff)` is the `mnu` counterpart of
   `SetNeutrinoHierarchy`. It first configures the hierarchy from the old
   non-relativistic proxy — `SetNeutrinoHierarchyBase` inverts that proxy exactly,
   so the mass splitting it produces is the one requested — then replaces each
   active eigenstate density with `omnuh2_from_mass` of its physical mass and
   renormalizes `Nu_mass_fractions`.

3. `standard_neutrino_neff` is now passed to Fortran for both entry points, as an
   optional argument on `SetNeutrinoHierarchy` defaulting to `default_nnu`.
   Previously the Fortran hierarchy code always used `default_nnu = 3.044` for the
   eigenstate degeneracies and mass splitting while Python used
   `standard_neutrino_neff` for the density conversion, so a non-default value was
   applied inconsistently (and, for `omnuh2_active`, not at all).

`omnuh2_active` keeps its meaning: it sets `Omega_nu h^2` for the active neutrinos
directly. `meffsterile` also keeps its Planck-paper definition,
`Omega_sterile h^2 = meffsterile/94.07`, so the sterile component is unchanged and
only the active neutrinos use the new mapping.

## Result changes at fixed mnu

`Omega_nu h^2` increases by roughly `15*zeta5/(2*zeta3) / am**2` per eigenstate,
i.e. it grows as the eigenstate masses fall:

| case | old `omnuh2` | new `omnuh2` | shift |
| --- | --- | --- | --- |
| 1 massive, `mnu=0.06` | 6.448666e-4 | 6.448994e-4 | +5.1e-5 |
| 3 degenerate, `mnu=0.06` | 6.448666e-4 | 6.451617e-4 | +4.6e-4 |
| 3 degenerate, `mnu=0.12` | 1.289733e-3 | 1.289881e-3 | +1.1e-4 |
| normal, `mnu=0.12` | 1.289733e-3 | 1.289896e-3 | +1.3e-4 |
| inverted, `mnu=0.11` | 1.182255e-3 | 1.182544e-3 | +2.4e-4 |
| 3 degenerate, `mnu=1.0` | 1.074778e-2 | 1.074779e-2 | +1.7e-6 |
| 3 degenerate, `mnu=0.01` | 1.074778e-4 | 1.092166e-4 | +1.6e-2 |

For minimal-mass models the shift is equivalent to changing `mnu` by under
`3e-5 eV`, which is far below any current or planned sensitivity. The two-eigenstate
mass fractions change slightly because the lighter eigenstate gets the larger
correction: for `mnu=0.11` inverted, `nu_mass_fractions[0]` moves from `0.915197`
to `0.915040`.

## Very light neutrinos

The old proxy could ask for a density below the massless value, which
`find_nu_mass_for_rho` correctly refuses to invert: `mnu` below roughly

| `num_massive_neutrinos` | relativistic floor `omnuh2` | old proxy `mnu` |
| --- | --- | --- |
| 1 | 5.6987e-6 | 5.30e-4 eV |
| 2 | 1.1397e-5 | 1.06e-3 eV |
| 3 | 1.7096e-5 | 1.59e-3 eV |

used to invert to zero mass, with `Omega_nu` then set to zero. Such masses now give
real eigenstate masses and the correct (relativistic) density.

Note that `omnuh2` therefore does not tend to zero as `mnu -> 0+`; it tends to the
relativistic floor above, because those neutrinos are still being counted as massive
eigenstates. Only `mnu = 0` exactly moves them into `Num_Nu_Massless`. The total
density is the same either way — `Omega_de` is continuous to 1e-12 across the
transition — but the split between `get_Omega("nu")` and `get_Omega("neutrino")`
changes.

## Round-trip accuracy

`mnu -> Omega_nu h^2 -> nu_masses` reproduces the input mass sum to better than
`1.3e-6 eV` absolute over `mnu` from `1e-5` to `10 eV`, for all hierarchies and
`num_massive_neutrinos` 1-3. The worst relative error, `4.9e-3` near
`mnu = 2.6e-4 eV`, is set by the low-mass branches of `find_nu_mass_for_rho`; in
absolute terms it is ~1e-6 eV and has no cosmological impact, so those branches were
left alone. `testNuMassRoundTrip` checks the round trip to `1e-5 eV` absolute across
the branch structure, checks the massive-eigenstate degeneracies against
`num_nu_massive * standard_neutrino_neff / 3`, and checks that feeding the resulting
`omnuh2` back in via `omnuh2_active` gives the same heating (and, for degenerate
masses, the identical model).

## Related fixes in the same series

- `find_nu_mass_for_rho` inversion accuracy (`ebac155a`): the near-non-relativistic
  branch now uses the analytic asymptotic inversion plus one Newton step with
  `drho`, instead of a finite-difference derivative from a second `rho` call; a
  low-mass inverse series was added for `rho <= 1.005`; and `rho <= 1` (unphysical,
  below the massless value) is separated out. This also fixed an aliasing bug where
  `nu_mass` was passed to `brentq` as both the `intent(in)` bracket end and the
  `intent(out)` root.
- Inverted hierarchy inversion hole (`e35458c4`): the two-eigenstate solve now
  solves for the lightest mass `m3` from zero using `sum_mnu_for_m3` and sets
  `m1 = sqrt(m3**2 + delta_mnu31)`, instead of Newton-solving for `m1` from a
  bracket starting at `sqrt(delta_mnu31)`, which failed just above the minimum
  inverted mass sum. Covered by `testInvertedNuMassThreshold`.

## Compatibility

- Results at fixed `mnu` change by the amounts tabulated above. To reproduce old
  numbers exactly, pass the old proxy density as `omnuh2_active` (with `mnu=None`)
  instead of `mnu`.
- The Fortran `SetNeutrinoHierarchy` interface gained a trailing optional
  `standard_neutrino_neff`; existing positional calls are unaffected.
