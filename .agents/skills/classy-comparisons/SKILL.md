---
name: classy-comparisons
description: "Use when comparing CAMB outputs against a local CLASS/classy install on Unix-like systems. Documents how to install the latest released CLASS version, report that version, and run the CLASS binary or Python wrapper."
argument-hint: "Describe what to compare, such as background quantities, CMB Cls, transfer functions, or matter power, and include the parameter set if known."
user-invocable: true
---

# CLASSy Comparisons

Use this skill when you need to compare CAMB against a local install of the
latest released CLASS/classy version available from the upstream `vX.Y.Z` tags.
CLASS/classy does not run on Windows, so these commands assume a Unix-like
shell.

## Installed Location

- Work directory: `${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}`
- Tarball: `${CLASS_WORKDIR}/class_public-${CLASS_VERSION}.tar.gz`
- Source tree: `${CLASS_WORKDIR}/class_public-${CLASS_VERSION}`
- CLASS binary: `${CLASS_WORKDIR}/class_public-${CLASS_VERSION}/class`
- Python wrapper target: `${CLASS_WORKDIR}/classy-local`

This install is local to the current machine. Reinstall with the latest
released CLASS tag and report the version:

```bash
CLASS_WORKDIR=${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}
export CLASS_WORKDIR
mkdir -p "${CLASS_WORKDIR}"
cd "${CLASS_WORKDIR}"
CLASS_TAG=$(
  python -c 'import re, subprocess
tags = subprocess.check_output(
    ["git", "ls-remote", "--tags", "--refs", "https://github.com/lesgourg/class_public.git", "v*"],
    text=True,
)
versions = []
for line in tags.splitlines():
    tag = line.rsplit("/", 1)[-1]
    if re.fullmatch(r"v\d+\.\d+\.\d+", tag):
        versions.append((tuple(map(int, tag[1:].split("."))), tag))
print(max(versions)[1])'
)
CLASS_VERSION=${CLASS_TAG#v}
CLASS_ROOT=${CLASS_WORKDIR}/class_public-${CLASS_VERSION}
CLASS_TARBALL=${CLASS_WORKDIR}/class_public-${CLASS_VERSION}.tar.gz
CLASS_URL=https://github.com/lesgourg/class_public/archive/refs/tags/${CLASS_TAG}.tar.gz
echo "Installing CLASS ${CLASS_VERSION} from ${CLASS_URL}"
rm -rf "${CLASS_ROOT}"
curl -L -o "${CLASS_TARBALL}" "${CLASS_URL}"
mkdir -p "${CLASS_ROOT}"
tar -xzf "${CLASS_TARBALL}" -C "${CLASS_ROOT}" --strip-components=1
cd "${CLASS_ROOT}"
make -j"$(python -c 'import os; print(os.cpu_count() or 1)')"
rm -rf "${CLASS_WORKDIR}/classy-local"
python -m pip install --target "${CLASS_WORKDIR}/classy-local" .
echo "Installed CLASS ${CLASS_VERSION} in ${CLASS_ROOT}"
```

## Quick Checks

Run the CLASS binary with the reference ini file:

```bash
CLASS_WORKDIR=${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}
export CLASS_WORKDIR
CLASS_ROOT=$(
  python -c 'import os, pathlib
roots = sorted(pathlib.Path(os.environ["CLASS_WORKDIR"]).glob("class_public-[0-9]*"))
print(roots[-1])'
)
echo "Using CLASS ${CLASS_ROOT##*/class_public-} from ${CLASS_ROOT}"
cd "${CLASS_ROOT}"
./class explanatory.ini
```

Check that the Python wrapper imports:

```bash
CLASS_WORKDIR=${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}
export CLASS_WORKDIR
PYTHONPATH="${CLASS_WORKDIR}/classy-local${PYTHONPATH:+:$PYTHONPATH}" python -c "from classy import Class; import os, pathlib; root = sorted(pathlib.Path(os.environ['CLASS_WORKDIR']).glob('class_public-[0-9]*'))[-1]; print(f'Using CLASS {root.name.removeprefix(\"class_public-\")}'); print(Class)"
```

## How to Run Comparisons

There is no dedicated CAMB-vs-CLASS driver in CAMB at the moment, so
comparisons should be run from the root of the CAMB checkout with the local
`classy` module on `PYTHONPATH`.

1. Start from the root of the CAMB checkout:

```bash
cd /path/to/camb
```

2. Make the local `classy` install visible to Python:

```bash
CLASS_WORKDIR=${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}
export CLASS_WORKDIR
export PYTHONPATH="${CLASS_WORKDIR}/classy-local${PYTHONPATH:+:$PYTHONPATH}"
```

3. Run a comparison driver from the repo root so `import camb` resolves to the
CAMB checkout:

```bash
python your_comparison_driver.py
```

## Matched Conventions

Always match recombination, helium, CMB temperature, lensing, units, and output
normalization explicitly. For RECFAST comparisons:

```python
pars.set_classes(recombination_model="Recfast")
class_params["recombination"] = "RECFAST"
```

Set the helium mass fraction explicitly on both sides; do not rely on a BBN
default or a code-specific inferred value:

```python
pars.set_cosmology(..., YHe=0.245, TCMB=2.7255)
class_params["YHe"] = 0.245
class_params["T_cmb"] = 2.7255
```

For Planck-era RECFAST matching in CAMB, use
`pars.Recomb.set_params(recfast_approx_model="planck")`.  CLASS/RecFastCLASS
defaults are close to these Planck-era RECFAST knobs, but robust comparison
drivers should still set them explicitly when the exact RECFAST approximation is
part of the test.

## CAMB Neutrino Mapping

Map CAMB massive neutrinos after `pars.set_cosmology(...)` as one CLASS
`ncdm` species per CAMB mass eigenstate:

```python
T_std = (4 / 11) ** (1 / 3)
n = pars.nu_mass_eigenstates
class_params.update({
    "N_ur": pars.num_nu_massless,
    "N_ncdm": n,
    "deg_ncdm": ",".join(str(pars.nu_mass_numbers[i]) for i in range(n)),
    "T_ncdm": ",".join(
        str(T_std * (pars.nu_mass_degeneracies[i] / pars.nu_mass_numbers[i]) ** 0.25)
        for i in range(n)
    ),
    "omega_ncdm": ",".join(str(pars.omnuh2 * pars.nu_mass_fractions[i]) for i in range(n)),
})
```

Use `omega_ncdm`, not `m_ncdm`, when trying to reproduce CAMB exactly: CAMB's
`mnu` is first converted to `omnuh2`, then CAMB solves the thermal mass that
gives that density, including a small velocity correction.

Example: CAMB `mnu=0.21, nnu=3.044, num_massive_neutrinos=3,
neutrino_hierarchy="degenerate"` becomes one CLASS species with
`N_ur=0`, `N_ncdm=1`, `deg_ncdm=3`, `T_ncdm=(4/11)**(1/3)*(3.044/3)**0.25`,
and `omega_ncdm=pars.omnuh2`.

For a hand-specified CAMB model with general eigenstate degeneracies, use the
same rule. For example, if CAMB has two eigenstates with
`nu_mass_numbers=[2, 1]`, `nu_mass_degeneracies=[2.1, 0.944]`, and
`nu_mass_fractions=[0.75, 0.25]`, then CLASS should use
`deg_ncdm="2,1"`,
`T_ncdm=f"{T_std*(2.1/2)**0.25},{T_std*(0.944/1)**0.25}"`, and
`omega_ncdm=f"{pars.omnuh2*0.75},{pars.omnuh2*0.25}"`.

CLASS neutrino-related parameters to check:

- Physics mapping: `N_ur`, `N_ncdm`, `omega_ncdm`/`Omega_ncdm`/`m_ncdm`,
  `T_ncdm`, `deg_ncdm`, and `ksi_ncdm`. For CAMB standard thermal neutrinos
  use `ksi_ncdm=0` or leave it unset.
- Nonthermal distributions: `use_ncdm_psd_files`, `ncdm_psd_filenames`, and
  `ncdm_psd_parameters`. Leave unset for CAMB thermal Fermi-Dirac neutrinos.
- Numerical neutrino controls: `tol_ncdm_bg`, `tol_ncdm`, `tol_M_ncdm`,
  `ncdm_N_momentum_bins`, `ncdm_quadrature_strategy`, `ncdm_maximum_q`,
  `l_max_ncdm`, `ncdm_fluid_approximation`, and
  `ncdm_fluid_trigger_tau_over_tau_k`. Record these when changed.
- Massless-neutrino perturbation parameters: `ceff2_ur`, `cvis2_ur`,
  `ur_fluid_approximation`, and `ur_fluid_trigger_tau_over_tau_k`. Leave at
  defaults for standard CAMB comparisons unless deliberately testing altered
  massless-neutrino physics.

## Comparison Rules

- Match cosmological parameters explicitly. Do not assume CAMB and CLASS
  defaults use the same neutrino, helium, or recombination conventions.
- Compare like with like: unlensed with unlensed, lensed with lensed, raw `C_l`
  with raw `C_l`, and the same unit normalization.
- Put both codes on a common `ell`, `k`, or redshift grid before taking
  differences.
- Record any non-default CLASS precision file or precision overrides used in the
  comparison.
- When debugging a discrepancy, first confirm that CLASS runs standalone before
  changing CAMB or CLASS precision settings.

## Lensed CMB Accuracy Checks

- Lensed spectra have high-ell convolution edge effects. For a plotted range
  ending at `lmax_plot`, compute CLASS lensed spectra beyond the plotted range
  and truncate only for the comparison. As a starting point, use
  `l_max_scalars = lmax_plot + 2000`
- For purely numerical comparisons, match the effective lensing/source cutoff,
  not just the named accuracy setting. If CAMB's reference has `max_eta_k`, first
  convert it to a physical source cutoff using CAMB's conformal age:

```python
camb_kmax_1Mpc = camb_max_eta_k / camb_eta0
```

- For CLASS 3.2 and later, the default CMB lensing-potential calculation uses
  the full-Limber scheme (`want_lcmb_full_limber = yes`) for `lCl`. Match the
  CAMB source reach by setting

```python
class_params["k_max_limber_over_l_max_scalars"] = camb_kmax_1Mpc / class_l_max_scalars
```

  Here `class_l_max_scalars` is the CLASS compute value, including any
  high-ell margin used before truncating to the plotted range. If you change the
  CLASS compute `l_max_scalars`, recompute this ratio so the physical
  `k_max [1/Mpc]` remains fixed.

- Do not use `k_max_tau0_over_l_max` to match the high-k CMB lensing reach for
  CLASS 3.2+ full-Limber `lCl` comparisons. That parameter controls the primary
  CMB transfer/source cutoff and can make flat lensed-CMB runs much slower
  without matching the full-Limber lensing cutoff being tested. Use it only when
  deliberately matching primary-CMB transfer cutoffs, or explicitly setting
  `want_lcmb_full_limber = no`.

- `P_k_max_h/Mpc` and `P_k_max_1/Mpc` control output matter-power calculations.
  They are not the CMB-lensing source cutoff and should not be used to match
  CAMB lensed-CMB source reach unless the comparison is explicitly about
  exported matter power.


## Timing Hygiene

- Set thread-control environment variable OMP_NUM_THREADS before importing CLASS if specific thread number is needed. Disgard first run which may be doing
one off cached calculations (esp. for non-flat).


## CLASS-native Commands

Useful CLASS commands from the extracted source tree:

```bash
CLASS_WORKDIR=${CLASS_WORKDIR:-${TMPDIR:-/tmp}/classy-comparisons}
export CLASS_WORKDIR
CLASS_ROOT=$(
  python -c 'import os, pathlib
roots = sorted(pathlib.Path(os.environ["CLASS_WORKDIR"]).glob("class_public-[0-9]*"))
print(roots[-1])'
)
echo "Using CLASS ${CLASS_ROOT##*/class_public-} from ${CLASS_ROOT}"
cd "${CLASS_ROOT}"
./class default.ini
./class explanatory.ini cl_permille.pre
python CPU.py -h
```
