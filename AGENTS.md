# Project Overview

This is the CAMB Python cosmology code, wrapping Fortran 2003 for numerics.

---

# Coding Standards

## Python
- Use Python 3.10+ with modern type-hint syntax.
- Use double-quotes for strings; maximum line length is 120.
- Code is auto-formatted by **ruff**, including import reordering.
- Use the **pre-commit** command line to run ruff check/format (e.g., `pre-commit run --all-files`); do not invoke `ruff` directly.
- Ruff/pyupgrade hooks and their configuration live in `.pre-commit-config.yaml`
- All code must be Python or Fortran only—no bash scripts except for devcontainer setup and GitHub Actions.

## Fortran
- Use modern Fortran 2003, compiled with **gfortran** or **ifort**.
- `/forutils` is a git submodule containing shared Fortran utilities and classes.
- For guidance on modifying Fortran or wrapped Fortran code, see `/docs/modifying_code.rst`.

---

# Development Workflow

## Installation & Build
- Use `pip install -e .` for the initial editable install.
- After the first install, use `python setup.py make` to quickly rebuild the Fortran for use with Python after any fortran changes.
- If `python setup.py make` produces stale `.mod` type-mismatch errors, run `python setup.py clean` first (not normally needed).

## Testing
- **Default test:** `python -m unittest camb.tests.camb_test`.
  - Do **not** run `camb.tests.hmcode_test` unless specifically relevant.
- **Long test** against precomputed results: `python fortran/tests/run_tests.py`, which calls `CAMB_test_files.py`. Only run this if explicitly asked.

## Configuration
- See `CONTRIBUTING.md` for installation, cloning, modification guidelines, contribution instructions, pre-commit setup, and VS Code configuration.

---

# Maintenance & Versioning
- When updating the version, update **both**:
  - the Python `__version__`, and
  - the version defined near the top of `fortran/config.f90`.
- When adding new `.ini` file parameters, also add the corresponding `.ini` write logic in `camb/_ini.py`.

---

# Physics & Numerical Guidance
- When plotting fractional differences for cross-spectra, use e.g. $\Delta TE / \sqrt{TT \cdot EE}$.
- High-accuracy lensing requires `lens_potential_accuracy = 8` or so (this mainly increases `kmax`).
- Use `camb/check_accuracy.py` to assess numerical stability by comparing a default run to a high-accuracy run.
