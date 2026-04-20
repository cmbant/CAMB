# Project instructions
- This is the CAMB python cosmology code source, wrapping fortran 2003 for numerics
- Follow the rules in this file.
- If `AGENTS.local.md` exists, read it after this file and let it override local environment details.

# Project rules and guidance
- Python uses python 3.10+ with modern type hint syntax
- Use double-quotes for strings; line length 120
- Code is auto-formatted by ruff, including import reordering.
- Fortran uses modern Fortran 2003, compiled with gfortran or ifort
- Install and run pre-commit hooks as needed for checking
- /forutils is a git submodule with shared fortran utilities and classes
- For fortran/wrapped fortran code modifications see /docs/modifying_code.rst
- When updating version, update both python __version__ and version defined near the top of fortran/config.f90
- Use pip -e for editable first install. After that use "python setup.py make" to quickly rebuild fortran for use with python.
- Default test: python -m unittest camb.tests.camb_test (do not run hmcode_test unless specifically relevant)
- Installation/clone/modification/contributing/pre-commit/vscode config: see /CONTRIBUTING.md
- When plotting fractional differences for cross-spectra, use e.g. Delta TE/sqrt(TT*EE)
- High accuracy lensing needs lens_potential_accuracy = 8 or so (mainly increases kmax)
