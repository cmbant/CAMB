# CAMB devcontainer

This devcontainer uses Python 3.14 and caps the container and default build runtime to 8 cores.

The Python environment is created inside the container filesystem, not in the workspace checkout, so it does not fight with host-mounted file permissions or filesystem performance.

When the container is created it:

- installs the CAMB build toolchain, including `gfortran`, `g++`, `make`, and `libgsl-dev`
- installs `uv` for locked, repeatable Python environment syncs
- downloads CosmoRec into `external/CosmoRec`
- clones HYREC-2 into `external/HYREC-2`
- patches the CosmoRec makefile to add `-fPIC`, which CAMB's Python wrapper requires
- installs the small Python dependency/tool set into a container-local virtual environment, rebuilds `camb` with `recfast`, `cosmorec`, and `hyrec` enabled, and adds the workspace to the environment with a `.pth` file
- configures `git core.hooksPath` to use the tracked `.githooks/` directory instead of mutating `.git/hooks`

The bootstrap intentionally avoids `pip install -e .` during container creation:

- `uv pip install ...` to install the small runtime and dev tool set directly
- `python setup.py make` to rebuild the linked Fortran library with the optional recombination backends enabled
- a small `.pth` file in the container-local environment so the workspace checkout is importable immediately without waiting on slow editable metadata generation

The container no longer exports a global `MAKEFLAGS`. CAMB's first clean Fortran build can fail under parallel make because module dependencies are not fully serialized, so the editable build step forces `MAKEFLAGS=` for reliability.

The `external/` directory stays in the workspace, so the downloaded sources are visible from Windows as normal files.

Useful commands inside the container:

```bash
uv pip install --python /home/vscode/.local/share/camb-devcontainer/.venv/bin/python3 numpy scipy sympy packaging ruff pre-commit
uv run python setup.py make
uv run python -m unittest camb.tests.camb_test
```

If you explicitly want a standards-based editable install after the container is up, you can still run `python -m pip install --no-build-isolation --no-deps -e .`, but that setuptools metadata phase is the part that remains slow on this bind-mounted checkout.

If you want Git hooks installed from inside the container, set `CAMB_INSTALL_PRECOMMIT=1` in the devcontainer environment and rebuild. On some host-mounted workspaces this still fails because `.git/hooks` does not allow the chmod operations that `pre-commit` performs.

The repository-local `.githooks/pre-commit` wrapper avoids that `.git/hooks` permission problem and works from either the container-local environment, a workspace `.venv`, or any `pre-commit` available on `PATH`.
