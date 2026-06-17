# CAMB devcontainer

This devcontainer uses Python 3.14 and caps the container and default build runtime to 8 cores.

The Python environment is created inside the container filesystem, not in the workspace checkout, so it does not fight with host-mounted file permissions or filesystem performance.

The repository-level VS Code settings stay host-friendly and only help VS Code discover a workspace `./.venv`.
The devcontainer sets its own interpreter path separately, so a Windows host virtual environment and the Linux
container virtual environment can coexist cleanly.

The devcontainer is designed to open either the main CAMB checkout or a linked Git worktree below that checkout, such
as `testtree`. It mounts the parent directory of the opened folder so the container can see the repository `.git`
metadata, and the startup script normalizes linked-worktree metadata to portable relative paths so the same checkout
works from both Windows and the Linux container.

Git inside the container is also told to ignore executable-bit changes for the mounted workspace. On Windows-backed
bind mounts, Linux often reports regular files as executable even when the Git index records them as non-executable,
which would otherwise produce spurious `100644 => 100755` diffs.

When the container is created it:

- installs the CAMB build toolchain, including `gfortran`, `g++`, `make`, and `libgsl-dev`
- installs `uv` for locked, repeatable Python environment syncs
- creates the container-local Python virtual environment and pre-installs the small Python dependency/tool set
- downloads CosmoRec into `external/CosmoRec`
- clones HYREC-2 into `external/HYREC-2`
- patches the CosmoRec makefile to add `-fPIC`, which CAMB's Python wrapper requires
- rebuilds `camb` with `recfast`, `cosmorec`, and `hyrec` enabled, and adds the workspace to the environment with a `.pth` file, but leaves recfast default for new builds.
- configures `core.hooksPath` to use the tracked `.githooks/pre-commit` wrapper, which can use `pre-commit` from the container venv or from a host `.venv`
- pre-installs the `pre-commit` hook environments into a container-local cache so later hook runs can work offline

The devcontainer also installs `findent` and `fortls` into the container-local virtual environment so the Modern Fortran
extension can format and index code without depending on host-side tools.

Codex state is stored in the Docker volume `camb-codex`, so sessions and container-side caches survive devcontainer
rebuilds and are shared across CAMB worktrees using this devcontainer. Only the host `~/.codex/auth.json` file is bind
mounted into the volume, which gives the container the host login without mounting host-specific Codex logs, caches,
SQLite state, or plugin files.

The container is configured as a trusted AI-agent sandbox: Codex starts with `approval_policy = "never"` and
`default_permissions = ":danger-full-access"` in its container `config.toml`, and Claude Code starts new sessions in
`bypassPermissions` mode.

The image and bootstrap intentionally avoid `pip install -e .`:

- the Dockerfile uses `uv pip install ...` to install the small runtime and dev tool set directly into the container-local virtual environment
- `post-create.sh` runs `python setup.py make` to rebuild the linked Fortran library with the optional recombination backends enabled
- `post-create.sh` writes a small `.pth` file in the container-local environment so the workspace checkout is importable immediately without waiting on slow editable metadata generation

Clean-tree CAMB builds automatically do the first Fortran sub-build with `-j1` until compiler-generated `.d` files
exist, after which `make -j` and `python setup.py make` can safely respect any `MAKEFLAGS` you set.

The first container bootstrap still needs network access to fetch the hook repositories listed in
`.pre-commit-config.yaml`, but once bootstrap succeeds the tracked Git hook can run later without internet access.

The `external/` directory stays in the workspace, so the downloaded sources are visible from Windows as normal files.

Useful commands inside the container:

```bash
uv pip install --python /home/vscode/.local/share/camb-devcontainer/.venv/bin/python3 numpy scipy sympy packaging ruff pre-commit
uv run python setup.py make
uv run python -m unittest camb.tests.camb_test
```

On the host, keep using a normal workspace `.venv` if you want. The tracked Git hook resolves the right
`pre-commit` executable at runtime instead of baking in one interpreter path during installation.

Because the workspace may still contain a host-side `.venv` that is not usable inside Linux, tools that do not honor the
devcontainer's selected interpreter may need to be pointed explicitly at
`/home/vscode/.local/share/camb-devcontainer/.venv/bin/python3`.

If you explicitly want a standards-based editable install after the container is up, you can still run `python -m pip install --no-build-isolation --no-deps -e .`.


The container includes builds of CosmoRec and HyRec. If needed, you can build with these by doing `export RECOMBINATION_FILES="recfast cosmorec hyrec"`
