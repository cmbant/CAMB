#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
venv_dir="${CAMB_DEVCONTAINER_VENV:-${UV_PROJECT_ENVIRONMENT:-/home/vscode/.local/share/camb-devcontainer/.venv}}"
venv_python="${venv_dir}/bin/python3"

if ! command -v uv >/dev/null 2>&1; then
    echo "uv is required for the CAMB devcontainer bootstrap." >&2
    exit 1
fi

cd "${workspace_dir}"
bash .devcontainer/scripts/install-recombination-deps.sh

mkdir -p "$(dirname "${venv_dir}")"
export UV_PROJECT_ENVIRONMENT="${venv_dir}"

if [[ ! -x "${venv_python}" ]]; then
    uv venv "${venv_dir}" --python /usr/local/bin/python3
fi

# Install the small declared tool/runtime set directly, then do the editable CAMB build separately.
uv pip install \
    --python "${venv_python}" \
    pip \
    setuptools \
    wheel \
    numpy \
    scipy \
    sympy \
    packaging \
    ruff \
    pre-commit

site_packages_dir="$("${venv_python}" - <<'PY'
import site

print(site.getsitepackages()[0])
PY
)"

# CAMB's clean Fortran build is not dependency-safe under parallel make, so force a serial build here.
MAKEFLAGS= "${venv_python}" setup.py make

# Modern editable installs spend a long time in setuptools metadata generation on this bind mount.
# Add the workspace directly to site-packages so imports resolve immediately in the devcontainer.
printf '%s\n' "${workspace_dir}" > "${site_packages_dir}/camb-devcontainer.pth"

if [[ -d .git ]]; then
    if ! git config core.hooksPath .githooks; then
        echo "Failed to configure repository-local Git hooks." >&2
    fi

    if [[ "${CAMB_INSTALL_PRECOMMIT:-0}" == "1" ]]; then
        if ! uv run pre-commit install-hooks; then
            echo "pre-commit hook environment installation failed; run 'uv run pre-commit install-hooks' manually if needed." >&2
        fi
    fi
fi

echo "CAMB devcontainer bootstrap finished."
