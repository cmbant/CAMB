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

# Keep any caller-provided MAKEFLAGS; CAMB falls back to a one-time serial sub-build
# automatically when compiler-generated dependency files do not exist yet.
"${venv_python}" setup.py make

# Modern editable installs spend a long time in setuptools metadata generation on this bind mount.
# Add the workspace directly to site-packages so imports resolve immediately in the devcontainer.
printf '%s\n' "${workspace_dir}" > "${site_packages_dir}/camb-devcontainer.pth"

if [[ -d .git ]]; then
    if ! git config core.hooksPath .githooks; then
        echo "Failed to configure repository-local Git hooks." >&2
    fi

    if [[ -f .githooks/pre-commit ]]; then
        chmod +x .githooks/pre-commit || true
    fi
fi

echo "CAMB devcontainer bootstrap finished."
