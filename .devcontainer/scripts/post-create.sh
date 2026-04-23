#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
venv_dir="${CAMB_DEVCONTAINER_VENV:-${UV_PROJECT_ENVIRONMENT:-/home/vscode/.local/share/camb-devcontainer/.venv}}"
venv_python="${venv_dir}/bin/python3"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
has_git_worktree=0

repair_writable_path() {
    local path="$1"

    if [[ ! -e "${path}" || -w "${path}" ]]; then
        return 0
    fi

    if ! command -v sudo >/dev/null 2>&1; then
        return 0
    fi

    sudo chown -R "$(id -u):$(id -g)" "${path}" >/dev/null 2>&1 || true
    sudo chmod -R u+rwX "${path}" >/dev/null 2>&1 || true
}

if ! command -v uv >/dev/null 2>&1; then
    echo "uv is required for the CAMB devcontainer bootstrap." >&2
    exit 1
fi

cd "${workspace_dir}"

bash "${script_dir}/post-start.sh"

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    has_git_worktree=1
fi

if [[ ! -f "${FORUTILSPATH:-${workspace_dir}/forutils}/Makefile" ]] && [[ "${has_git_worktree}" == "1" ]]; then
    git submodule update --init --recursive forutils
fi

repair_writable_path "${workspace_dir}/external"
repair_writable_path "${workspace_dir}/fortran/Releaselib"
repair_writable_path "${workspace_dir}/fortran/Debuglib"
repair_writable_path "${workspace_dir}/forutils/Releaselib"
repair_writable_path "${workspace_dir}/forutils/Debuglib"

bash .devcontainer/scripts/install-recombination-deps.sh

mkdir -p "$(dirname "${venv_dir}")"
export UV_PROJECT_ENVIRONMENT="${venv_dir}"

if [[ ! -x "${venv_python}" ]]; then
    uv venv "${venv_dir}" --python /usr/local/bin/python3
fi

mkdir -p "${PRE_COMMIT_HOME:-/home/vscode/.cache/pre-commit}"

# Install the small declared tool/runtime set directly, then do the editable CAMB build separately.
uv pip install \
    --upgrade \
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

if [[ "${has_git_worktree}" == "1" ]]; then
    if ! git config core.hooksPath .githooks; then
        echo "Failed to configure repository-local Git hooks." >&2
    fi

    if [[ -f .githooks/pre-commit ]]; then
        chmod +x .githooks/pre-commit >/dev/null 2>&1 || true
    fi

    if ! "${venv_python}" -m pre_commit install-hooks; then
        echo "pre-commit hook environment installation failed; offline hook runs may not work until you rerun '${venv_python} -m pre_commit install-hooks' with network access." >&2
    fi
fi

echo "CAMB devcontainer bootstrap finished."
