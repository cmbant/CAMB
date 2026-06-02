#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
venv_dir="${CAMB_DEVCONTAINER_VENV:-${UV_PROJECT_ENVIRONMENT:-/home/vscode/.local/share/camb-devcontainer/.venv}}"
venv_python="${venv_dir}/bin/python3"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
has_git_worktree=0
workspace_git_dir="${GIT_DIR:-}"
workspace_git_work_tree="${GIT_WORK_TREE:-}"

workspace_git() {
    local -a git_args

    git_args=()
    if [[ -n "${workspace_git_dir}" ]]; then
        git_args+=("--git-dir=${workspace_git_dir}")
    fi
    if [[ -n "${workspace_git_work_tree}" ]]; then
        git_args+=("--work-tree=${workspace_git_work_tree}")
    fi

    env -u GIT_DIR -u GIT_WORK_TREE -u GIT_COMMON_DIR -u GIT_INDEX_FILE git "${git_args[@]}" "$@"
}

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

ensure_codex_cli() {
    local npm_bin

    if command -v codex >/dev/null 2>&1; then
        return 0
    fi

    npm_bin="$(command -v npm || true)"
    if [[ -z "${npm_bin}" ]]; then
        echo "npm is required to install the Codex CLI." >&2
        exit 1
    fi

    sudo env "PATH=${PATH}" "${npm_bin}" install -g @openai/codex
}

cd "${workspace_dir}"

bash "${script_dir}/post-start.sh"

ensure_codex_cli

if workspace_git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    has_git_worktree=1
fi

unset GIT_DIR GIT_WORK_TREE GIT_COMMON_DIR GIT_INDEX_FILE

if [[ ! -f "${FORUTILSPATH:-${workspace_dir}/forutils}/Makefile" ]] && [[ "${has_git_worktree}" == "1" ]]; then
    workspace_git submodule update --init --recursive forutils
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
    if ! workspace_git config core.hooksPath .githooks; then
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
