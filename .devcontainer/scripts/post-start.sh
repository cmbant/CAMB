#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
host_gitconfig="/tmp/devcontainer-host.gitconfig"

safe_git() {
    env -u GIT_DIR -u GIT_WORK_TREE -u GIT_COMMON_DIR -u GIT_INDEX_FILE git "$@"
}

normalize_worktree_metadata() {
    local workspace_name git_file admin_dir gitdir_file

    git_file="${workspace_dir}/.git"
    if [[ ! -f "${git_file}" ]] || [[ -d "${git_file}" ]]; then
        return 0
    fi

    workspace_name="$(basename "${workspace_dir}")"
    admin_dir="$(cd "${workspace_dir}/.." && pwd)/.git/worktrees/${workspace_name}"
    gitdir_file="${admin_dir}/gitdir"

    if [[ ! -d "${admin_dir}" ]] || [[ ! -f "${gitdir_file}" ]]; then
        return 0
    fi

    printf 'gitdir: ../.git/worktrees/%s\n' "${workspace_name}" > "${git_file}"
    printf '../../%s/.git\n' "${workspace_name}" > "${gitdir_file}"
}

normalize_worktree_metadata

ensure_safe_directory() {
    local repo_dir="$1"
    local existing_dir

    if [[ ! -e "${repo_dir}/.git" ]]; then
        return 0
    fi

    while IFS= read -r existing_dir; do
        if [[ "${existing_dir}" == "${repo_dir}" ]]; then
            return 0
        fi
    done < <(safe_git config --global --get-all safe.directory || true)

    safe_git config --global --add safe.directory "${repo_dir}"
}

ensure_safe_directory "${workspace_dir}"
ensure_safe_directory "${workspace_dir}/forutils"

if [[ ! -f "${host_gitconfig}" ]]; then
    exit 0
fi

cd "${HOME}"

host_user_name="$(safe_git config --file "${host_gitconfig}" --get user.name || true)"
host_user_email="$(safe_git config --file "${host_gitconfig}" --get user.email || true)"

if [[ -n "${host_user_name}" ]] && [[ -z "$(safe_git config --global --get user.name || true)" ]]; then
    safe_git config --global user.name "${host_user_name}"
fi

if [[ -n "${host_user_email}" ]] && [[ -z "$(safe_git config --global --get user.email || true)" ]]; then
    safe_git config --global user.email "${host_user_email}"
fi
