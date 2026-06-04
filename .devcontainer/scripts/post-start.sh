#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
host_gitconfig="/tmp/devcontainer-host.gitconfig"

safe_git() {
    env -u GIT_DIR -u GIT_WORK_TREE -u GIT_COMMON_DIR -u GIT_INDEX_FILE git "$@"
}

ensure_writable_dir() {
    local path="$1"

    mkdir -p "${path}"

    if [[ -w "${path}" ]]; then
        return 0
    fi

    if command -v sudo >/dev/null 2>&1; then
        sudo chown "$(id -u):$(id -g)" "${path}" >/dev/null 2>&1 || true
    fi
}

ensure_writable_dir "${HOME}/.codex"

ensure_codex_full_access() {
    local config_dir config_file tmp_file

    config_dir="${HOME}/.codex"
    config_file="${config_dir}/config.toml"

    mkdir -p "${config_dir}"
    tmp_file="$(mktemp)"

    if [[ -f "${config_file}" ]]; then
        awk '
            BEGIN { in_table = 0 }
            /^[[:space:]]*\[/ { in_table = 1 }
            !in_table && /^[[:space:]]*(approval_policy|default_permissions|sandbox_mode)[[:space:]]*=/ { next }
            { print }
        ' "${config_file}" > "${tmp_file}"
    fi

    {
        printf 'approval_policy = "never"\n'
        printf 'default_permissions = ":danger-full-access"\n'
        if [[ -s "${tmp_file}" ]]; then
            printf '\n'
            cat "${tmp_file}"
        fi
    } > "${config_file}"

    rm -f "${tmp_file}"
}

ensure_codex_full_access

normalize_worktree_metadata() {
    local admin_dir candidate git_file gitdir_file gitdir_path gitdir_path_normalized search_dir worktree_id
    local admin_gitdir_rel workspace_gitdir_rel

    git_file="${workspace_dir}/.git"
    if [[ ! -f "${git_file}" ]] || [[ -d "${git_file}" ]]; then
        return 0
    fi

    gitdir_path="$(sed -n 's/^[[:space:]]*gitdir:[[:space:]]*//p' "${git_file}" | head -n 1)"
    if [[ -z "${gitdir_path}" ]]; then
        return 0
    fi

    gitdir_path_normalized="${gitdir_path//\\//}"
    case "${gitdir_path_normalized}" in
        /*)
            if [[ -d "${gitdir_path_normalized}" ]]; then
                admin_dir="$(realpath --no-symlinks "${gitdir_path_normalized}")"
            fi
            ;;
        *:/*)
            ;;
        *)
            candidate="${workspace_dir}/${gitdir_path_normalized}"
            if [[ -d "${candidate}" ]]; then
                admin_dir="$(realpath --no-symlinks "${candidate}")"
            fi
            ;;
    esac

    if [[ -z "${admin_dir:-}" ]]; then
        case "${gitdir_path_normalized}" in
            *"/worktrees/"*)
                worktree_id="${gitdir_path_normalized##*/worktrees/}"
                search_dir="${workspace_dir}"
                while [[ "${search_dir}" != "/" ]]; do
                    candidate="${search_dir}/.git/worktrees/${worktree_id}"
                    if [[ -d "${candidate}" ]]; then
                        admin_dir="$(realpath --no-symlinks "${candidate}")"
                        break
                    fi
                    search_dir="$(dirname "${search_dir}")"
                done
                ;;
        esac
    fi

    gitdir_file="${admin_dir:-}/gitdir"
    if [[ -z "${admin_dir:-}" ]] || [[ ! -f "${gitdir_file}" ]]; then
        return 0
    fi

    workspace_gitdir_rel="$(realpath --no-symlinks --relative-to="${workspace_dir}" "${admin_dir}")"
    admin_gitdir_rel="$(realpath --no-symlinks --relative-to="${admin_dir}" "${git_file}")"

    printf 'gitdir: %s\n' "${workspace_gitdir_rel}" > "${git_file}"
    printf '%s\n' "${admin_gitdir_rel}" > "${gitdir_file}"
}

normalize_worktree_metadata

if [[ -f "${workspace_dir}/.githooks/pre-commit" ]]; then
    chmod +x "${workspace_dir}/.githooks/pre-commit" >/dev/null 2>&1 || true
fi

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



setup_workspace_claude_link() {
    local claude_link="${workspace_dir}/.claude"
    local agents_dir="${workspace_dir}/.agents"

    if [[ ! -d "${agents_dir}" ]]; then
        return 0
    fi

    if [[ -L "${claude_link}" || -e "${claude_link}" ]]; then
        return 0
    fi

    ln -s ".agents" "${claude_link}"
}

setup_workspace_claude_link

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
