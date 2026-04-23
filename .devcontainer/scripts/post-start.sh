#!/usr/bin/env bash

set -euo pipefail

host_gitconfig="/tmp/devcontainer-host.gitconfig"

if [[ ! -f "${host_gitconfig}" ]]; then
    exit 0
fi

cd "${HOME}"

host_user_name="$(git config --file "${host_gitconfig}" --get user.name || true)"
host_user_email="$(git config --file "${host_gitconfig}" --get user.email || true)"

if [[ -n "${host_user_name}" ]] && [[ -z "$(git config --global --get user.name || true)" ]]; then
    git config --global user.name "${host_user_name}"
fi

if [[ -n "${host_user_email}" ]] && [[ -z "$(git config --global --get user.email || true)" ]]; then
    git config --global user.email "${host_user_email}"
fi
