#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
external_dir="${workspace_dir}/external"
cosmorec_dir="${COSMOREC_PATH:-${external_dir}/CosmoRec}"
hyrec_dir="${HYREC_PATH:-${external_dir}/HYREC-2}"
cosmorec_url="${COSMOREC_URL:-https://www.cita.utoronto.ca/~jchluba/Recombination/_Downloads_/CosmoRec.v2.0.3b.tar.gz}"
hyrec_repo="${HYREC_REPO:-https://github.com/nanoomlee/HYREC-2.git}"

safe_git() {
    env -u GIT_DIR -u GIT_WORK_TREE -u GIT_COMMON_DIR -u GIT_INDEX_FILE git "$@"
}

cosmorec_dir="${cosmorec_dir%/}"
hyrec_dir="${hyrec_dir%/}"
patched_cosmorec_flags=0

mkdir -p "${external_dir}"

if [[ ! -f "${cosmorec_dir}/Makefile" ]]; then
    archive_path="${external_dir}/CosmoRec.v2.0.3b.tar.gz"
    temp_dir="$(mktemp -d)"
    trap 'rm -rf "${temp_dir}"' EXIT

    if [[ ! -f "${archive_path}" ]]; then
        curl --fail --show-error --location --retry 3 --retry-delay 2 "${cosmorec_url}" -o "${archive_path}"
    fi

    mkdir -p "${temp_dir}/extract"
    tar --no-same-owner --no-same-permissions -xzf "${archive_path}" -C "${temp_dir}/extract"

    extracted_root="$(find "${temp_dir}/extract" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
    if [[ -z "${extracted_root}" ]]; then
        echo "Failed to locate extracted CosmoRec sources" >&2
        exit 1
    fi

    rm -rf "${cosmorec_dir}"
    mv "${extracted_root}" "${cosmorec_dir}"
fi

chmod -R u+rwX "${cosmorec_dir}" >/dev/null 2>&1 || true

cosmorec_makefile_in="${cosmorec_dir}/Makefile.in"
if [[ -f "${cosmorec_makefile_in}" ]] && ! grep -Eq '^[[:space:]]*CXXFLAGS.*-fPIC' "${cosmorec_makefile_in}"; then
    sed -i -E '/^[[:space:]]*CXXFLAGS[[:space:]]*=/ s/$/ -fPIC/' "${cosmorec_makefile_in}"
    patched_cosmorec_flags=1
fi

if [[ -f "${cosmorec_makefile_in}" ]]; then
    if grep -Eq '^[[:space:]]*CXXFLAGSLOC[[:space:]]*=.*RECFASTPPPATH' "${cosmorec_makefile_in}"; then
        sed -i -E '/^[[:space:]]*CXXFLAGSLOC[[:space:]]*=/ s|[[:space:]]+-D RECFASTPPPATH=.*$||' "${cosmorec_makefile_in}"
        patched_cosmorec_flags=1
    fi

    if ! grep -Eq '^[[:space:]]+-D RECFASTPPPATH=' "${cosmorec_makefile_in}"; then
        sed -i '/-D COSMORECPATH=/i\              -D RECFASTPPPATH=\\"$(PWD)/$(DEV_DIR)/Cosmology/Recfast++/\\" \\' "${cosmorec_makefile_in}"
        patched_cosmorec_flags=1
    fi
fi

cosmorec_makefile="${cosmorec_dir}/Makefile"
if [[ -f "${cosmorec_makefile}" ]] && ! grep -Eq '^[[:space:]]*CXXFLAGS.*-fPIC' "${cosmorec_makefile}"; then
    printf '\nCXXFLAGS += -fPIC\n' >> "${cosmorec_makefile}"
    patched_cosmorec_flags=1
fi

if [[ "${patched_cosmorec_flags}" == "1" ]]; then
    make -C "${cosmorec_dir}" tidy >/dev/null 2>&1 || true
fi

if [[ ! -d "${hyrec_dir}/.git" ]]; then
    rm -rf "${hyrec_dir}"
    safe_git clone --depth 1 --filter=blob:none "${hyrec_repo}" "${hyrec_dir}"
fi

chmod -R u+rwX "${hyrec_dir}" >/dev/null 2>&1 || true

if [[ ! -f "${hyrec_dir}/Makefile" ]]; then
    echo "HYREC-2 checkout at ${hyrec_dir} does not contain a Makefile" >&2
    exit 1
fi
