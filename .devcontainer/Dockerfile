ARG VARIANT="3.14-bookworm"
FROM mcr.microsoft.com/devcontainers/python:${VARIANT}

ENV CAMB_DEVCONTAINER_VENV=/home/vscode/.local/share/camb-devcontainer/.venv \
    PRE_COMMIT_HOME=/home/vscode/.cache/pre-commit \
    UV_CACHE_DIR=/home/vscode/.cache/uv \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/home/vscode/.local/share/camb-devcontainer/.venv \
    VIRTUAL_ENV=/home/vscode/.local/share/camb-devcontainer/.venv \
    PATH=/home/vscode/.local/share/camb-devcontainer/.venv/bin:${PATH}

RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    g++ \
    gfortran \
    git \
    libgsl-dev \
    make \
    pkg-config \
    ripgrep \
    tar \
    unzip \
    xz-utils \
    && git config --system --add safe.directory '*' \
    && curl -LsSf https://astral.sh/uv/install.sh -o /tmp/install-uv.sh \
    && env UV_INSTALL_DIR=/usr/local/bin sh /tmp/install-uv.sh \
    && rm -f /tmp/install-uv.sh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p \
    /home/vscode/.local/share/camb-devcontainer \
    /home/vscode/.codex \
    "${PRE_COMMIT_HOME}" \
    "${UV_CACHE_DIR}" \
    && chown -R vscode:vscode /home/vscode/.local /home/vscode/.cache /home/vscode/.codex \
    && su vscode -c 'uv venv /home/vscode/.local/share/camb-devcontainer/.venv --python /usr/local/bin/python3' \
    && su vscode -c 'uv pip install \
        --upgrade \
        --python /home/vscode/.local/share/camb-devcontainer/.venv/bin/python3 \
        pip \
        setuptools \
        wheel \
        findent \
        fortls \
        matplotlib \
        numpy \
        scipy \
        sympy \
        packaging \
        ruff \
        pre-commit'
