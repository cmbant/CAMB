ARG VARIANT="3.14-bookworm"
FROM mcr.microsoft.com/devcontainers/python:${VARIANT}

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
    tar \
    unzip \
    xz-utils \
    && curl -LsSf https://astral.sh/uv/install.sh -o /tmp/install-uv.sh \
    && env UV_INSTALL_DIR=/usr/local/bin sh /tmp/install-uv.sh \
    && rm -f /tmp/install-uv.sh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
