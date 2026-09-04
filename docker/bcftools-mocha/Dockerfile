FROM eclipse-temurin:25-jdk

# ---------------------------------------------------------------
# Versions
# ---------------------------------------------------------------
ENV BCFTOOLS_VERSION=1.20
ENV MOCHA_COMMIT=95686b7b65f53a490513be76bb120e5fc20a8bcf

ENV DEBIAN_FRONTEND=noninteractive
ENV BCFTOOLS_PLUGINS=/opt/bcftools/plugins
ENV PATH=/opt/bcftools/bin:$PATH

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# ---------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        make \
        wget \
        ca-certificates \
        git \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libgsl-dev \
        libncurses-dev \
        libdeflate-dev; \
    rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------
# Download BCFtools, fetch MoChA at an immutable commit, build, install
# ---------------------------------------------------------------
RUN set -eux; \
    cd /tmp; \
    wget -q "http://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"; \
    tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2"; \
    \
    git clone https://github.com/freeseek/mocha.git /tmp/mocha; \
    cd /tmp/mocha; \
    git checkout --detach "${MOCHA_COMMIT}"; \
    cp \
        mocha.c \
        mocha.h \
        mochatools.c \
        extendFMT.c \
        beta_binom.h \
        genome_rules.h \
        "/tmp/bcftools-${BCFTOOLS_VERSION}/plugins/"; \
    \
    cd "/tmp/bcftools-${BCFTOOLS_VERSION}"; \
    make -j"$(nproc)"; \
    \
    # bcftools' release tarball bundles htslib as a subdirectory; build it \
    # separately to get tabix (and bgzip/htsfile as a bonus). \
    cd htslib-${BCFTOOLS_VERSION}; \
    make -j"$(nproc)" tabix bgzip htsfile; \
    cd ..; \
    \
    mkdir -p /opt/bcftools/bin /opt/bcftools/plugins; \
    cp bcftools /opt/bcftools/bin/; \
    cp htslib-${BCFTOOLS_VERSION}/tabix /opt/bcftools/bin/; \
    cp htslib-${BCFTOOLS_VERSION}/bgzip /opt/bcftools/bin/; \
    cp htslib-${BCFTOOLS_VERSION}/htsfile /opt/bcftools/bin/; \
    cp plugins/*.so /opt/bcftools/plugins/; \
    \
    printf 'BCFTOOLS_VERSION=%s\nMOCHA_COMMIT=%s\n' \
        "${BCFTOOLS_VERSION}" "${MOCHA_COMMIT}" \
        > /opt/bcftools/BUILD_INFO; \
    /opt/bcftools/bin/bcftools --version > /opt/bcftools/BCFTOOLS_VERSION; \
    \
    cd /; \
    rm -rf \
        /tmp/mocha \
        "/tmp/bcftools-${BCFTOOLS_VERSION}" \
        "/tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2"



RUN set -eux; echo "=== BCFtools ==="; bcftools --version

RUN set -eux; echo "=== tabix ==="; tabix --version

RUN set -eux; echo "=== Java ==="; java --version

RUN set -eux; echo "=== Installed plugins ==="; bcftools plugin -l

RUN set -eux; \
    echo "=== Checking MoChA plugins ==="; \
    test -f "${BCFTOOLS_PLUGINS}/mocha.so"; \
    test -f "${BCFTOOLS_PLUGINS}/mochatools.so"; \
    test -f "${BCFTOOLS_PLUGINS}/extendFMT.so"

RUN set -eux; \
    bcftools plugin -l | grep -qx "mocha"; \
    bcftools plugin -l | grep -qx "mochatools"; \
    bcftools plugin -l | grep -qx "extendFMT"

RUN set -eux; \
    echo "=== MoChA command ==="; \
    (bcftools +mocha 2>&1 || true ) | grep -q "Usage:.*bcftools +mocha"

RUN set -eux; \
    echo "=== Build information ==="; \
    cat /opt/bcftools/BUILD_INFO

RUN echo "MoChA installation test passed."