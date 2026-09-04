FROM eclipse-temurin:25-jdk-noble

# ---------------------------------------------------------------
# Versions
# ---------------------------------------------------------------
ENV BCFTOOLS_VERSION=1.20
ENV MOCHA_COMMIT=95686b7b65f53a490513be76bb120e5fc20a8bcf

ENV DEBIAN_FRONTEND=noninteractive

# BCFtools installs plugins here with:
#   make install prefix=/opt/bcftools
ENV BCFTOOLS_PLUGINS=/opt/bcftools/libexec/bcftools
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
# Download BCFtools, fetch MoChA at immutable commit,
# build and install
# ---------------------------------------------------------------
RUN set -eux; \
    cd /tmp; \
    \
    wget -q \
        "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"; \
    \
    tar -xjf \
        "bcftools-${BCFTOOLS_VERSION}.tar.bz2"; \
    \
    git clone \
        https://github.com/freeseek/mocha.git \
        /tmp/mocha; \
    \
    cd /tmp/mocha; \
    git checkout --detach "${MOCHA_COMMIT}"; \
    \
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
    \
    # Build BCFtools and all plugins using the same installation prefix
    make -j"$(nproc)" prefix=/opt/bcftools; \
    make install prefix=/opt/bcftools; \
    \
    # -----------------------------------------------------------
    # Build HTSlib utilities bundled with the BCFtools release.
    # -----------------------------------------------------------
    cd "htslib-${BCFTOOLS_VERSION}"; \
    make -j"$(nproc)" tabix bgzip htsfile; \
    cp tabix bgzip htsfile /opt/bcftools/bin/; \
    \
    # -----------------------------------------------------------
    # Record build information
    # -----------------------------------------------------------
    printf 'BCFTOOLS_VERSION=%s\nMOCHA_COMMIT=%s\n' \
        "${BCFTOOLS_VERSION}" \
        "${MOCHA_COMMIT}" \
        > /opt/bcftools/BUILD_INFO; \
    \
    /opt/bcftools/bin/bcftools --version \
        > /opt/bcftools/BCFTOOLS_VERSION; \
    \
    # -----------------------------------------------------------
    # Cleanup
    # -----------------------------------------------------------
    cd /; \
    rm -rf \
        /tmp/mocha \
        "/tmp/bcftools-${BCFTOOLS_VERSION}" \
        "/tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# ---------------------------------------------------------------
# Basic executable tests
# ---------------------------------------------------------------
RUN set -eux; \
    echo "=== BCFtools ==="; \
    bcftools --version; \
    \
    echo "=== tabix ==="; \
    tabix --version; \
    \
    echo "=== Java ==="; \
    java --version

# ---------------------------------------------------------------
# Check all installed plugins
# ---------------------------------------------------------------
RUN set -eux; \
    echo "=== BCFTOOLS_PLUGINS ==="; \
    echo "${BCFTOOLS_PLUGINS}"; \
    \
    echo "=== Installed plugins ==="; \
    ls -la "${BCFTOOLS_PLUGINS}"; \
    \
    bcftools plugin -l

# ---------------------------------------------------------------
# Check that required MoChA / BCFtools plugins exist
# ---------------------------------------------------------------
RUN set -eux; \
    echo "=== Checking required plugins ==="; \
    \
    test -f "${BCFTOOLS_PLUGINS}/mocha.so"; \
    test -f "${BCFTOOLS_PLUGINS}/mochatools.so"; \
    test -f "${BCFTOOLS_PLUGINS}/extendFMT.so"; \
    test -f "${BCFTOOLS_PLUGINS}/fixploidy.so"; \
    \
    echo "Required plugin files found."

# ---------------------------------------------------------------
# Check that BCFtools can discover the plugins
# ---------------------------------------------------------------
RUN set -eux; \
    echo "=== Checking plugin discovery ==="; \
    \
    bcftools plugin -l | grep -qx "mocha"; \
    bcftools plugin -l | grep -qx "mochatools"; \
    bcftools plugin -l | grep -qx "extendFMT"; \
    bcftools plugin -l | grep -qx "fixploidy"; \
    \
    echo "Required plugins successfully discovered by BCFtools."

# ---------------------------------------------------------------
# Build information
# ---------------------------------------------------------------
RUN set -eux; \
    echo "=== Build information ==="; \
    cat /opt/bcftools/BUILD_INFO

# ---------------------------------------------------------------
# Final confirmation
# ---------------------------------------------------------------
RUN echo "MoChA installation test passed."