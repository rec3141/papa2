# syntax=docker/dockerfile:1

# ============================================================
# Stage 1 — build: compile libpapa2.so from source
# ============================================================
FROM ubuntu:22.04 AS build

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        make \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

COPY src/ src/
COPY Makefile .

RUN make libpapa2.so

# ============================================================
# Stage 2 — runtime: Python 3.12-slim + papa2 package
# ============================================================
FROM python:3.12-slim AS runtime

# Runtime zlib — libpapa2.so links against libz at dlopen time
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        zlib1g \
        libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# Install numpy (required by papa2 at import time)
RUN pip install --no-cache-dir numpy

# ── Locate site-packages and install the papa2 package ────────────────────
# We copy files directly rather than running `pip install .` to avoid
# triggering setup.py's `make libpapa2.so` hook (no build toolchain in this
# stage).  The package_data entry in pyproject.toml declares *.so, so a
# direct copy is equivalent to what `pip install` would do after the build.

# Determine site-packages path at build time and deposit files there.
# _cdada.py resolution order:
#   1. dirname(dirname(abspath(__file__))) + "/libpapa2.so"
#      → <site-packages>/libpapa2.so   (preferred)
#   2. "libpapa2.so"  (CWD fallback)
# We satisfy (1) by placing libpapa2.so one level above the papa2/ package dir.

RUN python3 -c "import site; print(site.getsitepackages()[0])" > /site_packages.txt

COPY papa2/ /tmp/papa2_pkg/
COPY --from=build /build/libpapa2.so /tmp/libpapa2.so

RUN SITE=$(cat /site_packages.txt) && \
    cp -r /tmp/papa2_pkg "$SITE/papa2" && \
    cp /tmp/libpapa2.so "$SITE/libpapa2.so" && \
    # Also register a minimal dist-info so `pip show papa2` works
    mkdir -p "$SITE/papa2-0.1.0.dist-info" && \
    printf 'Metadata-Version: 2.1\nName: papa2\nVersion: 0.1.0\n' \
        > "$SITE/papa2-0.1.0.dist-info/METADATA" && \
    printf 'papa2\npapa2-0.1.0.dist-info\n' \
        > "$SITE/papa2-0.1.0.dist-info/RECORD"

# Place libpapa2.so in the system library path as well so that any dynamic
# linker lookup (LD_LIBRARY_PATH / ldconfig) also resolves it
RUN cp /tmp/libpapa2.so /usr/local/lib/libpapa2.so && ldconfig

# Smoke-test: verify import and version before finalising the image
RUN python3 -c "import papa2; print('papa2 version:', papa2.__version__)"

WORKDIR /data

CMD ["python3"]
