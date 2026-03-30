# papa2

Python-first amplicon denoising derived from the standalone CPU-only port of DADA2.

This repo contains:
- a pure Python package in `papa2/`
- a standalone C/C++ core built as `libdada2.so`
- pytest-based regression coverage for the Python wrapper and denoising parity checks

This repo does not include the original R package structure.

## Status

The current branch is CPU-only. The previous hybrid GPU path was profiled and shelved because its exact mode was dominated by CPU fallback work.

## Development Setup

```bash
mamba env create -f environment.yml
conda activate dada2-dev
make clean libdada2.so
pytest
```

## Quick Smoke Test

```bash
python - <<'PY'
import papa2
print(papa2.__version__)
PY
```

## Project Layout

- `papa2/`: Python package
- `src/`: standalone native core used by the Python bindings
- `tests/`: pytest suite, including bundled small FASTQ fixtures

## Provenance

The denoising algorithm and much of the native core logic originate from DADA2:
- https://github.com/benjjneb/dada2
