# papa2

Python-first amplicon denoising derived from the standalone CPU-only port of DADA2.

This repo contains:
- a pure Python package in `papa2/`
- a standalone C/C++ core built as `libdada2.so`
- pytest-based regression coverage for the Python wrapper and denoising parity checks

This repo does not include the original R package structure.

## Status

The current project is CPU-only.

## Provenance

`papa2` was split out of the development work that happened in:
- fork: [rec3141/dada2](https://github.com/rec3141/dada2)
- branch: [gpu-python](https://github.com/rec3141/dada2/tree/gpu-python)

That line of work branched from upstream DADA2 at:
- upstream repo: [benjjneb/dada2](https://github.com/benjjneb/dada2)
- branch-point commit: [72da7700b58290e40cdce4b0856314aecf2b9dc4](https://github.com/benjjneb/dada2/commit/72da7700b58290e40cdce4b0856314aecf2b9dc4)
- commit date: `2026-02-13`

The initial `papa2` import is a Python-only extraction of that line of work.

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

## Upstream

The denoising algorithm and much of the native core logic originate from DADA2:
- [benjjneb/dada2](https://github.com/benjjneb/dada2)
