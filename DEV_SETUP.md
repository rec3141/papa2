# Dev Setup

## Conda Environment

Create the development environment from the repo root:

```bash
mamba env create -f environment.yml
conda activate papa2-dev
```

## CPU-Only Build

Build the standalone shared library:

```bash
make clean libpapa2.so
```

## Quick Checks

Verify the Python package can load:

```bash
python - <<'PY'
import papa2
print("papa2 import ok:", papa2.__version__)
PY
```

Run Python tests:

```bash
pytest
```

## Notes

- This branch is CPU-only.
