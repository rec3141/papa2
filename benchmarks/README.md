# Benchmarks

This directory contains a reproducible benchmark harness for comparing:
- the Python package in this repo (`papa2`)
- the R `dada2` package or a local R package source tree

For the current benchmark writeup and real-data results, see
[BENCHMARKS.md](./BENCHMARKS.md).

## Scope

The harness is designed to sweep:
- total input size across 5 orders of magnitude
- sample count across 3 orders of magnitude

The default synthetic matrix uses:
- total bases: `1e4`, `1e5`, `1e6`, `1e7`, `1e8`, `1e9`
- sample counts: `1`, `10`, `100`, `1000`

Each case is generated deterministically from the bundled FASTQ seed files in
`tests/data/`.

## Outputs

Each benchmark run writes:
- `stdout.log`
- `stderr.log`
- `/usr/bin/time -v` output
- sampled telemetry (`telemetry.csv`)
- parsed summary metrics (`summary.json`)

Suite runs also write a collated `results.csv`.

## Main entrypoints

- `benchmarks/run_suite.py`
- `benchmarks/plot_results.py`

## Example

```bash
python benchmarks/run_suite.py \
  --workspace .benchmarks \
  --engines python r \
  --matrix small \
  --r-repo /data/dada2_gpu
```

Then generate plots:

```bash
python benchmarks/plot_results.py \
  --results .benchmarks/results.csv \
  --out-dir .benchmarks/plots
```
