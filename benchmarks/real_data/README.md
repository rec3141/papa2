# Real-data benchmark harness

Times every pipeline step (filter, dereplicate, learn errors, denoise,
merge, sequence table, chimera removal, SILVA taxonomy) for papa2 vs
R dada2 on real BioProjects, with replicates in randomized order and
resource tracking (peak process-tree RSS, CPU time).

## Datasets used for the v1.0 announcement

| Label | BioProject | Layout | Samples | truncLen / maxEE |
|-------|-----------|--------|---------|------------------|
| Mock | PRJEB6244 | 2x250 | 24 | (230,210) / (2,2) |
| Gut | PRJEB27564 | 2x300 | 24 | (300,260) / (2,5) |
| MOSAiC | PRJNA895866 | 2x151 | 120 | (145,145) / (2,2) |

Raw reads come from ENA (`https://ftp.sra.ebi.ac.uk/vol1/fastq/...`);
put each dataset's files in a directory as `<ACC>_1.fastq.gz` /
`<ACC>_2.fastq.gz` plus an `accs.txt` listing the accessions. The SILVA
138.1 train set is
`https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz`.

## Running

One engine, one dataset, resource-tracked:

```bash
RESMON_OUT=run.json python resmon.py \
  Rscript bench_r.R RAW_DIR OUT_DIR 230 210 2 2 silva_train.fa.gz
RESMON_OUT=run.json python resmon.py \
  python bench_py.py RAW_DIR OUT_DIR 230 210 2 2 silva_train.fa.gz
```

Each writes `TIME <step> <seconds>` lines to stdout and full outputs
(filtered FASTQs, error matrices, per-sample denoised sequences,
sequence tables, taxonomy) to OUT_DIR for parity comparison between the
engines. Run every (dataset, engine, replicate) combination in a single
shuffled sequence on an otherwise idle machine.

`make_tables.py RUNS_DIR` aggregates `<ds>_<engine>_rep<N>.log` and
`.resmon.json` files into the LaTeX table used in the announcement.
