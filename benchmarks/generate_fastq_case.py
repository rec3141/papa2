from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from matrix import BenchmarkCase, build_matrix


def read_fastq_records(path: Path) -> list[tuple[str, str, str, str]]:
    records: list[tuple[str, str, str, str]] = []
    with gzip.open(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            records.append((header, seq, plus, qual))
    if not records:
        raise ValueError(f"No FASTQ records found in {path}")
    return records


def write_case(case: BenchmarkCase, seeds: list[Path], out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = case.to_dict()
    files: list[str] = []

    seed_records = [read_fastq_records(seed) for seed in seeds]

    for sample_idx in range(case.sample_count):
        sample_path = out_dir / f"{case.case_id}_sample_{sample_idx:04d}.fastq.gz"
        records = seed_records[sample_idx % len(seed_records)]
        target_bases = case.per_sample_bases
        bases_written = 0
        rec_idx = 0

        with gzip.open(sample_path, "wt") as out:
            while bases_written < target_bases:
                header, seq, plus, qual = records[rec_idx % len(records)]
                tag = f"{sample_idx}_{rec_idx}"
                out.write(f"{header.rstrip()}_{tag}\n")
                out.write(seq)
                out.write(plus)
                out.write(qual)
                bases_written += len(seq.strip())
                rec_idx += 1
        files.append(str(sample_path))

    manifest["files"] = files
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default="smoke")
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--seed",
        action="append",
        default=[],
        help="FASTQ seed file. Defaults to tests/data/sam1F.fastq.gz and sam2F.fastq.gz.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    seeds = [Path(s) for s in args.seed] or [
        root / "tests" / "data" / "sam1F.fastq.gz",
        root / "tests" / "data" / "sam2F.fastq.gz",
    ]

    case = next(c for c in build_matrix(args.matrix) if c.case_id == args.case_id)
    manifest = write_case(case, seeds, Path(args.out_dir))
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
