from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from matrix import build_matrix


def run_checked(cmd: list[str], cwd: Path) -> None:
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd), check=True)


def ensure_python_build(root: Path) -> None:
    if (root / "libpapa2.so").exists():
        return
    run_checked(["make", "libpapa2.so"], cwd=root)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--matrix", default="small")
    parser.add_argument("--engines", nargs="+", default=["python", "r"])
    parser.add_argument("--r-repo", default="")
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--rscript-bin", default="Rscript")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    bench_dir = root / "benchmarks"
    workspace = Path(args.workspace)
    workspace.mkdir(parents=True, exist_ok=True)

    if "python" in args.engines:
        ensure_python_build(root)

    cases = build_matrix(args.matrix)
    rows: list[dict] = []

    for case in cases:
        case_dir = workspace / case.case_id
        data_dir = case_dir / "data"
        run_checked(
            [
                args.python_bin,
                str(bench_dir / "generate_fastq_case.py"),
                "--matrix",
                args.matrix,
                "--case-id",
                case.case_id,
                "--out-dir",
                str(data_dir),
            ],
            cwd=root,
        )
        manifest_path = data_dir / "manifest.json"

        for engine in args.engines:
            engine_dir = case_dir / engine
            if engine == "python":
                command = [
                    args.python_bin,
                    str(bench_dir / "run_python_case.py"),
                    "--manifest",
                    str(manifest_path),
                ]
            elif engine == "r":
                command = [
                    args.rscript_bin,
                    str(bench_dir / "run_r_case.R"),
                    "--manifest",
                    str(manifest_path),
                ]
                if args.r_repo:
                    command.extend(["--r-repo", args.r_repo])
            else:
                raise ValueError(f"Unknown engine: {engine}")

            run_checked(
                [
                    args.python_bin,
                    str(bench_dir / "monitor_command.py"),
                    "--label",
                    f"{case.case_id}:{engine}",
                    "--out-dir",
                    str(engine_dir),
                    "--cwd",
                    str(root),
                    "--",
                    *command,
                ],
                cwd=root,
            )

            summary = json.loads((engine_dir / "summary.json").read_text())
            summary.update(case.to_dict())
            summary["engine"] = engine
            rows.append(summary)

    results_path = workspace / "results.csv"
    if rows:
        fieldnames = sorted({k for row in rows for k in row.keys()})
        with open(results_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    print(results_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
