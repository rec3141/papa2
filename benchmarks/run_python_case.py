from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import papa2


def _normalize_results(results):
    if isinstance(results, dict):
        return [results]
    return list(results)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    files = manifest["files"]
    learn_nbases = int(manifest["learn_nbases"])

    err = papa2.learn_errors(files, nbases=learn_nbases, verbose=False)
    results = _normalize_results(papa2.dada(files, err=err, verbose=False))

    summary = {
        "engine": "python",
        "case_id": manifest["case_id"],
        "sample_count": len(files),
        "learn_nbases": learn_nbases,
        "n_asv_total": int(sum(len(r["denoised"]) for r in results)),
        "n_reads_assigned": int(sum(sum(r["denoised"].values()) for r in results)),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
