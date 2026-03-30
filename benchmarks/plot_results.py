from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_metric(df: pd.DataFrame, metric: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for engine, g in df.groupby("engine"):
        g = g.sort_values("total_bases")
        ax.plot(g["total_bases"], g[metric], marker="o", label=engine)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total bases")
    ax.set_ylabel(metric)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.results)

    for metric in ["wall_s", "user_s", "max_rss_kb", "sampled_peak_rss_kb"]:
        if metric in df.columns:
            plot_metric(df, metric, out_dir / f"{metric}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
