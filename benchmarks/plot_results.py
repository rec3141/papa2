from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


METRICS = [
    ("wall_s", "Wall Time (s)", "log"),
    ("monitor_wall_s", "Monitor Wall Time (s)", "log"),
    ("user_s", "User CPU Time (s)", "log"),
    ("sys_s", "System CPU Time (s)", "log"),
    ("max_rss_kb", "Max RSS (kB)", "log"),
    ("sampled_peak_rss_kb", "Sampled Peak RSS (kB)", "log"),
    ("time_cpu_percent", "CPU Utilization from /usr/bin/time (%)", "linear"),
    ("sampled_avg_cpu_percent", "Average CPU Utilization (%)", "linear"),
    ("sampled_peak_cpu_percent", "Peak CPU Utilization (%)", "linear"),
]

X_AXES = [
    ("total_bases", "Total Bases", "log"),
    ("sample_count", "Sample Count", "log"),
]


def plot_metric(df: pd.DataFrame, metric: str, metric_label: str,
                metric_scale: str, x_col: str, x_label: str,
                x_scale: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for engine, group in df.groupby("engine"):
        group = group.sort_values([x_col, "total_bases", "sample_count"])
        ax.plot(group[x_col], group[metric], marker="o", linewidth=2, label=engine)
    ax.set_xscale(x_scale)
    if metric_scale == "log":
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(metric_label)
    ax.set_title(f"{metric_label} vs {x_label}")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.results)

    for metric, metric_label, metric_scale in METRICS:
        if metric not in df.columns:
            continue
        metric_dir = out_dir / metric
        metric_dir.mkdir(parents=True, exist_ok=True)
        for x_col, x_label, x_scale in X_AXES:
            if x_col not in df.columns:
                continue
            plot_metric(
                df,
                metric,
                metric_label,
                metric_scale,
                x_col,
                x_label,
                x_scale,
                metric_dir / f"{x_col}.png",
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
