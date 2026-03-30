from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import time
from pathlib import Path

import psutil


TIME_FIELDS = {
    "User time (seconds)": "user_s",
    "System time (seconds)": "sys_s",
    "Percent of CPU this job got": "time_cpu_percent",
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "elapsed_hms",
    "Maximum resident set size (kbytes)": "max_rss_kb",
    "Exit status": "exit_status",
}


def parse_hms(text: str) -> float:
    parts = text.strip().split(":")
    if len(parts) == 2:
        mins, secs = parts
        return int(mins) * 60 + float(secs)
    if len(parts) == 3:
        hours, mins, secs = parts
        return int(hours) * 3600 + int(mins) * 60 + float(secs)
    return float(text)


def parse_time_file(path: Path) -> dict:
    metrics: dict[str, object] = {}
    for line in path.read_text().splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        out_key = TIME_FIELDS.get(key)
        if out_key is None:
            continue
        metrics[out_key] = value
    if "elapsed_hms" in metrics:
        metrics["wall_s"] = parse_hms(str(metrics["elapsed_hms"]))
    for key in ("user_s", "sys_s"):
        if key in metrics:
            metrics[key] = float(metrics[key])
    if "max_rss_kb" in metrics:
        metrics["max_rss_kb"] = int(metrics["max_rss_kb"])
    if "time_cpu_percent" in metrics:
        metrics["time_cpu_percent"] = float(str(metrics["time_cpu_percent"]).rstrip("%"))
    if "exit_status" in metrics:
        metrics["exit_status"] = int(metrics["exit_status"])
    return metrics


def process_tree(proc: psutil.Process) -> list[psutil.Process]:
    try:
        return [proc] + proc.children(recursive=True)
    except psutil.Error:
        return [proc]


def sample_tree(proc: psutil.Process) -> tuple[float, int]:
    cpu_total = 0.0
    rss_total = 0
    for child in process_tree(proc):
        try:
            times = child.cpu_times()
            cpu_total += times.user + times.system
            rss_total += child.memory_info().rss
        except psutil.Error:
            continue
    return cpu_total, rss_total


def append_sample(writer, now: float, start: float, cpu_now: float, rss_now: int,
                  last_wall: float, last_cpu: float) -> tuple[float, float]:
    dt_wall = max(now - last_wall, 1e-9)
    dt_cpu = max(cpu_now - last_cpu, 0.0)
    writer.writerow(
        {
            "t_s": now - start,
            "rss_bytes": rss_now,
            "cpu_s": cpu_now,
            "cpu_pct": 100.0 * dt_cpu / dt_wall,
        }
    )
    return now, cpu_now


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--cwd", default=".")
    parser.add_argument("--sample-interval", type=float, default=0.2)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        raise SystemExit("No command provided")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = out_dir / "stdout.log"
    stderr_path = out_dir / "stderr.log"
    time_path = out_dir / "time.txt"
    telemetry_path = out_dir / "telemetry.csv"

    wrapped = ["/usr/bin/time", "-v", "-o", str(time_path), *args.command]
    start = time.time()
    with open(stdout_path, "w") as stdout, open(stderr_path, "w") as stderr:
        proc = subprocess.Popen(wrapped, cwd=args.cwd, stdout=stdout, stderr=stderr)

    ps_proc = psutil.Process(proc.pid)
    last_wall = time.time()
    last_cpu, last_rss = sample_tree(ps_proc)
    peak_sampled_rss = last_rss

    with open(telemetry_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["t_s", "rss_bytes", "cpu_s", "cpu_pct"])
        writer.writeheader()

        while proc.poll() is None:
            time.sleep(args.sample_interval)
            now = time.time()
            cpu_now, rss_now = sample_tree(ps_proc)
            peak_sampled_rss = max(peak_sampled_rss, rss_now)
            last_wall, last_cpu = append_sample(
                writer, now, start, cpu_now, rss_now, last_wall, last_cpu
            )
            last_rss = rss_now

        now = time.time()
        cpu_now, rss_now = sample_tree(ps_proc)
        peak_sampled_rss = max(peak_sampled_rss, rss_now)
        last_wall, last_cpu = append_sample(
            writer, now, start, cpu_now, rss_now, last_wall, last_cpu
        )
        last_rss = rss_now

    end = time.time()
    time_metrics = parse_time_file(time_path)
    telemetry_rows = telemetry_path.read_text().splitlines()[1:]
    cpu_pcts = []
    for row in telemetry_rows:
        parts = row.split(",")
        if len(parts) == 4:
            cpu_pcts.append(float(parts[3]))

    fallback_cpu = False
    sampled_avg_cpu = (sum(cpu_pcts) / len(cpu_pcts)) if cpu_pcts else 0.0
    sampled_peak_cpu = max(cpu_pcts) if cpu_pcts else 0.0
    if sampled_peak_cpu == 0.0 and "time_cpu_percent" in time_metrics:
        sampled_avg_cpu = float(time_metrics["time_cpu_percent"])
        sampled_peak_cpu = float(time_metrics["time_cpu_percent"])
        fallback_cpu = True

    summary = {
        "label": args.label,
        "command": args.command,
        "cwd": os.path.abspath(args.cwd),
        "sample_interval_s": args.sample_interval,
        "sampled_peak_rss_kb": int(peak_sampled_rss / 1024),
        "sampled_avg_cpu_percent": sampled_avg_cpu,
        "sampled_peak_cpu_percent": sampled_peak_cpu,
        "sampled_cpu_fallback": fallback_cpu,
        "monitor_wall_s": end - start,
    }
    summary.update(time_metrics)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return int(summary.get("exit_status", 0))


if __name__ == "__main__":
    raise SystemExit(main())
