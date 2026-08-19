#!/usr/bin/env python3
"""Run a command, tracking wall time, peak RSS summed over the process
tree (sampled at 4 Hz), and total CPU time (user+sys, incl. reaped
children via wait4 rusage).  Writes a JSON line to the file in $RESMON_OUT.
"""
import json, os, subprocess, sys, threading, time
import psutil

cmd = sys.argv[1:]
out_path = os.environ.get("RESMON_OUT", "resmon.json")

t0 = time.monotonic()
proc = subprocess.Popen(cmd)
peak = 0
stop = False

def sampler():
    global peak
    try:
        p = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        return
    while not stop:
        total = 0
        try:
            procs = [p] + p.children(recursive=True)
            for q in procs:
                try:
                    total += q.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        except psutil.NoSuchProcess:
            break
        peak = max(peak, total)
        time.sleep(0.25)

th = threading.Thread(target=sampler, daemon=True)
th.start()
_, status, ru = os.wait4(proc.pid, 0)
wall = time.monotonic() - t0
stop = True
th.join(timeout=2)

rec = {
    "cmd": " ".join(cmd[:3]),
    "wall_s": round(wall, 1),
    "peak_rss_gb": round(peak / 1e9, 2),
    "cpu_user_s": round(ru.ru_utime, 1),
    "cpu_sys_s": round(ru.ru_stime, 1),
    "exit": os.waitstatus_to_exitcode(status),
}
with open(out_path, "w") as f:
    json.dump(rec, f)
print("[resmon]", json.dumps(rec), file=sys.stderr)
sys.exit(max(0, rec["exit"]))
