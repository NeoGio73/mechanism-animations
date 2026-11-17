"""Restart one after-path single-point lane between jobs, so it loads new code, without stopping any ORCA job
(written 2026-11-17 for the tertiary mirror-frame change, J. Picione: the running lane has the old
pipeline/steps.py loaded).

The lane is stopped while the core ledger is locked and only at a moment when it holds no cores (both workers
waiting): with the lock held no job can start, and a lane without cores runs no ORCA job. Its waiting requests
disappear with its processes (pipeline/cores._read drops dead waiters). The new lane is started at once through
pipeline/run_singlepoints.sh with no window and reuses every finished job (pipeline.orca.run).
Stop the restart guard first (it would count the stop) and start it again afterwards.
Usage (started hidden or from a shell): python -X utf8 -m pipeline.lane_restart <reaction> [<why>]
"""
from __future__ import annotations

import base64
import subprocess
import sys
import time
from pathlib import Path

import psutil

from . import cores

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "docs/sp_switch.log"
PYW = r"C:\Users\giova\AppData\Local\Python\pythoncore-3.15-64\pythonw.exe"
BASH = r"C:\Program Files\Git\usr\bin\bash.exe"


def log(msg: str) -> None:
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")


def lane_roots(rx: str) -> list[psutil.Process]:
    """The lane's outermost processes: the `bash -c` started by launch_hidden (its tree holds everything else)."""
    out = []
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        cl = " ".join(p.info["cmdline"] or [])
        if p.info["name"] == "bash.exe" and " -c " in f" {cl} " and f"run_singlepoints.sh' {rx}" in cl:
            out.append(p)
    return out


def lane_pids(rx: str) -> list[int]:
    pids = []
    for p in psutil.process_iter(["pid", "cmdline"]):
        cl = " ".join(p.info["cmdline"] or [])
        if f"run_singlepoints.sh' {rx}" in cl or f"run_singlepoints.sh {rx}" in cl \
                or f"-m pipeline reactions/{rx} " in cl:
            pids.append(p.pid)
    return pids


def main(rx: str, why: str) -> None:
    roots = lane_roots(rx)
    if len(roots) != 1:
        raise SystemExit(f"{rx}: expected one running lane, found {len(roots)}")
    log(f"{rx}: restart between jobs requested ({why}); waiting until the lane holds no cores")
    t0 = time.time()
    while True:
        with cores._locked():
            s = cores._read()
            held = [h for h in s["held"] if f"reactions/{rx} " in h.get("lane", "")]
            if not held:
                before = lane_pids(rx)
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(roots[0].pid)], capture_output=True)
                time.sleep(1)
                left = [p for p in before if psutil.pid_exists(p)]
                for p in left:
                    subprocess.run(["taskkill", "/F", "/T", "/PID", str(p)], capture_output=True)
                break
        time.sleep(0.5)
    log(f"{rx}: old lane stopped while it held no cores (waited {time.time() - t0:.0f} s; pids {before}; "
        "no ORCA job was running in it)")
    cmd = (f'"{BASH}" -c "export PATH=/usr/bin:/bin:$PATH; bash \'/c/Users/giova/Mechanism Animations/pipeline/'
           f'run_singlepoints.sh\' {rx}"')
    subprocess.run([PYW, str(ROOT / "pipeline/launch_hidden.py"), base64.b64encode(cmd.encode()).decode(),
                    str(ROOT)])
    log(f"{rx}: new lane started (pipeline/run_singlepoints.sh, new code; finished jobs reused)")
    with (ROOT / "docs/end_check_all.log").open("a", encoding="utf-8") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} note: {rx} lane restarted between jobs by "
                f"pipeline/lane_restart.py ({why}); the next 'start' line is the new lane, no job was stopped\n")


if __name__ == "__main__":
    main(sys.argv[1], " ".join(sys.argv[2:]) or "new code")
