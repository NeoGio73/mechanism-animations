"""Switch the running after-path single-point lanes from low priority to normal priority on all 8 cores
(J. Picione, 2026-11-16 evening), without stopping any ORCA job.

A running lane has the old queue code loaded, so each lane is restarted with pipeline/run_singlepoints.sh. First a
gate goes into the core ledger: a waiting entry on a folder held by this helper (0 cores). Low-priority lanes start a
job only when nothing waits, so the gate stops every old lane from starting a new job; normal lanes skip it (its
folder is busy). A lane that holds no cores is stopped and restarted at once; a lane running a job is stopped right
after that job ends and releases its cores, then restarted. Finished jobs are reused (pipeline.orca.run). When every
lane is switched the gate is removed. Log: docs/level_test_queue/allocation.log and docs/sp_switch.log.
Usage (started hidden): python -X utf8 -m pipeline.sp_switch <reaction> [<reaction> ...]
"""
from __future__ import annotations

import base64
import os
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
GATE = "GATE sp_switch"


def log(msg: str) -> None:
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")


def gate(on: bool) -> None:
    pid = os.getpid()
    with cores._locked():
        s = cores._read()
        s["held"] = [h for h in s["held"] if h.get("folder") != GATE]
        s["waiting"] = [w for w in s["waiting"] if w.get("folder") != GATE]
        if on:
            e = {"pid": pid, "lane": "sp_switch (gate)", "what": "gate: no new low-priority start", "cores": 0,
                 "mem_mb": 0, "folder": GATE, "since": time.strftime("%H:%M:%S"), "beat": time.time(), "front": False}
            s["held"].append({**e, "id": f"{pid}-gate-h"})
            s["waiting"].append({**e, "id": f"{pid}-{time.time_ns()}"})
        cores._write(s)
    cores._log(f"{'gate on ' if on else 'gate off'} sp_switch (pid {pid}): low-priority lanes start no new job")


def old_lane_procs(rx: str) -> list[psutil.Process]:
    out = []
    for p in psutil.process_iter(["pid", "cmdline"]):
        cl = " ".join(p.info["cmdline"] or [])
        if (f"run_after_path_low.sh' {rx}" in cl or f"run_after_path_low.sh {rx}" in cl
                or (f"-m pipeline reactions/{rx} from_stationary" in cl and "CORES_FIT" not in cl)):
            out.append(p)
    return out


def old_holds(rx: str) -> list[dict]:
    s = cores._read()
    return [h for h in s["held"] if h.get("low") and f"reactions/{rx} " in h.get("lane", "")]


def stop_old(rx: str) -> None:
    procs = old_lane_procs(rx)
    for p in procs:
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)], capture_output=True)
    log(f"{rx}: old low-priority lane stopped between jobs (pids {[p.pid for p in procs]})")


def start_new(rx: str) -> None:
    cmd = (f'"{BASH}" -c "export PATH=/usr/bin:/bin:$PATH; bash \'/c/Users/giova/Mechanism Animations/pipeline/'
           f'run_singlepoints.sh\' {rx}"')
    b64 = base64.b64encode(cmd.encode("utf-8")).decode()
    subprocess.run([PYW, str(ROOT / "pipeline/launch_hidden.py"), b64, str(ROOT)])
    log(f"{rx}: new lane started (normal priority, all 8 cores; pipeline/run_singlepoints.sh)")


def main(rxs: list[str]) -> None:
    log(f"switch start: {rxs}")
    gate(True)
    todo = list(rxs)
    beat = time.time()
    while todo:
        for rx in list(todo):
            if not old_holds(rx):
                time.sleep(1)                       # let a just-released lane reach the gate
                if old_holds(rx):
                    continue
                stop_old(rx)
                start_new(rx)
                todo.remove(rx)
        if time.time() - beat > 60:                 # keep the gate's waiting entry fresh (cores.LOST_S)
            with cores._locked():
                s = cores._read()
                for w in s["waiting"]:
                    if w.get("folder") == GATE:
                        w["beat"] = time.time()
                cores._write(s)
            beat = time.time()
        time.sleep(1)
    gate(False)
    log("switch done: every lane on normal priority")


if __name__ == "__main__":
    main(sys.argv[1:])
