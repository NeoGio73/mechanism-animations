"""ORCA processes that are not this project's (J. Picione, 2026-11-21).

His own projects may run ORCA on this machine during his machine block (2026-11-21: single points from a folder
in Downloads). Those jobs are not ours: never stop, rerun or count them, and never treat them as a stall. Before
the block is released (docs/level_test_queue/window_release), and so before the lanes restart, this check must
show none of them still running; if one is, wait for it to finish.

Usage:  python -X utf8 -m pipeline.foreign_orca            prints them; exit code 1 if any is running
A process counts as foreign when it is an ORCA/MPI process (orca*, mpiexec, smpd) whose working folder or
command line is not inside this repository. Read-only: it never touches a process.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[1]
PREFIXES = ("orca", "mpiexec", "smpd")


def foreign() -> list[dict]:
    root = str(ROOT).lower()
    out = []
    for p in psutil.process_iter(["pid", "name", "create_time"]):
        name = (p.info["name"] or "").lower()
        if not name.startswith(PREFIXES):
            continue
        try:
            cwd = p.cwd()
        except psutil.Error:
            cwd = ""
        try:
            cmd = " ".join(p.cmdline())
        except psutil.Error:
            cmd = ""
        if root in cwd.lower() or root in cmd.lower() or root.replace("\\", "/") in cmd.lower():
            continue
        out.append({"pid": p.info["pid"], "name": p.info["name"], "cwd": cwd, "cmd": cmd[:200],
                    "started": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.info["create_time"] or 0))})
    return out


def main() -> None:
    f = foreign()
    if not f:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} no ORCA process from outside this project is running")
        sys.exit(0)
    for r in f:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} running: {r['name']} pid {r['pid']} since {r['started']}, "
              f"folder {r['cwd'] or '(not readable)'}")
    sys.exit(1)


if __name__ == "__main__":
    main()
