"""Point 3b (J. Picione, 2026-11-11): reproduce the step-test crashes on purpose.

Two short IRC jobs (SN2 B3LYP-D3(BJ)/ma-def2-SVP SMD water TS, forward and backward, the 0.10 settings
of the step test, MaxIter 60, 2 cores each) are run
  shared   both in one scratch folder at the same time (twice: rounds 1 and 2)
  separate each in its own folder at the same time (once)
A watcher records every change of redirect.tmp in each folder (time, size, first line that differs)
and copies the per-evaluation property files every 5 s, so they survive normal termination. The D3
values every evaluation used are then checked against D3 recomputed alone (as in check 1, part A).

The core ledger holds 4 cores, 12000 MB for each round; ORCA is started here directly, because the
queue itself never puts two ORCA jobs in one folder.
Usage: python -m pipeline.redirect_repro run
Results: reactions/sn2_hydroxide_bromomethane_masvp/02_irc/attempts/redirect_repro/results.json
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

from . import cores, orca
from .geometry import read_xyz
from .steps import Reaction, irc_init_scale

ROOT = Path(__file__).resolve().parents[1]
RXD = ROOT / "reactions/sn2_hydroxide_bromomethane_masvp"
BASE = RXD / "02_irc/attempts/redirect_repro"
SD, STEP, MAXITER = 0.02359, 0.10, 60


def inp(rx, d: str) -> str:
    sc = irc_init_scale(RXD / "01_ts/ts/ts.hess", STEP)
    blk = (f'%irc\n  MaxIter {MAXITER}\n  Direction {d}\n  InitHess read\n  Hess_Filename "ts.hess"\n'
           f'  Scale_Displ_SD {SD:.5f}\n  Adapt_Scale_Displ false\n  Init_Displ length\n'
           f'  Scale_Init_Displ {sc:.5f}\n  PrintLevel 1\nend')
    return orca.make_input(rx.keywords("IRC"), read_xyz(RXD / "01_ts/ts/ts.xyz"), rx.cfg["charge"],
                           rx.cfg["mult"], blocks=rx.solvent_block + "\n" + blk, nprocs=2,
                           maxcore_mb=rx.cfg["maxcore_mb"])


def watch(folders: list[Path], stop: threading.Event, log: list):
    last = {}
    while not stop.is_set():
        for f in folders:
            r = f / "redirect.tmp"
            try:
                st = r.stat()
                key = (st.st_mtime_ns, st.st_size)
                if last.get(r) != key:
                    txt = r.read_text(errors="replace")
                    c6 = next((ln.strip() for ln in txt.splitlines() if "C6(AA)" in ln), "")
                    ed = next((ln.strip() for ln in txt.splitlines() if "Edisp" in ln), "")
                    log.append({"t": time.strftime("%H:%M:%S"), "folder": f.name, "size": st.st_size,
                                "md5": hashlib.md5(txt.encode()).hexdigest()[:8], "c6": c6, "edisp": ed})
                    last[r] = key
            except (FileNotFoundError, PermissionError):
                pass
            for p in f.glob("irc_*_IRC_*.property.txt"):
                try:
                    shutil.copy(p, f / ("kept_" + p.name))
                except OSError:
                    pass
        time.sleep(0.5)


def round_(rx, name: str, shared: bool) -> dict:
    wd = BASE / name
    if wd.exists():
        shutil.rmtree(wd)
    dirs = {d: (wd if shared else wd / d) for d in ("forward", "backward")}
    for d, p in dirs.items():
        p.mkdir(parents=True, exist_ok=True)
        shutil.copy(RXD / "01_ts/ts/ts.hess", p / "ts.hess")
        (p / f"irc_{d[0]}.inp").write_text(inp(rx, d))
    rid = cores.acquire(4, 12000, f"redirect repro {name} ({'one folder' if shared else 'two folders'})",
                        lane="redirect_repro")
    env = os.environ.copy()
    env["PATH"] = f"{orca.MSMPI_BIN};{orca.ORCA_DIR};" + env.get("PATH", "")
    stop, wlog = threading.Event(), []
    th = threading.Thread(target=watch, args=(sorted(set(dirs.values())), stop, wlog), daemon=True)
    th.start()
    t0 = time.time()
    procs = {}
    try:
        for d, p in dirs.items():
            fh = open(p / f"irc_{d[0]}.out", "w")
            procs[d] = (subprocess.Popen([str(orca.ORCA_DIR / "orca.exe"), f"irc_{d[0]}.inp"], cwd=p,
                                         stdout=fh, stderr=subprocess.STDOUT, env=env), fh, time.time())
        ends = {}
        while len(ends) < len(procs):
            for d, (pr, fh, ts) in procs.items():
                if d not in ends and pr.poll() is not None:
                    ends[d] = {"exit_code": pr.returncode, "wall_s": round(time.time() - ts),
                               "ended": time.strftime("%H:%M:%S")}
                    fh.close()
            time.sleep(1)
    finally:
        stop.set()
        th.join(5)
        cores.release(rid)
    res = {"folders": "one" if shared else "two", "started": time.strftime("%H:%M:%S", time.localtime(t0))}
    for d, p in dirs.items():
        out = (p / f"irc_{d[0]}.out").read_text(errors="replace")
        n = sum(1 for ln in out.splitlines() if ln.strip()[:1].isdigit() and "-2689." in ln)
        res[d] = {**ends[d], "terminated_normally": "ORCA TERMINATED NORMALLY" in out,
                  "irc_iterations": n, "error_lines": [ln.strip() for ln in out.splitlines()
                                                         if "error" in ln.lower() or "abort" in ln.lower()][-5:]}
    res["redirect_writes"] = wlog
    res["n_redirect_writes"] = len(wlog)
    return res


def d3_check(rx, res_all: dict) -> dict:
    """D3 every evaluation used (kept property files) vs D3 recomputed alone."""
    from . import shared_folder_check as sfc
    sfc.CK = BASE / "d3check"
    sfc.CK.mkdir(parents=True, exist_ok=True)
    sfc.LOG = orca.RunLog(sfc.CK / "run_log.md")
    out = {}
    for p in sorted(BASE.glob("*/**/kept_irc_*_IRC_*.property.txt")) + sorted(BASE.glob("*/kept_irc_*_IRC_*.property.txt")):
        tag = f"{p.parent.relative_to(BASE)}_{p.stem}".replace("\\", "_").replace("/", "_")
        if tag in out:
            continue
        ev = sfc.read_property(p)
        d3 = sfc.d3_batch(tag, ev)
        dd = [abs(a - e["vdw"]) for a, e in zip(d3, ev)]
        out[tag] = {"n_evaluations": len(ev), "max_abs_Eh": max(dd), "n_over_1e-9": sum(x > 1e-9 for x in dd)}
        print(tag, out[tag], flush=True)
    return out


if __name__ == "__main__":
    rx = Reaction.load(RXD)
    BASE.mkdir(parents=True, exist_ok=True)
    R = {}
    for name, shared in (("shared_round1", True), ("separate", False), ("shared_round2", True)):
        R[name] = round_(rx, name, shared)
        print(name, json.dumps({k: v for k, v in R[name].items() if k != "redirect_writes"}), flush=True)
        (BASE / "results.json").write_text(json.dumps(R, indent=1), encoding="utf-8")
    R["d3_check"] = d3_check(rx, R)
    (BASE / "results.json").write_text(json.dumps(R, indent=1), encoding="utf-8")
