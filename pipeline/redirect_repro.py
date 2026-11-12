"""Point 3b (J. Picione, 2026-11-11): reproduce the step-test crashes on purpose.

Two short IRC jobs (SN2 B3LYP-D3(BJ)/ma-def2-SVP SMD water TS, forward and backward, the 0.10 settings
of the step test, MaxIter 60, 2 cores each) are run
  shared   both in one scratch folder at the same time (twice: rounds 1 and 2)
  separate each in its own folder at the same time (once)
  overload each in its own folder, beside two load workers that keep starting fresh 4-core ORCA jobs:
           12 cores on 8 physical cores, as at 16:07 (once)
A watcher records every change of redirect.tmp in each folder (time, size, C6 and Edisp lines) and
copies the per-evaluation property files every 0.5 s, so they survive normal termination. The D3
values every evaluation used are then checked against D3 recomputed alone (as in check 1, part A).

Low priority, never beside a production job (J. Picione, 2026-11-11): it runs after the step-test
verdict and the E2 lane, and every round takes the whole budget (8 cores, 28000 MB) from the core
ledger and waits until no other ORCA process is running. ORCA is started here directly, because the
queue itself never puts two ORCA jobs in one folder. If no round crashes, the cause is recorded as
unknown.
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


LOAD_XYZ = ROOT / "reactions/e2_ethoxide_2_bromopropane/01_ts/level_test/M06-2X/2-bromopropane/2-bromopropane.xyz"


def load_worker(k: int, stop: threading.Event, env: dict, log: list) -> None:
    """Overload round: back-to-back 4-core ORCA single points (M06-2X/ma-def2-TZVP, 2-bromopropane,
    maxcore 1000), each in a new folder, until stop is set; every start is a fresh ORCA job, as at 16:06:39."""
    i = 0
    while not stop.is_set():
        wd = BASE / "overload" / f"load{k}" / f"run{i:03d}"
        wd.mkdir(parents=True, exist_ok=True)
        g = read_xyz(LOAD_XYZ)
        (wd / "sp.inp").write_text(orca.make_input("M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF DefGrid3",
                                                   g, 0, 1, nprocs=4, maxcore_mb=1000))
        with open(wd / "sp.out", "w") as fh:
            pr = subprocess.Popen([str(orca.ORCA_DIR / "orca.exe"), "sp.inp"], cwd=wd, stdout=fh,
                                  stderr=subprocess.STDOUT, env=env)
            log.append({"t": time.strftime("%H:%M:%S"), "load": k, "run": i, "event": "start"})
            while pr.poll() is None:
                if stop.is_set():
                    subprocess.run(["taskkill", "/F", "/T", "/PID", str(pr.pid)], capture_output=True)
                    pr.wait()
                time.sleep(1)
        for tmp in wd.glob("*.tmp*"):
            tmp.unlink(missing_ok=True)
        for big in wd.glob("*.gbw"):
            big.unlink(missing_ok=True)
        i += 1


def _foreign_orca() -> list[int]:
    import psutil
    return [p.pid for p in psutil.process_iter(["name"]) if (p.info["name"] or "").lower() == "orca.exe"]


def round_(rx, name: str, shared: bool, overload: bool = False) -> dict:
    wd = BASE / name
    if wd.exists():
        shutil.rmtree(wd)
    dirs = {d: (wd if shared else wd / d) for d in ("forward", "backward")}
    for d, p in dirs.items():
        p.mkdir(parents=True, exist_ok=True)
        shutil.copy(RXD / "01_ts/ts/ts.hess", p / "ts.hess")
        (p / f"irc_{d[0]}.inp").write_text(inp(rx, d))
    # J. Picione: never beside a production job. The round takes the whole budget, so it starts only
    # when nothing else holds cores and nothing can start beside it; no ORCA outside the ledger either.
    # J. Picione, 2026-11-11: no reproduction round during or near a reserved window (from 3 h before
    # the reservation starts until it is released)
    while cores.window_active(3.0):
        time.sleep(60)
    rid = cores.acquire(8, 28000, f"redirect repro {name} (whole machine; "
                        f"{'one folder' if shared else 'two folders'}{', 12 cores of load' if overload else ''})",
                        lane="redirect_repro")
    while _foreign_orca():
        time.sleep(30)
    env = os.environ.copy()
    env["PATH"] = f"{orca.MSMPI_BIN};{orca.ORCA_DIR};" + env.get("PATH", "")
    stop, wlog = threading.Event(), []
    th = threading.Thread(target=watch, args=(sorted(set(dirs.values())), stop, wlog), daemon=True)
    th.start()
    loads, llog = [], []
    if overload:
        for k in (1, 2):
            t = threading.Thread(target=load_worker, args=(k, stop, env, llog), daemon=True)
            t.start()
            loads.append(t)
        time.sleep(15)                       # the load is running before the IRCs start
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
        for t in loads:
            t.join(120)
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
    if overload:
        runs = sorted((BASE / name).glob("load*/run*/sp.out"))
        res["load"] = {"orca_starts": len(llog), "starts": llog,
                       "ended_normally": sum("ORCA TERMINATED NORMALLY" in r.read_text(errors="replace") for r in runs),
                       "note": "2 load workers × 4 cores beside the 2 × 2-core IRC jobs: 12 cores on 8; the "
                               "last load job of each worker is ended when the IRC jobs have ended"}
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
    # J. Picione, 2026-11-11: the shared-folder reproduction as planned, then one overload round
    # (separate folders, 12 cores on 8); every round has the machine to itself
    for name, shared, ovl in (("shared_round1", True, False), ("separate", False, False),
                              ("shared_round2", True, False), ("overload", False, True)):
        if name in R or ((BASE / "results.json").exists()
                         and name in json.loads((BASE / "results.json").read_text(encoding="utf-8"))):
            continue
        R[name] = round_(rx, name, shared, overload=ovl)
        print(name, json.dumps({k: v for k, v in R[name].items() if k != "redirect_writes"}), flush=True)
        (BASE / "results.json").write_text(json.dumps(R, indent=1), encoding="utf-8")
    R["d3_check"] = d3_check(rx, R)
    (BASE / "results.json").write_text(json.dumps(R, indent=1), encoding="utf-8")
