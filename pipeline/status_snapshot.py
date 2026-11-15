"""Automatic status snapshot, written without the Claude app (2026-11-15).

The Claude desktop app can be closed by an automatic Microsoft Store update (it was at 00:45 on 2026-11-15) and
does not reopen by itself; while it is closed, no 07:00 / 17:00 report is written. This script runs outside the
app (started through WMI) and writes a plain, mechanical status file at 06:45 and 16:45:
docs/reports/<date>_<0645|1645>_auto.md. It reads files only; it never starts, stops or changes a calculation,
and it gives no interpretation. The full report is still written by Claude when the app is open.

Usage:  python -m pipeline.status_snapshot            (write one snapshot now)
        python -m pipeline.status_snapshot loop        (06:45 and 16:45 each day until 2026-11-20 12:00)
"""
from __future__ import annotations

import json
import re
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "docs/reports"
TIMES = [(6, 45), (16, 45)]
UNTIL = datetime(2026, 11, 20, 12, 0)
H = 2625.499639


def _tail(p: Path, n: int) -> list[str]:
    try:
        return p.read_text(encoding="utf-8", errors="replace").splitlines()[-n:]
    except OSError:
        return []


def _lanes() -> list[str]:
    out = []
    for rx in sorted((ROOT / "reactions").glob("*_m062x")):
        for d in ("backward", "forward"):
            wd = rx / "02_irc/true_irc" / d
            if not wd.exists():
                continue
            n = len(list(wd.glob("g*")))
            pj = wd / "path.json"
            if pj.exists():
                js = json.loads(pj.read_text(encoding="utf-8"))
                state = f"finished: {js['ended']}, {js['n_points']} points, {js['n_engrad']} gradient jobs"
            else:
                con = rx / f"02_irc/console_true_irc_{d}.txt"
                done = [ln for ln in _tail(con, 400) if ln.startswith("[done]")]
                last = done[-1] if done else ""
                m = re.search(r"E=(-?\d+\.\d+)", last)
                age = (time.time() - con.stat().st_mtime) / 60 if con.exists() else None
                state = (f"running: {n} gradient jobs so far" + (f", last energy {m[1]} Eh" if m else "")
                         + (f", console last written {age:.0f} min ago" if age is not None else ""))
            out.append(f"| {rx.name} | {d} | {state} |")
        if (rx / "02_irc/irc_path.json").exists():
            out.append(f"| {rx.name} | joined | 02_irc/irc_path.json present |")
    return out


def _single_points() -> list[str]:
    out = []
    for rx in sorted((ROOT / "reactions").glob("*_m062x")):
        base = rx / "02_irc/path_sp"
        if not base.exists():
            continue
        for lv in sorted(p for p in base.iterdir() if p.is_dir()):
            outs = list(lv.glob("f*/f.out"))
            ok = sum("ORCA TERMINATED NORMALLY" in o.read_text(errors="replace") for o in outs)
            out.append(f"| {rx.name} | {lv.name} | {ok} finished of {len(outs)} started |")
    return out


def snapshot(now: datetime | None = None) -> Path:
    now = now or datetime.now()
    stall = ROOT / "docs/STALL_ACTIVE"
    L = [f"# Automatic status {now:%Y-%m-%d %H:%M} (written without Claude; facts only, no interpretation)", "",
         "**Stall check:** " + (f"STALL ACTIVE: {stall.read_text(encoding='utf-8', errors='replace').strip()}"
                                if stall.exists() else "none (docs/STALL_ACTIVE absent)"),
         "", "## Production paths", "", "| Reaction | Side | State |", "|---|---|---|", *_lanes(),
         "", "## Single points along finished paths", "", "| Reaction | Level | Jobs |", "|---|---|---|",
         *_single_points(), "", "## Last lines of docs/production_true_irc.log", "", "```",
         *_tail(ROOT / "docs/production_true_irc.log", 6), "```", "",
         "## Last lines of the core queue log", "", "```",
         *[ln[:200] for ln in _tail(ROOT / "docs/level_test_queue/allocation.log", 4)], "```", "",
         "If this file exists but no full report from Claude follows it, the Claude app was probably closed "
         "(for example by an automatic update); the calculations themselves run outside the app."]
    REPORTS.mkdir(parents=True, exist_ok=True)
    p = REPORTS / f"{now:%Y-%m-%d_%H%M}_auto.md"
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return p


def _next(now: datetime) -> datetime:
    cands = [now.replace(hour=h, minute=m, second=0, microsecond=0) + timedelta(days=k)
             for k in (0, 1) for h, m in TIMES]
    return min(c for c in cands if c > now)


def loop() -> None:
    while True:
        t = _next(datetime.now())
        if t > UNTIL:
            return
        while datetime.now() < t:
            time.sleep(min(300, max(1, (t - datetime.now()).total_seconds())))
        try:
            snapshot(t)
        except Exception as e:          # a snapshot must never stop the loop
            (REPORTS / f"{t:%Y-%m-%d_%H%M}_auto_FAILED.txt").write_text(repr(e), encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "loop":
        loop()
    else:
        print(snapshot())
