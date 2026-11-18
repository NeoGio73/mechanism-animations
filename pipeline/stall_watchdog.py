"""Stall watchdog (J. Picione, 2026-11-13), separate from the queue's own code (pipeline/cores.py is
not imported). Every minute it reads the queue state file and the running processes directly:

  stall = jobs are waiting in the queue, the free cores would fit at least one of them, and no ORCA
          job has started for 10 minutes (newest orca.exe process start, or newest "start" line
          in the queue log, whichever is later).

On a stall it appends one line to docs/stall_log.md and creates docs/STALL_ACTIVE (read at the top of
every report and before every time estimate). When jobs start again it appends a "cleared" line and
removes the marker. Usage: python -m pipeline.stall_watchdog [--once]

Planned pause (J. Picione, 2026-11-15): from the start of a reserved window's reservation
(docs/level_test_queue/window.json) until it is released, waiting jobs are expected, so no stall is raised;
the pause and its end go to the stall log. After the release the 10-minute limit counts from the release.
"""
from __future__ import annotations

import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[1]
Q = ROOT / "docs/level_test_queue"
STATE, QLOG, WINDOW = Q / "allocation.json", Q / "allocation.log", Q / "window.json"
LOG, MARK = ROOT / "docs/stall_log.md", ROOT / "docs/STALL_ACTIVE"
LIMIT_S = 600
CORES = 8


def last_queue_start() -> float:
    try:
        with open(QLOG, "rb") as fh:
            fh.seek(0, 2)
            fh.seek(max(0, fh.tell() - 400_000))
            tail = fh.read().decode("utf-8", "replace").splitlines()
    except OSError:
        return 0.0
    # 'end' counts too (fix 2026-11-18): the idle time runs from the moment cores came free, not from the last
    # start. Before, a long job on all 8 cores followed by a few seconds' gap between one job's end and the next
    # start read as a stall (false alarms 09:06 and 12:50 on 2026-11-18).
    for ln in reversed(tail):
        m = re.match(r"(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)  (start|end) ", ln)
        if m:
            return datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S").timestamp()
    return 0.0


def last_orca_start() -> float:
    t = 0.0
    for p in psutil.process_iter(["name", "create_time"]):
        if (p.info["name"] or "").lower() == "orca.exe":
            t = max(t, p.info["create_time"] or 0.0)
    return t


def planned_pause() -> tuple[bool, float]:
    """(in a planned pause now, release time of the last window or 0)."""
    try:
        st = json.loads(WINDOW.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False, 0.0
    rel = st.get("released_at")
    rel_t = datetime.strptime(rel, "%Y-%m-%d %H:%M:%S").timestamp() if rel else 0.0
    if st.get("status") == "released":
        return False, rel_t
    t0 = datetime.strptime(st["reservation_from"], "%Y-%m-%d %H:%M").timestamp()
    return (st.get("status") != "planned" or time.time() >= t0), rel_t


def check() -> dict:
    try:
        s = json.loads(STATE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"ok": None, "why": "queue state unreadable this minute"}
    alive = lambda pid: psutil.pid_exists(pid)
    held = [h for h in s.get("held", []) if alive(h["pid"])]
    waiting = [w for w in s.get("waiting", []) if alive(w["pid"])]
    free = CORES - sum(h["cores"] for h in held)
    fits = [w for w in waiting if w["cores"] <= free]
    pause, rel_t = planned_pause()
    last = max(last_queue_start(), last_orca_start(), rel_t)
    idle = time.time() - last if last else float("inf")
    stalled = bool(waiting) and bool(fits) and idle > LIMIT_S and not pause
    return {"ok": not stalled, "pause": pause, "waiting": len(waiting), "free": free, "idle_s": idle, "last": last,
            "head": waiting[0]["what"] if waiting else None, "lane": waiting[0]["lane"] if waiting else None}


def note(line: str) -> None:
    new = not LOG.exists()
    with open(LOG, "a", encoding="utf-8") as fh:
        if new:
            fh.write("# Stall log\n\nWritten by pipeline/stall_watchdog.py (separate from the queue code).\n\n")
        fh.write(line + "\n")


def main(once: bool = False) -> None:
    active = MARK.exists()
    in_pause = None
    while True:
        try:
            r = check()
            now = time.strftime("%Y-%m-%d %H:%M")
            if r.get("pause") is not None and r["pause"] != in_pause:
                if r["pause"]:
                    note(f"- {now}  planned pause (reserved window): no stall alarm until the window is released")
                elif in_pause:
                    note(f"- {now}  planned pause over (window released): stall check active again, "
                         f"10-minute limit counted from the release")
                in_pause = r["pause"]
            if r["ok"] is False and not active:
                last = time.strftime("%H:%M", time.localtime(r["last"])) if r["last"] else "never"
                line = (f"- {now}  STALL: {r['waiting']} job(s) queued, {r['free']} of {CORES} cores free, "
                        f"no ORCA job started since {last} ({r['idle_s'] / 60:.0f} min). Head of queue: "
                        f"{r['lane']}: {r['head']}")
                note(line)
                MARK.write_text(line + "\n", encoding="utf-8")
                active = True
            elif r["ok"] is True and active:
                note(f"- {now}  cleared: jobs starting again")
                MARK.unlink(missing_ok=True)
                active = False
        except Exception as e:                       # the watchdog itself must keep running
            note(f"- {time.strftime('%Y-%m-%d %H:%M')}  watchdog error: {type(e).__name__}: {e}")
        if once:
            print(json.dumps(check(), default=str))
            return
        time.sleep(60)


if __name__ == "__main__":
    main(once="--once" in sys.argv)
