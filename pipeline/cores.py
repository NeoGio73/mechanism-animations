"""Core and memory budget shared by every lane (J. Picione, 2026-11-11).

The machine has 8 physical cores and 31 GB RAM. All ORCA jobs together may use at most 8 cores and
about 28 GB of declared memory (nprocs × %maxcore). Every ORCA job started through pipeline.orca.run
asks this ledger first and waits until its cores and memory are free; jobs are served first come,
first served. A lane is one running script (for example "level_test e2_arbitration").

  docs/level_test_queue/allocation.json   what is held and what is waiting (live state)
  docs/level_test_queue/allocation.log    one line per lane start, job start and job end

Holds belong to a process id; a hold whose process has ended is dropped automatically, so a crashed
lane cannot keep cores. Lanes started before this ledger existed are registered by hand with
`hold_for_pid` (cores reserved until that process ends).

CLI: python -m pipeline.cores status
     python -m pipeline.cores hold <pid> <cores> <mem_mb> <lane> <what>
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import psutil

CORES = 8
MEM_MB = 28000
Q = Path(__file__).resolve().parents[1] / "docs/level_test_queue"
STATE = Q / "allocation.json"
LOG = Q / "allocation.log"
LOCK = Q / "allocation.lock"
_lane_logged: set[int] = set()


def lane_name() -> str:
    if os.environ.get("MA_LANE"):
        return os.environ["MA_LANE"]
    a = sys.argv
    return " ".join([Path(a[0]).stem if a else "python"] + a[1:3])


@contextmanager
def _locked():
    Q.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    while True:
        try:
            fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except FileExistsError:
            try:
                if time.time() - LOCK.stat().st_mtime > 60:      # a lock left by a killed process
                    LOCK.unlink(missing_ok=True)
                    continue
            except FileNotFoundError:
                continue
            if time.time() - t0 > 600:
                raise RuntimeError("allocation ledger lock held for 10 min")
            time.sleep(0.2)
    try:
        yield
    finally:
        os.close(fd)
        LOCK.unlink(missing_ok=True)


def _alive(pid: int) -> bool:
    try:
        return psutil.Process(pid).status() != psutil.STATUS_ZOMBIE
    except psutil.Error:
        return False


def _read() -> dict:
    try:
        s = json.loads(STATE.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        s = {}
    s.setdefault("budget", {"cores": CORES, "mem_mb": MEM_MB})
    s.setdefault("held", [])
    s.setdefault("waiting", [])
    s.setdefault("paused", [])
    gone = [h for h in s["held"] if not _alive(h["pid"])]
    for h in gone:
        _log(f"release (process ended)  {h['lane']}: {h['what']}  {h['cores']} cores, {h['mem_mb']} MB")
    s["held"] = [h for h in s["held"] if _alive(h["pid"])]
    s["waiting"] = [w for w in s["waiting"] if _alive(w["pid"])]
    return s


def _write(s: dict) -> None:
    s["in_use"] = {"cores": sum(h["cores"] for h in s["held"]), "mem_mb": sum(h["mem_mb"] for h in s["held"])}
    tmp = STATE.with_suffix(".tmp")
    tmp.write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp, STATE)


def _log(msg: str) -> None:
    Q.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}\n")


def _summary(s: dict) -> str:
    c = sum(h["cores"] for h in s["held"])
    m = sum(h["mem_mb"] for h in s["held"])
    lanes: dict[str, list[int]] = {}
    for h in s["held"]:
        lanes.setdefault(h["lane"], [0, 0])
        lanes[h["lane"]][0] += h["cores"]
        lanes[h["lane"]][1] += h["mem_mb"]
    per = "; ".join(f"{k} {v[0]} cores {v[1]} MB" for k, v in lanes.items()) or "nothing"
    return f"in use {c}/{CORES} cores, {m}/{MEM_MB} MB ({per})"


def request(inp: str) -> tuple[int, int]:
    """Cores and declared memory (MB) of an ORCA input: nprocs × maxcore."""
    n = re.search(r"nprocs\s+(\d+)", inp, re.I)
    mc = re.search(r"%maxcore\s+(\d+)", inp, re.I)
    cores = int(n.group(1)) if n else 1
    return cores, cores * (int(mc.group(1)) if mc else 1000)


def acquire(cores: int, mem_mb: int, what: str, pid: int | None = None, lane: str | None = None,
            folder: str | None = None, front: bool = False) -> str:
    """Wait until the job fits in the budget and no older request is waiting; then hold it.
    folder: the job folder; two ORCA jobs never run in the same folder at once (ORCA writes fixed-name
    files such as redirect.tmp there), so a job waits while its folder is busy, without holding up
    the jobs behind it. front: the request goes to the head of the queue (reserved windows)."""
    if cores > CORES or mem_mb > MEM_MB:
        raise RuntimeError(f"{what}: {cores} cores / {mem_mb} MB exceeds the budget {CORES} / {MEM_MB} MB")
    pid = pid or os.getpid()
    lane = lane or lane_name()
    rid = f"{pid}-{time.time_ns()}"
    with _locked():
        s = _read()
        if pid not in _lane_logged and not any(h["pid"] == pid for h in s["held"]):
            _log(f"lane start  {lane} (pid {pid}); {_summary(s)}")
            _lane_logged.add(pid)
        entry = {"id": rid, "pid": pid, "lane": lane, "what": what, "cores": cores,
                 "mem_mb": mem_mb, "folder": folder, "since": time.strftime("%H:%M:%S")}
        if front:
            s["waiting"].insert(0, entry)
        else:
            s["waiting"].append(entry)
        _write(s)
    logged_wait = False
    while True:
        with _locked():
            s = _read()
            c = sum(h["cores"] for h in s["held"])
            m = sum(h["mem_mb"] for h in s["held"])
            busy = {h.get("folder") for h in s["held"]} - {None}
            ready = [w for w in s["waiting"] if w.get("folder") not in busy]
            first = bool(ready) and ready[0]["id"] == rid
            if first and c + cores <= CORES and m + mem_mb <= MEM_MB:
                me = next(w for w in s["waiting"] if w["id"] == rid)
                s["waiting"] = [w for w in s["waiting"] if w["id"] != rid]
                me["since"] = time.strftime("%H:%M:%S")
                s["held"].append(me)
                _write(s)
                _log(f"start  {lane}: {what}  {cores} cores, {mem_mb} MB -> {_summary(s)}")
                return rid
            if not logged_wait:
                _log(f"wait   {lane}: {what}  needs {cores} cores, {mem_mb} MB; {_summary(s)}")
                logged_wait = True
        time.sleep(10)


def release(rid: str, note: str = "") -> None:
    with _locked():
        s = _read()
        h = next((x for x in s["held"] if x["id"] == rid), None)
        s["held"] = [x for x in s["held"] if x["id"] != rid]
        _write(s)
        if h:
            _log(f"end    {h['lane']}: {h['what']}  {h['cores']} cores, {h['mem_mb']} MB{note} -> {_summary(s)}")


def hold_for_pid(pid: int, cores: int, mem_mb: int, lane: str, what: str) -> str:
    """Reserve cores for a process that does not use the ledger itself (a lane started before it)."""
    return acquire(cores, mem_mb, what, pid=pid, lane=lane)


def set_paused(entry: dict | None, lane: str) -> None:
    if entry:
        entry["lane"] = lane
    with _locked():
        s = _read()
        s["paused"] = [p for p in s["paused"] if p["lane"] != lane] + ([entry] if entry else [])
        _write(s)
        _log(f"{'pause ' if entry else 'resume'} {lane}; {_summary(s)}")


def resume_after(marker: str, paused_json: str, n_cores: int, mem_mb: int) -> None:
    """Queue a paused lane behind another: wait for the marker file, then for n_cores and mem_mb in
    the budget (held for the paused lane's process until it ends), then resume its processes."""
    m, pj = Path(marker), Path(paused_json)
    e = json.loads(pj.read_text(encoding="utf-8"))
    lane = "level_test " + e["lane"] if not e["lane"].startswith("level_test") else e["lane"]
    _log(f"queued {lane}: resumes after {m.name} exists, with {n_cores} cores, {mem_mb} MB")
    while not m.exists():
        time.sleep(60)
    hold_for_pid(e["root_pid"], n_cores, mem_mb, lane,
                 f"paused lane resumed: {n_cores // 2} chains × 2 cores, maxcore 2000")
    for pid in e["pids"]:
        try:
            psutil.Process(pid).resume()
        except psutil.Error:
            pass
    set_paused(None, lane)


WINDOW = Q / "window.json"


def window(start: str, cores_n: int, mem_mb: int, who: str) -> None:
    """A reserved window for other work on the machine (J. Picione, 2026-11-11: Unity from about 04:00,
    up to 1 h). From `start` (YYYY-mm-dd HH:MM, two hours before the window) a reservation of
    cores_n / mem_mb goes to the head of the queue. No job starts while it waits; once it is held,
    everything else together stays within (CORES − cores_n) / (MEM_MB − mem_mb). Running jobs are
    never stopped or paused. The reservation is held by this process until the release file
    window_release appears (J. Picione says he is done); then the full budget returns.
    The reservation is data in the ledger, so lanes started before this code respect it too."""
    t0 = time.mktime(time.strptime(start, "%Y-%m-%d %H:%M"))
    st = {"who": who, "reservation_from": start, "cores": cores_n, "mem_mb": mem_mb,
          "left_for_jobs": {"cores": CORES - cores_n, "mem_mb": MEM_MB - mem_mb},
          "status": "planned", "pid": os.getpid()}
    WINDOW.write_text(json.dumps(st, indent=1), encoding="utf-8")
    _log(f"window planned: {who}; from {start} no job starts unless all jobs together stay within "
         f"{CORES - cores_n} cores, {MEM_MB - mem_mb} MB")
    early = Q / "window_start_now"            # written when the verdict comes first: E2 then runs
    while time.time() < t0 and not early.exists():      # within the window budget from the start
        time.sleep(30)
    if early.exists():
        st["started_early"] = early.read_text(encoding="utf-8").strip()
        early.unlink(missing_ok=True)
    st["status"] = "waiting for running jobs to fall within the window budget"
    st["waiting_since"] = time.strftime("%Y-%m-%d %H:%M:%S")
    WINDOW.write_text(json.dumps(st, indent=1), encoding="utf-8")
    rid = acquire(cores_n, mem_mb, f"reserved window: {who}", lane="window", front=True)
    with _locked():
        s = _read()
        others = [h for h in s["held"] if h["id"] != rid]
    st["status"] = "granted"
    st["granted_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    st["others_at_grant"] = {"cores": sum(h["cores"] for h in others), "mem_mb": sum(h["mem_mb"] for h in others),
                             "jobs": [f"{h['lane']}: {h['what']}" for h in others]}
    st["orca_ram_mb_at_grant"] = round(sum(p.memory_info().private for p in psutil.process_iter(["name"])
                                           if (p.info["name"] or "").lower().startswith(("orca", "mpiexec", "smpd")))
                                       / 2**20)
    WINDOW.write_text(json.dumps(st, indent=1), encoding="utf-8")
    _log(f"window granted: {who}; other jobs now {st['others_at_grant']['cores']} cores, "
         f"{st['others_at_grant']['mem_mb']} MB declared, ORCA processes {st['orca_ram_mb_at_grant']} MB in RAM")
    rel = Q / "window_release"
    while not rel.exists():
        time.sleep(10)
    release(rid, "  (window over)")
    st["status"] = "released"
    st["released_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    WINDOW.write_text(json.dumps(st, indent=1), encoding="utf-8")
    _log(f"window released: {who}; full budget {CORES} cores, {MEM_MB} MB again")
    rel.unlink(missing_ok=True)


def guard(protect: str, until_marker: str, suspend_pids: list[int]) -> None:
    """J. Picione, 2026-11-11: until the step-size verdict, start a new ORCA job only when no protected
    job (held entry whose job path contains `protect`) is running. While one runs: a blocker stands at
    the head of the queue (cores 9, so it never fits and nothing behind it starts), and the driver
    processes in suspend_pids (lanes started before the ledger, which launch jobs without asking it)
    are suspended; their running ORCA jobs are not touched. Lifted when no protected job runs; ends
    when until_marker exists. The blocker is data in the ledger, so older lanes respect it too."""
    me = os.getpid()
    bid = f"guard-{me}"
    entry = {"id": bid, "pid": me, "lane": "guard", "folder": None, "cores": CORES + 1, "mem_mb": 0,
             "what": f"no new ORCA job while a {protect} job runs (J. Picione, until the verdict)",
             "since": time.strftime("%H:%M:%S")}
    state = None
    _log(f"guard on: no new ORCA job while a '{protect}' job runs; drivers held between jobs: {suspend_pids}")
    while not Path(until_marker).exists():
        with _locked():
            s = _read()
            prot = [h for h in s["held"] if protect in h["what"]]
            have = any(w["id"] == bid for w in s["waiting"])
            if prot and not have:
                entry["since"] = time.strftime("%H:%M:%S")
                s["waiting"].insert(0, entry)
            elif not prot and have:
                s["waiting"] = [w for w in s["waiting"] if w["id"] != bid]
            _write(s)
        new = "block" if prot else "open"
        for pid in suspend_pids:
            try:
                p = psutil.Process(pid)
                p.suspend() if prot else p.resume()
            except psutil.Error:
                pass
        if new != state:
            _log(f"guard {new}: " + ("; ".join(h["what"] for h in prot) if prot else "no protected job running"))
            state = new
        time.sleep(15)
    with _locked():
        s = _read()
        s["waiting"] = [w for w in s["waiting"] if w["id"] != bid]
        _write(s)
    for pid in suspend_pids:
        try:
            psutil.Process(pid).resume()
        except psutil.Error:
            pass
    _log("guard off: verdict given")


def window_active(margin_h: float = 3.0) -> bool:
    """True from margin_h hours before the reservation starts until the window is released: the
    crash reproduction (whole machine) does not start a round then."""
    if not WINDOW.exists():
        return False
    st = json.loads(WINDOW.read_text(encoding="utf-8"))
    if st.get("status") == "released":
        return False
    t0 = time.mktime(time.strptime(st["reservation_from"], "%Y-%m-%d %H:%M"))
    return time.time() >= t0 - margin_h * 3600


def status() -> str:
    with _locked():
        s = _read()
        _write(s)
    lines = [_summary(s)]
    lines += [f"  held    {h['lane']}: {h['what']}  {h['cores']} cores {h['mem_mb']} MB since {h['since']}"
              for h in s["held"]]
    lines += [f"  waiting {w['lane']}: {w['what']}  {w['cores']} cores {w['mem_mb']} MB since {w['since']}"
              for w in s["waiting"]]
    lines += [f"  paused  {p['lane']} since {p.get('paused_at')}" for p in s["paused"]]
    return "\n".join(lines)


if __name__ == "__main__":
    if sys.argv[1] == "status":
        print(status())
    elif sys.argv[1] == "hold":
        pid, c, m = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        print(hold_for_pid(pid, c, m, sys.argv[5], sys.argv[6]))
    elif sys.argv[1] == "guard":
        guard(sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4].split(",") if x])
    elif sys.argv[1] == "window":
        window(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif sys.argv[1] == "resume_after":
        resume_after(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
