"""Tests for two queue-safety fixes (run: .venv/Scripts/python.exe -X utf8 tests/test_queue_safety.py).

1. The verdict guard (pipeline/cores.guard). On Windows every suspend adds to a count and every resume
   removes one; the old guard suspended a driver every 15 s and resumed it once, which left the
   bridged-ion driver frozen for 1.5 h on 2026-11-12. These tests use a dummy process and a scratch
   copy of the ledger, never the live queue or a real job.
2. The clean-folder rule in pipeline/orca.run: a rerun after a crash never starts among the crashed
   run's files.
"""
import json
import subprocess
import sys
import tempfile
import threading
import time
import types
from pathlib import Path

import psutil

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import cores, orca  # noqa: E402

_real_sleep = time.sleep


def _dummy(counter: Path) -> subprocess.Popen:
    """A process that writes a rising counter, so 'frozen' and 'running' can be seen from outside."""
    code = (f"import time\nn=0\nwhile True:\n    n+=1\n    open(r'{counter}','w').write(str(n))\n"
            f"    time.sleep(0.05)\n")
    # the venv python.exe is a launcher that starts the real interpreter as a child; use the real one
    return subprocess.Popen([getattr(sys, "_base_executable", sys.executable), "-c", code])


def _advancing(counter: Path) -> bool:
    a = int(counter.read_text() or 0)
    _real_sleep(0.6)
    return int(counter.read_text() or 0) > a


def _scratch_ledger(d: Path) -> None:
    cores.Q, cores.STATE, cores.LOG, cores.LOCK = d, d / "allocation.json", d / "allocation.log", d / "allocation.lock"
    # guard waits 15 s per round; run it 75 times faster here
    cores.time = types.SimpleNamespace(**{k: getattr(time, k) for k in dir(time) if not k.startswith("_")})
    cores.time.sleep = lambda s: _real_sleep(s / 75)


def test_windows_counts_suspends(tmp: Path) -> None:
    """The cause: 15 suspends and one resume leave the process frozen; _resume_fully wakes it."""
    c = tmp / "count_a.txt"
    p = _dummy(c)
    _real_sleep(0.5)
    try:
        q = psutil.Process(p.pid)
        for _ in range(15):
            q.suspend()
        q.resume()
        assert q.status() == psutil.STATUS_STOPPED and not _advancing(c), "one resume should not wake it"
        cores._resume_fully(p.pid)
        assert q.status() != psutil.STATUS_STOPPED and _advancing(c)
    finally:
        p.kill()
    print("ok  cause reproduced: 15 suspends + 1 resume = still frozen; _resume_fully wakes it")


def test_guard(tmp: Path) -> None:
    """Guard blocks for many rounds, opens, blocks again for many rounds, then ends: the driver must
    run whenever no protected job runs, and run after the guard ends."""
    c = tmp / "count_b.txt"
    p = _dummy(c)
    _real_sleep(0.5)
    marker = tmp / "VERDICT_GIVEN"
    try:
        rid = cores.acquire(1, 100, "step_convergence/s005/forward/irc_f (fake)", pid=psutil.Process().pid)
        t = threading.Thread(target=cores.guard, args=("step_convergence/s0", str(marker), [p.pid]), daemon=True)
        t.start()
        _real_sleep(4)                      # about 20 guard rounds while 'protected job' runs
        assert psutil.Process(p.pid).status() == psutil.STATUS_STOPPED and not _advancing(c), "should be held"
        print("ok  driver held while a protected job runs (about 20 guard rounds)")
        cores.release(rid)
        _real_sleep(1)
        assert psutil.Process(p.pid).status() != psutil.STATUS_STOPPED and _advancing(c), "should run"
        print("ok  driver runs as soon as no protected job runs")
        rid = cores.acquire(1, 100, "step_convergence/s010/backward/irc_b (fake)", pid=psutil.Process().pid)
        _real_sleep(4)
        assert psutil.Process(p.pid).status() == psutil.STATUS_STOPPED
        marker.write_text("FAIL")           # verdict given while still blocked
        t.join(timeout=10)
        assert not t.is_alive()
        assert psutil.Process(p.pid).status() != psutil.STATUS_STOPPED and _advancing(c), "should run after guard"
        cores.release(rid)
        print("ok  driver runs after the guard ends while blocked (the 2026-11-12 case)")
        st = json.loads(cores.STATE.read_text(encoding="utf-8"))
        assert not st["waiting"], "guard blocker left in the ledger"
        print("ok  guard blocker removed from the ledger")
    finally:
        p.kill()


def test_crashed_run_moved(tmp: Path) -> None:
    d = tmp / "chain"
    d.mkdir()
    files = {"plus.inp": "x", "plus.out": "...error termination", "plus.gbw": "x", "plus_trj.xyz": "x",
             "redirect.tmp": "x", "plus_min.inp": "x", "plus_min.gbw": "x", "p00.inp": "x",
             "p00.out": "****ORCA TERMINATED NORMALLY****"}
    for k, v in files.items():
        (d / k).write_text(v)
    orca._move_crashed_run(d, "plus", d / "plus.out", d / "plus.stopped.txt")
    moved = [x for x in d.iterdir() if x.is_dir()]
    assert len(moved) == 1 and moved[0].name.startswith("crashed_") and moved[0].name.endswith("_plus")
    assert sorted(x.name for x in moved[0].iterdir()) == ["plus.gbw", "plus.inp", "plus.out", "plus_trj.xyz", "redirect.tmp"]
    assert sorted(x.name for x in d.iterdir() if x.is_file()) == ["p00.inp", "p00.out", "plus_min.gbw", "plus_min.inp"]
    orca._move_crashed_run(d, "p00", d / "p00.out", d / "p00.stopped.txt")   # finished normally: untouched
    assert (d / "p00.out").exists()
    print("ok  crashed run moved whole to crashed_<time>_plus/; other chain jobs' files and finished jobs untouched")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_crashed_run_moved(tmp)
        _scratch_ledger(tmp)
        test_windows_counts_suspends(tmp)
        test_guard(tmp)
    print("all passed")
