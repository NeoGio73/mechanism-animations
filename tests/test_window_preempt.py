"""Test of the window pre-emption (J. Picione, 2026-11-15: a job still running 5 minutes before a reserved
window is stopped, queued again, and rerun after the window in a clean folder; noted in the queue log).
Scratch ledger and a stand-in for ORCA (a Python process); never the live queue or a real job.
Run: .venv/Scripts/python.exe -X utf8 tests/test_window_preempt.py
"""
import json
import subprocess
import sys
import tempfile
import threading
import time
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import cores, orca  # noqa: E402

PY = getattr(sys, "_base_executable", sys.executable)


def main(tmp: Path) -> None:
    q = tmp / "queue"
    q.mkdir()
    cores.Q, cores.STATE, cores.LOG, cores.LOCK = q, q / "allocation.json", q / "allocation.log", q / "allocation.lock"
    cores.WINDOW = q / "window.json"
    cores.POLL_S = 0.2
    cores.ORCA_PROC_PREFIXES = ("orca", "mpiexec", "smpd", "python")   # the stand-in runs as python.exe
    cores.PREEMPT_AWARE_SINCE = 0.0
    calls = []
    real_popen = subprocess.Popen

    def fake_popen(args, cwd=None, stdout=None, stderr=None, env=None):
        calls.append(time.time())
        if len(calls) == 1:     # first run: a long job, still running when the window comes
            code = "import time; print('running', flush=True); time.sleep(600)"
        else:                   # the rerun after the window: finishes normally
            code = "print('FINAL SINGLE POINT ENERGY      -1.000000000000'); print('****ORCA TERMINATED NORMALLY****')"
        return real_popen([PY, "-c", code], cwd=cwd, stdout=stdout, stderr=stderr)

    # only orca.run's ORCA launch is replaced; taskkill and everything else use the real subprocess module
    orca.subprocess = types.SimpleNamespace(Popen=fake_popen, run=subprocess.run, STDOUT=subprocess.STDOUT,
                                            TimeoutExpired=subprocess.TimeoutExpired)
    jd = tmp / "reactions" / "x" / "g0001"
    log = orca.RunLog(tmp / "run_log.md")
    inp = "! fake\n%pal nprocs 4 end\n%maxcore 1000\n\n* xyz 0 1\nH 0 0 0\n*\n"
    out = {}
    t_job = threading.Thread(target=lambda: out.update(res=orca.run(jd, "j", inp, log, note="test job")), daemon=True)
    t_job.start()
    time.sleep(2)                                      # the job is running
    now = time.time()
    start = time.strftime("%Y-%m-%d %H:%M", time.localtime(now - 60))
    clear_by = time.strftime("%Y-%m-%d %H:%M", time.localtime(now + 300 + 60))   # pre-emption within ~1 min
    t_win = threading.Thread(target=cores.window, args=(start, 8, 8000, "test window", clear_by), daemon=True)
    t_win.start()
    for _ in range(240):
        st = json.loads(cores.WINDOW.read_text(encoding="utf-8")) if cores.WINDOW.exists() else {}
        if st.get("status") == "granted":
            break
        time.sleep(0.5)
    assert st.get("status") == "granted", f"window not granted: {st}"
    assert st["others_at_grant"]["cores"] == 0, st["others_at_grant"]
    assert st.get("preempted") and st["preempted"][0]["orca_pids"], st.get("preempted")
    assert len(calls) == 1 and not out, "the job must wait behind the window, not rerun during it"
    moved = [p for p in jd.iterdir() if p.is_dir() and p.name.startswith("preempted_")]
    assert len(moved) == 1 and (moved[0] / "j.out").exists(), list(jd.iterdir())
    (q / "window_release").write_text("done")
    t_job.join(60)
    assert out.get("res") and out["res"].ok and len(calls) == 2, (out, calls)
    qlog = cores.LOG.read_text(encoding="utf-8")
    assert "preempted " in qlog and "window granted" in qlog and "window released" in qlog
    i_pre, i_grant, i_rel = qlog.index("preempted "), qlog.index("window granted"), qlog.index("window released")
    i_rerun = qlog.rindex("start  ")
    assert i_pre < i_grant < i_rel < i_rerun, "order: pre-empt, grant, release, rerun"
    assert "stopped before the reserved window" in (tmp / "run_log.md").read_text(encoding="utf-8")
    print("window pre-emption: passed (stopped, queued, granted with 0 other cores, rerun after release in a clean folder)")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        main(Path(td))
