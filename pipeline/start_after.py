"""Start one command when the processes it waits for have ended (J. Picione, 2026-11-25: start the explicit-water
test as soon as the CPU is free after the SMD-scheme check, without waiting for his go).

Waits until no process command line contains the marker text AND every console file named ends with its
'finished' line; then runs the command once and exits with its return code. If a console never says
'finished' (a lane died), it does not start: it writes the reason and exits, so nothing starts on a failed check.

Usage: python -X utf8 -m pipeline.start_after "<marker>" <console1>[,<console2>...] -- <command ...>
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[1]


def say(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}", flush=True)


def alive(marker: str) -> bool:
    for p in psutil.process_iter(["cmdline"]):
        try:
            cl = " ".join(p.info["cmdline"] or [])
        except (psutil.Error, TypeError):
            continue
        if marker in cl and "pipeline.start_after" not in cl:
            return True
    return False


def main() -> None:
    i = sys.argv.index("--")
    marker, consoles, cmd = sys.argv[1], [ROOT / c for c in sys.argv[2].split(",")], sys.argv[i + 1:]
    say(f"waiting for '{marker}' to end; then: {' '.join(cmd)}")
    while alive(marker):
        time.sleep(60)
    time.sleep(5)
    bad = [c.name for c in consoles
           if not c.exists() or "finished" not in (c.read_text(encoding="utf-8", errors="replace").strip()
                                                     .splitlines() or [""])[-1]]
    if bad:
        say(f"NOT STARTED: {bad} did not end with 'finished'")
        sys.exit(1)
    say("start")
    r = subprocess.run(cmd, cwd=ROOT)
    say(f"ended, return code {r.returncode}")
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
