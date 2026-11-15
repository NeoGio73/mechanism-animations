"""Steps taken by every running production IRC side, for the reports (J. Picione, 2026-11-15): a side past
FLAG_STEPS steps (about twice the longest side so far, 188) is flagged, so a runaway is caught within about half
a day; the hard backstop stays at 3000 steps (pipeline.true_irc_prod.MAX_CYCLES).
The step count is the last cycle line pysisyphus prints in 02_irc/console_true_irc_<side>.txt.
Usage: python -m pipeline.side_steps   (prints one line per running side; exit code 1 if any is flagged)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLAG_STEPS = 400
CYCLE = re.compile(r"^\s+(\d+)\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*$")


def running_sides() -> list[dict]:
    out = []
    for side in sorted(ROOT.glob("reactions/*_m062x/02_irc/true_irc/*")):
        if not side.is_dir() or side.name not in ("forward", "backward") or (side / "path.json").exists():
            continue
        con = side.parents[1] / f"console_true_irc_{side.name}.txt"
        steps = 0
        if con.exists():
            for ln in con.read_text(encoding="utf-8", errors="replace").splitlines():
                m = CYCLE.match(ln)
                if m:
                    steps = int(m[1])
        out.append({"reaction": side.parents[2].name, "side": side.name, "steps": steps,
                    "flag": steps > FLAG_STEPS})
    return out


if __name__ == "__main__":
    rows = running_sides()
    for r in rows:
        print(f"{r['reaction']} {r['side']}: {r['steps']} steps" + (f"  FLAG: past {FLAG_STEPS} steps" if r["flag"] else ""))
    if not rows:
        print("no production side running")
    sys.exit(1 if any(r["flag"] for r in rows) else 0)
