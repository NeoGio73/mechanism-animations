"""IBO single points on the final M06-2X paths, started by themselves after E2's last path single point
(J. Picione's post-E2 order, 2026-11-23: Section 3.2 numbers, IBO tracking, renders, SMD-scheme check, water test).

Start condition (J. Picione, 2026-11-25; pipeline/ibo_check.gate): all 507 E2 DLPNO-CCSD(T) frames end with ORCA
TERMINATED NORMALLY, E2's level file is written, no E2 lane process is alive and no E2 job holds cores. If the E2
lane ends without completing, the lanes hold and docs/IBO_HOLD.md says why.

One lane runs pipeline/orbital_steps.orbital_sps (one working-level single point with IAO/IBO localization per IRC
frame, every frame, the tertiary mirror half included) for each reaction in its list, then pipeline/ibo_check.finish
for that reaction: the output check, and only if it passes, tracking, cube files and (tertiary) the mirror
symmetry check. A failed check stops that reaction (not tracked) and the lane goes on with the next one. Two lanes
run side by side, 4 cores each, through the ordinary core queue (2 × 4 × 3000 MB = 24 000 MB declared). Finished
frames are reused, so a lane stopped by an outage restarts where it was.

Usage: python -X utf8 -m pipeline.ibo_lanes <reaction>[,<reaction>...]
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

from .ibo_check import finish, gate
from .orbital_steps import orbital_sps
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]


def say(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}", flush=True)


def main(names: list[str]) -> None:
    say("waiting for E2 to be complete (all 507 frames normal, no E2 process, no E2 cores held)")
    gate()
    for n in names:
        t0 = time.time()
        say(f"{n}: IBO single points start")
        try:
            orbital_sps(Reaction.load(ROOT / "reactions" / n))
        except Exception as e:                      # a failed frame: not a result; stop this reaction, report it
            say(f"{n}: IBO single points STOPPED: {e}")
            st = ROOT / "reactions" / n / "03_orbitals/status.json"
            st.parent.mkdir(exist_ok=True)
            st.write_text(json.dumps({"reaction": n, "stage": f"STOPPED: single point failed ({e})"}, indent=1),
                          encoding="utf-8")
            continue
        say(f"{n}: IBO single points done ({(time.time() - t0) / 3600:.2f} h); output check, tracking, cubes")
        try:
            st = finish(n)
            say(f"{n}: {st.get('stage')}; check {st.get('check')}"
                + (f"; mirror worst overlap {st['mirror_check']['worst_overlap']}" if "mirror_check" in st else ""))
        except Exception:
            say(f"{n}: finish step FAILED (reaction stopped, not tracked further):\n{traceback.format_exc()}")
    say("lane finished")


if __name__ == "__main__":
    main(sys.argv[1].split(","))
