"""IBO single points on the final M06-2X paths, started by themselves after E2's last path single point
(J. Picione's post-E2 order, 2026-11-23: Section 3.2 numbers, IBO tracking, renders, SMD-scheme check, water test).

One lane runs pipeline/orbital_steps.orbital_sps (one working-level single point with IAO/IBO localization per IRC
frame, every frame, the tertiary mirror half included, as he asked) for each reaction in its list, one after the
other. Two lanes run side by side, 4 cores each, through the ordinary core queue. Tracking and cube files are cheap
and are run afterwards by hand, after the outputs are checked. Finished frames are reused (pipeline.orca.run), so a
lane stopped by an outage restarts where it was.

Usage: python -X utf8 -m pipeline.ibo_lanes <reaction>[,<reaction>...]
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from .orbital_steps import orbital_sps
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
GATE = ROOT / "reactions/e2_ethoxide_2_bromopropane_m062x/02_irc/energies_DLPNO-CCSD(T)_ma-def2-TZVP.json"


def main(names: list[str]) -> None:
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} waiting for E2's last path single point ({GATE.name})", flush=True)
    while not GATE.exists():
        time.sleep(60)
    for n in names:
        t0 = time.time()
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {n}: IBO single points start", flush=True)
        orbital_sps(Reaction.load(ROOT / "reactions" / n))
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {n}: IBO single points done ({(time.time() - t0) / 3600:.2f} h)",
              flush=True)
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} lane finished", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","))
