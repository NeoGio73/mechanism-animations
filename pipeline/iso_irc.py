"""Isotope fairness check (J. Picione, 2026-11-12, item 2): ORCA IRC on the SN2 TS with deuterium
masses on all four hydrogens, the settings of the ORCA 0.10 step-test path (s010: Scale_Displ_SD
0.02359, Adapt_Scale_Displ false, Init_Displ length, Scale_Init_Displ 0.02819, TS Hessian read), 30
steps per side. A true IRC changes with the masses; a Cartesian steepest-descent path changes only
through its first step (the transition-state mode it starts along).

Masses are set on the coordinate lines (ORCA "M = <mass>"); the output's mass table and the direction
of the first displacement (compared with the normal-mass and deuterium-mass mass-weighted TS modes)
show whether ORCA used them.
Folder: <step_convergence>/iso_d4/<forward|backward>/; comparison: pipeline/iso_compare.py
Usage: python -m pipeline.iso_irc [hess]   (hess: deuterium masses also in the TS Hessian file, folder iso_d4_hess)
"""
from __future__ import annotations

import re
import shutil
from concurrent.futures import ThreadPoolExecutor

from . import orca
from .irc_angle import st

D_MASS = 2.01410178
OUT = st.OUT / "iso_d4"
OUT_H = st.OUT / "iso_d4_hess"   # variant: deuterium masses also in the TS Hessian file


def inp_for(direction: str) -> str:
    src = (st.jd("s010", direction) / f"irc_{direction[0]}.inp").read_text()
    src = re.sub(r"MaxIter \d+", "MaxIter 30", src)
    lines = []
    for ln in src.splitlines():
        p = ln.split()
        if len(p) == 4 and p[0] == "H":
            ln = f"{ln}  M = {D_MASS}"
        lines.append(ln)
    return "\n".join(lines) + "\n"


def d4_hess(src, dst):
    """The TS Hessian file with deuterium masses on H in its $atoms block (the Hessian itself unchanged)."""
    L = src.read_text().splitlines()
    i = L.index("$atoms")
    n = int(L[i + 1])
    for j in range(i + 2, i + 2 + n):
        p = L[j].split()
        if p[0] == "H":
            L[j] = f" H     {D_MASS:.5f}   " + "   ".join(p[2:])
    dst.write_text("\n".join(L) + "\n")


def one(direction: str, hess_d4: bool = False):
    wd = (OUT_H if hess_d4 else OUT) / direction
    wd.mkdir(parents=True, exist_ok=True)
    if hess_d4:
        d4_hess(st.jd("s010", direction) / "ts.hess", wd / "ts.hess")
    else:
        shutil.copy(st.jd("s010", direction) / "ts.hess", wd / "ts.hess")
    r = orca.run(wd, f"irc_{direction[0]}", inp_for(direction), orca.RunLog(wd.parent / "run_log_iso.md"),
                 note=f"isotope check: ORCA IRC {direction}, D on all H, s010 settings, 30 steps"
                      + (", D masses also in the Hessian file" if hess_d4 else ""))
    return direction, r.ok


def main(hess_d4: bool = False):
    (OUT_H if hess_d4 else OUT).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(2) as ex:
        print(list(ex.map(lambda d: one(d, hess_d4), ("forward", "backward"))))


if __name__ == "__main__":
    import sys
    main(len(sys.argv) > 1 and sys.argv[1] == "hess")
