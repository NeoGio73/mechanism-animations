"""Step test, near the TS (J. Picione, 2026-11-12): the three pairs of report4 restricted to |s| <= 1.0
bohr·amu½, the region that sets the barrier and holds the coupled-cluster maximum. Same matching as
report4 (each frame of the coarser path vs the finer path interpolated at the same C–O / C–Br).
Also the B3LYP energy profile of each path near the TS on a common C–Br grid (product side), and the
largest gradient along that stretch, for the coupled-cluster-maximum question.
Result: <step_convergence>/near_ts.json
Usage: python -m pipeline.step_near_ts
"""
from __future__ import annotations

import json

import numpy as np

from .irc_angle import st

PAIRS = [("0.10 vs new 0.05", "s010", "s005f10"), ("new 0.05 vs 0.025", "s005f10", "s0025"),
         ("old 0.05 vs new 0.05", "s005", "s005f10")]


def main() -> dict:
    ts = st.read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    ends = st.adaptive_ends()
    P = {n: st.load(n, ts, ends) for n in st.R4}
    out = {"window_abs_s": 1.0, "pairs": {}}
    for title, a, b in PAIRS:
        c = st.compare(P[a], P[b], span=(-1.0, 1.0))
        row = {}
        for side, v in c.items():
            if not v.get("n"):
                row[side] = {"n": 0}
                continue
            rows = v["rows"]
            iE = max(range(len(rows)), key=lambda i: abs(rows[i][2]))
            iB = max(range(len(rows)), key=lambda i: abs(rows[i][3]))
            row[side] = {"match_by": v["match_by"], "other": v["other"], "n": len(rows),
                         "max_dE_kJ": rows[iE][2], "at_s": rows[iE][0], "at_dist_A": rows[iE][1],
                         "max_dbond_A": rows[iB][3], "bond_at_s": rows[iB][0], "bond_at_dist_A": rows[iB][1],
                         "pass": abs(rows[iE][2]) <= st.TOL_E and abs(rows[iB][3]) <= st.TOL_B}
        out["pairs"][title] = row
    # frame spacing near the TS on each path, and the B3LYP energy on a common C–Br grid (product side)
    grid = np.round(np.arange(2.39, 2.53, 0.02), 3)
    out["spacing_near_TS"] = {n: float(np.median([r["step"] for r in P[n] if 0 < abs(r["s"]) <= 1.0])) for n in st.R4}
    prof = {}
    for n in st.R4:
        R = sorted([r for r in P[n] if r["side"] in ("product", "TS")], key=lambda r: r["CBr"])
        x = np.array([r["CBr"] for r in R])
        prof[n] = {f"{g:.2f}": float(np.interp(g, x, [r["E"] for r in R])) for g in grid}
    out["B3LYP_E_on_CBr_grid_kJ"] = prof
    (st.OUT / "near_ts.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
