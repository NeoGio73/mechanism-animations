"""Coupled-cluster maximum near the TS on the true (mass-weighted) 0.10 IRC (J. Picione, 2026-11-12,
item 5), the same procedure as pipeline/cc_near_ts.py on the ORCA paths.

DLPNO-CCSD(T) TightPNO/ma-def2-TZVP, SMD water, 4 cores, each job alone in its own folder, at every
frame with -1.2 <= s <= +0.4 bohr·amu½ of a pipeline/true_irc2.py run (native 0.10 step), plus the
TS (the same TS job as cc_near_ts, reused). Highest frame and the parabola through it and its two
neighbours; energies relative to the TS geometry at the same level; s from the TS, + toward the
reactant complex. Compared with the ORCA paths' result (cc_near_ts.json).
Result: <step_convergence>/cc_true_irc_<run>.json and .md
Usage: python -m pipeline.cc_true_irc <run>
"""
from __future__ import annotations

import json
import sys

import numpy as np

from .cc_near_ts import LEVEL, W, peak, sp
from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, st

T2 = st.OUT / "true_irc2"
H = 2625.499639


def frames(run):
    ts = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    sym = list(ts.symbols)
    m = np.repeat([MASS[x] for x in sym], 3)
    x0 = np.asarray(ts.coords)
    out = [dict(s=0.0, CBr=st.dist(x0, st.CBR), CO=st.dist(x0, st.CO), xyz=x0, side="TS", k=0)]
    for d in ("forward", "backward"):
        X = [np.asarray(f.coords) for f in read_xyz_frames(T2 / run / d / "path.xyz")]
        X = [x0] + [c for c in X if np.max(np.abs(c - x0)) > 1e-6]
        side = "reactant" if st.dist(X[-1], st.CO) > st.dist(X[-1], st.CBR) else "product"
        s = 0.0
        for k in range(1, len(X)):
            s += float(np.linalg.norm(((X[k] - X[k - 1]) / BOHR).ravel() * np.sqrt(m)))
            out.append(dict(s=s if side == "reactant" else -s, CBr=st.dist(X[k], st.CBR), CO=st.dist(X[k], st.CO),
                            xyz=X[k], side=side, k=k))
    return sym, sorted(out, key=lambda r: r["s"])


def main(run: str) -> dict:
    from concurrent.futures import ThreadPoolExecutor
    sym, F = frames(run)
    E_ts = sp("TS", F[[r["side"] for r in F].index("TS")], sym)
    sel = [f for f in F if f["side"] != "TS" and -1.2 <= f["s"] <= 0.4]
    with ThreadPoolExecutor(2) as ex:
        E = list(ex.map(lambda f: sp(f"true_irc2_{run}", f, sym), sel))
    rows = [dict(f, E=0.0) for f in F if f["side"] == "TS"] + [dict(f, E=(e - E_ts) * H) for f, e in zip(sel, E)]
    rows.sort(key=lambda r: r["s"])
    pk = peak(rows)
    orca_ref = json.loads((st.OUT / "cc_near_ts.json").read_text(encoding="utf-8"))["paths"]
    res = {"level": LEVEL, "run": run, "E_TS_Eh": E_ts, "peak": pk,
           "computed": [{k: r[k] for k in ("side", "k", "s", "CBr", "CO", "E")} for r in rows],
           "orca_paths": {n: v["peak"] for n, v in orca_ref.items()}}
    (st.OUT / f"cc_true_irc_{run}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    L = [f"# Coupled-cluster maximum near the TS: true IRC ({run}) against the ORCA paths", "",
         "DLPNO-CCSD(T) TightPNO/ma-def2-TZVP, SMD water; kJ/mol relative to the TS geometry.", "",
         "| Path | Highest frame s | C–Br / C–O (Å) | E | Parabola s | E | C–Br (Å) |", "|---|---|---|---|---|---|---|"]
    for lab, v in [(f"true IRC {run}", pk)] + [(n, p) for n, p in res["orca_paths"].items()]:
        f, p = v["frame"], v["parabola"]
        L.append(f"| {lab} | {f['s']:+.3f} | {f['CBr']:.4f} / {f['CO']:.4f} | {f['E']:+.3f} | "
                 + (f"{p['s']:+.3f} | {p['E_kJ']:+.3f} ({p['E_kJ'] / 4.184:+.4f} kcal/mol) | {p['CBr_A']:.4f} |" if p else "— | — | — |"))
    (st.OUT / f"cc_true_irc_{run}.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return res


if __name__ == "__main__":
    main(sys.argv[1])
