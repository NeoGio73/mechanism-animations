"""Completeness of the pysisyphus IRC sides (J. Picione, 2026-11-12, item 3: "each side must reach a
real end point"): optimize from the last frame of each 0.05 side to the nearest minimum, with
frequencies, at the same level as the paths. The distance from each path end (0.10 and 0.05, and the
ORCA 0.10 / 0.05 ends) to that minimum, in kJ/mol and Å (RMSD after alignment), says how far short of
the minimum each side stopped. There is no optimized complex at this level (the endpoints job was
stopped before it ran, 2026-11-11), so this is the only reference.

Level: B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX, SMD water, VeryTightSCF, TightOpt Freq, 2 cores.
Folders: <step_convergence>/true_irc/ends/<reactant|product>/; result: true_irc/ends.json
Usage: python -m pipeline.true_irc_ends
"""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from . import orca
from .geometry import read_xyz, read_xyz_frames
from .irc_angle import st

T = st.OUT / "true_irc"
LEVEL = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX VeryTightSCF CPCM(water)"
SMD = '%cpcm smd true SMDsolvent "WATER" end'


def rmsd(A, B):
    A = A - A.mean(0)
    B = B - B.mean(0)
    U, _, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((A @ R - B) ** 2).sum(1).mean()))


def last(path):
    f = read_xyz_frames(path)[-1]
    return np.asarray(f.coords), float(f.comment.split()[-1])


def opt(side: str) -> dict:
    d = {"reactant": "forward", "product": "backward"}[side]
    x, _ = last(T / "s0050" / d / f"{side}_trj.xyz")
    sym = list(read_xyz(st.ROOT / "01_ts/ts/ts.xyz").symbols)
    geo = "\n".join(f"{a:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}" for a, p in zip(sym, x))
    inp = (f"! {LEVEL} TightOpt Freq\n%pal nprocs 2 end\n%maxcore 3000\n{SMD}\n%geom MaxIter 200 end\n\n"
           f"* xyz -1 1\n{geo}\n*\n")
    wd = T / "ends" / side
    r = orca.run(wd, side, inp, orca.RunLog(T / "ends" / "run_log_ends.md"),
                 note=f"minimum nearest the end of the pysisyphus 0.05 {side} side (completeness check)", front=True)
    xm = np.asarray(read_xyz(wd / f"{side}.xyz").coords) if r.ok else None
    res = {"ok": bool(r.ok), "opt_converged": r.opt_converged, "E_eh": r.energy_eh,
           "imaginary_cm": r.imaginary if r.ok else None, "ends": {}}
    if not (r.ok and r.opt_converged):
        return res
    for name, p in ((f"pysisyphus 0.10", T / "s0100" / d / f"{side}_trj.xyz"),
                    (f"pysisyphus 0.05", T / "s0050" / d / f"{side}_trj.xyz"),
                    ("ORCA 0.10", st.trj("s010", d)), ("ORCA 0.05", st.trj("s005", d))):
        xe, e = last(p)
        res["ends"][name] = {"above_min_kJ": (e - r.energy_eh) * st.H, "rmsd_A": rmsd(xe, xm),
                             "CO": st.dist(xe, st.CO), "CBr": st.dist(xe, st.CBR)}
    res["min_CO"], res["min_CBr"] = st.dist(xm, st.CO), st.dist(xm, st.CBR)
    return res


def main() -> dict:
    (T / "ends").mkdir(exist_ok=True)
    with ThreadPoolExecutor(2) as ex:
        out = dict(zip(("reactant", "product"), ex.map(opt, ("reactant", "product"))))
    (T / "ends.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
