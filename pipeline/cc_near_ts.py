"""Coupled-cluster maximum near the TS on the four step-test paths (J. Picione, 2026-11-12).

DLPNO-CCSD(T) TightPNO/ma-def2-TZVP, SMD water (the per-frame path level of the SN2 pilot), 4 cores,
each job alone in its own folder, at frames of the B3LYP-D3(BJ)/ma-def2-SVP step-test paths:
  stage 1  the TS and every frame of the 0.10 path with -1.2 <= s <= +0.4 bohr·amu½
  stage 2  on the other three paths, every frame within 0.25 bohr·amu½ of the 0.10 path's maximum
           (located on each path by the same C–Br), plus the matched near-TS points of the gradient
           check (|s| <= 1) on all four paths
Per path: the highest frame (s, C–Br, C–O, energy) and a parabola through it and its two neighbours
(position and height of the maximum between frames). Energies relative to the TS geometry at the same
level. s: bohr·amu½ from the TS along each path, + toward the reactant complex, - toward products.
Result: <step_convergence>/cc_near_ts.json and cc_near_ts.md
Usage: python -m pipeline.cc_near_ts
"""
from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from . import orca
from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, TARGETS, st

LEVEL = "DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)"
SMD = '%cpcm smd true SMDsolvent "WATER" end'
W = st.OUT / "cc_near_ts"
W.mkdir(exist_ok=True)
LOG = orca.RunLog(W / "run_log_cc.md")
H = 2625.499639


def frames(n):
    """TS-ordered frames of one path: list of dict(s, CBr, CO, xyz, side, k)."""
    ts = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    sym = list(ts.symbols)
    m = np.repeat([MASS[x] for x in sym], 3)
    x0 = np.asarray(ts.coords)
    out = [dict(s=0.0, CBr=st.dist(x0, st.CBR), CO=st.dist(x0, st.CO), xyz=x0, side="TS", k=0)]
    for side, (d, _, _) in TARGETS.items():
        X = [x0] + [np.asarray(f.coords) for f in read_xyz_frames(st.trj(n, d))]
        X = [X[0]] + [c for a, c in zip(X, X[1:]) if np.max(np.abs(c - a)) > 1e-7]
        s = 0.0
        for k in range(1, len(X)):
            s += float(np.sqrt((m * ((X[k] - X[k - 1]) / BOHR).ravel() ** 2).sum()))
            out.append(dict(s=s if side == "reactant" else -s, CBr=st.dist(X[k], st.CBR), CO=st.dist(X[k], st.CO),
                            xyz=X[k], side=side, k=k))
    return sym, sorted(out, key=lambda r: r["s"])


def sp(n, f, sym):
    wd = W / ("TS" if f["side"] == "TS" else f"{n}/{f['side']}_{f['k']:04d}")
    o = wd / "f.out"
    if not (o.exists() and "ORCA TERMINATED NORMALLY" in o.read_text(errors="replace")):
        wd.mkdir(parents=True, exist_ok=True)
        geo = "\n".join(f"{a:2s} {x:16.8f} {y:16.8f} {z:16.8f}" for a, (x, y, z) in zip(sym, f["xyz"]))
        inp = f"! {LEVEL}\n%pal nprocs 4 end\n%maxcore 3000\n{SMD}\n\n* xyz -1 1\n{geo}\n*\n"
        r = orca.run(wd, "f", inp, LOG, note=f"CC near TS: {n} {f['side']} frame {f['k']}")
        if not r.ok:
            raise orca.OrcaError(f"{wd} did not end normally")
        for t in wd.glob("*.tmp*"):
            t.unlink()
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", o.read_text(errors="replace"))[-1])


def run(jobs):
    with ThreadPoolExecutor(2) as ex:
        return list(ex.map(lambda a: sp(*a), jobs))


def peak(rows):
    i = int(np.argmax([r["E"] for r in rows]))
    top = rows[i]
    fit = None
    if 0 < i < len(rows) - 1:
        s = np.array([rows[j]["s"] for j in (i - 1, i, i + 1)])
        e = np.array([rows[j]["E"] for j in (i - 1, i, i + 1)])
        c = np.polyfit(s, e, 2)
        if c[0] < 0:
            sm = -c[1] / (2 * c[0])
            fit = {"s": float(sm), "E_kJ": float(np.polyval(c, sm)),
                   "CBr_A": float(np.interp(sm, s, [rows[j]["CBr"] for j in (i - 1, i, i + 1)])),
                   "CO_A": float(np.interp(sm, s, [rows[j]["CO"] for j in (i - 1, i, i + 1)]))}
    return {"frame": {k: top[k] for k in ("side", "k", "s", "CBr", "CO", "E")}, "parabola": fit,
            "neighbours_computed": 0 < i < len(rows) - 1}


def main():
    F = {n: frames(n) for n in st.R4}
    sym = F["s010"][0]
    E_ts = run([("TS", F["s010"][1][[r["side"] for r in F["s010"][1]].index("TS")], sym)])[0]
    done = {n: {} for n in st.R4}

    def compute(n, sel):
        todo = [f for f in sel if f["side"] != "TS" and (f["side"], f["k"]) not in done[n]]
        for f, e in zip(todo, run([(n, f, sym) for f in todo])):
            done[n][(f["side"], f["k"])] = (e - E_ts) * H

    # stage 1
    compute("s010", [f for f in F["s010"][1] if -1.2 <= f["s"] <= 0.4])
    rows10 = [dict(f, E=done["s010"].get((f["side"], f["k"]), 0.0 if f["side"] == "TS" else None))
              for f in F["s010"][1]]
    rows10 = [r for r in rows10 if r["E"] is not None]
    cbr_star = peak(rows10)["parabola"]["CBr_A"] if peak(rows10)["parabola"] else peak(rows10)["frame"]["CBr"]
    # stage 2 (and the matched near-TS points of the gradient check on all four paths)
    for n in st.R4:
        fr = F[n][1]
        prod = sorted([f for f in fr if f["side"] in ("product", "TS")], key=lambda f: f["CBr"])
        s_star = float(np.interp(cbr_star, [f["CBr"] for f in prod], [f["s"] for f in prod]))
        sel = [f for f in fr if abs(f["s"] - s_star) <= 0.25]
        for side, (_, bond, targets) in TARGETS.items():
            sf = [f for f in fr if f["side"] == side]
            for t in targets:
                f = min(sf, key=lambda f: abs((f["CO"] if side == "reactant" else f["CBr"]) - t))
                if abs(f["s"]) <= 1.0:
                    sel.append(f)
        compute(n, sel)
    res = {"level": LEVEL, "E_TS_Eh": E_ts, "CBr_of_max_on_0.10_path_A": cbr_star, "paths": {}}
    for n in st.R4:
        rows = [dict(f, E=0.0 if f["side"] == "TS" else done[n].get((f["side"], f["k"])))
                for f in F[n][1]]
        rows = [r for r in rows if r["E"] is not None]
        res["paths"][n] = {"label": st.LABEL[n], "peak": peak(rows),
                           "computed": [{k: r[k] for k in ("side", "k", "s", "CBr", "CO", "E")} for r in rows]}
    (st.OUT / "cc_near_ts.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    (st.OUT / "cc_near_ts.md").write_text(md(res), encoding="utf-8")
    print(md(res), flush=True)
    return res


def md(o):
    L = ["# Coupled-cluster maximum near the TS, four step-test paths", "",
         f"DLPNO-CCSD(T) TightPNO/ma-def2-TZVP, SMD water, at B3LYP-D3(BJ)/ma-def2-SVP path frames. "
         "Energies relative to the TS geometry at the same level.", "",
         "| Path | Highest frame: s (bohr·amu½) | C–Br / C–O (Å) | E (kJ/mol) | Parabola: s | E (kJ/mol) | C–Br (Å) |",
         "|---|---|---|---|---|---|---|"]
    for n, v in o["paths"].items():
        f, p = v["peak"]["frame"], v["peak"]["parabola"]
        L.append(f"| {v['label']} | {f['s']:+.3f} | {f['CBr']:.4f} / {f['CO']:.4f} | {f['E']:+.3f} ({f['E'] / 4.184:+.4f} kcal/mol) | "
                 + (f"{p['s']:+.3f} | {p['E_kJ']:+.3f} ({p['E_kJ'] / 4.184:+.4f} kcal/mol) | {p['CBr_A']:.4f} |" if p else "— | — | — |"))
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
