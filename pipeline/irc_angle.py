"""Gradient check on the four step-test paths (J. Picione, 2026-11-12, "Option B").

At matched points (the frame nearest each target C–O on the reactant side and C–Br on the product side)
compute energy and gradient at the IRC's own level (B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX, TightSCF, SMD
water), alone in its own folder, and the angle between the downhill gradient and the path direction:
  Cartesian      angle(-g, t)               ORCA 6.1 follows steepest descent in Cartesian coordinates
                                            ("default and only option", manual 4.4.2): 0° = on its own path
  mass-weighted  angle(-M^-1/2 g, M^1/2 t)  the IRC as defined (Fukui): 0° = on the true IRC
t = x(i+1) - x(i-1), the central chord through the frame (neighbouring frames on the same path; the TS
for the first frame). Also the forward chord x(i+1) - x(i).
Result: <step_convergence>/angle_check.json and angle_check.md
Usage: python -m pipeline.irc_angle
"""
from __future__ import annotations

import importlib.util
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca
from .geometry import read_xyz, read_xyz_frames

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "reactions/sn2_hydroxide_bromomethane_masvp/00_inputs/step_convergence_test.py"
_s = importlib.util.spec_from_file_location("sct", F)
st = importlib.util.module_from_spec(_s)
_s.loader.exec_module(st)
OUT = st.OUT
W = OUT / "angle_check"
W.mkdir(exist_ok=True)
LOG = orca.RunLog(W / "run_log_angle.md")
LEVEL = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water)"
SMD = '%cpcm smd true SMDsolvent "WATER" end'
MASS = {"H": 1.00794, "C": 12.0107, "O": 15.9994, "Br": 79.904}
TARGETS = {"reactant": ("forward", st.CO, [2.20, 2.25, 2.35, 2.50, 2.70, 2.90]),
           "product": ("backward", st.CBR, [2.40, 2.45, 2.50, 2.55, 2.65, 2.75, 2.85, 2.95, 3.10])}
BOHR = 0.529177210903


def ang(a, b):
    c = float(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, c)))))


def engrad(n, side, k, sym, xyz):
    wd = W / n / f"{side}_{k:04d}"
    pf = wd / "e.property.txt"
    if not (pf.exists() and "ORCA TERMINATED NORMALLY" in (wd / "e.out").read_text(errors="replace")):
        wd.mkdir(parents=True, exist_ok=True)
        geo = "\n".join(f"{s:2s} {x:16.8f} {y:16.8f} {z:16.8f}" for s, (x, y, z) in zip(sym, xyz))
        inp = f"! {LEVEL} EnGrad NoPop\n%pal nprocs 2 end\n%maxcore 3000\n{SMD}\n\n* xyz -1 1\n{geo}\n*\n"
        r = orca.run(wd, "e", inp, LOG, note=f"angle check: {n} {side} frame {k}")
        if not r.ok:
            raise orca.OrcaError(f"{n} {side} {k} did not end normally")
        for t in wd.glob("*.tmp*"):
            t.unlink()
    t = pf.read_text()
    E = float(re.search(r"&FinalEnergy \[&Type \"Double\"\]\s+(\S+)", t).group(1))
    gb = re.search(r"\$SCF_Nuc_Gradient(.*?)\$End", t, re.S).group(1)
    g = np.array([float(v) for v in re.findall(r"^\s*\d+\s+(-?\d\.\d+e[-+]\d+)\s*$", gb, re.M)]).reshape(-1, 3)
    return E, g


def points():
    ts = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    sym = list(ts.symbols)
    m = np.repeat([MASS[s] for s in sym], 3)
    x0 = np.asarray(ts.coords)
    todo = []
    for n in st.R4:
        for side, (d, bond, targets) in TARGETS.items():
            fr = [x0] + [np.asarray(f.coords) for f in read_xyz_frames(st.trj(n, d))]
            # drop exact repeats (ORCA writes the last frame twice when it converges)
            X = [fr[0]] + [c for a, c in zip(fr, fr[1:]) if np.max(np.abs(c - a)) > 1e-7]
            b = np.array([st.dist(c, bond) for c in X])
            s = np.concatenate([[0.0], np.cumsum([np.sqrt((m * ((X[i + 1] - X[i]) / BOHR).ravel() ** 2).sum())
                                                 for i in range(len(X) - 1)])])
            for tgt in targets:
                k = int(np.argmin(np.abs(b[1:-1] - tgt))) + 1      # a frame with neighbours on both sides
                todo.append(dict(path=n, side=side, target=tgt, k=k, bond=float(b[k]),
                                 s=float(s[k]) * (1 if side == "reactant" else -1), sym=sym, m=m,
                                 prev=X[k - 1], here=X[k], nxt=X[k + 1]))
    return todo


def main():
    todo = points()
    with ThreadPoolExecutor(2) as ex:
        eg = list(ex.map(lambda p: engrad(p["path"], p["side"], p["k"], p["sym"], p["here"]), todo))
    rows = []
    for p, (E, g) in zip(todo, eg):
        g = g.ravel()
        m = p["m"]
        tc = ((p["nxt"] - p["prev"]) / BOHR).ravel()
        tf = ((p["nxt"] - p["here"]) / BOHR).ravel()
        rows.append({"path": p["path"], "label": st.LABEL[p["path"]], "side": p["side"], "target_A": p["target"],
                     "frame": p["k"], "bond_A": round(p["bond"], 4), "s": round(p["s"], 3),
                     "E_rel_TS_kJ": (E - st.E_TS) * st.H, "max_grad": float(np.max(np.abs(g))),
                     "cart_deg": ang(-g, tc), "mw_deg": ang(-g / np.sqrt(m), tc * np.sqrt(m)),
                     "cart_fwd_deg": ang(-g, tf), "mw_fwd_deg": ang(-g / np.sqrt(m), tf * np.sqrt(m))})
    summ = {}
    for n in st.R4:
        R = [r for r in rows if r["path"] == n]
        summ[n] = {k: {"median": float(np.median([r[k] for r in R])), "max": float(max(r[k] for r in R))}
                   for k in ("cart_deg", "mw_deg", "cart_fwd_deg", "mw_fwd_deg")}
        near = [r for r in R if abs(r["s"]) <= 1.0]
        summ[n]["near_TS_|s|<=1"] = {k: float(max(r[k] for r in near)) if near else None for k in ("cart_deg", "mw_deg")}
    out = {"rows": rows, "summary": summ}
    (OUT / "angle_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    (OUT / "angle_check.md").write_text(md(out), encoding="utf-8")
    print(md(out), flush=True)
    return out


def md(o):
    L = ["# Gradient check on the four step-test paths (SN2, B3LYP-D3(BJ)/ma-def2-SVP, SMD water)", "",
         "Angle between the downhill gradient and the path direction (central chord through the frame), in degrees.",
         "Cartesian: the path ORCA 6.1 follows (its only coordinate option). Mass-weighted: the IRC as defined.", "",
         "| Path | Cartesian median | Cartesian largest | Mass-weighted median | Mass-weighted largest |",
         "|---|---|---|---|---|"]
    for n, v in o["summary"].items():
        L.append(f"| {st.LABEL[n]} | {v['cart_deg']['median']:.2f} | {v['cart_deg']['max']:.2f} | "
                 f"{v['mw_deg']['median']:.2f} | {v['mw_deg']['max']:.2f} |")
    for n in st.R4:
        L += ["", f"## {st.LABEL[n]}", "",
              "| Side | Matched bond (Å) | s (bohr·amu½) | E rel. TS (kJ/mol) | Cartesian (°) | Mass-weighted (°) |",
              "|---|---|---|---|---|---|"]
        for r in [r for r in o["rows"] if r["path"] == n]:
            bn = "C–O" if r["side"] == "reactant" else "C–Br"
            L.append(f"| {r['side']} | {bn} {r['bond_A']:.4f} | {r['s']:+.3f} | {r['E_rel_TS_kJ']:.3f} | "
                     f"{r['cart_deg']:.2f} | {r['mw_deg']:.2f} |")
    return "\n".join(L) + "\n"


def oh_check() -> dict:
    """Where the gradient leaves the path (no new jobs; the gradients of main()): at every matched point
    the part of the downhill gradient perpendicular to the path, per atom; the O–H length and the force
    stretching it (Eh/bohr, + = the bond wants to lengthen); the kink between the two chords through the
    frame. Also the O–H length along each path's product side (every frame, C–Br 2.38–3.15 Å).
    Result: <step_convergence>/oh_check.json"""
    rows = []
    for p in points():
        E, g = engrad(p["path"], p["side"], p["k"], p["sym"], p["here"])
        x = p["here"]
        t1 = ((x - p["prev"]) / BOHR).ravel()
        t2 = ((p["nxt"] - x) / BOHR).ravel()
        tc = (t1 + t2) / np.linalg.norm(t1 + t2)
        gu = -g.ravel()
        perp = gu - np.dot(gu, tc) * tc
        u = x[3] - x[2]
        r = float(np.linalg.norm(u))
        rows.append({"path": p["path"], "side": p["side"], "target_A": p["target"], "frame": p["k"],
                     "kink_deg": ang(t1, t2), "cart_deg": ang(gu, tc),
                     "g_along": float(np.dot(gu, tc)), "g_perp": float(np.linalg.norm(perp)),
                     "g_perp_by_atom": [float(v) for v in np.linalg.norm(perp.reshape(-1, 3), axis=1)],
                     "OH_A": r, "F_OH": float(-np.dot(g[3] - g[2], u / r) / 2)})
    prof = {}
    from .cc_near_ts import frames
    for n in st.R4:
        _, fr = frames(n)
        P = [f for f in fr if f["side"] == "product" and 2.38 <= f["CBr"] <= 3.15]
        oh = np.array([np.linalg.norm(f["xyz"][3] - f["xyz"][2]) for f in P])
        d = np.diff(oh)
        prof[n] = {"frames": len(P), "OH_min_A": float(oh.min()), "OH_max_A": float(oh.max()),
                   "direction_changes": int(np.sum(np.sign(d[1:]) != np.sign(d[:-1]))),
                   "CBr_OH": [[round(f["CBr"], 4), round(float(o), 5)] for f, o in zip(P, oh)]}
    ts = np.asarray(read_xyz(st.ROOT / "01_ts/ts/ts.xyz").coords)
    out = {"OH_at_TS_A": float(np.linalg.norm(ts[3] - ts[2])), "points": rows, "OH_along_product_side": prof}
    (OUT / "oh_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    import sys
    oh_check() if sys.argv[1:] == ["oh"] else main()
