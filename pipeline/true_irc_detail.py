"""Details behind the mass-weighted IRC test report (J. Picione, 2026-11-12, items 1–3), on the
pysisyphus paths of pipeline/true_irc.py, after pipeline/true_irc_report.py has run.

1. Measured mass-weighted step lengths per side (first, near the TS |s| <= 1, beyond).
2. Mass-weighted angle by frame: percentiles, inside / beyond the bond-changing region (5–95 % of the
   key-bond change, ends from the adaptive IRC), and at frames above 2°: the step-by-step angles
   (incoming step, outgoing step, turn between them) to tell a chord artefact from a path that runs
   off the gradient, and which atoms carry the off-path part of the mass-weighted gradient.
3. End of each side: energy, key bonds, per-step energy change over the last 20 steps, and the
   distance to the other paths' ends.
4. O–H along both sides of every path (pysisyphus with its relaxed value; ORCA 0.10 and 0.05 raw).
Result: <step_convergence>/true_irc/detail.json
Usage: python -m pipeline.true_irc_detail
"""
from __future__ import annotations

import json

import numpy as np

from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, ang, st
from .oh_relaxed import main as oh_main
from .true_irc_ends import rmsd
from .true_irc_report import TAGS, T, chords, load_side

WIN = {"reactant": ("CO", st.CO, 2.20, 3.22), "product": ("CBr", st.CBR, 2.38, 3.51)}


def side_detail(sd, x0, sym, m3, ends, k):
    sq = np.sqrt(m3)
    X, G, E = sd["X"], sd["G"], sd["E"]
    P = [x0] + X
    steps = [float(np.linalg.norm(((P[i + 1] - P[i]) / BOHR).ravel() * sq)) for i in range(len(X))]
    s = np.cumsum(steps)
    (co0, cb0), (co1, cb1) = ends
    rows = []
    for i, ((x, ch), g) in enumerate(zip(chords(x0, X), G[:len(X) - 1])):
        t = (ch / BOHR).ravel() * sq
        gm = -g.ravel() / sq
        a = ((P[i + 1] - P[i]) / BOHR).ravel() * sq
        b = ((P[i + 2] - P[i + 1]) / BOHR).ravel() * sq
        co, cb = st.dist(x, st.CO), st.dist(x, st.CBR)
        inreg = any(0.05 <= p <= 0.95 for p in ((co - co0) / (co1 - co0), (cb - cb0) / (cb1 - cb0)))
        row = {"frame": i + 1, "s": float(s[i]), "in_region": inreg, "plain": ang((ch / BOHR).ravel(), -g.ravel()),
               "mw": ang(t, gm), "mw_in_step": ang(a, gm), "mw_out_step": ang(b, gm), "turn": ang(a, b),
               "rms_grad": float(np.sqrt(np.mean(g ** 2))), "CO": co, "CBr": cb}
        if row["mw"] > 2.0:
            tu = t / np.linalg.norm(t)
            perp = (gm - np.dot(gm, tu) * tu).reshape(-1, 3)
            fr = (perp ** 2).sum(1) / (perp ** 2).sum()
            row["off_path_by_atom"] = {f"{sym[j]}{j}": round(float(f), 2) for j, f in enumerate(fr)}
        rows.append(row)
    mw = np.array([r["mw"] for r in rows])

    def zone(rr):
        if not rr:
            return None
        v = np.array([r["mw"] for r in rr])
        return {"n": len(rr), "s": [rr[0]["s"], rr[-1]["s"]], "mw_median": float(np.median(v)),
                "mw_max": float(v.max()), "plain_median": float(np.median([r["plain"] for r in rr]))}
    nm, ij, lo, hi = WIN[sd["side"]]
    oh = []
    for x, g in zip(X, G):
        if lo <= st.dist(x, ij) <= hi:
            u = x[3] - x[2]
            r = float(np.linalg.norm(u))
            F = float(-np.dot(g[3] - g[2], u / r) / 2)
            oh.append((r, abs(F / k * BOHR)))
    return {
        "steps_mw": {"first": steps[0], "median_|s|<=1": float(np.median([a for a, b in zip(steps, s) if b <= 1])),
                     "median_|s|>1": float(np.median([a for a, b in zip(steps, s) if b > 1])),
                     "min": float(min(steps[1:])), "max": float(max(steps[1:]))},
        "s_end": float(s[-1]),
        "mw_percentiles": {p: float(np.percentile(mw, p)) for p in (50, 90, 95, 99)},
        "first_frame_mw": rows[0]["mw"],
        "in_region_excluding_first": zone([r for r in rows[1:] if r["in_region"]]),
        "beyond_region": zone([r for r in rows if not r["in_region"]]),
        "max_turn": float(max(r["turn"] for r in rows)),
        "frames_over_2deg": [r for r in rows if r["mw"] > 2.0],
        "end": {"E_rel_TS_kJ": (E[-1] - sd["E_ts"]) * st.H, "CO": st.dist(X[-1], st.CO), "CBr": st.dist(X[-1], st.CBR),
                "last20_dE_kJ": [(E[j] - E[j - 1]) * st.H for j in range(len(X) - 20, len(X))]},
        "OH": {"window": f"{nm} {lo}–{hi} Å", "min": min(o[0] for o in oh), "max": max(o[0] for o in oh),
               "max_dev_from_relaxed": max(o[1] for o in oh)},
    }


def main() -> dict:
    ts = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    sym = list(ts.symbols)
    m3 = np.repeat([MASS[a] for a in sym], 3)
    x0 = np.asarray(ts.coords)
    ends = st.adaptive_ends()
    k = oh_main()["k_OH_Eh_per_bohr2"]
    out, last = {}, {}
    for tag in TAGS:
        for d in ("forward", "backward"):
            side, sd = load_side(tag, d, ts)
            sd["side"] = side
            key = f"pysisyphus {TAGS[tag]:.2f} {side}"
            out[key] = side_detail(sd, x0, sym, m3, ends, k)
            last[key] = sd["X"][-1]
    for name in ("s010", "s005"):
        for d, side in (("forward", "reactant"), ("backward", "product")):
            fr = read_xyz_frames(st.trj(name, d))
            c = np.asarray(fr[-1].coords)
            nm, ij, lo, hi = WIN[side]
            r = [float(np.linalg.norm(np.asarray(f.coords)[3] - np.asarray(f.coords)[2]))
                 for f in fr if lo <= st.dist(np.asarray(f.coords), ij) <= hi]
            key = f"ORCA {name} {side}"
            out[key] = {"end": {"E_rel_TS_kJ": (float(fr[-1].comment.split()[-1]) - st.E_TS) * st.H,
                                "CO": st.dist(c, st.CO), "CBr": st.dist(c, st.CBR)},
                        "OH": {"window": f"{nm} {lo}–{hi} Å", "min": min(r), "max": max(r)}}
            last[key] = c
    for side in ("reactant", "product"):
        ref = f"pysisyphus 0.05 {side}"
        for key in last:
            if key.endswith(side):
                out[key]["end"]["rmsd_to_pysis005_end_A"] = rmsd(last[key], last[ref])
                out[key]["end"]["E_minus_pysis005_end_kJ"] = out[key]["end"]["E_rel_TS_kJ"] - out[ref]["end"]["E_rel_TS_kJ"]
    (T / "detail.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for key, v in out.items():
        print(key, json.dumps({a: b for a, b in v.items() if a != "frames_over_2deg"}, default=str)[:900])
    return out


if __name__ == "__main__":
    main()
