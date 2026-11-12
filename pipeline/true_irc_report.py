"""Report on the mass-weighted IRC test (pysisyphus EulerPC) against the ORCA step-test paths
(J. Picione, 2026-11-12).

1. Angle check, both consistent forms, at every interior frame (central chord through the frame):
     plain          angle(Δx, -g)                       0° on a Cartesian steepest-descent path
     mass-weighted  angle(M^1/2 Δx, -M^-1/2 g)          0° on a mass-weighted IRC
     Cartesian form angle(Δx, -M^-1 g)                  0° on a mass-weighted IRC (same condition,
                                                         plain metric; a different number when not 0)
   ORCA paths: the 60 matched points of pipeline/irc_angle.py (gradients from its jobs).
   pysisyphus paths: every frame, with the gradient the IRC computed at that frame.
2. Completeness of each pysisyphus side: how it ended, steps, drop from the TS, RMS gradient at the
   last point, and the project's stopping rule (pipeline/steps.irc_stop_index).
3. Pair tables (same matching as report4): pysisyphus 0.10 vs 0.05, and pysisyphus 0.10 vs ORCA 0.10,
   whole bond-changing region and |s| <= 1.
4. O–H along every path and its relaxed value (O–H + F/k, k from the TS Hessian; pipeline/oh_relaxed.py).
5. Cost: ORCA calls per step and seconds per call (pysisyphus); seconds per step (ORCA IRC).
Result: <step_convergence>/true_irc/report.json and report.md
Usage: python -m pipeline.true_irc_report
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, ang, engrad, points, st
from .oh_relaxed import main as oh_main

T = st.OUT / "true_irc"
TAGS = {"s0100": 0.10, "s0050": 0.05}
H = st.H


def load_side(tag, direction, ts):
    d = T / tag / direction
    js = json.loads((d / "path.json").read_text(encoding="utf-8"))
    fr = read_xyz_frames(d / "path.xyz")
    X = [np.asarray(f.coords) for f in fr]
    G = [np.asarray(g).reshape(-1, 3) for g in js["gradients"]]
    E = [float(f.comment.split()[-1]) for f in fr]
    x_ts = np.asarray(ts.coords)
    keep = [i for i, x in enumerate(X) if np.max(np.abs(x - x_ts)) > 1e-6]
    X, G, E = [X[i] for i in keep], [G[i] for i in keep], [E[i] for i in keep]
    side = "product" if st.dist(X[-1], st.CO) < st.dist(X[-1], st.CBR) else "reactant"
    e0 = float(re.search(r"E = (-?\d+\.\d+) Eh", (d / "run_log_true_irc.md").read_text(encoding="utf-8")).group(1))
    con = (T / f"console_{direction}_{TAGS[tag]:.2f}.txt").read_text(encoding="utf-8", errors="replace")
    why = next((m for m in ("rms(grad) converged!", "Energy increased!", "Energy converged!", "IRC steps exceeded")
                if m in con), "not stated")
    return side, {"X": X, "G": G, "E": E, "E_ts": e0, "json": js, "ended": why, "dir": d}


def chords(X0, X):
    """Central chords through every interior frame of one side (TS prepended)."""
    P = [X0] + X
    return [(P[i], P[i + 1] - P[i - 1]) for i in range(1, len(P) - 1)]


def angles_at(x_chord, g, m):
    dx = (x_chord / BOHR).ravel()
    g = g.ravel()
    return {"plain": ang(dx, -g), "mw": ang(dx * np.sqrt(m), -g / np.sqrt(m)), "cart_Minv": ang(dx, -g / m)}


def write_trj(side_data, path):
    with open(path, "w") as fh:
        for x, e in zip(side_data["X"], side_data["E"]):
            fh.write(f"{len(x)}\ntrue IRC E {e:.12f}\n" + "".join(
                f"{a:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n" for a, p in zip(SYM, x)))


def main() -> dict:
    global SYM
    ts = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
    SYM = list(ts.symbols)
    m = np.repeat([MASS[s] for s in SYM], 3)
    x0 = np.asarray(ts.coords)
    from .steps import Reaction, irc_stop_index
    rx = Reaction.load(st.ROOT)
    k = oh_main()["k_OH_Eh_per_bohr2"]
    out = {"angles": {}, "completeness": {}, "pairs": {}, "oh": {}, "cost": {}}
    P = {}
    for tag in TAGS:
        sides = dict(load_side(tag, d, ts) for d in ("forward", "backward"))
        # angles at every interior frame, both consistent forms
        rows = []
        for side, sd in sides.items():
            sgn = 1 if side == "reactant" else -1
            s = 0.0
            prev = x0
            for (x, ch), g in zip(chords(x0, sd["X"]), sd["G"][:len(sd["X"]) - 1]):
                s += float(np.sqrt((m * ((x - prev) / BOHR).ravel() ** 2).sum()))
                prev = x
                rows.append({"side": side, "s": sgn * s, **angles_at(ch, g, m)})
        out["angles"][tag] = summarize(rows)
        # completeness
        comp = {}
        for side, sd in sides.items():
            g_last = sd["G"][-1].ravel()
            trj = sd["dir"] / f"{side}_trj.xyz"
            write_trj(sd, trj)
            hit = irc_stop_index(rx, read_xyz_frames(trj))
            drop = (sd["E_ts"] - sd["E"][-1]) * H
            rms = float(np.sqrt(np.mean(g_last ** 2)))
            ok = (sd["ended"] == "rms(grad) converged!" and len(sd["X"]) >= 10 and drop >= 1.0) or hit is not None
            comp[side] = {"ended": sd["ended"], "steps": len(sd["X"]), "drop_from_TS_kJ": drop,
                          "rms_grad_last": rms, "stopping_rule_met_at_step": None if hit is None else int(hit[0]),
                          "verdict": "complete" if ok else "INCOMPLETE"}
        out["completeness"][tag] = comp
        # path in the step-test format (energies relative to this run's own TS energy)
        e_ts = sides["reactant"]["E_ts"]
        p = st.path(sides["reactant"]["dir"] / "reactant_trj.xyz", sides["product"]["dir"] / "product_trj.xyz", ts)
        P[tag] = st.region_flags(st.with_energy(p, e_ts), st.adaptive_ends())
        # O–H along the product side, with its relaxed value
        sd = sides["product"]
        oh = []
        for x, g in zip(sd["X"], sd["G"]):
            u = x[3] - x[2]
            r = float(np.linalg.norm(u))
            F = float(-np.dot(g[3] - g[2], u / r) / 2)
            oh.append([round(st.dist(x, st.CBR), 4), round(r, 5), round(r + F / k * BOHR, 5)])
        oh = [o for o in oh if 2.38 <= o[0] <= 3.15]
        out["oh"][tag] = {"CBr_OH_relaxed": oh, "max_abs_dev_A": max(abs(o[1] - o[2]) for o in oh),
                          "OH_range_A": [min(o[1] for o in oh), max(o[1] for o in oh)]}
        # cost
        walls = [w for d in ("forward", "backward") for w in sides_json(tag, d)["wall_s_per_call"] if w]
        steps = sum(len(sd["X"]) for sd in sides.values())
        out["cost"][tag] = {"orca_calls": len(walls), "steps": steps, "calls_per_step": len(walls) / steps,
                            "median_s_per_call": float(np.median(walls)), "s_per_step": float(np.sum(walls) / steps)}
    # ORCA 0.10 path for comparison, and the Cartesian-form angle at its 60 matched points
    P["orca_s010"] = st.load("s010", ts, st.adaptive_ends())
    orca_rows = []
    for pt in points():
        _, g = engrad(pt["path"], pt["side"], pt["k"], pt["sym"], pt["here"])
        orca_rows.append({"path": pt["path"], "side": pt["side"], "s": pt["s"], **angles_at(pt["nxt"] - pt["prev"], g, m)})
    for n in st.R4:
        out["angles"][f"orca_{n}"] = summarize([r for r in orca_rows if r["path"] == n])
    for title, a, b in (("pysisyphus 0.10 vs 0.05", "s0100", "s0050"),
                        ("ORCA 0.10 vs pysisyphus 0.10", "orca_s010", "s0100")):
        res = {}
        for zone, span in (("whole region", None), ("|s| <= 1", (-1.0, 1.0))):
            c = st.compare(P[a], P[b], span=span)
            res[zone] = {sd: {k2: v for k2, v in vv.items() if k2 != "rows"} for sd, vv in c.items()}
            for sd, vv in c.items():
                rows = vv.get("rows", [])
                bad = [r for r in rows if abs(r[2]) > st.TOL_E or abs(r[3]) > st.TOL_B]
                res[zone][sd]["n_failing"] = len(bad)
        out["pairs"][title] = res
    # ORCA cost per step (run log walls)
    out["cost"]["orca_irc"] = orca_cost()
    (T / "report.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    (T / "report.md").write_text(md(out), encoding="utf-8")
    print(md(out))
    return out


def sides_json(tag, d):
    return json.loads((T / tag / d / "path.json").read_text(encoding="utf-8"))


def summarize(rows):
    o = {}
    for zone, sel in (("near_TS_|s|<=1", [r for r in rows if abs(r["s"]) <= 1.0]),
                      ("far_|s|>1", [r for r in rows if abs(r["s"]) > 1.0]), ("all", rows)):
        o[zone] = {"n": len(sel), **{f"{k}_{f}": (float(getattr(np, f)([r[k] for r in sel])) if sel else None)
                                      for k in ("plain", "mw", "cart_Minv") for f in ("median", "max")}}
    return o


def orca_cost():
    t = (st.ROOT / "run_log.md").read_text(encoding="utf-8")
    comp = json.loads((st.OUT / "report4.json").read_text(encoding="utf-8"))["completeness"]
    res = {}
    for n in ("s010", "s005f10"):
        walls = [int(w) for w in re.findall(rf"step_convergence/{n}/(?:forward|backward)/irc_\w` \| `[^`]*` \| (\d+) s", t)][-2:]
        steps = sum(v["steps"] for v in comp[n].values())
        res[n] = {"walls_s_2_cores": walls, "steps": steps, "s_per_step": sum(walls) / steps}
    return res


def md(o):
    L = ["# Mass-weighted IRC test (pysisyphus 1.0.0, EulerPC) on SN2, against the ORCA step-test paths", "",
         "Angles in degrees between the path direction (central chord) and the downhill gradient, median / largest.",
         "plain: Δx vs −g. mass-weighted: M^½Δx vs −M^−½g. Cartesian form: Δx vs −M^−1g.", "",
         "| Path | Zone | n | plain | mass-weighted | Cartesian form |", "|---|---|---|---|---|---|"]
    for p, z in o["angles"].items():
        for zone, v in z.items():
            if zone == "all" or not v["n"]:
                continue
            f = lambda k: f"{v[k + '_median']:.2f} / {v[k + '_max']:.2f}"
            L.append(f"| {p} | {zone} | {v['n']} | {f('plain')} | {f('mw')} | {f('cart_Minv')} |")
    L += ["", "## Completeness (pysisyphus sides)", "", "| Path | Side | Ended | Steps | Drop from TS (kJ/mol) | RMS grad last | Stopping rule met | Verdict |",
          "|---|---|---|---|---|---|---|---|"]
    for t, c in o["completeness"].items():
        for sd, v in c.items():
            L.append(f"| {t} | {sd} | {v['ended']} | {v['steps']} | {v['drop_from_TS_kJ']:.2f} | {v['rms_grad_last']:.1e} | "
                     f"{v['stopping_rule_met_at_step']} | {v['verdict']} |")
    for title, z in o["pairs"].items():
        L += ["", f"## {title}", "", "| Zone | Side | Points | Largest energy difference (kJ/mol) | at | Largest other-bond difference (Å) | Within limits |",
              "|---|---|---|---|---|---|---|"]
        for zone, sides in z.items():
            for sd, v in sides.items():
                if not v.get("n"):
                    L.append(f"| {zone} | {sd} | 0 | — | — | — | — |")
                    continue
                L.append(f"| {zone} | {sd} ({v['match_by']}) | {v['n']} | {v['max_dE_kJ']:+.3f} ({v['max_dE_kJ'] / 4.184:+.4f} kcal/mol) | "
                         f"{v['match_by']} {v['at_dist_A']:.3f} Å | {v['other']} {v['max_dbond_A']:+.4f} | "
                         f"{'yes' if v['pass'] else 'no'} ({v['n_failing']} outside) |")
    L += ["", "## O–H on the product side (Å)", ""]
    for t, v in o["oh"].items():
        L.append(f"- {t}: O–H {v['OH_range_A'][0]:.4f}–{v['OH_range_A'][1]:.4f}, largest distance from the relaxed value {v['max_abs_dev_A']:.4f}")
    L += ["", "## Cost", ""]
    for t, v in o["cost"].items():
        L.append(f"- {t}: {json.dumps(v)}")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
