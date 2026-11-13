"""Reports for the pipeline/true_irc2.py runs (J. Picione, 2026-11-12, items 2–4).

tails    item 3: for each tail run (downhill from the same frame of the 0.05 path) and, for reference,
         the same stretch of the earlier pysisyphus 0.10 and 0.05 paths: mass-weighted angle at every
         interior frame (central chord, M^1/2 Δx vs −M^−1/2 g), median / 90th percentile / largest,
         frames above 10°, largest turn between steps; ORCA calls (gradients, Hessians), seconds;
         how the run ended; energy at the end above the optimized minimum (true_irc/ends.json); and
         the largest energy difference from the 0.05 path at the same key bond (kJ/mol).
isotopes item 2: deuterium against normal masses at the same settings, for ORCA (iso_d4 vs s010)
         and pysisyphus (d4_eulerpc vs h_eulerpc): at every deuterium frame, the normal-mass path
         interpolated to the same key bond (C–O reactant side, C–Br product side); RMSD after
         alignment and largest atom displacement (Å), other key bond difference (Å); and the first
         displacement against the normal- and deuterium-mass TS modes.
final    item 4: a final run (both sides from the TS): steps, how each side ended, the stopping rule,
         angles, end against the optimized minima, cost.
Result: <step_convergence>/true_irc2/<what>.json and .md
Usage: python -m pipeline.true_irc2_report tails|isotopes|final <run>
"""
from __future__ import annotations

import json
import re
import sys

import numpy as np

from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, ang, st
from .mass_weight_evidence import read_hess_matrix
from .true_irc_ends import rmsd

T2 = st.OUT / "true_irc2"
T1 = st.OUT / "true_irc"
TS = read_xyz(st.ROOT / "01_ts/ts/ts.xyz")
SYM = list(TS.symbols)
X0 = np.asarray(TS.coords)
D_MASS = 2.01410178
TAIL_START = {"reactant": 80, "product": 170}   # as pipeline/true_irc2.py (that module needs pysisyphus)


def m3(masses=None):
    m = masses if masses is not None else [MASS[a] for a in SYM]
    return np.repeat(np.asarray(m, float), 3)


def load_run(run, d):
    js = json.loads((T2 / run / d / "path.json").read_text(encoding="utf-8"))
    X = [np.asarray(f.coords) for f in read_xyz_frames(T2 / run / d / "path.xyz")]
    G = [np.asarray(g).reshape(-1, 3) for g in js["gradients"]]
    return X, G, js


def angle_rows(X, G, mm, prepend=None):
    sq = np.sqrt(mm)
    P = ([prepend] if prepend is not None else []) + X
    off = 1 if prepend is not None else 0
    rows, s = [], 0.0
    for i in range(1, len(P) - 1):
        a = ((P[i] - P[i - 1]) / BOHR).ravel() * sq
        b = ((P[i + 1] - P[i]) / BOHR).ravel() * sq
        s += float(np.linalg.norm(a))
        g = G[i - off].ravel()
        rows.append({"s_from_start": s, "mw": ang(a + b, -g / sq), "turn": ang(a, b),
                     "CO": st.dist(P[i], st.CO), "CBr": st.dist(P[i], st.CBR), "rms_grad": float(np.sqrt(np.mean(g ** 2)))})
    return rows


def stats(rows):
    v = np.array([r["mw"] for r in rows])
    return {"n": len(rows), "mw_median": float(np.median(v)), "mw_p90": float(np.percentile(v, 90)),
            "mw_max": float(v.max()), "n_over_10deg": int((v > 10).sum()), "max_turn": float(max(r["turn"] for r in rows))}


def ref_curve(side):
    d = "forward" if side == "reactant" else "backward"
    fr = read_xyz_frames(T1 / "s0050" / d / f"{side}_trj.xyz")
    key = st.CO if side == "reactant" else st.CBR
    b = np.array([st.dist(np.asarray(f.coords), key) for f in fr])
    e = np.array([float(f.comment.split()[-1]) for f in fr])
    o = np.argsort(b)
    return b[o], e[o]


def tails():
    ends = json.loads((T1 / "ends.json").read_text(encoding="utf-8"))
    out = {}
    for side, d in (("reactant", "forward"), ("product", "backward")):
        key = st.CO if side == "reactant" else st.CBR
        bref, eref = ref_curve(side)
        e_min = ends[side]["E_eh"]
        # earlier paths over the same stretch (from the start frame's key bond onwards)
        b_start = st.dist(np.asarray(read_xyz_frames(T1 / "s0050" / d / f"{side}_trj.xyz")[TAIL_START[side]].coords), key)
        for tag in ("s0100", "s0050"):
            js = json.loads((T1 / tag / d / "path.json").read_text(encoding="utf-8"))
            X = [np.asarray(f.coords) for f in read_xyz_frames(T1 / tag / d / f"{side}_trj.xyz")]
            G = [np.asarray(g).reshape(-1, 3) for g in js["gradients"]]
            keep = [i for i, x in enumerate(read_xyz_frames(T1 / tag / d / "path.xyz")) if np.max(np.abs(np.asarray(x.coords) - X0)) > 1e-6]
            G = [G[i] for i in keep]
            k0 = next(i for i, x in enumerate(X) if st.dist(x, key) >= b_start)
            rows = angle_rows(X[k0:], G[k0:], m3())
            out[f"{side}: earlier pysisyphus {tag} (same stretch)"] = {**stats(rows)}
        for integ in ("eulerpc", "eulerpc_r5", "gs"):
            run = f"tail_{integ}_{side}"
            if not (T2 / run / "downhill" / "path.json").exists():
                out[f"{side}: {run}"] = None
                continue
            X, G, js = load_run(run, "downhill")
            rows = angle_rows(X, G, m3())
            E = js["energies"]
            bb = np.array([st.dist(x, key) for x in X])
            inside = (bb >= bref[0]) & (bb <= bref[-1])
            dE = (np.array(E)[inside] - np.interp(bb[inside], bref, eref)) * st.H
            walls = [c["wall_s"] for c in js["calls"] if c["wall_s"]]
            out[f"{side}: {run}"] = {**stats(rows), "ended": js["ended"], "steps": js["n_points"] - 1,
                                      "n_engrad": js["n_engrad"], "n_hessian": js["n_hessian"],
                                      "engrad_per_step": js["n_engrad"] / max(1, js["n_points"] - 1),
                                      "s_engrad": float(np.median([c["wall_s"] for c in js["calls"] if c["kind"] == "EnGrad"])),
                                      "s_hessian": (float(np.median([c["wall_s"] for c in js["calls"] if c["kind"] == "Freq"]))
                                                    if js["n_hessian"] else None),
                                      "s_per_step": float(sum(walls)) / max(1, js["n_points"] - 1),
                                      "end_above_min_kJ": (E[-1] - e_min) * st.H,
                                      "end_rmsd_to_min_A": rmsd(X[-1], np.asarray(read_xyz(T1 / "ends" / side / f"{ends[side]['job']}.xyz").coords)),
                                      "max_abs_dE_vs_0.05_path_kJ": float(np.max(np.abs(dE))) if len(dE) else None}
    (T2 / "tails.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = ["# Tail test: mass-weighted angle along the flat tails (downhill from the same frame of the 0.05 path)", "",
         "| Run | n | mw median | mw 90 % | mw largest | frames > 10° | largest turn | ended | gradients / step | Hessians | s per step | end above min (kJ/mol) | ΔE vs 0.05 path (kJ/mol) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, v in out.items():
        if v is None:
            L.append(f"| {k} | not run |")
            continue
        L.append(f"| {k} | {v['n']} | {v['mw_median']:.2f} | {v['mw_p90']:.2f} | {v['mw_max']:.2f} | {v['n_over_10deg']} | "
                 f"{v['max_turn']:.1f} | {v.get('ended', '—')} | {v.get('engrad_per_step', float('nan')):.2f} | {v.get('n_hessian', '—')} | "
                 f"{v.get('s_per_step', float('nan')):.1f} | {v.get('end_above_min_kJ', float('nan')):.3f} | "
                 f"{(v.get('max_abs_dE_vs_0.05_path_kJ') or float('nan')):.3f} |")
    (T2 / "tails.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return out


def interp_path(X, key):
    b = np.array([st.dist(x, key) for x in X])
    o = np.argsort(b)
    b, XX = b[o], np.array(X)[o]

    def at(v):
        if not (b[0] <= v <= b[-1]):
            return None
        j = int(np.clip(np.searchsorted(b, v), 1, len(b) - 1))
        w = (v - b[j - 1]) / (b[j] - b[j - 1]) if b[j] > b[j - 1] else 0.0
        return XX[j - 1] * (1 - w) + XX[j] * w
    return at


def compare_paths(Xd, Xh, side):
    key, other = (st.CO, st.CBR) if side == "reactant" else (st.CBR, st.CO)
    at = interp_path(Xh, key)
    rows = []
    for i, x in enumerate(Xd):
        y = at(st.dist(x, key))
        if y is None:
            continue
        rows.append({"frame": i + 1, "key_bond": st.dist(x, key), "rmsd": rmsd(x, y),
                     "max_atom": float(np.max(np.linalg.norm(x - y, axis=1))), "other_bond_diff": st.dist(x, other) - st.dist(y, other)})
    return rows


def first_step_modes(x1):
    """Angle (degrees) of the first displacement with the normal- and deuterium-mass mass-weighted TS modes
    (both expressed as Cartesian displacements)."""
    H = read_hess_matrix(st.ROOT / "01_ts/ts/ts.hess")
    res = {}
    for lab, masses in (("normal", [MASS[a] for a in SYM]), ("D4", [D_MASS if a == "H" else MASS[a] for a in SYM])):
        mm = m3(masses)
        Hm = H / np.sqrt(np.outer(mm, mm))
        w, V = np.linalg.eigh(Hm)
        v = V[:, 0] / np.sqrt(mm)                       # Cartesian displacement along the mw mode
        dx = (x1 - X0).ravel()
        res[lab] = min(ang(dx, v), ang(dx, -v))
    return res


def isotopes():
    out = {}
    for lab, (dpath, hpath) in {
            "ORCA (Cartesian steepest descent)": (lambda d: st.OUT / "iso_d4" / d / f"irc_{d[0]}_IRC_{d[0].upper()}_trj.xyz",
                                                  lambda d: st.trj("s010", d)),
            "ORCA, D masses also in the Hessian file": (lambda d: st.OUT / "iso_d4_hess" / d / f"irc_{d[0]}_IRC_{d[0].upper()}_trj.xyz",
                                                        lambda d: st.trj("s010", d)),
            "pysisyphus EulerPC (mass-weighted IRC)": (lambda d: T2 / "d4_eulerpc" / d / "path.xyz",
                                                       lambda d: T2 / "h_eulerpc" / d / "path.xyz")}.items():
        for d in ("forward", "backward"):
            if not (dpath(d).exists() and hpath(d).exists()):
                out[f"{lab} {d}"] = None
                continue
            Xd = [x for x in (np.asarray(f.coords) for f in read_xyz_frames(dpath(d))) if np.max(np.abs(x - X0)) > 1e-6]
            Xh = [X0] + [x for x in (np.asarray(f.coords) for f in read_xyz_frames(hpath(d))) if np.max(np.abs(x - X0)) > 1e-6]
            Xd = [Xd[0]] + [c for a, c in zip(Xd, Xd[1:]) if np.max(np.abs(c - a)) > 1e-7]
            side = "reactant" if st.dist(Xd[-1], st.CO) > st.dist(Xd[-1], st.CBR) else "product"
            rows = compare_paths(Xd[1:], Xh, side)       # beyond the first step
            hd = [x for x in (np.asarray(f.coords) for f in read_xyz_frames(hpath(d))) if np.max(np.abs(x - X0)) > 1e-6][0]
            out[f"{lab} {d}"] = {"side": side, "n_frames_D": len(Xd), "n_compared": len(rows),
                                  "first_step_angle_to_modes_deg": {"D path": first_step_modes(Xd[0]), "normal path": first_step_modes(hd)},
                                  "first_point_rmsd_D_vs_normal_A": rmsd(Xd[0], hd),
                                  "rmsd_A": {"median": float(np.median([r["rmsd"] for r in rows])), "max": float(max(r["rmsd"] for r in rows)),
                                             "last": rows[-1]["rmsd"]} if rows else None,
                                  "max_atom_A": float(max(r["max_atom"] for r in rows)) if rows else None,
                                  "other_bond_max_A": float(max(abs(r["other_bond_diff"]) for r in rows)) if rows else None,
                                  "key_bond_range_A": [rows[0]["key_bond"], rows[-1]["key_bond"]] if rows else None,
                                  "rows": rows}
    (T2 / "isotopes.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = ["# Isotope check: deuterium on all four H against normal masses, same settings", "",
         "Beyond the first step, each deuterium frame against the normal-mass path at the same key bond (Å).", "",
         "| Path | Side | frames | RMSD median | RMSD largest | RMSD last | largest atom shift | other key bond | first step vs normal / D4 TS mode (°) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for k, v in out.items():
        if not v or not v["rmsd_A"]:
            L.append(f"| {k} | not available |")
            continue
        f = v["first_step_angle_to_modes_deg"]["D path"]
        L.append(f"| {k} | {v['side']} | {v['n_compared']} | {v['rmsd_A']['median']:.4f} | {v['rmsd_A']['max']:.4f} | {v['rmsd_A']['last']:.4f} | "
                 f"{v['max_atom_A']:.4f} | {v['other_bond_max_A']:.4f} | {f['normal']:.2f} / {f['D4']:.2f} |")
    (T2 / "isotopes.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return out


def final(run):
    ends = json.loads((T1 / "ends.json").read_text(encoding="utf-8"))
    out = {}
    for d in ("forward", "backward"):
        X, G, js = load_run(run, d)
        keep = [i for i, x in enumerate(X) if np.max(np.abs(x - X0)) > 1e-6]
        X, G, E = [X[i] for i in keep], [G[i] for i in keep], [js["energies"][i] for i in keep]
        side = "reactant" if st.dist(X[-1], st.CO) > st.dist(X[-1], st.CBR) else "product"
        rows = angle_rows(X, G, m3(), prepend=X0)
        near = [r for r in rows if r["s_from_start"] <= 1.0]
        far = [r for r in rows if r["s_from_start"] > 1.0]
        steps = [float(np.linalg.norm(((b - a) / BOHR).ravel() * np.sqrt(m3()))) for a, b in zip([X0] + X[:-1], X)]
        xm = np.asarray(read_xyz(T1 / "ends" / side / f"{ends[side]['job']}.xyz").coords)
        walls = [c["wall_s"] for c in js["calls"] if c["wall_s"]]
        out[side] = {"ended": js["ended"], "own_stop": js["own_stop"], "steps": len(X),
                     "step_mw": {"first": steps[0], "median": float(np.median(steps[1:])), "min": float(min(steps[1:])), "max": float(max(steps[1:]))},
                     "s_end": float(sum(steps)), "angles_|s|<=1": stats(near), "angles_|s|>1": stats(far) if far else None,
                     "E_end_rel_TS_kJ": (E[-1] - E_TS()) * st.H, "end_above_min_kJ": (E[-1] - ends[side]["E_eh"]) * st.H,
                     "end_rmsd_to_min_A": rmsd(X[-1], xm), "end_CO": st.dist(X[-1], st.CO), "end_CBr": st.dist(X[-1], st.CBR),
                     "n_engrad": js["n_engrad"], "n_hessian": js["n_hessian"], "s_per_step": float(sum(walls)) / len(X)}
    (T2 / f"final_{run}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return out


def E_TS():
    return st.E_TS


if __name__ == "__main__":
    {"tails": lambda: tails(), "isotopes": lambda: isotopes(), "final": lambda: final(sys.argv[2])}[sys.argv[1]]()
