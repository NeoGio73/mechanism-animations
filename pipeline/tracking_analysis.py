"""Tracking numbers for the manuscript, from files already on disk (J. Picione, 2026-11-25 evening). No ORCA run.

  table          all six reactions: number of tracked orbitals; for each active (changing) orbital the lowest
                 frame-to-frame overlap with its frame and s; swaps, i.e. frames where an orbital's largest overlap
                 with the next frame is not the orbital it was matched to, and near-ties (second-largest overlap
                 above 0.5)  -> docs/results/tracking_table.md/.json
  da_rotation    Diels–Alder, frames whose active-orbital overlap is below 0.9: the 3×3 frame-to-frame overlap
                 matrix of the three active orbitals, its singular values, and the share of each orbital's overlap
                 that stays inside the three (sum of squares over the 3 vs over all 17)
                 -> reactions/da_butadiene_ethylene_m062x/03_orbitals/active_rotation.json
  centroids <rx> basis-independent measure: centroid of each active orbital at every frame from its cube file
                 (r weighted by φ², normalized on the grid), reduced to one coordinate f: the projection of the
                 centroid on the line from the start site to the end site in that frame's geometry (site = mean
                 position of the atoms listed in reaction.toml start/end; 0 = at the start site, 1 = at the end
                 site). Fastest change of f and the frame where f is halfway between its two end values
                 -> <rx>/03_orbitals/centroids.json

Mulliken populations elsewhere in the pipeline (tracking.json, orbital_order.py) are Mulliken populations of the
IBOs in the full ma-def2-SVP basis (C·SC summed over each atom's basis functions), not IAO charges.

Usage: python -X utf8 -m pipeline.tracking_analysis table | da_rotation | centroids <reaction folder name>
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

from .orbital_steps import _fdir, _frames
from .orbitals import load_npz, overlap
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
BOHR = 0.529177210903
SIX = ["sn2_hydroxide_bromomethane_m062x", "e2_ethoxide_2_bromopropane_m062x", "nucadd_cyanide_acetone_m062x",
       "da_butadiene_ethylene_m062x", "hshift_3_methyl_2_butyl_m062x", "hshift_2_3_dimethyl_2_butyl_m062x"]


def _load(name):
    rx = Reaction.load(ROOT / "reactions" / name)
    tr = json.loads((rx.root / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    s = [float(r["s_bohr_amu05"]) for r in csv.DictReader(open(rx.root / "02_irc/energies.csv", encoding="utf-8"))]
    return rx, tr, np.array(s)


def _active(tr):
    ec = tr.get("expected_changes") or {}
    ids = [it["tracked_id"] for it in ec.get("items", []) if it.get("found")]
    lab = {it["tracked_id"]: it["label"] for it in ec.get("items", []) if it.get("found")}
    for t in tr["tracked"]:
        if t["changes"] and t["id"] not in ids:
            ids.append(t["id"])
            lab[t["id"]] = t.get("label") or f"{t['label_start']} → {t['label_end']}"
    return ids, lab


def _pair_matrix(rx, tr, syms, k, ids_a, ids_b):
    """|<orbital i at frame k-1 | orbital j at frame k>| for tracked ids (signs as tracked)."""
    A = load_npz(_fdir(rx, k - 1) / "f_loc.npz", syms)
    B = load_npz(_fdir(rx, k) / "f_loc.npz", syms)
    pa, pb = tr["per_frame"][k - 1], tr["per_frame"][k]
    ia = [pa["orbital_index"][i] for i in ids_a]
    ib = [pb["orbital_index"][j] for j in ids_b]
    M = overlap(A, B, ia, ib)
    sa = np.array([pa["sign"][i] for i in ids_a])[:, None]
    sb = np.array([pb["sign"][j] for j in ids_b])[None, :]
    return M * sa * sb


def table() -> dict:
    out = {}
    for name in SIX:
        rx, tr, s = _load(name)
        syms = _frames(rx)[0].symbols
        m = len(tr["tracked"])
        act, lab = _active(tr)
        rows = []
        for i in act:
            ov = [(tr["per_frame"][k]["overlap_with_previous"][i], k) for k in range(1, len(tr["per_frame"]))]
            v, k = min(ov)
            rows.append({"orbital": lab[i], "tracked_id": i, "lowest_overlap": round(v, 4), "frame": k,
                         "s": round(float(s[k]), 3)})
        swaps, ties = [], []
        allids = list(range(m))
        for k in range(1, len(tr["per_frame"])):
            M = np.abs(_pair_matrix(rx, tr, syms, k, allids, allids))
            for i in allids:
                order = np.argsort(-M[i])
                if order[0] != i:
                    swaps.append({"frame": k, "s": round(float(s[k]), 3), "tracked_id": i,
                                  "matched": round(float(M[i, i]), 4), "best_other": int(order[0]),
                                  "best_other_overlap": round(float(M[i, order[0]]), 4)})
                sec = M[i, order[1]] if order[0] == i else M[i, order[0]]
                if sec > 0.5:
                    ties.append({"frame": k, "tracked_id": i, "second_overlap": round(float(sec), 4)})
        lowest_all = min((min(pf["overlap_with_previous"]), pf["frame"]) for pf in tr["per_frame"][1:])
        out[name] = {"tracked": m, "frames": len(tr["per_frame"]), "active": rows, "swaps": swaps,
                     "near_ties_above_0.5": ties,
                     "lowest_overlap_all_tracked": {"value": round(lowest_all[0], 4), "frame": lowest_all[1],
                                                    "s": round(float(s[lowest_all[1]]), 3)}}
        print(name, m, rows, "swaps", len(swaps), "ties", len(ties), flush=True)
    (ROOT / "docs/results/tracking_table.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                          encoding="utf-8")
    md = ["# IBO tracking along the six final paths", "",
          "Frame-to-frame overlap = |<orbital at frame k−1 | matched orbital at frame k>| (03_orbitals/tracking.json). "
          "A swap is a frame where an orbital's largest overlap with the next frame is with an orbital other than the "
          "one it was matched to. s in bohr·amu½.", "",
          "| Reaction | Tracked orbitals | Active orbital | Lowest frame-to-frame overlap | Frame | s | Swaps |",
          "|---|---|---|---|---|---|---|"]
    for name, r in out.items():
        for j, a in enumerate(r["active"]):
            md.append(f"| {name if j == 0 else ''} | {r['tracked'] if j == 0 else ''} | {a['orbital']} | "
                      f"{a['lowest_overlap']:.3f} | {a['frame']} | {a['s']:+.2f} | "
                      f"{len(r['swaps']) if j == 0 else ''} |")
    md += ["", "Lowest overlap over all tracked orbitals: " + "; ".join(
        f"{n}: {r['lowest_overlap_all_tracked']['value']:.3f} (frame {r['lowest_overlap_all_tracked']['frame']})"
        for n, r in out.items()) + ".",
        "", "Near-ties (another orbital overlapping above 0.5): " + "; ".join(
        f"{n}: {len(r['near_ties_above_0.5'])}" for n, r in out.items()) + "."]
    (ROOT / "docs/results/tracking_table.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return out


def da_rotation() -> dict:
    name = "da_butadiene_ethylene_m062x"
    rx, tr, s = _load(name)
    syms = _frames(rx)[0].symbols
    act, lab = _active(tr)
    m = len(tr["tracked"])
    rows = []
    for k in range(1, len(tr["per_frame"])):
        if min(tr["per_frame"][k]["overlap_with_previous"][i] for i in act) >= 0.9:
            continue
        M3 = _pair_matrix(rx, tr, syms, k, act, act)
        Mall = _pair_matrix(rx, tr, syms, k, act, list(range(m)))
        sv = np.linalg.svd(M3, compute_uv=False)
        rows.append({"frame": k, "s": round(float(s[k]), 3), "orbitals": [lab[i] for i in act],
                     "M": np.round(M3, 4).tolist(), "singular_values": np.round(sv, 5).tolist(),
                     "det": round(float(np.linalg.det(M3)), 5),
                     "share_inside_three": np.round((M3 ** 2).sum(1) / (Mall ** 2).sum(1), 5).tolist()})
        print(rows[-1], flush=True)
    res = {"reaction": name, "rows": rows}
    (rx.root / "03_orbitals/active_rotation.json").write_text(json.dumps(res, indent=1, ensure_ascii=False),
                                                             encoding="utf-8")
    return res


def _cube(p: Path):
    lines = gzip.open(p, "rt").read().splitlines()
    nat = abs(int(lines[2].split()[0]))
    origin = np.array([float(v) for v in lines[2].split()[1:4]])
    dims, axes = [], []
    for i in range(3):
        q = lines[3 + i].split()
        dims.append(int(q[0]))
        axes.append(float(q[1 + i]))
    start = 6 + nat + (1 if int(lines[2].split()[0]) < 0 else 0)
    v = np.array(" ".join(lines[start:]).split(), dtype=np.float64).reshape(dims)
    return v, origin, np.array(axes)


def centroids(name: str) -> dict:
    rx, tr, s = _load(name)
    frames = _frames(rx)
    items = [it for it in tr["expected_changes"]["items"] if it["found"]]
    res = {"reaction": name, "coordinate": "f = projection of the orbital centroid on the line start site → end "
           "site (sites = mean position of the reaction.toml start/end atoms, in that frame's geometry); "
           "0 at the start site, 1 at the end site", "orbitals": {}}
    F = {it["label"]: [] for it in items}
    C = {it["label"]: [] for it in items}
    for k, fr in enumerate(frames):
        X = np.asarray(fr.coords)
        for it in items:
            v, o, a = _cube(rx.root / f"03_orbitals/cubes/f{k:03d}/ibo{it['tracked_id']}.cube.gz")
            w = v * v
            g = [o[d] + a[d] * np.arange(v.shape[d]) for d in range(3)]
            tot = w.sum()
            c = np.array([(w.sum(axis=tuple(x for x in range(3) if x != d)) * g[d]).sum() / tot
                          for d in range(3)]) * BOHR
            ps, pe = X[it["start"]].mean(0), X[it["end"]].mean(0)
            if set(it["start"]) == set(it["end"]):
                raise ValueError("start and end sites coincide")
            f = float((c - ps) @ (pe - ps) / ((pe - ps) @ (pe - ps)))
            F[it["label"]].append(f)
            C[it["label"]].append(c.round(4).tolist())
        if k % 50 == 0:
            print(f"frame {k}", flush=True)
    for lab, f in F.items():
        f = np.array(f)
        rate = np.abs(np.gradient(f, s))
        k = int(rate.argmax())
        w = np.where(rate > 0.5 * rate[k])[0]
        half = 0.5 * (f[0] + f[-1])
        kh = int(np.argmax((f - half) * np.sign(f[-1] - f[0]) >= 0))
        res["orbitals"][lab] = {"f_start": round(float(f[0]), 3), "f_end": round(float(f[-1]), 3),
                                "fastest_frame": k, "fastest_s": round(float(s[k]), 3),
                                "max_rate_per_bohr_amu05": round(float(rate[k]), 3),
                                "half_change_frame": kh, "half_change_s": round(float(s[kh]), 3),
                                "rate_above_half_max_s": [round(float(s[w.min()]), 2), round(float(s[w.max()]), 2)],
                                "f": np.round(f, 4).tolist(), "centroid_A": C[lab]}
        print(lab, {x: y for x, y in res["orbitals"][lab].items() if x not in ("f", "centroid_A")}, flush=True)
    (rx.root / "03_orbitals/centroids.json").write_text(json.dumps(res, indent=1, ensure_ascii=False),
                                                       encoding="utf-8")
    return res


if __name__ == "__main__":
    {"table": lambda: table(), "da_rotation": lambda: da_rotation(),
     "centroids": lambda: centroids(sys.argv[2])}[sys.argv[1]]()
