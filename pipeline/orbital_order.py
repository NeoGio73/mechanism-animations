"""Where along the path each expected orbital change happens fastest, against the bond distances
(J. Picione, 2026-11-25, E2 point 2: does the orbital order match "proton first, double bond next, bromide last"?).

Orbitals: Mulliken atomic populations of each tracked IBO of reaction.toml [[expected_orbital_changes]], recomputed
from the localized orbitals at every frame (03_orbitals/frames/fNNN/f_loc.npz) with no cut-off (tracking.json keeps
only populations above a cut-off, which makes false jumps where an atom crosses it). Rate of change at a frame =
sum over atoms of |dp/ds| (central differences in s), in e per bohr·amu½. Reported: the frame of the largest rate,
the frame where half of the total end-to-end change is reached, and the s range where the rate is above half its
maximum. Distances: the frame where |d(distance)/ds| is largest, for each distance column of 02_irc/energies.csv.

Usage: python -X utf8 -m pipeline.orbital_order <reaction folder name>   (writes 03_orbitals/orbital_order.json)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

from .orbital_steps import _fdir, _frames
from .orbitals import load_npz
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]


def run(name: str) -> dict:
    rx = Reaction.load(ROOT / "reactions" / name)
    R = rx.root
    tr = json.loads((R / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    rows = list(csv.DictReader(open(R / "02_irc/energies.csv", encoding="utf-8")))
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    syms = _frames(rx)[0].symbols
    items = [it for it in tr["expected_changes"]["items"] if it["found"]]
    P = {it["label"]: [] for it in items}
    for k, pf in enumerate(tr["per_frame"]):
        F = load_npz(_fdir(rx, k) / "f_loc.npz", syms)
        for it in items:
            P[it["label"]].append(F.populations([pf["orbital_index"][it["tracked_id"]]])[0])
    out = {"reaction": name, "ts_frame": int(json.loads((R / "02_irc/irc_meta.json").read_text())["ts_frame"]),
           "units": "s in bohr·amu½; rate in e per bohr·amu½", "orbitals": {}, "distances": {}}
    for lab, p in P.items():
        p = np.array(p)
        rate = np.abs(np.gradient(p, s, axis=0)).sum(1)
        k = int(rate.argmax())
        w = np.where(rate > 0.5 * rate[k])[0]
        tot = np.abs(p - p[0]).sum(1)
        half = int(np.argmax(tot >= 0.5 * tot[-1]))
        main = [int(j) for j in np.argsort(-np.abs(p[-1] - p[0]))[:2]]
        out["orbitals"][lab] = {
            "fastest_frame": k, "fastest_s": round(float(s[k]), 3), "max_rate": round(float(rate[k]), 3),
            "half_change_frame": half, "half_change_s": round(float(s[half]), 3),
            "rate_above_half_max_s": [round(float(s[w.min()]), 2), round(float(s[w.max()]), 2)],
            "main_atoms": {f"{syms[j]}{j}": [round(float(p[0, j]), 3), round(float(p[k, j]), 3),
                                              round(float(p[-1, j]), 3)] for j in main},
            "main_atoms_note": "population at reactant end, at the fastest frame, at product end"}
    for c in [c for c in rows[0] if c.startswith("d_")]:
        d = np.array([float(r[c]) for r in rows])
        g = np.abs(np.gradient(d, s))
        k = int(g.argmax())
        w = np.where(g > 0.5 * g[k])[0]
        out["distances"][c] = {"fastest_frame": k, "fastest_s": round(float(s[k]), 3),
                               "distance_there_A": round(float(d[k]), 3),
                               "rate_above_half_max_s": [round(float(s[w.min()]), 2), round(float(s[w.max()]), 2)]}
    (R / "03_orbitals/orbital_order.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


if __name__ == "__main__":
    run(sys.argv[1])
