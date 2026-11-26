"""Hyperconjugation into the cationic centres along the two hydride-shift paths, from the IBOs already on disk
(J. Picione, 2026-11-26, items 6 and 7 and the addition on the new cationic carbon). No ORCA job.

Share of an orbital on an atom = Mulliken population of that atom in the normalized IBO (sums to 1 over atoms), in the
full ma-def2-SVP basis, recomputed from 03_orbitals/frames/fNNN/f_loc.npz with no cut-off, using the frame-to-frame
orbital matching already stored in 03_orbitals/tracking.json. Same caveat as for E2 (docs/results/e2_orbital_measures.md):
these are Mulliken populations, not IAO charges of single orbitals, which ORCA's orca_loc cannot print; they can be
slightly negative. Tracking calls an atom a carrier of an orbital when its share exceeds POP_CUT (0.15), and calls an
orbital changed when its set of carriers differs between the two ends of the path.

Runs at idle priority on one thread so it cannot slow the explicit-water test.
Orbital names in the JSON come from the two atoms with the largest shares at the reactant end, in chemical numbering.
Usage: python -X utf8 -m pipeline.hshift_hyperconj   (writes 03_orbitals/hyperconjugation.json for both reactions;
the reading is in docs/results/hshift_hyperconjugation.md)
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import csv
import ctypes
import json
from pathlib import Path

import numpy as np

from .orbital_steps import POP_CUT, _fdir, _frames
from .orbitals import load_npz
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
SEC = "hshift_3_methyl_2_butyl_m062x"
TER = "hshift_2_3_dimethyl_2_butyl_m062x"


def _idle() -> None:
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)  # IDLE
    except Exception:
        pass


def _pops(name: str):
    rx = Reaction.load(ROOT / "reactions" / name)
    R = rx.root
    tr = json.loads((R / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    syms = _frames(rx)[0].symbols
    P = []
    for k, pf in enumerate(tr["per_frame"]):
        F = load_npz(_fdir(rx, k) / "f_loc.npz", syms)
        P.append(F.populations(pf["orbital_index"]))
    rows = list(csv.DictReader(open(R / "02_irc/energies.csv", encoding="utf-8")))
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    meta = json.loads((R / "02_irc/irc_meta.json").read_text(encoding="utf-8"))
    return rx, tr, syms, np.array(P), s, rows, int(meta["ts_frame"])


# Chemical numbering of the 3-methyl-2-butyl cation, CH3(C1)–C2H+–C3H(CH3)–C4H3, kept along the whole path
# (atom indices 0-based in the xyz; C2 = 1, C3 = 3, migrating H = 9 from reaction.toml [key_atoms]).
# C4 and the 3-methyl carbon are the two methyls on C3; the 3-methyl is the one that leans over the empty p orbital
# of C2 in the cation (card_notes.toml note 3), atom 5.
SEC_C = {0: "C1", 1: "C2", 3: "C3", 4: "C4", 5: "CH₃ on C3"}
# 2,3-dimethyl-2-butyl cation, (CH3)2C2+–C3H(CH3)2: 0 C2, 1 C3, 2 migrating H, 3–4 methyls on C2, 5–6 methyls on C3
TER_C = {0: "C2", 1: "C3", 3: "CH₃ on C2", 4: "CH₃ on C2", 5: "CH₃ on C3", 6: "CH₃ on C3"}


def _name(row: np.ndarray, syms, cmap, bonded_h) -> str:
    """σ label from the two atoms with the largest shares, in chemical numbering."""
    a, b = [int(x) for x in np.argsort(-row)[:2]]
    if syms[a] == "H":
        a, b = b, a
    if syms[b] == "H":
        if a in cmap and cmap[a].startswith("CH₃"):
            return f"σ(C–H) of the {cmap[a]}"
        return f"σ({cmap.get(a, syms[a] + str(a))}–H)" + (" (migrating H)" if b == bonded_h else "")
    na, nb = cmap.get(a, f"C{a}"), cmap.get(b, f"C{b}")
    if na.startswith("CH₃"):
        na, nb = nb, na
    if nb.startswith("CH₃ on "):
        own = nb.removeprefix("CH₃ on ")
        return f"σ({own}–CH₃)" if own == na else f"σ({na}–CH₃ on {own})"
    na, nb = sorted((na, nb))
    return f"σ({na}–{nb})"


def _fast(p: np.ndarray, s: np.ndarray) -> dict:
    g = np.abs(np.gradient(p, s))
    k = int(g.argmax())
    w = np.where(g > 0.5 * g[k])[0]
    tot = abs(p[-1] - p[0])
    half = int(np.argmax(np.abs(p - p[0]) >= 0.5 * tot)) if tot > 0 else None
    return {"fastest_frame": k, "fastest_s": round(float(s[k]), 2), "rate_per_bohr_amu05": round(float(g[k]), 3),
            "fast_stretch_s": [round(float(s[w.min()]), 2), round(float(s[w.max()]), 2)],
            "half_change_s": round(float(s[half]), 2) if half is not None else None}


def secondary() -> dict:
    rx, tr, syms, P, s, rows, ts = _pops(SEC)
    n, m, _ = P.shape
    C2, C3, H = 1, 3, 9
    end = n - 1
    pts = {"reactant end": 0, "TS (C2–C3 rotation TS)": ts, "product end": end}
    names = {j: _name(P[0, j], syms, SEC_C, H) for j in range(m)}
    deloc = {}
    for j, lab in ((6, "σ(C1–H) hyperconjugating into C2"), (14, "σ(C3–CH₃) hyperconjugating into C2")):
        deloc[lab] = {"tracked_id": j,
                      "share_on_C2": {k: round(float(P[i, j, C2]), 3) for k, i in pts.items()},
                      "share_on_C3": {k: round(float(P[i, j, C3]), 3) for k, i in pts.items()},
                      "change_of_share_on_C2": _fast(P[:, j, C2], s),
                      "carriers_reactant_end": tr["tracked"][j]["start"], "carriers_product_end": tr["tracked"][j]["end"]}
    dist = {}
    for c in [c for c in rows[0] if c.startswith("d_")]:
        d = np.array([float(r[c]) for r in rows])
        dist[c] = {**_fast(d, s), "at_reactant_end_A": round(float(d[0]), 3), "at_TS_A": round(float(d[ts]), 3),
                   "at_product_end_A": round(float(d[-1]), 3)}
    # the new cationic carbon C3: orbitals not carried by C3 at either end (bonds next to it, not on it)
    on_c3 = [j for j in range(m) if P[0, j, C3] > POP_CUT or P[-1, j, C3] > POP_CUT]
    gain = []
    for j in range(m):
        if j in on_c3:
            continue
        gain.append({"tracked_id": j, "orbital": names[j],
                     "share_on_C3_reactant_end": round(float(P[0, j, C3]), 3),
                     "share_on_C3_TS": round(float(P[ts, j, C3]), 3),
                     "share_on_C3_product_end": round(float(P[-1, j, C3]), 3),
                     "largest_share_on_C3_anywhere": round(float(P[:, j, C3].max()), 3),
                     "frame_of_largest": int(P[:, j, C3].argmax())})
    gain.sort(key=lambda g: -(g["share_on_C3_product_end"] - g["share_on_C3_reactant_end"]))
    # same view for C2 at the reactant end, for comparison: the orbitals next to C2 there
    on_c2 = [j for j in range(m) if j not in (6, 14) and (P[0, j, C2] > POP_CUT or P[-1, j, C2] > POP_CUT)]
    c2_side = sorted(({"tracked_id": j, "orbital": names[j],
                       "share_on_C2_reactant_end": round(float(P[0, j, C2]), 3),
                       "share_on_C2_product_end": round(float(P[-1, j, C2]), 3)}
                      for j in range(m) if j not in on_c2), key=lambda g: -g["share_on_C2_reactant_end"])
    # orbitals other than the three that change whose set of carriers differs from the reactant-end set at some frame
    # (the comparison's "flagged frames"): which atom crosses the threshold, where, and how far
    at = lambda r: frozenset(int(a) for a in np.where(r > POP_CUT)[0])
    transient = []
    for j in range(m):
        if j in (6, 7, 14):
            continue
        ks = [k for k in range(n) if at(P[k, j]) != at(P[0, j])]
        for a in sorted(set().union(*[at(P[k, j]) ^ at(P[0, j]) for k in ks]) if ks else []):
            kk = int(np.argmax(P[:, j, a]))
            transient.append({"tracked_id": j, "orbital": names[j], "atom": SEC_C.get(a, f"{syms[a]}{a}"),
                              "frames": [ks[0], ks[-1]], "s": [round(float(s[ks[0]]), 2), round(float(s[ks[-1]]), 2)],
                              "largest_share": round(float(P[kk, j, a]), 3), "at_frame": kk,
                              "at_s": round(float(s[kk]), 2), "share_reactant_end": round(float(P[0, j, a]), 3),
                              "share_product_end": round(float(P[-1, j, a]), 3)})
    sigma_ch_c2 = {k: round(float(P[i, 7, C2]), 3) for k, i in pts.items()}
    out = {"reaction": SEC, "ts_frame": ts, "n_frames": n, "threshold_share": POP_CUT,
           "transient_threshold_crossings": transient,
           "migrating_CH_share_on_C2": sigma_ch_c2,
           "energy_kJ": {k: round(float(rows[i]["E_kJ_per_mol"]), 2) for k, i in pts.items()},
           "numbering": "3-methyl-2-butyl cation, CH3(C1)–C2H+–C3H(CH3)–C4H3, kept along the whole path; at the "
                        "product end the cationic carbon is C3 (C2 of 2-methyl-2-butyl in that cation's own name)",
           "points": {k: {"frame": i, "s": round(float(s[i]), 2)} for k, i in pts.items()},
           "delocalized_orbitals": deloc, "distances": dist,
           "next_to_new_cation_C3": gain, "orbitals_on_C3": [names[j] for j in on_c3],
           "next_to_old_cation_C2_reactant_end": c2_side,
           "orbital_names": names}
    (rx.root / "03_orbitals/hyperconjugation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                               encoding="utf-8")
    return out


def tertiary() -> dict:
    rx, tr, syms, P, s, rows, ts = _pops(TER)
    n, m, _ = P.shape
    C2, C3, H = 0, 1, 2
    exp_ids = {e["tracked_id"] for e in tr["expected_changes"]["items"]}
    names = {j: _name(P[0, j], syms, TER_C, H) for j in range(m)}
    rest = []
    for j in range(m):
        if j in exp_ids:
            continue
        carriers0 = set(np.where(P[0, j] > POP_CUT)[0])
        for c, lab in ((C2, "C2"), (C3, "C3")):
            if c in carriers0:
                continue
            rest.append({"tracked_id": j, "orbital": names[j], "atom": lab,
                         "share_reactant_end": round(float(P[0, j, c]), 3), "share_TS": round(float(P[ts, j, c]), 3),
                         "share_product_end": round(float(P[-1, j, c]), 3),
                         "largest_anywhere": round(float(P[:, j, c].max()), 3),
                         "frame_of_largest": int(P[:, j, c].argmax())})
    rest.sort(key=lambda r: -r["largest_anywhere"])
    carriers_change = [j for j in range(m) if j not in exp_ids and any(
        set(np.where(P[k, j] > POP_CUT)[0]) != set(np.where(P[0, j] > POP_CUT)[0]) for k in range(n))]
    out = {"reaction": TER, "ts_frame": ts, "n_frames": n, "threshold_share": POP_CUT,
           "numbering": "2,3-dimethyl-2-butyl cation, (CH3)2C2+–C3H(CH3)2; C2 cationic at the reactant end, C3 at the "
                        "product end (mirror half)",
           "other_orbitals_changing_carriers_at_any_frame": carriers_change,
           "largest_shares_on_C2_or_C3_of_other_orbitals": rest[:8], "orbital_names": names}
    (rx.root / "03_orbitals/hyperconjugation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                               encoding="utf-8")
    return out


if __name__ == "__main__":
    _idle()
    a = secondary()
    b = tertiary()
    print(json.dumps({"secondary": a, "tertiary": b}, indent=1, ensure_ascii=False))
