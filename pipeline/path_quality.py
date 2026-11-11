"""Path-quality rule (decided by J. Picione, 2026-11-11): within the bond-changing region, flag any
IRC step longer than 1.5 times the median step of the path; the card states every flagged step.

Steps are differences of the mass-weighted path length s between consecutive IRC frames
(02_irc/energies.csv). The bond-changing region is the set of frames at which at least one
forming or breaking key bond has covered more than 5 % and less than 95 % of its total change
from the first to the last frame; a step counts if either of its frames is in the region.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

FACTOR = 1.5
LO, HI = 0.05, 0.95


def quality(rx) -> dict:
    rows = list(csv.DictReader(open(rx.root / "02_irc/energies.csv", encoding="utf-8")))
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    frames = [int(r["frame"]) for r in rows]
    if "key_bonds" in rx.cfg:
        kb = [b for b in rx.key_bonds() if b["role"] in ("forms", "breaks")]
        pairs = [tuple(sorted(b["atoms"])) for b in kb]
    else:
        ka = rx.cfg["key_atoms"]
        pairs = [tuple(sorted((ka["center"], ka["nucleophile"]))), tuple(sorted((ka["center"], ka["leaving_group"])))]
    cols = []
    for c in rows[0]:
        if c.startswith("d_") and c.endswith("_A"):
            parts = c[:-2].split("_")
            if tuple(sorted((int(parts[-2]), int(parts[-1])))) in pairs:
                cols.append(c)
    D = np.array([[float(r[c]) for c in cols] for r in rows])
    prog = (D - D[0]) / np.where(np.abs(D[-1] - D[0]) > 1e-9, D[-1] - D[0], 1.0)
    inreg = ((prog > LO) & (prog < HI)).any(axis=1)
    ds = np.diff(s)
    med = float(np.median(ds))
    flagged = []
    for k, d in enumerate(ds):
        if (inreg[k] or inreg[k + 1]) and d > FACTOR * med:
            flagged.append({"from_frame": frames[k], "to_frame": frames[k + 1], "step_s": float(d),
                            "ratio_to_median": float(d / med), "s_from": float(s[k]), "s_to": float(s[k + 1])})
    out = {"median_step_s": med, "factor": FACTOR, "region_frames": [frames[i] for i in np.where(inreg)[0][[0, -1]]]
           if inreg.any() else None, "bonds": cols, "flagged": flagged,
           "first_steps_at_TS": None}
    ts = int(np.argmin(np.abs(s)))
    if 0 < ts < len(s) - 1:
        out["first_steps_at_TS"] = {"backward_s": float(s[ts] - s[ts - 1]), "forward_s": float(s[ts + 1] - s[ts])}
    (rx.root / "02_irc/path_quality.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def card_text(q: dict) -> str:
    """Card sentence: count, the largest step, the two steps at the TS; full list in path_quality.json."""
    base = (f"median step {q['median_step_s']:.3f} bohr·amu½; steps longer than {q['factor']:g} × median in "
            f"the bond-changing region (IRC frames {q['region_frames'][0]}–{q['region_frames'][1]}): ")
    fl = q["flagged"]
    if not fl:
        return base + "none"
    one = lambda f: (f"frames {f['from_frame']}→{f['to_frame']} (s {f['s_from']:+.2f} → {f['s_to']:+.2f}), "
                     f"{f['step_s']:.3f} bohr·amu½ = {f['ratio_to_median']:.1f} × median")
    if len(fl) <= 4:
        return base + "; ".join(one(f) for f in fl)
    big = max(fl, key=lambda f: f["step_s"])
    at_ts = [f for f in fl if abs(f["s_from"]) < 1e-6 or abs(f["s_to"]) < 1e-6]
    txt = base + f"{len(fl)} steps, longest {one(big)}"
    rest = [f for f in at_ts if f is not big]
    if rest:
        txt += "; at the TS also " + "; ".join(one(f) for f in rest)
    return txt + " (full list in `02_irc/path_quality.json`)"
