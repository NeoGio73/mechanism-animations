"""End check for production paths (J. Picione, 2026-11-14).

A side that ends by the stopping rule more than 0.5 kJ/mol above its nearest optimized minimum
(02_irc/true_irc/angles.json) is accepted, with a note on the card, only if an optimization from its
last frame reaches that minimum. If it reaches a different structure, that is reported as a decision.

Optimization: ORCA Opt Freq at the working level (reaction.toml method, SMD, VeryTightSCF), as the
M06-2X minima were made (pipeline/stationary_m062x.py), MaxIter 200, one job in
02_irc/true_irc/end_check/<side>/. "Reaches that minimum": energy within SAME_E_KJ of it and RMSD after
alignment within SAME_RMSD_A, with no imaginary frequency.

Writes 02_irc/true_irc/end_check.json and .md and a run-log note. Called by true_irc_prod.join after the
angle check; also
Usage: python -m pipeline.end_check <reaction folder>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from . import orca
from .geometry import Geometry as PGeom
from .geometry import read_xyz
from .path_angles import HARTREE_KJ, ROOT, _minima, _xyz_frames, rmsd
from .steps import Reaction

LIMIT_KJ = 0.5                       # J. Picione, 2026-11-14
SAME_E_KJ, SAME_RMSD_A = 0.1, 0.1    # "reaches that minimum"


def check(rx_root: Path) -> dict:
    rx = Reaction.load(rx_root)
    wd = rx.root / "02_irc/true_irc"
    ang = json.loads((wd / "angles.json").read_text(encoding="utf-8"))
    mins = {m["name"]: m for m in _minima(rx.root)}
    out = {"reaction": rx.root.name, "limit_kJ": LIMIT_KJ,
           "same_minimum": f"|dE| <= {SAME_E_KJ} kJ/mol and RMSD <= {SAME_RMSD_A} Å, no imaginary frequency",
           "sides": {}}
    for d, v in ang["sides"].items():
        n = (v.get("end") or {}).get("nearest")
        if n is None or v["ended"] != "own stopping rule":
            out["sides"][d] = {"needed": False, "why": "no nearest minimum or not ended by the stopping rule"}
            continue
        if n["end_above_min_kJ"] <= LIMIT_KJ:
            out["sides"][d] = {"needed": False, "end_above_min_kJ": n["end_above_min_kJ"], "minimum": n["minimum"]}
            continue
        sym, x = _xyz_frames(wd / d / "path.xyz")[-1]
        jd = wd / "end_check" / d
        jd.mkdir(parents=True, exist_ok=True)
        r = orca.run(jd, "opt", rx.inp("Opt Freq", PGeom(sym, x, f"last frame of the {d} side"),
                                       blocks="%geom MaxIter 200 end"), rx.log,
                     note=f"end check: Opt Freq from the last frame of the {d} side")
        if not (r.ok and r.opt_converged):
            res = {"needed": True, "minimum": n["minimum"], "end_above_min_kJ": n["end_above_min_kJ"],
                   "end_rmsd_A": n["rmsd_A"], "outcome": "INCOMPLETE: optimization did not finish normally",
                   "verdict": "decision"}
        else:
            xo = np.asarray(read_xyz(jd / "opt.xyz").coords)
            m = mins[n["minimum"]]
            dE = (r.energy_eh - m["E_eh"]) * HARTREE_KJ
            dR = rmsd(xo, m["xyz"])
            others = sorted(((k, rmsd(xo, mm["xyz"])) for k, mm in mins.items()), key=lambda t: t[1])
            same = abs(dE) <= SAME_E_KJ and dR <= SAME_RMSD_A and not r.imaginary
            res = {"needed": True, "minimum": n["minimum"], "end_above_min_kJ": n["end_above_min_kJ"],
                   "end_rmsd_A": n["rmsd_A"], "opt_E_minus_min_kJ": dE, "opt_rmsd_to_min_A": dR,
                   "opt_imaginary_cm": r.imaginary, "nearest_minima_after_opt": others,
                   "outcome": "reaches the same minimum" if same else "reaches a DIFFERENT structure",
                   "verdict": "accepted with a note" if same else "decision"}
        out["sides"][d] = res
    out["verdict"] = ("decision" if any(s.get("verdict") == "decision" for s in out["sides"].values())
                      else "ok")
    (wd / "end_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = [f"# End check, production path: {rx.root.name}", "",
         f"A side ending more than {LIMIT_KJ} kJ/mol above its nearest optimized minimum is accepted with a note only if "
         f"Opt Freq from its last frame reaches that minimum ({out['same_minimum']}). Verdict: **{out['verdict']}**", "",
         "| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | RMSD (Å) | imaginary | outcome |",
         "|---|---|---|---|---|---|---|---|---|"]
    for d, s in out["sides"].items():
        if not s.get("needed"):
            e = s.get("end_above_min_kJ")
            L.append(f"| {d} | {s.get('minimum', '—')} | {e:+.3f} | — | no | — | — | — | within the limit |" if e is not None
                     else f"| {d} | — | — | — | no | — | — | — | {s.get('why')} |")
        else:
            oe = f"{s['opt_E_minus_min_kJ']:+.3f}" if "opt_E_minus_min_kJ" in s else "—"
            orr = f"{s['opt_rmsd_to_min_A']:.3f}" if "opt_rmsd_to_min_A" in s else "—"
            L.append(f"| {d} | {s['minimum']} | {s['end_above_min_kJ']:+.3f} | {s['end_rmsd_A']:.3f} | yes | {oe} | {orr} | "
                     f"{len(s.get('opt_imaginary_cm') or [])} | {s['outcome']} |")
    (wd / "end_check.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    rx.log.note(f"`02_irc/true_irc/end_check.md`: end check {out['verdict']}: " + "; ".join(
        f"{d}: {s.get('outcome', 'not needed')}" for d, s in out["sides"].items()))
    print("\n".join(L))
    return out


if __name__ == "__main__":
    p = Path(sys.argv[1])
    check(p if p.is_absolute() else ROOT / "reactions" / p)
