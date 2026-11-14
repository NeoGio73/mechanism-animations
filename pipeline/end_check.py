"""End check for production paths (J. Picione, 2026-11-14).

A side that ends by the stopping rule more than 0.5 kJ/mol above its nearest optimized minimum
(02_irc/true_irc/angles.json) is accepted, with a note on the card, only if an optimization from its
last frame reaches that minimum. If it reaches a different structure, that is reported as a decision.

Optimization (J. Picione, 2026-11-14: the same convergence settings as the reference minimum): Opt Freq
or TightOpt Freq as in the reference minimum's first job (01_stationary/<name>), at the working level
(reaction.toml method, SMD, VeryTightSCF), MaxIter 200 (a cycle cap, not a convergence setting). If it ends
with an imaginary mode, the procedure used for every M06-2X minimum (pipeline/stationary_m062x.py) follows:
displace 0.1 Å along the lowest mode and re-optimize once with TightOpt Freq. One job per folder in
02_irc/true_irc/end_check/<side>/.
"Reaches that minimum" (J. Picione, 2026-11-14, refined the same day): energy within SAME_E_KJ of it; no imaginary
frequency; RMSD over heavy atoms only (after aligning the heavy atoms) within SAME_RMSD_HEAVY_A; every key
distance of the reaction ([[key_bonds]]) within SAME_KEY_A of the reference minimum. The all-atom RMSD is reported
for information. A side that fails any criterion is reported as a decision with the numbers.

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
SAME_E_KJ = 0.1             # kJ/mol
SAME_RMSD_HEAVY_A = 0.2      # Å, heavy atoms only
SAME_KEY_A = 0.05            # Å, each key distance


def _heavy_rmsd(sym, a, b) -> float:
    h = [i for i, s in enumerate(sym) if s != "H"]
    return rmsd(np.asarray(a)[h], np.asarray(b)[h])


def _ref_opt_keyword(rx_root: Path, name: str) -> str:
    """Opt keyword of the reference minimum's first job (its .inp in 01_stationary/<name>)."""
    inp = rx_root / "01_stationary" / name / f"{name}.inp"
    first = inp.read_text(encoding="utf-8", errors="replace").splitlines()[0].split()
    for k in ("VeryTightOpt", "TightOpt", "LooseOpt", "Opt"):
        if k in first:
            return k
    raise RuntimeError(f"no Opt keyword in {inp}")


def check(rx_root: Path) -> dict:
    rx = Reaction.load(rx_root)
    wd = rx.root / "02_irc/true_irc"
    ang = json.loads((wd / "angles.json").read_text(encoding="utf-8"))
    mins = {m["name"]: m for m in _minima(rx.root)}
    out = {"reaction": rx.root.name, "limit_kJ": LIMIT_KJ,
           "same_minimum": (f"|dE| <= {SAME_E_KJ} kJ/mol, no imaginary frequency, heavy-atom RMSD <= {SAME_RMSD_HEAVY_A} Å, "
                            f"each key distance within {SAME_KEY_A} Å (all-atom RMSD reported for information)"),
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
        kw = _ref_opt_keyword(rx.root, n["minimum"])
        r = orca.run(jd, "opt", rx.inp(f"{kw} Freq", PGeom(sym, x, f"last frame of the {d} side"),
                                       blocks="%geom MaxIter 200 end"), rx.log,
                     note=f"end check: {kw} Freq from the last frame of the {d} side (as the reference minimum)")
        jobs = [{"job": "opt", "keyword": kw, "imaginary_cm": r.imaginary if r.ok else None}]
        last = "opt"
        if r.ok and r.opt_converged and r.imaginary:
            h = orca.read_hess(jd / "opt.hess")
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            gp = read_xyz(jd / "opt.xyz")
            gp = PGeom(gp.symbols, np.asarray(gp.coords) + 0.1 * v / np.linalg.norm(v, axis=1).max(),
                       "displaced 0.1 Å along the lowest imaginary mode")
            r = orca.run(jd, "opt_tight", rx.inp("TightOpt Freq", gp, blocks="%geom MaxIter 200 end"), rx.log,
                         note=f"end check: opt had imaginary mode(s) {[round(f, 1) for f in jobs[0]['imaginary_cm']]} cm-1; "
                              "displaced 0.1 Å along the lowest, TightOpt Freq (as for every M06-2X minimum)")
            jobs.append({"job": "opt_tight", "keyword": "TightOpt", "imaginary_cm": r.imaginary if r.ok else None})
            last = "opt_tight"
        if not (r.ok and r.opt_converged):
            res = {"needed": True, "minimum": n["minimum"], "end_above_min_kJ": n["end_above_min_kJ"],
                   "end_rmsd_A": n["rmsd_A"], "outcome": "INCOMPLETE: optimization did not finish normally",
                   "verdict": "decision"}
        else:
            go = read_xyz(jd / f"{last}.xyz")
            xo = np.asarray(go.coords)
            m = mins[n["minimum"]]
            dE = (r.energy_eh - m["E_eh"]) * HARTREE_KJ
            dR = rmsd(xo, m["xyz"])
            dH = _heavy_rmsd(go.symbols, xo, m["xyz"])
            gm = PGeom(go.symbols, m["xyz"], "reference minimum")
            keys = {b["label"]: {"optimized_A": go.distance(*b["atoms"]), "minimum_A": gm.distance(*b["atoms"])}
                    for b in rx.key_bonds(go.symbols)}
            for k_ in keys.values():
                k_["diff_A"] = k_["optimized_A"] - k_["minimum_A"]
            others = sorted(((k, rmsd(xo, mm["xyz"])) for k, mm in mins.items()), key=lambda t: t[1])
            ok = {"energy": abs(dE) <= SAME_E_KJ, "no imaginary frequency": not r.imaginary,
                  "heavy-atom RMSD": dH <= SAME_RMSD_HEAVY_A,
                  "key distances": all(abs(k_["diff_A"]) <= SAME_KEY_A for k_ in keys.values())}
            failed = [k for k, v in ok.items() if not v]
            outcome = "reaches the same minimum" if not failed else "fails: " + ", ".join(failed)
            verdict = "accepted with a note" if not failed else "decision"
            res = {"needed": True, "minimum": n["minimum"], "end_above_min_kJ": n["end_above_min_kJ"],
                   "end_rmsd_A": n["rmsd_A"], "jobs": jobs, "final_job": last, "opt_E_minus_min_kJ": dE,
                   "opt_rmsd_heavy_A": dH, "opt_rmsd_all_atoms_A": dR, "key_distances": keys,
                   "opt_imaginary_cm": r.imaginary, "nearest_minima_after_opt": others, "criteria": ok,
                   "outcome": outcome, "verdict": verdict}
        out["sides"][d] = res
    out["verdict"] = ("decision" if any(s.get("verdict") == "decision" for s in out["sides"].values())
                      else "ok")
    (wd / "end_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = [f"# End check, production path: {rx.root.name}", "",
         f"A side ending more than {LIMIT_KJ} kJ/mol above its nearest optimized minimum is accepted with a note only if "
         f"Opt Freq from its last frame reaches that minimum ({out['same_minimum']}). Verdict: **{out['verdict']}**", "",
         "| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | "
         "heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for d, s in out["sides"].items():
        if not s.get("needed"):
            e = s.get("end_above_min_kJ")
            L.append(f"| {d} | {s.get('minimum', '—')} | {e:+.3f} | — | no | — | — | — | — | — | within the limit |" if e is not None
                     else f"| {d} | — | — | — | no | — | — | — | — | — | {s.get('why')} |")
        else:
            oe = f"{s['opt_E_minus_min_kJ']:+.3f}" if "opt_E_minus_min_kJ" in s else "—"
            oh = f"{s['opt_rmsd_heavy_A']:.3f}" if "opt_rmsd_heavy_A" in s else "—"
            oa = f"{s['opt_rmsd_all_atoms_A']:.3f}" if "opt_rmsd_all_atoms_A" in s else "—"
            kd = (max(((k, v["diff_A"]) for k, v in s["key_distances"].items()), key=lambda t: abs(t[1]))
                  if s.get("key_distances") else None)
            ks = f"{kd[1]:+.3f} ({kd[0]})" if kd else "—"
            L.append(f"| {d} | {s['minimum']} | {s['end_above_min_kJ']:+.3f} | {s['end_rmsd_A']:.3f} | yes | {oe} | {oh} | {oa} | "
                     f"{ks} | {len(s.get('opt_imaginary_cm') or [])} | {s['outcome']} |")
    (wd / "end_check.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    rx.log.note(f"`02_irc/true_irc/end_check.md`: end check {out['verdict']}: " + "; ".join(
        f"{d}: {s.get('outcome', 'not needed')}" for d, s in out["sides"].items()))
    print("\n".join(L))
    return out


if __name__ == "__main__":
    p = Path(sys.argv[1])
    check(p if p.is_absolute() else ROOT / "reactions" / p)
