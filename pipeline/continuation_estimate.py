"""How long would each side take to get within 0.5 kJ/mol of its end minimum if continued from its last frame?
(J. Picione's checkpoint, 2026-11-16: estimate from the tertiary continuation test before starting any continuation;
any side over about a day is reported to him before it starts.)

Model, per side, from its own last frames:
  gap0  = E(last frame) - E(end minimum)
  d0    = energy drop per IRC point over the side's last 10 points (kJ/mol per point)
  harmonic descent along a soft mode at fixed arc-length steps: the drop per point shrinks as sqrt(gap), so
      n_model = 2 * (gap0 - sqrt(gap0 * 0.5)) / d0   points to reach 0.5 kJ/mol above the minimum,
  plus 20 points for the flat steps that end the side.
  Calibration: the same model applied to the tertiary open side as it stood before the test (stopped 4.1 kJ/mol
  short) is compared with the points the continuation actually needed; the ratio multiplies every estimate.
  Before the test has ended the ratio is not known and the estimates are uncalibrated (marked so).
Time: points x the side's own wall time per point (wall_total_s / n_points, Gonzalez–Schlegel micro-iterations
included) + one exact Hessian at the start frame (scaled from the tertiary test: Hessian wall / gradient wall).
Usage: python -m pipeline.continuation_estimate
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HARTREE_KJ = 2625.4996394799
NEAR = 0.5
FLAT = 20

SIDES = [  # (reaction, side, folder under 02_irc, end minimum)
    ("hshift_3_methyl_2_butyl_m062x", "backward", "true_irc/backward", "hshift_methylbutyl/secondary_eclip"),
    ("nucadd_cyanide_acetone_m062x", "backward", "true_irc/backward", "cyanide/reactant_complex"),
    ("e2_ethoxide_2_bromopropane_m062x", "backward", "true_irc/backward", "e2/reactant_complex_irc"),
    ("e2_ethoxide_2_bromopropane_m062x", "forward", "true_irc/forward", "e2/product_complex_irc"),
]
CAL = ("hshift_2_3_dimethyl_2_butyl_m062x", "backward", "attempts/open_side_stopped_short_1218/backward",
       "hshift_tertiary/open_cation")
TEST = "attempts/continuation_test/backward"


def _st():
    return json.loads((ROOT / "docs/stationary_m062x.json").read_text(encoding="utf-8"))


def side_model(rx, folder, minimum, st=None):
    st = st or _st()
    j = json.loads((ROOT / "reactions" / rx / "02_irc" / folder / "path.json").read_text(encoding="utf-8"))
    E = [e * HARTREE_KJ for e in j["energies"]]
    em = st[minimum]["E_eh"] * HARTREE_KJ
    gap0 = E[-1] - em
    k = min(10, len(E) - 1)
    d0 = (E[-1 - k] - E[-1]) / k
    n_model = 2 * (gap0 - math.sqrt(gap0 * NEAR)) / d0 if gap0 > NEAR and d0 > 0 else (0.0 if gap0 <= NEAR else math.inf)
    return {"reaction": rx, "folder": folder, "minimum": minimum, "points": j["n_points"], "gap0_kJ": gap0,
            "drop_per_point_kJ": d0, "n_model": n_model, "wall_per_point_s": j["wall_total_s"] / j["n_points"],
            "wall_per_engrad_s": j["wall_total_s"] / max(1, j["n_engrad"])}


def calibration(st=None):
    """Ratio actual / model points on the tertiary open side, or None while the test runs."""
    rx, _, folder, minimum = CAL
    m = side_model(rx, folder, minimum, st)
    p = ROOT / "reactions" / rx / "02_irc" / TEST / "path.json"
    if not p.exists():
        return m, None, None
    j = json.loads(p.read_text(encoding="utf-8"))
    st = st or _st()
    em = st[minimum]["E_eh"] * HARTREE_KJ
    c = j["continuation"]
    E = [e * HARTREE_KJ for e in j["energies"]]
    new = E[c["earlier_points"] - 1:]
    hit = next((i for i, e in enumerate(new) if e - em <= NEAR), None)
    if hit is None:
        return m, None, {"ended_above_kJ": new[-1] - em, "new_points": len(new) - 1, "test": j}
    ratio = hit / m["n_model"] if m["n_model"] > 0 else None
    return m, ratio, {"points_to_near": hit, "new_points": len(new) - 1, "test": j}


def estimate() -> dict:
    st = _st()
    cal_model, ratio, info = calibration(st)
    hess_ratio = 357 / 40.8  # tertiary test: start-frame Hessian 357 s, one gradient about 41 s (4 cores)
    if info and info.get("test"):
        t = info["test"]
        h = (t["continuation"].get("hessian_at_start") or {})
        if isinstance(h, dict) and h.get("wall_s"):
            hess_ratio = h["wall_s"] / (t["wall_total_s"] / max(1, t["n_engrad"]))
    rows = []
    for rx, d, folder, minimum in SIDES:
        m = side_model(rx, folder, minimum, st)
        n = (m["n_model"] * (ratio if ratio else 1.0)) + FLAT
        h = (n * m["wall_per_point_s"] + hess_ratio * m["wall_per_engrad_s"]) / 3600
        rows.append({**m, "side": d, "points_est": n, "hours_est": h, "over_a_day": h > 24})
    out = {"calibrated": ratio is not None, "ratio_actual_over_model": ratio,
           "calibration_side": {**cal_model, **({k: v for k, v in (info or {}).items() if k != "test"})},
           "hessian_over_gradient": hess_ratio, "sides": rows}
    (ROOT / "docs/continuation_estimate.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    o = estimate()
    c = o["calibration_side"]
    print(f"calibration: tertiary open side gap {c['gap0_kJ']:.2f} kJ/mol, drop {c['drop_per_point_kJ']:.4f}/point, "
          f"model {c['n_model']:.0f} points; ratio {o['ratio_actual_over_model']} ({'calibrated' if o['calibrated'] else 'NOT calibrated yet'})")
    for r in o["sides"]:
        print(f"{r['reaction']} {r['side']}: {r['gap0_kJ']:.2f} kJ/mol above {r['minimum']}, drop {r['drop_per_point_kJ']:.4f} kJ/mol/point, "
              f"{r['wall_per_point_s']:.0f} s/point -> {r['points_est']:.0f} points, {r['hours_est']:.1f} h"
              + ("  OVER A DAY" if r["over_a_day"] else ""))
    sys.exit(0)
