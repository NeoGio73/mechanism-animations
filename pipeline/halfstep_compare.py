"""Half-step check of a production side near its TS (J. Picione, 2026-11-15; first use: tertiary hydride shift,
open side, whose 13 frames next to the TS fail the angle check).

The check run recomputes the side from the TS at half the step, all other settings the same
(`python -m pipeline.true_irc_prod <rx> backward step=0.025 check_dir=attempts/halfstep_check max_cycles=80`).
This script compares the two at equal distance along the path: for each production frame (arc length s from
the TS, mass-weighted, bohr·amu½) the check frame with the nearest s (they differ by at most half a check step;
no interpolation), heavy-atom RMSD between the two frames and the angle check value of each. Nothing in the
production path is changed by this script.
Usage: python -m pipeline.halfstep_compare <reaction> <side> <check_dir> [s_max]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from .path_angles import ROOT, _xyz_frames, rmsd, side_angles


def _side(base: Path, d: str, x_ts, sym):
    js = json.loads((base / d / "path.json").read_text(encoding="utf-8"))
    fr = _xyz_frames(base / d / "path.xyz")
    keep = [k for k, (_, x) in enumerate(fr) if np.max(np.abs(x - x_ts)) > 1e-6]
    X = [fr[k][1] for k in keep]
    G = [np.asarray(js["gradients"][k]) for k in keep]
    return X, side_angles(X, G, x_ts, sym), js


def compare(rx: str, d: str, check_dir: str, s_max: float = 2.0) -> dict:
    R = ROOT / "reactions" / rx
    sym, x_ts = _xyz_frames(R / "01_stationary/ts/ts.xyz")[0]
    heavy = [i for i, s in enumerate(sym) if s != "H"]
    Xp, Ap, jp = _side(R / "02_irc/true_irc", d, x_ts, sym)
    Xc, Ac, jc = _side(R / "02_irc" / check_dir, d, x_ts, sym)
    sc = np.array([r["s_abs"] for r in Ac])
    rows = []
    for r in Ap:
        if r["s_abs"] > s_max:
            break
        j = int(np.argmin(np.abs(sc - r["s_abs"])))
        if abs(sc[j] - r["s_abs"]) > 0.75 * jc["step_bohr_amu05"]:       # beyond the check run's reach
            continue
        i = r["frame_from_ts"] - 1
        rows.append({"prod_frame": i, "s_prod": r["s_abs"], "check_frame": j, "s_check": float(sc[j]),
                     "heavy_rmsd_A": rmsd(Xp[i][heavy], Xc[j][heavy]),
                     "angle_prod_deg": r["mw_deg"], "angle_check_deg": Ac[j]["mw_deg"]})
    out = {"reaction": rx, "side": d, "check_dir": check_dir, "step_prod": jp["step_bohr_amu05"],
           "step_check": jc["step_bohr_amu05"], "s_max": s_max, "rows": rows,
           "max_heavy_rmsd_A": max(r["heavy_rmsd_A"] for r in rows) if rows else None}
    (R / "02_irc" / check_dir / f"compare_{d}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = [f"# Half-step check, {rx}, {d} side", "",
         f"Production step {out['step_prod']} vs check step {out['step_check']} bohr·amu½, compared at equal s "
         f"(nearest check frame, no interpolation) up to s = {s_max} bohr·amu½.", "",
         "| production frame | s (bohr·amu½) | check frame | s check | heavy-atom RMSD (Å) | angle, production | angle, check |",
         "|---|---|---|---|---|---|---|"]
    L += [f"| {r['prod_frame']} | {r['s_prod']:.3f} | {r['check_frame']} | {r['s_check']:.3f} | {r['heavy_rmsd_A']:.4f} | "
          f"{r['angle_prod_deg']:.1f}° | {r['angle_check_deg']:.1f}° |" for r in rows]
    (R / "02_irc" / check_dir / f"compare_{d}.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    o = compare(a[0], a[1], a[2], float(a[3]) if len(a) > 3 else 2.0)
    print(json.dumps({k: v for k, v in o.items() if k != "rows"}, indent=1))
