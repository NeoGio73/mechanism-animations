"""End the tertiary continuation test early (J. Picione, 2026-11-16 afternoon): once the rerun from the TS has at least
50 frames past the old side's last frame, compare it with the continuation (criteria fixed beforehand, per frame),
stop both runs and keep their frames.

Both runs started before the clean pause existed, so they are stopped by ending their processes right after
pysisyphus has written its data file (<side>/pysis/<kind>_irc_data.h5: coordinates, energies and gradients of every
IRC point so far). This step turns those files into the usual side records (path.json, path.xyz), so the comparison,
the join and the angle check read them as any other side:
  rerun         02_irc/true_irc/backward                 its points, TS-outward
  continuation  02_irc/attempts/continuation_test/backward  the old side's points, then the continuation's
Points pysisyphus had not yet written to the data file at the stop are left out (named in the record).
Usage: python -m pipeline.test_stop materialize
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

from .geometry import Geometry

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions/hshift_2_3_dimethyl_2_butyl_m062x"
BOHR = 0.529177210903
OLD = "attempts/open_side_stopped_short_1218/backward"
STOP_NOTE = ("stopped at the early comparison (J. Picione, 2026-11-16: compare once the rerun has 50 frames past the "
             "old side's last frame, then stop both runs and keep their frames)")


def _h5(p: Path):
    with h5py.File(p, "r") as h:
        return list(h["atoms"][:]), h["coords"][:], h["energies"][:], h["gradients"][:]


def _write(wd: Path, sym, C_bohr, E, G, rec: dict):
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C_bohr, E):
            fh.write(Geometry(sym, np.asarray(c).reshape(-1, 3) * BOHR, "").to_xyz(f"true IRC backward E {e:.12f}"))
    rec.update({"n_points": len(E), "energies": [float(e) for e in E], "gradients": [np.asarray(g).tolist() for g in G]})
    (wd / "path.json").write_text(json.dumps(rec), encoding="utf-8")


def _syms(raw):
    return [(a.decode() if isinstance(a, bytes) else str(a)).capitalize() for a in raw]


def _wall(wd: Path):
    g = sorted(wd.glob("g*/j.out")) + sorted(wd.glob("h*/j.out"))
    t = [p.stat().st_mtime for p in g]
    return (max(t) - min(t)) if t else None, len(list(wd.glob("g*"))), len(list(wd.glob("h*")))


def materialize() -> dict:
    out = {}
    tso = json.loads((ROOT / "docs/stationary_m062x.json").read_text(encoding="utf-8"))["hshift_tertiary/open_cation"]
    # rerun from the TS
    wd = RX / "02_irc/true_irc/backward"
    sym, C, E, G = _h5(wd / "pysis/backward_irc_data.h5")
    sym = _syms(sym)
    trj_n = (wd / "pysis/backward_irc.trj").read_text().count("backward IRC, step")
    wall, ng, nh = _wall(wd)
    rec = {"reaction": RX.name, "direction": "backward", "method": "as 02_irc/true_irc/forward (production settings, "
           "step 0.05 bohr·amu½), stopping rule (c) with hshift_tertiary/open_cation", "ended": STOP_NOTE,
           "check_run": None, "continuation": None, "own_stop": None,
           "end_min": {"minimum": "hshift_tertiary/open_cation", "E_eh": tso["E_eh"], "near_min_kJ": 0.5},
           "n_engrad": ng, "n_hessian": nh, "wall_total_s": wall, "points_in_trajectory_at_stop": trj_n,
           "points_not_kept": trj_n - len(E), "source": "pysis/backward_irc_data.h5"}
    _write(wd, sym, C, E, G, rec)
    out["rerun"] = {"points": len(E), "E_last_above_min_kJ": (E[-1] - tso["E_eh"]) * 2625.4996394799, "not_kept": trj_n - len(E)}
    # continuation
    wd = RX / "02_irc/attempts/continuation_test/backward"
    od = RX / "02_irc" / OLD
    oj = json.loads((od / "path.json").read_text(encoding="utf-8"))
    from .path_angles import _xyz_frames
    ofr = _xyz_frames(od / "path.xyz")
    oC = [x.ravel() / BOHR for _, x in ofr]
    sym2, C2, E2, G2 = _h5(wd / "pysis/downhill_irc_data.h5")
    trj_n = (wd / "pysis/downhill_irc.trj").read_text().count("IRC, step")
    wall, ng, nh = _wall(wd)
    cont = {"from": OLD, "earlier_points": len(oC), "new_points": len(E2) - 1,
            "start_frame_energy_recomputed_minus_earlier_kJ": (E2[0] - oj["energies"][-1]) * 2625.499639,
            "start_frame_max_coordinate_difference_bohr": float(np.max(np.abs(C2[0] - oC[-1]))),
            "hessian_at_start": "exact (ORCA Freq at the start frame, h0000), then Bofill updates",
            "earlier_n_engrad": oj.get("n_engrad"), "earlier_wall_s": oj.get("wall_total_s")}
    rec = {"reaction": RX.name, "direction": "backward", "method": "continuation downhill from the last frame of " + OLD
           + " (production settings, step 0.05 bohr·amu½, exact Hessian at the start frame)", "ended": STOP_NOTE,
           "check_run": "attempts/continuation_test", "continuation": cont, "own_stop": None,
           "end_min": {"minimum": "hshift_tertiary/open_cation", "E_eh": tso["E_eh"], "near_min_kJ": 0.5},
           "n_engrad": ng, "n_hessian": nh, "wall_total_s": wall, "points_in_trajectory_at_stop": trj_n,
           "points_not_kept": trj_n - len(E2), "source": "pysis/downhill_irc_data.h5"}
    _write(wd, _syms(sym2), list(oC[:-1]) + list(C2), list(oj["energies"][:-1]) + list(E2),
           [np.asarray(g) for g in oj["gradients"][:-1]] + list(G2), rec)
    out["continuation"] = {"points": len(oC) - 1 + len(E2), "new_points": len(E2) - 1,
                           "E_last_above_min_kJ": (E2[-1] - tso["E_eh"]) * 2625.4996394799, "not_kept": trj_n - len(E2)}
    (RX / "02_irc/attempts/continuation_test/stop.json").write_text(
        json.dumps({"time": time.strftime("%Y-%m-%d %H:%M:%S"), **out}, indent=1), encoding="utf-8")
    return out




def install(agree: bool) -> dict:
    """J. Picione, 2026-11-16 afternoon, item 4: if the test agrees, the continuation's frames (up to the stop) are the
    tertiary open side, ended there; otherwise the rerun's frames are. The side not used is kept as a record in
    02_irc/attempts/. Either way the side is accepted as it ended (accepted.json): the 0.5 kJ/mol target is dropped
    and the end is confirmed by optimization from its last frame at the join (end check)."""
    import shutil
    D = RX / "02_irc"
    stamp = time.strftime("%H%M")
    cont, rer = D / "attempts/continuation_test/backward", D / "true_irc/backward"
    tso = json.loads((ROOT / "docs/stationary_m062x.json").read_text(encoding="utf-8"))["hshift_tertiary/open_cation"]
    if agree:
        keep = D / f"attempts/rerun_from_ts_test_{stamp}"
        keep.mkdir(parents=True)
        shutil.move(str(rer), str(keep / "backward"))
        for c in ("console_true_irc_backward.txt",):
            if (D / c).exists():
                shutil.move(str(D / c), str(keep / c))
        (keep / "WHY.txt").write_text(
            "NOT the production open side. Rerun from the TS for the continuation test (J. Picione 2026-11-16), stopped "
            "at the early comparison once it had 50 frames past the old side's last frame. The test agreed, so the "
            "continuation's frames are the open side (02_irc/true_irc/backward); this rerun is kept as the test record "
            "(02_irc/attempts/continuation_test/compare_backward.md).\n", encoding="utf-8")
        shutil.copytree(cont, rer)
        used, other = "continuation from the last frame of the earlier side", keep.relative_to(RX).as_posix()
    else:
        used, other = "rerun from the TS", "02_irc/attempts/continuation_test"
    js = json.loads((rer / "path.json").read_text(encoding="utf-8"))
    above = (js["energies"][-1] - tso["E_eh"]) * 2625.4996394799
    acc = {"decision": "J. Picione, 2026-11-16 (afternoon): the 0.5 kJ/mol target is dropped; the tertiary open side "
                       f"ends where the test runs were stopped ({used}); the end is confirmed by optimization from its "
                       "last frame and the paper lists how far above its minimum it ended",
           "ended": js["ended"], "used": used, "other_kept_in": other, "agree": agree,
           "end_above_min_kJ": above, "points": js["n_points"],
           "card_note": (f"The open-cation side of the path ends {above:.2f} kJ/mol ({above / 4.184:.2f} kcal/mol) above the "
                         "optimized open cation, in the slow final approach along soft modes; "
                         + ("continued from the earlier side's last frame (checked against a rerun from the TS)" if agree
                            else "rerun from the TS") + " and stopped there by decision; the end is confirmed by "
                         "optimization from the last frame.")}
    (rer / "accepted.json").write_text(json.dumps(acc, indent=1), encoding="utf-8")
    return acc


if __name__ == "__main__":
    if sys.argv[1:] == ["materialize"]:
        print(json.dumps(materialize(), indent=1))
    elif sys.argv[1:2] == ["install"]:
        print(json.dumps(install(sys.argv[2] == "agree"), indent=1))
