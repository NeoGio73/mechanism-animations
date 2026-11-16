"""Production mass-weighted IRC at the working level (J. Picione's rule, 2026-11-13).

The method validated on SN2 (pipeline/true_irc2.py, run final_gs_s005): pysisyphus 1.0.0
Gonzalez–Schlegel with the 2026-11-13 fixes (Bofill update, translations and rotations projected out
in the current frame, bracketed λ search, a failed search is a failure), step 0.05 bohr·amu½ measured
natively in mass-weighted coordinates, exact TS Hessian at the start (01_stationary/ts/ts.hess), then
updates only; each side ended by the project's stopping rule (|ΔE| < 0.05 kJ/mol per step for 20
consecutive steps after the last bond event); pysisyphus's own gradient and energy stops switched off.
Energies and gradients: ORCA at the reaction's working level (reaction.toml "method", SMD solvent,
its nprocs and maxcore), one ORCA job per folder (02_irc/true_irc/<direction>/g####).

Result in 02_irc/, the same layout as pipeline.steps.irc: irc_IRC_Full_trj.xyz (backward reversed +
TS + forward), irc_IRC_F.xyz, irc_IRC_B.xyz, irc_path.json; plus true_irc/<direction>/path.json with
every call. A side that does not end by the stopping rule makes the path incomplete; nothing is joined.
Usage: .venv-pysis/Scripts/python.exe -X utf8 -m pipeline.true_irc_prod <reaction folder> <forward|backward|join>
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

from . import orca
from .geometry import Geometry as PGeom
from .geometry import read_xyz, read_xyz_frames
from .mass_weight_evidence import read_hess_matrix
from .steps import NEAR_MIN_KJ, Reaction, _write_json, end_min

ROOT = Path(__file__).resolve().parents[1]
BOHR = 0.529177210903
STEP = 0.05            # bohr·amu½ (validated: SN2, B3LYP-D3(BJ)/ma-def2-SVP, tail test and full path)
MAX_CYCLES = 3000
METHOD = ("pysisyphus 1.0.0 Gonzalez–Schlegel, mass-weighted, step 0.05 bohr·amu½, exact TS Hessian then "
          "Bofill updates, translations/rotations projected, project stopping rule (2026-11-13 fixes; amended by J. Picione "
          "2026-11-15: flat steps count only >= 1 kJ/mol below the TS, end-minimum test on 40 consecutive frames)")


def ts_paths(rx):
    st = rx.root / "01_stationary/ts/ts"
    return Path(f"{st}.xyz"), Path(f"{st}.hess"), Path(f"{st}.out")


def e_final(out: Path) -> float:
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out.read_text(errors="replace"))[-1])


class JobFailed(Exception):
    """An ORCA job of the path did not end normally (crash or kill). The side is incomplete, never a result."""


def run_side(rx_root: Path, direction: str, test: bool = False, step: float | None = None,
             check_dir: str | None = None, max_cycles: int | None = None, cont_from: str | None = None,
             out: str | None = None) -> dict:
    """test: 3 steps into 02_irc/smoketest/<direction> (checks the driver; never a result).
    step / check_dir / max_cycles: a check run (J. Picione, 2026-11-15: tertiary open side at half the step),
    in 02_irc/<check_dir>/<direction>, all other settings as in production; never a production path.
    cont_from (J. Picione, 2026-11-16): continue the side in 02_irc/<cont_from> downhill from its last frame
    (pysisyphus downhill run: no displacement, exact Hessian at that frame by an ORCA Freq job, then Bofill updates;
    step, patches and stopping rule as in production, the rule seeing the whole side). path.json/path.xyz in the new
    folder hold the whole side: the earlier frames, then the new ones (the start frame recomputed, its energy checked).
    out: the folder (relative to 02_irc) instead of the default."""
    STEP_ = step or STEP
    from pysisyphus.calculators.Calculator import Calculator
    from pysisyphus.Geometry import Geometry
    from .true_irc2 import GSmw, IntegratorFailure

    rx = Reaction.load(rx_root)
    tsx, tsh, tso = ts_paths(rx)
    ts = read_xyz(tsx)
    sym = list(ts.symbols)
    E_TS = e_final(tso)
    x_ts = np.asarray(ts.coords).ravel() / BOHR
    wd = rx.root / ("02_irc/smoketest" if test else f"02_irc/{check_dir}" if check_dir else "02_irc/true_irc") / direction
    if out:
        wd = rx.root / "02_irc" / out
    old = None
    if cont_from:
        od = rx.root / "02_irc" / cont_from
        oj = json.loads((od / "path.json").read_text(encoding="utf-8"))
        ofr = read_xyz_frames(od / "path.xyz")
        old = {"dir": cont_from, "json": oj, "C": [np.asarray(f.coords).ravel() / BOHR for f in ofr],
               "E": list(oj["energies"]), "G": [np.asarray(g) for g in oj["gradients"]]}
        assert len(old["C"]) == len(old["E"]) == len(old["G"]), "earlier side: frames, energies and gradients differ in number"
    wd.mkdir(parents=True, exist_ok=False)          # a clean folder; a rerun goes to a new one by hand
    (wd / "pysis").mkdir()
    os.chdir(wd / "pysis")

    class RxCalc(Calculator):
        def __init__(self, **kw):
            super().__init__(charge=rx.cfg["charge"], mult=rx.cfg["mult"], **kw)
            self.n, self.calls = 0, []

        def _job(self, atoms, coords, kind):
            jd = wd / f"{'g' if kind == 'EnGrad' else 'h'}{self.n:04d}"
            self.n += 1
            c = np.asarray(coords).reshape(-1, 3)
            if not np.all(np.isfinite(c)):
                raise IntegratorFailure("non-finite (NaN) coordinates from the integrator; no ORCA job sent")
            jd.mkdir(parents=True, exist_ok=False)
            g = PGeom(sym, c * BOHR, f"{kind} {self.n - 1}")
            r = orca.run(jd, "j", rx.inp(f"{kind} NoPop", g), rx.log,
                         note=("smoke test (not a result) " if test else "") + (f"check run {check_dir} " if check_dir else "") + f"true IRC production {direction} {kind} {self.n - 1} (step {STEP_})")
            if not r.ok:
                # not a RuntimeError: pysisyphus's GonzalezSchlegel.step catches RuntimeError and reports "constrained
                # search did not converge" (2026-11-15: a job killed by the app update was mislabelled that way)
                raise JobFailed(f"ORCA job {jd.name} did not end normally")
            self.calls.append({"folder": jd.name, "kind": kind, "wall_s": r.wall_s})
            for f in jd.glob("*.tmp*"):
                f.unlink()
            return jd, r

        def _engrad(self, atoms, coords):
            jd, r = self._job(atoms, coords, "EnGrad")
            v = [ln.strip() for ln in (jd / "j.engrad").read_text().splitlines() if ln.strip() and not ln.startswith("#")]
            nat = int(v[0])
            return float(v[1]), np.array([float(x) for x in v[2:2 + 3 * nat]])

        def get_energy(self, atoms, coords, **kw):
            return {"energy": self._engrad(atoms, coords)[0]}

        def get_forces(self, atoms, coords, **kw):
            E, g = self._engrad(atoms, coords)
            return {"energy": E, "forces": -g}

        def get_hessian(self, atoms, coords, **kw):
            if np.max(np.abs(np.asarray(coords) - x_ts)) < 1e-5:
                return {"energy": E_TS, "hessian": read_hess_matrix(tsh)}
            jd, r = self._job(atoms, coords, "Freq")
            return {"energy": r.energy_eh, "hessian": read_hess_matrix(jd / "j.hess")}

    geom = Geometry(sym, old["C"][-1] if old else x_ts)
    calc = RxCalc(out_dir=str(wd / "pysis"))
    geom.set_calculator(calc)
    class GSprod(GSmw):
        """pysisyphus's IRC loop always stops a side when one step raises the energy ("Energy increased!"); no
        setting turns it off. On a flat surface the per-step change is below the SCF scatter (0.034 kJ/mol), so a
        noise-level rise ended sides early (2026-11-14: tertiary hydride shift forward 0.015 kJ/mol below the TS;
        Diels–Alder forward at its minimum before the 20 flat steps). The protocol has pysisyphus's own stops
        off, so this one is off too: the loop reads energy_increased as False; every rise is recorded in
        energy_rises and reported (patch 0009)."""
        @property
        def energy_increased(self):
            return False

        @energy_increased.setter
        def energy_increased(self, value):
            if value and len(getattr(self, "irc_energies", [])) >= 2:
                self.energy_rises.append({"point": len(self.irc_energies) - 1,
                                          "rise_kJ": (self.irc_energies[-1] - self.irc_energies[-2]) * 2625.499639})

    irc = GSprod(geom, hessian_recalc=None, step_length=STEP_, max_cycles=3 if test else (max_cycles or MAX_CYCLES), rms_grad_thresh=1e-12,
               energy_thresh=-1.0, hessian_init="calc", out_dir=str(wd / "pysis"),
               forward=direction == "forward", backward=direction == "backward", displ="length", displ_length=STEP_,
               downhill=bool(old))
    irc.rx = rx
    em = None if (test or (check_dir and not old)) else end_min(rx, direction)      # J. Picione, 2026-11-16 (reaction.toml [irc_end_min])
    if old:
        irc.prefix_frames = [PGeom(sym, c.reshape(-1, 3) * BOHR, f"x E {e:.12f}") for c, e in zip(old["C"][:-1], old["E"][:-1])]
        irc.prefix_gradients = old["G"][:-1]
    irc.e_min_eh = em[1] if em else None
    irc.ts_frame = PGeom(sym, np.asarray(ts.coords), f"TS E {E_TS:.12f}")
    irc.own_stop = None
    irc.pause_file, irc.paused = wd / "PAUSE", None
    irc.fallbacks, irc.micro_log, irc.micro_max_hit, irc.energy_rises = [], [], [], []
    t0, failure = time.time(), None
    try:
        irc.run()
    except (IntegratorFailure, JobFailed) as e:
        failure = str(e)
    if failure is None:
        key = "downhill" if old else direction
        C = [np.asarray(c) for c in getattr(irc, f"{key}_coords")]
        E = list(getattr(irc, f"{key}_energies"))
        G = [np.asarray(g) for g in getattr(irc, f"{key}_gradients")]
    else:
        n = min(len(irc.irc_coords), len(irc.irc_energies), len(irc.irc_gradients))
        C, E, G = ([np.asarray(c) for c in irc.irc_coords[:n]], list(irc.irc_energies[:n]),
                   [np.asarray(g) for g in irc.irc_gradients[:n]])
    if C and np.linalg.norm(C[-1] - x_ts) < np.linalg.norm(C[0] - x_ts):      # TS-outward order
        C, E, G = C[::-1], E[::-1], G[::-1]
    ended = (f"FAILED: {failure}" if failure else f"paused ({irc.paused})" if irc.paused else
             "own stopping rule" if irc.own_stop else "max cycles")
    cont = None
    if old and C:
        cont = {"from": old["dir"], "earlier_points": len(old["C"]), "new_points": len(C) - 1,
                "start_frame_energy_recomputed_minus_earlier_kJ": (E[0] - old["E"][-1]) * 2625.499639,
                "start_frame_max_coordinate_difference_bohr": float(np.max(np.abs(C[0] - old["C"][-1]))),
                "hessian_at_start": "exact (ORCA Freq at the start frame, h0000), then Bofill updates",
                "earlier_n_engrad": old["json"].get("n_engrad"), "earlier_wall_s": old["json"].get("wall_total_s")}
        C, E, G = old["C"][:-1] + C, old["E"][:-1] + E, old["G"][:-1] + G
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C, E):
            fh.write(PGeom(sym, c.reshape(-1, 3) * BOHR, "").to_xyz(f"true IRC {direction} E {e:.12f}"))
    res = {"reaction": rx.root.name, "direction": direction, "method": METHOD.replace("step 0.05", f"step {STEP_}"), "level": rx.cfg["method"],
           "solvent": rx.cfg["solvent"], "nprocs": rx.cfg["nprocs"], "step_bohr_amu05": STEP_, "ended": ended,
           "check_run": check_dir, "continuation": cont,
           "own_stop": irc.own_stop, "end_min": ({"minimum": em[0], "E_eh": em[1], "near_min_kJ": NEAR_MIN_KJ} if em else None),
           "n_points": len(C), "E_TS_eh": E_TS,
           "n_engrad": sum(c["kind"] == "EnGrad" for c in calc.calls), "n_hessian": sum(c["kind"] == "Freq" for c in calc.calls),
           "calls": calc.calls, "wall_total_s": time.time() - t0, "energies": E, "gradients": [g.tolist() for g in G],
           "gs_micro_steps": len(irc.micro_log), "gs_micro_cycle_limit_hit_at": irc.micro_max_hit,
           "energy_rises": irc.energy_rises, "energy_increase_stop": "off (patch 0009, 2026-11-14)"}
    (wd / "path.json").write_text(json.dumps(res), encoding="utf-8")
    rx.log.note(("SMOKE TEST, not a result: " if test else f"CHECK RUN (step {STEP_}), not a production path: " if check_dir else "") + f"`{wd.relative_to(rx.root).as_posix()}`: production true IRC ({METHOD}; ORCA {rx.cfg['method']}, SMD "
                f"{rx.cfg['solvent']}, {rx.cfg['nprocs']} cores): "
                + (f"continued downhill from the last frame of `02_irc/{cont['from']}` ({cont['earlier_points']} earlier "
                   f"points, {cont['new_points']} new; start frame recomputed {cont['start_frame_energy_recomputed_minus_earlier_kJ']:+.4f} kJ/mol); "
                   if cont else "")
                + (f"ended by the stopping rule{f' (flat steps counted within {NEAR_MIN_KJ:g} kJ/mol of {em[0]}, J. Picione 2026-11-16)' if em else ''}, {len(C)} points" if ended == "own stopping rule"
                   else f"INCOMPLETE ({ended}), not a result")
                + f"; {res['n_engrad']} gradients, {res['n_hessian']} Hessians, {res['wall_total_s'] / 3600:.1f} h")
    print(json.dumps({k: v for k, v in res.items() if k not in ("gradients", "energies", "calls")}))
    return res


def _rerun_marked_sides(rx: Reaction) -> None:
    """J. Picione, 2026-11-16: a side marked for a redo (02_irc/true_irc/rerun_<side>.json, written by Claude for an
    approved redo) is redone here, at the start of the join, so the chain and the restart guard see one lane.
    The marker's "method":
      "pending"  - the method is not settled yet (the continuation test); the join waits, polling every minute;
      "rerun"    - the side's folder and console go to 02_irc/attempts/<archive_as>/ (not results, WHY.txt) and the
                   side runs from the TS in a clean folder with the settings in reaction.toml;
      "continue" - the side's folder becomes true_irc/<side>_part1 (kept: its frames are the first part of the side)
                   and the side is continued downhill from its last frame into a clean true_irc/<side> (run_side
                   cont_from), whose path.json holds the whole side.
    Either way the earlier join's outputs are copied to the attempts folder as a record and its end-check folders
    moved there (the new end check runs in clean folders); the marker moves there too.
    If this process is killed, the guard archives the unfinished side and the chain reruns it from the TS."""
    import shutil
    import subprocess
    wd = rx.root / "02_irc"
    for d in ("backward", "forward"):
        mk = wd / "true_irc" / f"rerun_{d}.json"
        if not mk.exists():
            continue
        m = json.loads(mk.read_text(encoding="utf-8"))
        if m.get("method", "rerun") == "pending":
            with open(ROOT / "docs/production_true_irc.log", "a", encoding="utf-8") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} production true IRC {rx.root.name}: join waits for the "
                         f"{d} side's method (continuation test, J. Picione 2026-11-16)\n")
            rx.log.note(f"join waits: the {d} side's redo method (continue from its last frame or rerun from the TS) is "
                        "settled by the continuation test on the tertiary open side (J. Picione 2026-11-16)")
            while m.get("method", "rerun") == "pending":
                time.sleep(60)
                m = json.loads(mk.read_text(encoding="utf-8"))
        method = m.get("method", "rerun")
        A = wd / "attempts" / m["archive_as"]
        A.mkdir(parents=True, exist_ok=False)
        for f in ("irc_IRC_Full_trj.xyz", "irc_IRC_B.xyz", "irc_IRC_F.xyz", "irc_path.json", "mirror_path.json",
                  "endpoints.json", "true_irc/angles.json", "true_irc/angles.md", "true_irc/angles_ts_region.md"):
            if (wd / f).exists():
                (A / "earlier_join").mkdir(exist_ok=True)
                shutil.copy2(wd / f, A / "earlier_join" / Path(f).name)
        moves = sorted((wd / "true_irc").glob("end_check*"))
        if method == "rerun":
            moves += [wd / "true_irc" / d, wd / f"console_true_irc_{d}.txt"]
        else:
            shutil.move(str(wd / "true_irc" / d), str(wd / "true_irc" / f"{d}_part1"))
            if (wd / f"console_true_irc_{d}.txt").exists():
                shutil.move(str(wd / f"console_true_irc_{d}.txt"), str(wd / f"console_true_irc_{d}_part1.txt"))
        for f in moves:
            if f.exists():
                shutil.move(str(f), str(A / f.name))
        (A / "WHY.txt").write_text(m["why"] + "\n", encoding="utf-8")
        shutil.move(str(mk), str(A / mk.name))
        t = time.strftime("%Y-%m-%d %H:%M:%S")
        what = ("rerun from the TS" if method == "rerun" else
                f"continued downhill from its last frame (first part kept as true_irc/{d}_part1)")
        with open(ROOT / "docs/production_true_irc.log", "a", encoding="utf-8") as fh:
            fh.write(f"{t} production true IRC {rx.root.name}: {d} side {what} at the start of the join "
                     f"({m['decision']}); records in 02_irc/attempts/{m['archive_as']}\n")
        rx.log.note(f"{d} side {what}: {m['decision']} (records in `02_irc/attempts/{m['archive_as']}`)")
        (wd / "logs" / d).mkdir(parents=True, exist_ok=True)
        con = wd / f"console_true_irc_{d}.txt"
        extra = [] if method == "rerun" else [f"continue_from=true_irc/{d}_part1", f"out=true_irc/{d}"]
        with open(con, "w", encoding="utf-8") as fh:
            rc = subprocess.run([sys.executable, "-X", "utf8", "-W", "ignore", "-u", "-m", "pipeline.true_irc_prod",
                                 str(rx.root), d] + extra, cwd=wd / "logs" / d, stdout=fh, stderr=subprocess.STDOUT,
                                env={**os.environ, "PYTHONPATH": str(ROOT)}).returncode
        with open(con, "a", encoding="utf-8") as fh:
            fh.write(f"end {time.strftime('%H:%M')} exit {rc}\n")
        if rc != 0:
            raise RuntimeError(f"{d} side ({method}) ended with exit {rc}: incomplete, nothing joined")


def join(rx_root: Path) -> dict:
    """Both sides ended by the stopping rule: write 02_irc in the layout of pipeline.steps.irc."""
    rx = Reaction.load(rx_root)
    _rerun_marked_sides(rx)
    tsx, _, tso = ts_paths(rx)
    ts = read_xyz(tsx)
    wd = rx.root / "02_irc"
    sides = {}
    for d, tag in (("forward", "F"), ("backward", "B")):
        js = json.loads((wd / "true_irc" / d / "path.json").read_text(encoding="utf-8"))
        acc = wd / "true_irc" / d / "accepted.json"       # a side accepted by J. Picione's decision, with its note
        js["accepted"] = json.loads(acc.read_text(encoding="utf-8")) if acc.exists() else None
        if js["ended"] != "own stopping rule" and not js["accepted"]:
            raise RuntimeError(f"{d} side ended '{js['ended']}': incomplete, nothing joined")
        fr = read_xyz_frames(wd / "true_irc" / d / "path.xyz")
        x0 = np.asarray(ts.coords)
        keep = [k for k, f in enumerate(fr) if np.max(np.abs(np.asarray(f.coords) - x0)) > 1e-6]
        sides[tag] = ([fr[k] for k in keep], [js["energies"][k] for k in keep], js)
    fb, eb, jb = sides["B"]
    ff, ef, jf = sides["F"]
    # TS energy from the path's own gradient job at the TS (g0000 of the forward side), so every frame's
    # energy comes from the same kind of job (the OptTS output differs by ~1e-5 Eh, SCF scatter)
    e_ts = e_final(wd / "true_irc/forward/g0000/j.out")
    frames, E = fb[::-1] + [ts] + ff, eb[::-1] + [e_ts] + ef
    (wd / "irc_IRC_Full_trj.xyz").write_text(
        "".join(g.to_xyz(f"IRC frame {k} E {e:.12f}") for k, (g, e) in enumerate(zip(frames, E))))
    fb[-1].write(wd / "irc_IRC_B.xyz", f"IRC B end E {eb[-1]:.12f}")
    ff[-1].write(wd / "irc_IRC_F.xyz", f"IRC F end E {ef[-1]:.12f}")
    stops = {t: ({"by": "rule", "last_step": s[2]["own_stop"]["stop_frame"],
                  "last_event_step": s[2]["own_stop"]["last_event_frame"]} if s[2]["own_stop"] else
                 {"by": f"accepted: {s[2]['accepted']['decision']}", "last_step": len(s[0]), "ended": s[2]["ended"],
                  "note": s[2]["accepted"]["card_note"]}) for t, s in sides.items()}
    _write_json(wd / "irc_path.json", {
        "energies_eh": E, "ts_index": len(fb), "n_backward": len(fb), "n_forward": len(ff), "stops": stops,
        "method": METHOD, "ts_source": "01_stationary/ts/ts",
        "note": "frames = backward (reversed) + TS + forward; energies from the true IRC (ORCA EnGrad at the working level)"})
    rx.log.note(f"`02_irc/irc_IRC_Full_trj.xyz`: true IRC joined, {len(frames)} frames (backward {len(fb)} + TS + "
                f"forward {len(ff)}), " + "; ".join(f"{t}: {v['by']}" for t, v in stops.items()))
    if rx.cfg.get("mirror_path"):
        # degenerate rearrangement (J. Picione, 2026-11-14): computed half to the bridged structure + its mirror image
        from .mirror_path import build
        m = build(rx, fb, eb, ff, ef, ts, e_ts)
        frames, E, info = m["frames"], m["energies"], m["info"]
        (wd / "irc_IRC_Full_trj.xyz").write_text(
            "".join(g.to_xyz(f"frame {k} E {e:.12f} {s}") for k, (g, e, s) in enumerate(zip(frames, E, m["frame_source"]))))
        frames[0].write(wd / "irc_IRC_B.xyz", f"path start (computed) E {E[0]:.12f}")
        frames[-1].write(wd / "irc_IRC_F.xyz", f"path end (mirror image) E {E[-1]:.12f}")
        _write_json(wd / "irc_path.json", {
            "energies_eh": E, "ts_index": info["ts_index"], "stops": stops, "method": METHOD,
            "ts_source": "01_stationary/ts/ts", "frame_source": m["frame_source"], "mirror_path": info,
            "note": ("computed IRC from the open cation through the TS to the bridged structure, then the mirror image "
                     "(C2⇄C3 swap) of those frames back out to the other open cation; mirror frames are not IRC points")})
        _write_json(wd / "mirror_path.json", info)
        rx.log.note(f"`02_irc/mirror_path.json`: path = computed half ({info['n_computed']} frames, TS to the bridged structure "
                    f"at frame {info['bridged_index']}) + mirror image ({info['n_mirror']} frames, C2⇄C3 swap). Bridged side: "
                    f"end {info['end_rmsd_to_reference_A']:.3f} Å from {info['reference']}, "
                    f"{'DRIFTED past it' if info['drifted_past_reference'] else 'no drift past it'} "
                    f"({info['frames_beyond_end_not_used']} frames beyond the end not used"
                    + ("; after the end the side slid back toward the open cation it started from"
                       if info.get("returned_toward_start_after_end") else "") + "); energy span TS to end "
                    f"{info['energy_span_TS_to_end_kJ']:.4f} kJ/mol against the SCF noise floor {info['noise_floor_kJ']} kJ/mol.")
    # angle check on every production path (J. Picione, 2026-11-13)
    from .path_angles import check
    a = check(rx.root)
    zs = "; ".join(f"{d}: {v['angles_near']['median_deg']:.2f}° / {v['angles_near']['p90_deg']:.2f}° within |s| ≤ 1, "
                   + (f"{v['angles_beyond']['median_deg']:.2f}° / {v['angles_beyond']['p90_deg']:.2f}° beyond"
                      if v["angles_beyond"] else "none beyond")
                   + (f", end nearest {v['end']['nearest']['minimum']} (RMSD {v['end']['nearest']['rmsd_A']:.3f} Å, "
                      f"{v['end']['nearest']['end_above_min_kJ']:+.3f} kJ/mol)" if v["end"] else "")
                   for d, v in a["sides"].items())
    rx.log.note(f"`02_irc/true_irc/angles.md`: angle check {'passes' if a['passes'] else 'FAILS (flagged)'} "
                f"(median / 90th percentile; rule median < 5°, 90 % < 15°): {zs}")
    # end check (J. Picione, 2026-11-14): a side ending > 0.5 kJ/mol above its minimum is optimized from its
    # last frame; accepted with a note only if that reaches the same minimum, otherwise a decision
    from .end_check import check as end_check
    ec = end_check(rx.root)
    _after_join(rx)
    return {"frames": len(frames), "ts_index": len(fb), "angles_pass": a["passes"], "end_check": ec["verdict"]}


def _after_join(rx: Reaction) -> None:
    """J. Picione, 2026-11-16 ("rerun it after the secondary side, and redo both joins"): 02_irc/after_join.json names
    joined reactions to withdraw once this join has finished. Their "join withdrawn" lines are written before the
    chain writes this reaction's "join exit" line, so the restart guard then continues the chain with them (it
    restarts the reactions of its list in a fixed order, so they cannot be withdrawn while this one still runs)."""
    mk = rx.root / "02_irc" / "after_join.json"
    if not mk.exists():
        return
    m = json.loads(mk.read_text(encoding="utf-8"))
    t = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(ROOT / "docs/production_true_irc.log", "a", encoding="utf-8") as fh:
        for w in m["withdraw"]:
            fh.write(f"{t} production true IRC {w['reaction']}: join withdrawn ({w['why']})\n")
    rx.log.note("after this join: " + "; ".join(f"join of {w['reaction']} withdrawn ({w['why']})" for w in m["withdraw"]))
    mk.rename(mk.with_name("after_join.done.json"))


if __name__ == "__main__":
    root = Path(sys.argv[1])
    root = root if root.is_absolute() else (ROOT / "reactions" / root)
    if sys.argv[2] == "join":
        print(join(root))
    else:
        kw = dict(a.split("=", 1) for a in sys.argv[3:] if "=" in a)      # step=0.025 check_dir=... max_cycles=80
        run_side(root, sys.argv[2], test="test" in sys.argv[3:], step=float(kw["step"]) if "step" in kw else None,
                 check_dir=kw.get("check_dir"), max_cycles=int(kw["max_cycles"]) if "max_cycles" in kw else None,
                 cont_from=kw.get("continue_from"), out=kw.get("out"))
