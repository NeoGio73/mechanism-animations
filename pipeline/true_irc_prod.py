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
from .steps import Reaction, _write_json

ROOT = Path(__file__).resolve().parents[1]
BOHR = 0.529177210903
STEP = 0.05            # bohr·amu½ (validated: SN2, B3LYP-D3(BJ)/ma-def2-SVP, tail test and full path)
MAX_CYCLES = 3000
METHOD = ("pysisyphus 1.0.0 Gonzalez–Schlegel, mass-weighted, step 0.05 bohr·amu½, exact TS Hessian then "
          "Bofill updates, translations/rotations projected, project stopping rule (2026-11-13 fixes)")


def ts_paths(rx):
    st = rx.root / "01_stationary/ts/ts"
    return Path(f"{st}.xyz"), Path(f"{st}.hess"), Path(f"{st}.out")


def e_final(out: Path) -> float:
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out.read_text(errors="replace"))[-1])


def run_side(rx_root: Path, direction: str, test: bool = False) -> dict:
    """test: 3 steps into 02_irc/smoketest/<direction> (checks the driver; never a result)."""
    from pysisyphus.calculators.Calculator import Calculator
    from pysisyphus.Geometry import Geometry
    from .true_irc2 import GSmw, IntegratorFailure

    rx = Reaction.load(rx_root)
    tsx, tsh, tso = ts_paths(rx)
    ts = read_xyz(tsx)
    sym = list(ts.symbols)
    E_TS = e_final(tso)
    x_ts = np.asarray(ts.coords).ravel() / BOHR
    wd = rx.root / ("02_irc/smoketest" if test else "02_irc/true_irc") / direction
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
                         note=("smoke test (not a result) " if test else "") + f"true IRC production {direction} {kind} {self.n - 1} (step {STEP})")
            if not r.ok:
                raise orca.OrcaError(f"{jd} did not end normally")
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

    geom = Geometry(sym, x_ts)
    calc = RxCalc(out_dir=str(wd / "pysis"))
    geom.set_calculator(calc)
    irc = GSmw(geom, hessian_recalc=None, step_length=STEP, max_cycles=3 if test else MAX_CYCLES, rms_grad_thresh=1e-12,
               energy_thresh=-1.0, hessian_init="calc", out_dir=str(wd / "pysis"),
               forward=direction == "forward", backward=direction == "backward", displ="length", displ_length=STEP)
    irc.rx = rx
    irc.ts_frame = PGeom(sym, np.asarray(ts.coords), f"TS E {E_TS:.12f}")
    irc.own_stop = None
    irc.fallbacks, irc.micro_log, irc.micro_max_hit = [], [], []
    t0, failure = time.time(), None
    try:
        irc.run()
    except IntegratorFailure as e:
        failure = str(e)
    if failure is None:
        C = [np.asarray(c) for c in getattr(irc, f"{direction}_coords")]
        E = list(getattr(irc, f"{direction}_energies"))
        G = [np.asarray(g) for g in getattr(irc, f"{direction}_gradients")]
    else:
        n = min(len(irc.irc_coords), len(irc.irc_energies), len(irc.irc_gradients))
        C, E, G = ([np.asarray(c) for c in irc.irc_coords[:n]], list(irc.irc_energies[:n]),
                   [np.asarray(g) for g in irc.irc_gradients[:n]])
    if C and np.linalg.norm(C[-1] - x_ts) < np.linalg.norm(C[0] - x_ts):      # TS-outward order
        C, E, G = C[::-1], E[::-1], G[::-1]
    ended = (f"FAILED: {failure}" if failure else "own stopping rule" if irc.own_stop else
             "energy increased" if getattr(irc, f"{direction}_energy_increased", False) else "max cycles")
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C, E):
            fh.write(PGeom(sym, c.reshape(-1, 3) * BOHR, "").to_xyz(f"true IRC {direction} E {e:.12f}"))
    res = {"reaction": rx.root.name, "direction": direction, "method": METHOD, "level": rx.cfg["method"],
           "solvent": rx.cfg["solvent"], "nprocs": rx.cfg["nprocs"], "step_bohr_amu05": STEP, "ended": ended,
           "own_stop": irc.own_stop, "n_points": len(C), "E_TS_eh": E_TS,
           "n_engrad": sum(c["kind"] == "EnGrad" for c in calc.calls), "n_hessian": sum(c["kind"] == "Freq" for c in calc.calls),
           "calls": calc.calls, "wall_total_s": time.time() - t0, "energies": E, "gradients": [g.tolist() for g in G],
           "gs_micro_steps": len(irc.micro_log), "gs_micro_cycle_limit_hit_at": irc.micro_max_hit}
    (wd / "path.json").write_text(json.dumps(res), encoding="utf-8")
    rx.log.note(("SMOKE TEST, not a result: " if test else "") + f"`{wd.relative_to(rx.root).as_posix()}`: production true IRC ({METHOD}; ORCA {rx.cfg['method']}, SMD "
                f"{rx.cfg['solvent']}, {rx.cfg['nprocs']} cores): "
                + (f"ended by the stopping rule, {len(C)} points" if ended == "own stopping rule"
                   else f"INCOMPLETE ({ended}), not a result")
                + f"; {res['n_engrad']} gradients, {res['n_hessian']} Hessians, {res['wall_total_s'] / 3600:.1f} h")
    print(json.dumps({k: v for k, v in res.items() if k not in ("gradients", "energies", "calls")}))
    return res


def join(rx_root: Path) -> dict:
    """Both sides ended by the stopping rule: write 02_irc in the layout of pipeline.steps.irc."""
    rx = Reaction.load(rx_root)
    tsx, _, tso = ts_paths(rx)
    ts = read_xyz(tsx)
    wd = rx.root / "02_irc"
    sides = {}
    for d, tag in (("forward", "F"), ("backward", "B")):
        js = json.loads((wd / "true_irc" / d / "path.json").read_text(encoding="utf-8"))
        if js["ended"] != "own stopping rule":
            raise RuntimeError(f"{d} side ended '{js['ended']}': incomplete, nothing joined")
        fr = read_xyz_frames(wd / "true_irc" / d / "path.xyz")
        x0 = np.asarray(ts.coords)
        keep = [k for k, f in enumerate(fr) if np.max(np.abs(np.asarray(f.coords) - x0)) > 1e-6]
        sides[tag] = ([fr[k] for k in keep], [js["energies"][k] for k in keep], js)
    fb, eb, jb = sides["B"]
    ff, ef, jf = sides["F"]
    e_ts = e_final(tso)
    frames, E = fb[::-1] + [ts] + ff, eb[::-1] + [e_ts] + ef
    (wd / "irc_IRC_Full_trj.xyz").write_text(
        "".join(g.to_xyz(f"IRC frame {k} E {e:.12f}") for k, (g, e) in enumerate(zip(frames, E))))
    fb[-1].write(wd / "irc_IRC_B.xyz", f"IRC B end E {eb[-1]:.12f}")
    ff[-1].write(wd / "irc_IRC_F.xyz", f"IRC F end E {ef[-1]:.12f}")
    stops = {t: {"by": "rule", "last_step": s[2]["own_stop"]["stop_frame"],
                 "last_event_step": s[2]["own_stop"]["last_event_frame"]} for t, s in sides.items()}
    _write_json(wd / "irc_path.json", {
        "energies_eh": E, "ts_index": len(fb), "n_backward": len(fb), "n_forward": len(ff), "stops": stops,
        "method": METHOD, "ts_source": "01_stationary/ts/ts",
        "note": "frames = backward (reversed) + TS + forward; energies from the true IRC (ORCA EnGrad at the working level)"})
    rx.log.note(f"`02_irc/irc_IRC_Full_trj.xyz`: true IRC joined, {len(frames)} frames (backward {len(fb)} + TS + "
                f"forward {len(ff)}), both sides ended by the stopping rule")
    return {"frames": len(frames), "ts_index": len(fb)}


if __name__ == "__main__":
    root = Path(sys.argv[1])
    root = root if root.is_absolute() else (ROOT / "reactions" / root)
    if sys.argv[2] == "join":
        print(join(root))
    else:
        run_side(root, sys.argv[2], test="test" in sys.argv[3:])
