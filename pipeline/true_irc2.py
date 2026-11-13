"""Mass-weighted IRC runs for J. Picione's five items of 2026-11-12 (evening): native mass-weighted step,
Hessian recalculation, a second integrator, isotopes, and the project's own stopping rule.
Runs in .venv-pysis (pysisyphus 1.0.0); every energy, gradient and Hessian is one ORCA 6.1.1 job,
alone in its own new folder, through the core ledger.

Integrators
  eulerpc   EulerPC (Hratchian–Frisch–Schlegel), Bofill updates, corrector "mbs", 5000 predictor
            substeps; step length measured natively in mass-weighted coordinates (bohr·amu½):
            pysisyphus measures it in plain bohr (EulerPC.get_integration_length), which is replaced
            here by the mass-weighted norm; the Euler substep is step/2500.
  gs        Gonzalez–Schlegel (J. Chem. Phys. 1989, 90, 2154): constrained optimization on a
            hypersphere of radius step/2 in mass-weighted coordinates (pysisyphus GonzalezSchlegel;
            its step is mass-weighted natively), BFGS updates, up to 20 micro-cycles, 1e-3 threshold.
  recalc N  exact Hessian every N steps (ORCA analytic Freq at the same level, 2 cores).
Stopping: the project's own rule (pipeline/steps.irc_stop_index: energy change below 0.05 kJ/mol per
step for 20 consecutive steps after the last bond event), or an energy rise (pysisyphus; always on),
or max_cycles. pysisyphus's RMS-gradient and energy-convergence stops are switched off.
Level: B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX, SMD water, VeryTightSCF, 2 cores (the step-test level).

Runs (name: settings), each direction in <step_convergence>/true_irc2/<name>/<direction>/:
  see RUNS below. "downhill" runs start at a frame of the pysisyphus 0.05 path (true_irc/s0050) and
  integrate downhill from there with an exact Hessian at the start (the tail test, item 3).
Usage: .venv-pysis/Scripts/python.exe -X utf8 -m pipeline.true_irc2 <run> <forward|backward|downhill>
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.irc import EulerPC, GonzalezSchlegel

from . import orca
from .geometry import Geometry as PGeom
from .geometry import read_xyz, read_xyz_frames
from .mass_weight_evidence import read_hess_matrix
from .steps import Reaction, irc_stop_index

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions/sn2_hydroxide_bromomethane_masvp"
SC = RX / "02_irc/attempts/step_convergence"
OUT = SC / "true_irc2"
LEVEL = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX VeryTightSCF CPCM(water)"
SMD = '%cpcm smd true SMDsolvent "WATER" end'
CHARGE, MULT, NPROCS = -1, 1, 2
BOHR = 0.529177210903
E_TS = -2689.43410109009
D4 = [(i, 2.01410178) for i in (3, 4, 5, 6)]
TAIL_START = {"reactant": 80, "product": 170}    # frames of true_irc/s0050: s 3.0 (C–O 2.74 Å) and 6.2 (C–Br 3.19 Å) bohr·amu½

RUNS = {
    # item 3: tails, downhill from a frame before the deviations
    "tail_eulerpc_reactant": dict(integ="eulerpc", recalc=None, start=("reactant", "forward"), max_cycles=120),
    "tail_eulerpc_product": dict(integ="eulerpc", recalc=None, start=("product", "backward"), max_cycles=120),
    "tail_eulerpc_r5_reactant": dict(integ="eulerpc", recalc=5, start=("reactant", "forward"), max_cycles=120),
    "tail_eulerpc_r5_product": dict(integ="eulerpc", recalc=5, start=("product", "backward"), max_cycles=120),
    "tail_gs_reactant": dict(integ="gs", recalc=None, start=("reactant", "forward"), max_cycles=120),
    "tail_gs_product": dict(integ="gs", recalc=None, start=("product", "backward"), max_cycles=120),
    # J. Picione's rule (2026-11-13): neither option passed at 0.10, so both are retested at 0.05 from the
    # same start frames (run with step 0.05; twice the cycles for the same path length)
    "tail_eulerpc_r5_s005_reactant": dict(integ="eulerpc", recalc=5, start=("reactant", "forward"), max_cycles=240),
    "tail_eulerpc_r5_s005_product": dict(integ="eulerpc", recalc=5, start=("product", "backward"), max_cycles=240),
    "tail_gs_s005_reactant": dict(integ="gs", recalc=None, start=("reactant", "forward"), max_cycles=240),
    "tail_gs_s005_product": dict(integ="gs", recalc=None, start=("product", "backward"), max_cycles=240),
    # item 2: isotope control on the true IRC (deuterium on all four H), 30 steps per side
    "d4_eulerpc": dict(integ="eulerpc", recalc=None, isotopes=D4, max_cycles=30),
    "h_eulerpc": dict(integ="eulerpc", recalc=None, max_cycles=30),
    # item 4: native 0.10 step, own stopping rule, both sides from the TS (integrator set when chosen)
    "final_eulerpc": dict(integ="eulerpc", recalc=None, max_cycles=1500),
    "final_eulerpc_r5": dict(integ="eulerpc", recalc=5, max_cycles=1500),
    "final_gs": dict(integ="gs", recalc=None, max_cycles=1500),
    # J. Picione's rule, 2026-11-13: GS passed the tail test at 0.05 (Hessian every 5 did not), so the
    # full SN2 validation path is GS at 0.05 (run with step 0.05)
    "final_gs_s005": dict(integ="gs", recalc=None, max_cycles=1500),
}


class IntegratorFailure(Exception):
    """Not a RuntimeError: pysisyphus's GonzalezSchlegel.step turns RuntimeError into 'converged'."""


class OrcaQueue(Calculator):
    """Energy + gradient (EnGrad) or Hessian (analytic Freq) from one ORCA job per call."""

    def __init__(self, wd: Path, front: bool, **kw):
        super().__init__(charge=CHARGE, mult=MULT, **kw)
        self.wd, self.front = wd, front
        self.n = 0
        self.log_ = orca.RunLog(wd / "run_log.md")
        self.calls = []
        self.ts = np.asarray(read_xyz(RX / "01_ts/ts/ts.xyz").coords).ravel() / BOHR

    def _job(self, atoms, coords, kind):
        wd = self.wd / f"{'g' if kind == 'EnGrad' else 'h'}{self.n:04d}"
        self.n += 1
        wd.mkdir(parents=True, exist_ok=False)
        c = np.asarray(coords).reshape(-1, 3)
        if not np.all(np.isfinite(c)):
            raise IntegratorFailure("the integrator produced non-finite (NaN) coordinates; no ORCA job sent")
        geo = "\n".join(f"{a:2s} {x:20.12f} {y:20.12f} {z:20.12f}" for a, (x, y, z) in zip(atoms, c))
        inp = (f"! {LEVEL} {kind} Bohrs NoPop\n%pal nprocs {NPROCS} end\n%maxcore 3000\n{SMD}\n\n"
               f"* xyz {CHARGE} {MULT}\n{geo}\n*\n")
        r = orca.run(wd, "j", inp, self.log_, note=f"true IRC ({self.wd.parent.name} {self.wd.name}) {kind} {self.n - 1}",
                     front=self.front)
        if not r.ok:
            raise orca.OrcaError(f"{wd} did not end normally")
        out = (wd / "j.out").read_text(errors="replace")
        m = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", out)
        wall = (int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]) + int(m[5]) / 1000) if m else None
        self.calls.append({"folder": wd.name, "kind": kind, "wall_s": wall})
        for f in wd.glob("*.tmp*"):
            f.unlink()
        return wd, r

    def _engrad(self, atoms, coords):
        wd, r = self._job(atoms, coords, "EnGrad")
        vals = [ln.strip() for ln in (wd / "j.engrad").read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        nat = int(vals[0])
        return float(vals[1]), np.array([float(v) for v in vals[2:2 + 3 * nat]])

    def get_energy(self, atoms, coords, **kw):
        return {"energy": self._engrad(atoms, coords)[0]}

    def get_forces(self, atoms, coords, **kw):
        E, g = self._engrad(atoms, coords)
        return {"energy": E, "forces": -g}

    def get_hessian(self, atoms, coords, **kw):
        if np.max(np.abs(np.asarray(coords) - self.ts)) < 1e-5:          # the TS: ORCA's TS Hessian
            return {"energy": E_TS, "hessian": read_hess_matrix(RX / "01_ts/ts/ts.hess")}
        wd, r = self._job(atoms, coords, "Freq")
        return {"energy": r.energy_eh, "hessian": read_hess_matrix(wd / "j.hess")}


class OwnStop:
    """The project's stopping rule, checked after every step (pysisyphus prints this line and stops
    with "Integrator indicated convergence!" when it is met)."""

    def get_additional_print(self):
        sym = [a.capitalize() for a in self.atoms]          # pysisyphus keeps symbols in lower case
        fr = [self.ts_frame] + [PGeom(sym, np.asarray(c).reshape(-1, 3) * BOHR, f"x E {e:.12f}")
                                for c, e in zip(self.irc_coords, self.irc_energies)]
        hit = irc_stop_index(self.rx, fr)
        if hit is not None:
            self.converged = True
            self.own_stop = {"direction": self.cur_direction, "last_event_frame": int(hit[0]), "stop_frame": int(hit[1])}
            return f"own stopping rule met (last bond event at frame {hit[0]}, stop at frame {hit[1]})"
        return ""


class EulerPCmw(OwnStop, EulerPC):
    """Fix of 2026-11-13 (both EulerPC tail runs ended in NaN): when pysisyphus's mBS corrector meets
    an oscillation on its interpolated surface it returns a point short of the step (0.027, then
    exactly 0 bohr·amu½ at C–O ≈ 3.1 Å); the next Hessian update divides by that zero step and fills
    the Hessian with NaN. Now a corrected step shorter than half the step length is replaced by the
    predictor point (full step length, Euler integration on the local quadratic model); every such step
    is counted (self.fallbacks) and reported. A zero step never reaches the Hessian update."""
    fallbacks: list

    def get_integration_length_func(self, init_mw_coords):
        def get_integration_length(cur_mw_coords):
            return np.linalg.norm(cur_mw_coords - init_mw_coords)      # bohr·amu½
        return get_integration_length

    def get_conv_fact(self, mw_grad, min_fact=2.0):
        return 2.0                                                      # Euler substep = step / 2500

    def corrector_step(self, init_mw_coords, step_length, dwi):
        predictor = self.mw_coords.copy()                # EulerPC.step sets the predictor point before this call
        x = super().corrector_step(init_mw_coords, step_length, dwi)
        got = float(np.linalg.norm(x - init_mw_coords))
        if got < 0.5 * step_length:
            self.fallbacks.append({"cycle": int(self.cur_cycle), "corrected_length": got})
            print(f"corrector oscillation at cycle {self.cur_cycle}: corrected step {got:.4f} < half the step; "
                  f"predictor point kept", flush=True)
            return predictor
        return x


def trans_rot_projector(cart_bohr: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """I − QQᵀ in mass-weighted Cartesians of the current frame, Q = orthonormal basis of the three
    translations and three rotations about the centre of mass. (pysisyphus's get_hessian_projector
    works in the standard orientation, not in the current frame, so it cannot be applied directly.)"""
    x = cart_bohr.reshape(-1, 3)
    m = np.asarray(masses)
    r = x - (m[:, None] * x).sum(0) / m.sum()
    sq = np.sqrt(m)[:, None]
    vecs = []
    for k in range(3):
        e = np.zeros(3)
        e[k] = 1.0
        vecs.append((sq * e).ravel())                     # translation
        vecs.append((sq * np.cross(e, r)).ravel())        # rotation
    U, s, _ = np.linalg.svd(np.array(vecs).T, full_matrices=False)
    Q = U[:, s > 1e-6 * s.max()]
    return np.eye(x.size) - Q @ Q.T


class GSmw(OwnStop, GonzalezSchlegel):
    """Fixes of 2026-11-13 (tail_gs_reactant "converged" after one step; replayed offline from its
    recorded Hessian and gradients): (1) the BFGS update after the first micro-step turned the lowest
    curvature from −0.015 to −0.62 au (BFGS is unstable where the curvature is negative, as in the flat
    tails), and (2) scipy's Newton search for λ then diverged; pysisyphus catches that and marks the run
    converged. Here: Bofill update (as EulerPC), translations and rotations projected out of the Hessian
    and the gradient, λ found by a bracketed search below the lowest eigenvalue (on that interval the
    constraint function rises monotonically from −r² to +∞, so the root is unique), and a failed search
    is a failure of the side, never convergence."""

    def micro_step(self, counter):
        from pysisyphus.optimizers.hessian_updates import bofill_update
        from scipy.optimize import brentq
        P = trans_rot_projector(self.geometry.cart_coords, self.geometry.masses)
        gradient = P @ self.mw_gradient
        coords_diff = P @ (self.mw_coords - self.prev_coords)
        gradient_diff = gradient - P @ self.prev_grad      # GS.step stores the unprojected start gradient
        self.prev_coords = self.mw_coords.copy()
        self.prev_grad = gradient.copy()
        if self.hessian_recalc and (self.micro_counter % self.hessian_recalc == 0):
            self.mw_hessian = self.geometry.mw_hessian
        elif np.linalg.norm(coords_diff) > 1e-10:
            dH, _ = bofill_update(self.mw_hessian, coords_diff, gradient_diff)
            self.mw_hessian += dH
        H = P @ self.mw_hessian @ P
        H = 0.5 * (H + H.T)
        w, V = np.linalg.eigh(H)
        big = np.abs(w) > 1e-8
        bw, bV = w[big], V[:, big]
        gs, ds = bV.T @ gradient, bV.T @ self.displacement
        r2 = (0.5 * self.step_length) ** 2

        def dx_of(lam):
            return -(gs - lam * ds) / (bw - lam)

        def f(lam):
            p = ds + dx_of(lam)
            return p @ p - r2

        b0 = bw[0]
        hi = b0 - 1e-10 * max(1.0, abs(b0))
        if not f(hi) > 0:
            raise IntegratorFailure(f"GS micro-step {counter}: no root below the lowest eigenvalue ({b0:.3e})")
        lo, k = b0 - 1.0, 0
        while f(lo) >= 0:
            lo, k = b0 - 2 * (b0 - lo), k + 1
            if k > 200:
                raise IntegratorFailure(f"GS micro-step {counter}: λ bracket not found")
        lam = brentq(f, lo, hi, xtol=1e-14, maxiter=500)
        dx = bV @ dx_of(lam)
        self.displacement += dx
        self.mw_coords += dx
        self.micro_counter += 1
        self.micro_log.append({"cycle": int(self.cur_cycle), "micro": int(counter), "lambda": float(lam),
                               "lowest_eigenvalue": float(b0), "dx_norm": float(np.linalg.norm(dx))})
        return dx, self.perp_component(gradient, self.displacement)

    def step(self):
        n0 = len(self.micro_log)
        super().step()
        if self.converged and not self.own_stop:   # pysisyphus's "Constrained search did not converge!"
            raise IntegratorFailure("GS constrained search did not converge (pysisyphus would call this converged)")
        if len(self.micro_log) - n0 >= self.max_micro_cycles:
            self.micro_max_hit.append(int(self.cur_cycle))


def main(name: str, direction: str, step: float = 0.10, front: bool = True) -> dict:
    cfg = RUNS[name]
    wd = OUT / name / direction
    wd.mkdir(parents=True, exist_ok=False)
    (wd / "pysis").mkdir()
    os.chdir(wd / "pysis")
    ts = read_xyz(RX / "01_ts/ts/ts.xyz")
    sym = list(ts.symbols)
    downhill = direction == "downhill"
    if downhill:
        side, d = cfg["start"]
        fr = read_xyz_frames(SC / "true_irc/s0050" / d / f"{side}_trj.xyz")
        start = np.asarray(fr[TAIL_START[side]].coords)
    else:
        start = np.asarray(ts.coords)
    geom = Geometry(sym, start.ravel() / BOHR, isotopes=cfg.get("isotopes"))
    calc = OrcaQueue(wd, front, out_dir=str(wd / "pysis"))
    geom.set_calculator(calc)
    kw = dict(step_length=step, max_cycles=cfg["max_cycles"], rms_grad_thresh=1e-12, energy_thresh=-1.0,
              hessian_init="calc", out_dir=str(wd / "pysis"))
    if downhill:
        kw.update(downhill=True, forward=False, backward=False)
    else:
        kw.update(forward=direction == "forward", backward=direction == "backward", displ="length", displ_length=step)
    if cfg["integ"] == "eulerpc":
        irc = EulerPCmw(geom, max_pred_steps=5000, hessian_recalc=cfg["recalc"], **kw)
    else:
        irc = GSmw(geom, hessian_recalc=cfg["recalc"], **kw)
    irc.rx = Reaction.load(RX)
    irc.ts_frame = PGeom(sym, np.asarray(ts.coords), f"TS E {E_TS:.12f}")
    irc.own_stop = None
    irc.fallbacks, irc.micro_log, irc.micro_max_hit = [], [], []
    t0 = time.time()
    failure = None
    try:
        irc.run()
    except IntegratorFailure as e:                 # keep the points before the failure; the side is incomplete
        failure = str(e)
    pre = direction
    if failure is None:
        C = [np.asarray(c) for c in getattr(irc, f"{pre}_coords")]
        E = list(getattr(irc, f"{pre}_energies"))
        G = [np.asarray(g) for g in getattr(irc, f"{pre}_gradients")]
    else:                                          # TS-outward order, as integrated
        n = min(len(irc.irc_coords), len(irc.irc_energies), len(irc.irc_gradients))
        C = [np.asarray(c) for c in irc.irc_coords[:n]]
        E = list(irc.irc_energies[:n])
        G = [np.asarray(g) for g in irc.irc_gradients[:n]]
    if not downhill:
        x0 = start.ravel() / BOHR
        if np.linalg.norm(C[-1] - x0) < np.linalg.norm(C[0] - x0):
            C, E, G = C[::-1], E[::-1], G[::-1]
    ended = (f"FAILED: {failure}" if failure else "own stopping rule" if irc.own_stop else
             "energy increased" if getattr(irc, f"{pre}_energy_increased", False) else
             "integrator indicated convergence" if getattr(irc, f"{pre}_is_converged", False) else "max cycles")
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C, E):
            fh.write(f"{len(sym)}\ntrue IRC {name} {direction} E {e:.12f}\n")
            fh.write("".join(f"{a:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n"
                             for a, (x, y, z) in zip(sym, c.reshape(-1, 3) * BOHR)))
    res = {"run": name, "direction": direction, "settings": {k: v for k, v in cfg.items()}, "step_bohr_amu05": step,
           "ended": ended, "own_stop": irc.own_stop, "n_points": len(C),
           "n_engrad": sum(c["kind"] == "EnGrad" for c in calc.calls), "n_hessian": sum(c["kind"] == "Freq" for c in calc.calls),
           "calls": calc.calls, "wall_total_s": time.time() - t0, "energies": E, "gradients": [g.tolist() for g in G],
           "masses": [float(m) for m in geom.masses],
           "code": "2026-11-13 fixes (EulerPC zero-step fallback; GS Bofill, projection, bracketed λ)",
           "corrector_fallbacks": irc.fallbacks, "gs_micro_steps": len(irc.micro_log),
           "gs_micro_cycle_limit_hit_at": irc.micro_max_hit, "gs_micro_log": irc.micro_log}
    (wd / "path.json").write_text(json.dumps(res), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k not in ("gradients", "energies", "calls")}))
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], float(a[2]) if len(a) > 2 else 0.10, front=not (len(a) > 3 and a[3] == "fifo"))
