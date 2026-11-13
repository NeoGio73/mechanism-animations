"""Cases for tests/test_pysisyphus_patches.py, run in a subprocess so that either the installed
(unmodified) pysisyphus 1.0.0 or a patched copy (put first on PYTHONPATH) is imported.

  python tests/pysis_patch_cases.py <case>     prints one JSON object

Cases
  toy_gs_override       pipeline/true_irc2.GSmw on the installed package (the code the runs use)
  toy_gs_patched        GonzalezSchlegel from the patched package, stop_func as the override's stop
  toy_epc_override      pipeline/true_irc2.EulerPCmw on the installed package
  toy_epc_patched       EulerPC from the patched package
  replay_original       recorded tail_gs_reactant first step (2026-11-13, before the fixes), pysisyphus 1.0.0
  replay_patched        the same recorded data through the patched micro-step
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(1, str(ROOT))          # pipeline.* (after a patched pysisyphus, if one is first)

from pysisyphus.calculators.Calculator import Calculator   # noqa: E402
from pysisyphus.Geometry import Geometry                     # noqa: E402

STOP_AFTER = 9         # toy runs: the "own stopping rule" stand-in ends a direction after this many points


# --------------------------------------------------------------------------------- toy potential
class Toy(Calculator):
    """Rotation- and translation-invariant model potential for four atoms (C, O, H, H): a Morse term
    for every pair with its own D, a, r0. Analytic gradient; Hessian by central differences of it."""
    P = {(0, 1): (0.30, 1.1, 2.3), (0, 2): (0.17, 0.9, 2.1), (0, 3): (0.15, 1.0, 2.0),
         (1, 2): (0.05, 0.8, 3.0), (1, 3): (0.04, 0.7, 3.4), (2, 3): (0.02, 0.6, 3.3)}

    def __init__(self, **kw):
        super().__init__(charge=0, mult=1, **kw)

    def eg(self, x):
        x = np.asarray(x, float).reshape(-1, 3)
        E, g = 0.0, np.zeros_like(x)
        for (i, j), (D, a, r0) in self.P.items():
            d = x[i] - x[j]
            r = np.linalg.norm(d)
            q = np.exp(-a * (r - r0))
            E += D * (1 - q) ** 2
            dEdr = 2 * D * (1 - q) * a * q
            g[i] += dEdr * d / r
            g[j] -= dEdr * d / r
        return E, g.ravel()

    def get_energy(self, atoms, coords, **kw):
        return {"energy": self.eg(coords)[0]}

    def get_forces(self, atoms, coords, **kw):
        E, g = self.eg(coords)
        return {"energy": E, "forces": -g}

    def get_hessian(self, atoms, coords, **kw):
        c = np.asarray(coords, float).ravel()
        h, H = 1e-5, np.zeros((c.size, c.size))
        for k in range(c.size):
            e = np.zeros(c.size)
            e[k] = h
            H[k] = (self.eg(c + e)[1] - self.eg(c - e)[1]) / (2 * h)
        return {"energy": self.eg(c)[0], "hessian": 0.5 * (H + H.T)}


START = np.array([[0.00, 0.00, 0.00], [2.95, 0.30, 0.10], [-0.80, 2.40, 0.20], [-0.60, -1.30, 1.70]])  # bohr
ATOMS = ["C", "O", "H", "H"]


def toy_geom():
    g = Geometry(ATOMS, START.ravel().copy())
    g.set_calculator(Toy())
    return g


def kw(out):
    return dict(step_length=0.05, max_cycles=40, rms_grad_thresh=1e-12, energy_thresh=-1.0,
                hessian_init="calc", downhill=True, out_dir=out)


def result(irc, extra=None):
    C = np.asarray(irc.downhill_coords)
    r = {"n": len(C), "coords": C.tolist(), "energies": list(map(float, irc.downhill_energies)),
         "converged": bool(irc.downhill_is_converged)}
    r.update(extra or {})
    return r


def stop_count(irc):
    return len(irc.irc_energies) >= STOP_AFTER


def toy_gs_override(out):
    from pipeline.true_irc2 import GSmw

    class T(GSmw):
        def get_additional_print(self):             # OwnStop with the stand-in rule
            if stop_count(self):
                self.converged = True
                self.own_stop = {"stand_in": True}
            return ""
    irc = T(toy_geom(), **kw(out))
    irc.own_stop, irc.micro_log, irc.micro_max_hit = None, [], []
    irc.run()
    return result(irc, {"micro_steps": len(irc.micro_log)})


def toy_gs_patched(out):
    from pysisyphus.irc.GonzalezSchlegel import GonzalezSchlegel
    irc = GonzalezSchlegel(toy_geom(), stop_func=stop_count, **kw(out))
    irc.run()
    return result(irc)


def toy_epc_override(out):
    from pipeline.true_irc2 import EulerPCmw

    class T(EulerPCmw):
        def get_additional_print(self):
            if stop_count(self):
                self.converged = True
                self.own_stop = {"stand_in": True}
            return ""
    irc = T(toy_geom(), max_pred_steps=5000, **kw(out))
    irc.own_stop, irc.fallbacks = None, []
    irc.run()
    return result(irc, {"fallbacks": irc.fallbacks})


def toy_epc_patched(out):
    from pysisyphus.irc.EulerPC import EulerPC
    irc = EulerPC(toy_geom(), max_pred_steps=5000, stop_func=stop_count, **kw(out))
    irc.run()
    return result(irc, {"fallbacks": irc.fallbacks})


# --------------------------------------------------------------------------------- recorded replay
REC = (ROOT / "reactions/sn2_hydroxide_bromomethane_masvp/02_irc/attempts/step_convergence/true_irc2/"
       "failed_before_fixes_1113/tail_gs_reactant/downhill")


def _job(folder):
    inp = (REC / folder / "j.inp").read_text().splitlines()
    i = next(k for k, ln in enumerate(inp) if ln.startswith("* xyz"))
    rows = [ln.split() for ln in inp[i + 1:] if ln.strip() and ln.strip() != "*"]
    sym = [r[0] for r in rows]
    x = np.array([[float(v) for v in r[1:4]] for r in rows]).ravel()           # Bohrs (keyword in the input)
    assert "Bohrs" in inp[0]
    return sym, x


def _engrad(folder):
    v = [ln.strip() for ln in (REC / folder / "j.engrad").read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    n = int(v[0])
    return float(v[1]), np.array([float(t) for t in v[2:2 + 3 * n]])


class Recorded(Calculator):
    """Energies and gradients of the recorded ORCA jobs, looked up by geometry; Hessian = h0001."""

    def __init__(self, **kw):
        super().__init__(charge=-1, mult=1, **kw)
        from pipeline.mass_weight_evidence import read_hess_matrix
        self.rec = []
        for f in ("g0000", "g0002", "g0003"):
            _, x = _job(f)
            E, g = _engrad(f)
            self.rec.append((f, x, E, g))
        self.H = read_hess_matrix(REC / "h0001/j.hess")
        self.x0 = self.rec[0][1]

    def _find(self, coords):
        c = np.asarray(coords).ravel()
        for f, x, E, g in self.rec:
            if np.max(np.abs(c - x)) < 1e-6:
                return f, E, g
        raise KeyError("no recorded job at this geometry")

    def get_energy(self, atoms, coords, **kw):
        return {"energy": self._find(coords)[1]}

    def get_forces(self, atoms, coords, **kw):
        _, E, g = self._find(coords)
        return {"energy": E, "forces": -g}

    def get_hessian(self, atoms, coords, **kw):
        assert np.max(np.abs(np.asarray(coords).ravel() - self.x0)) < 1e-6
        return {"energy": self.rec[0][2], "hessian": self.H}


def _replay_setup(out):
    from pysisyphus.irc.GonzalezSchlegel import GonzalezSchlegel
    sym, x0 = _job("g0000")
    geom = Geometry(sym, x0)
    geom.set_calculator(Recorded())
    gs = GonzalezSchlegel(geom, step_length=0.10, max_cycles=1, rms_grad_thresh=1e-12, energy_thresh=-1.0,
                          hessian_init="calc", downhill=True, out_dir=out)
    gs.init_hessian = geom.calculator.H            # as IRC.run (hessian_init="calc" at the start point)
    gs.prepare("downhill")
    return gs


def _replay_step_start(gs):
    """GonzalezSchlegel.step up to the micro-cycles (unchanged by the patches)."""
    grad0 = gs.mw_gradient
    gs.prev_grad, gs.prev_coords = grad0, gs.mw_coords
    pivot_step = 0.5 * gs.step_length * -grad0 / np.linalg.norm(grad0)
    pivot = gs.mw_coords + pivot_step
    gs.mw_coords = pivot + pivot_step
    gs.displacement = pivot_step
    return pivot


def _lowest(H, gs, P=None):
    if P is not None:
        H = P @ H @ P
    w = np.linalg.eigvalsh(0.5 * (H + H.T))
    return float(w[np.abs(w) > 1e-8][0])


def replay_original(out):
    """pysisyphus 1.0.0: micro-step 0 at the recorded first guess (g0002) must land on g0003; micro-step 1
    (gradient g0003) then updates with BFGS and searches λ with Newton."""
    gs = _replay_setup(out)
    pivot = _replay_step_start(gs)
    low0 = _lowest(gs.mw_hessian, gs)
    gs.micro_step(0)
    hit = gs.geometry.calculator._find(gs.coords)[0]            # replay reproduces the recorded run
    try:
        gs.micro_step(1)
        outcome = "λ found"
    except RuntimeError as e:
        outcome = f"RuntimeError: {str(e)[:120]}"
    return {"micro0_lands_on": hit, "lowest_start": low0, "lowest_after_update": _lowest(gs.mw_hessian, gs),
            "micro1": outcome}


def replay_patched(out):
    """Same inputs as micro-step 1 of the original run (state after the original micro-step 0, gradient
    g0003), through the patched micro-step: Bofill update, projection, bracketed λ."""
    import importlib
    M = importlib.import_module("pysisyphus.irc.GonzalezSchlegel")
    gs = _replay_setup(out)
    _replay_step_start(gs)
    _, x2 = _job("g0003")
    # state after the recorded micro-step 0 (pysisyphus 1.0.0): reproduce it with the original update
    from pysisyphus.optimizers.hessian_updates import bfgs_update
    g1 = gs.mw_gradient                                            # at the first guess (g0002)
    dH, _ = bfgs_update(gs.mw_hessian, gs.mw_coords - gs.prev_coords, g1 - gs.prev_grad)
    H_after0 = gs.mw_hessian + dH
    x2_mw = x2 * np.sqrt(np.repeat(gs.geometry.masses, 3))
    gs.displacement = gs.displacement + (x2_mw - gs.mw_coords)
    gs.prev_coords, gs.prev_grad = gs.mw_coords.copy(), g1.copy()
    gs.mw_coords = x2_mw
    gs.mw_hessian = H_after0
    P = M.trans_rot_projector(gs.geometry.cart_coords, gs.geometry.masses)
    low_before = _lowest(gs.mw_hessian, gs, P)
    dx, _ = gs.micro_step(1)
    p = gs.displacement
    # the constraint holds in the space without translations and rotations of the frame where the
    # micro-step was taken (P, computed there above)
    P2 = P
    return {"lowest_before_update": low_before, "lowest_after_update": _lowest(gs.mw_hessian, gs, P),
            "on_sphere_err": float(abs(np.linalg.norm(P2 @ p) - 0.5 * gs.step_length)),
            "full_norm_minus_r": float(np.linalg.norm(p) - 0.5 * gs.step_length), "dx_norm": float(np.linalg.norm(dx))}


if __name__ == "__main__":
    case = sys.argv[1]
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as out:
        os.chdir(out)
        print("RESULT " + json.dumps(globals()[case](out)))
