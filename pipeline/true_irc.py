"""Mass-weighted IRC test on SN2 with pysisyphus (J. Picione, 2026-11-12: "true IRC, feasibility first").

Runs in its own environment (.venv-pysis: Python 3.12 from uv's user folder, pysisyphus 1.0.0, GPL-3;
the pipeline's own .venv is Python 3.15, for which h5py, needed by pysisyphus, has no package). The
pipeline itself is imported unchanged: every energy and gradient is one ORCA 6.1.1 job, alone in its own
new folder, through the core ledger (pipeline/orca.py, pipeline/cores.py).

Integrator: EulerPC (Hratchian–Schlegel Euler predictor–corrector, Bofill Hessian updates) in
mass-weighted Cartesian coordinates. pysisyphus measures step_length in plain (un-mass-weighted)
Cartesian bohr (EulerPC.get_integration_length), like ORCA's Scale_Displ_SD; only the first displacement
(displ_length) is in bohr·amu½. The step target in bohr·amu½ (0.10, or 0.05 for the step test of this
method) is therefore converted with the same calibration as the ORCA step-test paths (0.10 -> 0.02359
bohr, 0.05 -> 0.01179 bohr), so both methods take the same Cartesian steps; the mass-weighted steps
are measured afterwards. First step 0.10 bohr·amu½ along the mass-weighted TS mode (displ "length").
max_pred_steps 5000 (pysisyphus default 500): with 500 Euler substeps the predictor did not reach the
step length in the flat tail and the path stalled (attempt 1). Start Hessian: ORCA's TS Hessian
(01_ts/ts/ts.hess, same level). End: RMS gradient below 5e-4 Eh/bohr (ORCA's TolRMSG) after the
inflection point, an energy rise, or 1500 steps. Each run writes its own pysisyphus log (irc.log in its
pysis folder).
Level: the B3LYP step-test level, B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX, SMD water, 2 cores, with VeryTightSCF
(decided for production 2026-11-12; pysisyphus ends a direction at the first energy rise, so SCF noise matters),
so the path compares directly with the four ORCA step-test paths.
Output: <step_convergence>/true_irc/s<step>/<direction>/: g####/ (one ORCA job each), path.xyz (TS-ordered
frames, comment = energy in Eh), path.json (gradients, timing).
Usage: .venv-pysis/Scripts/python.exe -X utf8 -m pipeline.true_irc forward|backward [step, default 0.10]
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.irc import EulerPC

from . import orca
from .geometry import read_xyz
from .mass_weight_evidence import read_hess_matrix

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions/sn2_hydroxide_bromomethane_masvp"
OUT = RX / "02_irc/attempts/step_convergence/true_irc"
LEVEL = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX VeryTightSCF CPCM(water)"   # VeryTightSCF: pysisyphus stops a direction at the first energy rise
SMD = '%cpcm smd true SMDsolvent "WATER" end'
CHARGE, MULT, NPROCS = -1, 1, 2
E_TS = -2689.43410109009


class OrcaQueue(Calculator):
    """Energy + gradient from one ORCA job per call, in a new folder, through the core ledger."""

    def __init__(self, wd: Path, hess: Path, **kw):
        super().__init__(charge=CHARGE, mult=MULT, **kw)
        self.wd, self.hess = wd, hess
        self.n = 0
        self.log_ = orca.RunLog(wd / "run_log_true_irc.md")
        self.calls = []

    def _engrad(self, atoms, coords):
        wd = self.wd / f"g{self.n:04d}"
        self.n += 1
        wd.mkdir(parents=True, exist_ok=False)
        c = np.asarray(coords).reshape(-1, 3)
        geo = "\n".join(f"{a:2s} {x:20.12f} {y:20.12f} {z:20.12f}" for a, (x, y, z) in zip(atoms, c))
        inp = (f"! {LEVEL} EnGrad Bohrs NoPop\n%pal nprocs {NPROCS} end\n%maxcore 3000\n{SMD}\n\n"
               f"* xyz {CHARGE} {MULT}\n{geo}\n*\n")
        r = orca.run(wd, "g", inp, self.log_, note=f"true IRC test (pysisyphus EulerPC), call {self.n - 1}",
                     front=True)
        if not r.ok:
            raise orca.OrcaError(f"{wd} did not end normally")
        t = (wd / "g.engrad").read_text().splitlines()
        vals = [ln.strip() for ln in t if ln.strip() and not ln.startswith("#")]
        nat = int(vals[0])
        E = float(vals[1])
        g = np.array([float(v) for v in vals[2:2 + 3 * nat]])
        out = (wd / "g.out").read_text(errors="replace")
        m = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", out)
        wall = (int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]) + int(m[5]) / 1000) if m else None
        self.calls.append({"folder": wd.name, "E": E, "wall_s": wall})
        for f in wd.glob("*.tmp*"):
            f.unlink()
        return E, g

    def get_energy(self, atoms, coords, **kw):
        E, _ = self._engrad(atoms, coords)
        return {"energy": E}

    def get_forces(self, atoms, coords, **kw):
        E, g = self._engrad(atoms, coords)
        return {"energy": E, "forces": -g}

    def get_hessian(self, atoms, coords, **kw):
        """ORCA's TS Hessian at the same level; only valid at the TS geometry."""
        ts = np.asarray(read_xyz(RX / "01_ts/ts/ts.xyz").coords).ravel() / 0.529177210903
        if np.max(np.abs(np.asarray(coords) - ts)) > 1e-5:
            raise RuntimeError("Hessian requested away from the TS; only the TS Hessian is available")
        return {"energy": E_TS, "hessian": read_hess_matrix(RX / "01_ts/ts/ts.hess")}


CART_PER_MW = {0.10: 0.02359, 0.05: 0.01179}   # bohr per target bohr·amu½: the ORCA step-test settings


def main(direction: str, step: float = 0.10) -> dict:
    wd = OUT / f"s{step:.3f}".replace(".", "") / direction
    wd.mkdir(parents=True, exist_ok=False)
    (wd / "pysis").mkdir()
    import os
    os.chdir(wd / "pysis")                          # pysisyphus writes irc.log to the working directory
    ts = read_xyz(RX / "01_ts/ts/ts.xyz")
    geom = Geometry(list(ts.symbols), np.asarray(ts.coords).ravel() / 0.529177210903)
    calc = OrcaQueue(wd, RX / "01_ts/ts/ts.hess", out_dir=str(wd / "pysis"))
    geom.set_calculator(calc)
    irc = EulerPC(geom, step_length=CART_PER_MW[step], max_cycles=1500, max_pred_steps=5000,
                  forward=direction == "forward",
                  backward=direction == "backward", displ="length", displ_length=0.10,
                  hessian_init="calc", rms_grad_thresh=5e-4, energy_thresh=1e-8, out_dir=str(wd / "pysis"))
    irc.run()
    C = [np.asarray(c) for c in getattr(irc, f"{direction}_coords")]
    E = list(getattr(irc, f"{direction}_energies"))
    G = [np.asarray(g) for g in getattr(irc, f"{direction}_gradients")]
    x0 = np.asarray(ts.coords).ravel() / 0.529177210903
    if np.linalg.norm(C[-1] - x0) < np.linalg.norm(C[0] - x0):    # order from the TS outwards
        C, E, G = C[::-1], E[::-1], G[::-1]
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C, E):
            xyz = c.reshape(-1, 3) * 0.529177210903
            fh.write(f"{len(ts.symbols)}\ntrue IRC {direction} E {e:.12f}\n")
            fh.write("".join(f"{a:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n" for a, (x, y, z) in zip(ts.symbols, xyz)))
    res = {"direction": direction, "integrator": "pysisyphus 1.0.0 EulerPC", "step_target_bohr_amu05": step,
           "step_length_bohr_cartesian": CART_PER_MW[step], "max_pred_steps": 5000,
           "first_step": 0.10, "n_points": len(C), "converged": bool(getattr(irc, f"{direction}_is_converged")),
           "n_orca_calls": calc.n, "wall_s_per_call": [c["wall_s"] for c in calc.calls],
           "gradients": [g.tolist() for g in G], "energies": E}
    (wd / "path.json").write_text(json.dumps(res), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k not in ("gradients", "energies", "wall_s_per_call")}))
    return res


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 0.10)
