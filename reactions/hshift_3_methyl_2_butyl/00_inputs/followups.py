"""Follow-ups on the 3-methyl-2-butyl cation (J. Picione, 2026-11-10), run after the fallback
TS and IRC are reported.

  method   Re-optimize the turned-away minima ("anti", "eclip"; gas phase) at M06-2X/ma-def2-SVP
           and RI-MP2/ma-def2-SVP (Opt Freq); report the stretched C3–CH₃ bond and the
           C2–C3–CH₃ angle next to B3LYP (does the partial methyl bridging survive?).
  scan     Relaxed scan of the C1–C2–C3–H dihedral from each turned-away minimum toward the
           aligned conformation (−90°), working level, DCM and gas phase, to locate a rotation
           barrier before the hydrogen shifts.
  ts       OptTS Freq from the highest scan point of one scan (argument: phase/start), then the IRC
           both ways if exactly one imaginary mode.

Atoms (0-based): 0 C1, 1 C2 (cationic), 2 H on C2, 3 C3, 4 and 5 methyl C on C3, 9 the H on C3.
Run: .venv/Scripts/python.exe -X utf8 reactions/hshift_3_methyl_2_butyl/00_inputs/followups.py <stage> [...]
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))

from pipeline import orca                     # noqa: E402
from pipeline.geometry import Geometry, read_xyz   # noqa: E402
from pipeline.steps import Reaction           # noqa: E402
from pipeline.units import kj                 # noqa: E402

CHECK = ROOT / "01_check"
METHODS = {"M06-2X": "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX TightSCF",
           "RI-MP2": "RI-MP2 ma-def2-SVP AutoAux TightSCF"}
H = 9


def inp(rx, kw, g, phase="gas", blocks=""):
    """Working-level style input; dcm = exactly the solvent setup of the step-1 jobs."""
    if phase == "dcm":
        return rx.inp(kw.replace(rx.cfg["method"], "").strip(), g, blocks=blocks)
    return orca.make_input(kw, g, rx.cfg["charge"], rx.cfg["mult"], blocks=blocks,
                           nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])


def measures(g):
    x = np.asarray(g.coords)
    d = lambda i, j: float(np.linalg.norm(x[i] - x[j]))
    def ang(i, j, k):
        a, b = x[i] - x[j], x[k] - x[j]
        return float(np.degrees(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))))
    me = max((4, 5), key=lambda m: d(3, m))           # the stretched C3–CH₃
    return {"stretched_methyl": me, "C3_CH3_A": round(d(3, me), 3),
            "other_C3_CH3_A": round(d(3, 9 - me), 3), "C2_C3_CH3_deg": round(ang(1, 3, me), 1),
            "C2_CH3_A": round(d(1, me), 3), "C3_H_A": round(d(3, H), 3), "C2_H_A": round(d(1, H), 3)}


def method_check(rx):
    out = {}
    for start in ("anti", "eclip"):
        g = read_xyz(CHECK / f"gas/{start}/{start}.xyz")
        out[f"B3LYP/{start}"] = measures(g)
        for lab, kw in METHODS.items():
            wd = CHECK / "method" / lab / start
            job = " Opt NumFreq" if "MP2" in lab else " Opt Freq"   # no analytic MP2 Hessian in ORCA
            r = orca.run(wd, start, inp(rx, kw + job, g), rx.log,
                         note=f"method check, gas phase, {lab}, from the B3LYP {start} minimum")
            name = start
            if r.ok and r.imaginary:
                # converged onto a saddle point: push along the imaginary mode (largest atom
                # 0.1 Å) and re-optimize once, same keywords (as the pipeline does for minima)
                h = orca.read_hess(wd / f"{start}.hess")
                v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
                gp = read_xyz(wd / f"{start}.xyz")
                gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
                name = f"{start}_min"
                r = orca.run(wd, name, inp(rx, kw + job, gp), rx.log,
                             note=f"rerun: {start} converged to a saddle point (imaginary "
                                  f"{[round(f, 1) for f in r.imaginary]} cm-1); displaced 0.1 Å along "
                                  "that mode, same keywords")
            row = {"ok": r.ok, "opt_converged": r.opt_converged, "E_eh": r.energy_eh,
                   "imaginary_cm": r.imaginary, "job": name}
            if (wd / f"{name}.xyz").exists():
                row.update(measures(read_xyz(wd / f"{name}.xyz")))
            out[f"{lab}/{start}"] = row
            print(lab, start, row, flush=True)
    (CHECK / "method_check.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")


def scan(rx):
    out = {}
    for phase in ("gas", "dcm"):
        for start in ("anti", "eclip"):
            g = read_xyz(CHECK / f"{phase}/{start}/{start}.xyz")
            x = np.asarray(g.coords)
            b0, b1, b2 = x[0] - x[1], x[3] - x[1], x[H] - x[3]
            b1 = b1 / np.linalg.norm(b1)
            v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
            phi0 = float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))
            target = -90.0
            n = int(round(abs(target - phi0) / 6.0)) + 1      # about 6° per step
            blocks = f"%geom Scan D 0 1 3 {H} = {phi0:.1f}, {target:.1f}, {n} end end"
            wd = CHECK / "scan" / phase / start
            r = orca.run(wd, "scan", inp(rx, rx.cfg["method"] + " Opt", g, phase, blocks), rx.log,
                         note=f"relaxed scan C1–C2–C3–H from {phi0:.0f}° to {target:.0f}°, {phase}, start {start}")
            dat = wd / "scan.relaxscanact.dat"
            pts = [list(map(float, ln.split())) for ln in dat.read_text().splitlines() if ln.strip()] if dat.exists() else []
            out[f"{phase}/{start}"] = {"ok": r.ok, "phi0_deg": phi0, "points": pts}
            if pts:
                E = np.array([p[1] for p in pts])
                print(phase, start, "dihedral, ΔE (kJ/mol):",
                      [(round(p[0]), round(kj(e - E[0]), 1)) for p, e in zip(pts, E)], flush=True)
    (CHECK / "scan.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    rx = Reaction.load(ROOT)
    {"method": method_check, "scan": scan}[sys.argv[1]](rx)
