"""Gas-phase comparison for the degenerate 1,2-hydride shift (J. Picione, 2026-11-10: SMD
dichloromethane, with gas phase run alongside). Same working level and same starting structures
as the solution run, no solvent model: TS (constrained pre-optimization, then OptTS Freq) and the
tertiary cation (Opt Freq). Results go to 01_ts/gas_phase_opt/ and gas_phase_opt.json.

Run: .venv/Scripts/python.exe -X utf8 reactions/hshift_2_3_dimethyl_2_butyl/00_inputs/gas_phase_opt.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))

from pipeline import orca                                     # noqa: E402
from pipeline.geometry import read_xyz                        # noqa: E402
from pipeline.steps import Reaction, _constraint_block, check_ts_mode   # noqa: E402
from pipeline.units import kj                                 # noqa: E402


def inp(rx, kw, g, blocks=""):
    return orca.make_input(f"{rx.cfg['method']} {kw}", g, rx.cfg["charge"], rx.cfg["mult"],
                           blocks=blocks, nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])


def main():
    rx = Reaction.load(ROOT)
    base = rx.d("01_ts/gas_phase_opt")
    out = {}

    g = read_xyz(ROOT / "00_inputs" / rx.cfg["reactants"][0]["start_xyz"])
    r = orca.run(base / "cation", "cation", inp(rx, "Opt Freq", g), rx.log,
                 note="gas phase: tertiary cation, Opt Freq")
    out["cation"] = {"ok": r.ok, "opt_converged": r.opt_converged, "E_eh": r.energy_eh,
                     "G_eh": r.gibbs_eh, "imaginary_cm": r.imaginary}

    g = read_xyz(ROOT / "00_inputs" / rx.cfg["ts_guess"])
    r = orca.run(base / "preopt", "preopt", inp(rx, "Opt", g, _constraint_block(rx)), rx.log,
                 note="gas phase: TS guess relaxed with forming/breaking bonds frozen")
    g = read_xyz(base / "preopt" / "preopt.xyz")
    r = orca.run(base / "ts", "ts", inp(rx, "OptTS Freq", g, "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end"),
                 rx.log, note="gas phase: TS search")
    gts = read_xyz(base / "ts" / "ts.xyz")
    out["ts"] = {"ok": r.ok, "opt_converged": r.opt_converged, "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
                 "imaginary_cm": r.imaginary,
                 "C2_H_A": round(gts.distance(0, 2), 3), "C3_H_A": round(gts.distance(1, 2), 3),
                 "C2_H_C3_deg": round(gts.angle(0, 2, 1), 1),
                 "mode_check": check_ts_mode(rx, base / "ts" / "ts.hess", gts)}
    if out["cation"]["E_eh"] and out["ts"]["E_eh"]:
        out["dE_barrier_kJ"] = round(kj(out["ts"]["E_eh"] - out["cation"]["E_eh"]), 2)
    if out["cation"]["G_eh"] and out["ts"]["G_eh"]:
        out["dG_barrier_298K_kJ"] = round(kj(out["ts"]["G_eh"] - out["cation"]["G_eh"]), 2)
    (ROOT / "01_ts" / "gas_phase_opt.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
