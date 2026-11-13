"""Step 1 for the hydride shift: is the 3-methyl-2-butyl cation a true minimum?

Three starting conformers, differing in how the C3–H bond (the hydrogen that would migrate) sits
relative to the empty p orbital on C2 (dihedral H–C3–C2 against the C2 plane):
  perp   C3–H parallel to the empty p orbital (best hyperconjugation, most prone to shift)
  eclip  C3–H in the C2 plane, eclipsing C2–H
  anti   C3–H in the C2 plane, anti to C2–H
Each is optimized with frequencies (Opt Freq, working level) in SMD dichloromethane and in the
gas phase, from a classical (unbridged) start. Each result is classified from the C2···H and
C3–H distances of the migrating hydrogen and the number of imaginary frequencies.

Run from the repository folder:
  .venv/Scripts/python.exe -X utf8 reactions/hshift_3_methyl_2_butyl/00_inputs/check_secondary_minimum.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))
from pipeline import orca                              # noqa: E402
from pipeline.geometry import Geometry, read_xyz         # noqa: E402
from pipeline.steps import Reaction                      # noqa: E402

rx = Reaction.load(ROOT)
CONFS = {"perp": 90.0, "eclip": 0.0, "anti": 180.0}      # dihedral O(=H new)–C2–C3–H, degrees


def start_geometry(dih: float) -> Geometry:
    """3-methylbutan-2-one with the C3–H dihedral set, then the carbonyl O replaced by H at
    1.09 Å: a planar sp2 C2 with CH3, H and C3, i.e. the classical secondary cation."""
    m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)C(C)C"))  # 0 C1, 1 C2, 2 O, 3 C3, 4 C4, 5 C5
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m)
    h3 = next(a.GetIdx() for a in m.GetAtomWithIdx(3).GetNeighbors() if a.GetSymbol() == "H")
    conf = m.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, 2, 1, 3, h3, dih)
    xyz = np.array(conf.GetPositions())
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    v = xyz[2] - xyz[1]
    xyz[2] = xyz[1] + 1.09 * v / np.linalg.norm(v)
    sym[2] = "H"
    # put the migrating H first among the H atoms is not needed; record its index instead
    return Geometry(sym, xyz), h3


def classify(g: Geometry, h: int) -> dict:
    d2, d3 = g.distance(1, h), g.distance(3, h)
    # bonded C–H ≈ 1.10 Å; a bridging H sits at roughly equal C···H distances of 1.2–1.5 Å
    if abs(d2 - d3) < 0.25 and max(d2, d3) < 1.6:
        kind = "hydrogen-bridged"
    elif d3 < 1.2 and d2 > 1.6:
        kind = "classical secondary (H still on C3)"
    elif d2 < 1.2 and d3 > 1.6:
        kind = "tertiary 2-methyl-2-butyl (H moved to C2)"
    else:
        kind = "partly bridged / unclear"
    return {"C2_H_A": round(d2, 3), "C3_H_A": round(d3, 3), "kind": kind,
            "C1_C2_C3_angle_deg": round(g.angle(0, 1, 3), 1)}


def main():
    out = {}
    for phase in ("dcm", "gas"):
        for name, dih in CONFS.items():
            g0, h = start_geometry(dih)
            wd = rx.d(f"01_check/{phase}/{name}")
            g0.write(wd / "start.xyz", f"3-methyl-2-butyl cation start, dihedral {dih:.0f}")
            if phase == "dcm":
                inp = rx.inp("Opt Freq", g0)
            else:
                inp = orca.make_input(f'{rx.cfg["method"]} Opt Freq', g0, 1, 1,
                                      nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
            r = orca.run(wd, name, inp, rx.log,
                         note=f"secondary-cation minimum check, {phase}, start {name}")
            if not r.ok:
                out[f"{phase}/{name}"] = {"ok": False}
                continue
            g = read_xyz(wd / f"{name}.xyz")
            out[f"{phase}/{name}"] = {"ok": True, "migrating_H": h, "E_eh": r.energy_eh,
                                      "G_eh": r.gibbs_eh, "opt_converged": r.opt_converged,
                                      "imaginary_cm": r.imaginary, **classify(g, h)}
            print(phase, name, out[f"{phase}/{name}"], flush=True)
    (ROOT / "01_check" / "secondary_minimum.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
