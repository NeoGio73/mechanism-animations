"""SN2: per-species solvation check (J. Picione, 2026-11-19). Small jobs, no path rerun.

Purpose: find where the error in the SN2 reaction free energy comes from (DLPNO-CCSD(T) -85.3 kJ/mol against
-100.8 kJ/mol measured): the gas-phase reaction energy, or the solvation of particular species.

Jobs (this script):
  1. Gas phase, working level: the same keywords and settings as the SMD jobs of 01_ts/fragments, without the
     solvent (no CPCM keyword, no %cpcm block). Opt Freq for hydroxide, bromomethane and methanol, started from the
     solution geometry; a single point for bromide (one atom).
  2. Gas-phase single points on those geometries at every [[diagnostic_levels]] entry, same inputs as
     01_ts/diagnostics without the solvent.
Output: 03_solvation/gas/<species>/ (working level), 03_solvation/gas/<species>/sp/<level>/sp (single points),
03_solvation/gas.json (energies). The analysis (solvation free energies and the error split) is
00_inputs/solvation_analysis.py.

Free energies: G(level) = E(level) + [G - E] from the working-level frequencies of the same phase and geometry, as
for every species in solution. Bromide: G - E from the standard formulas (pipeline/thermo.atom), the same in both
phases. Started with no window through pipeline/launch_hidden.py; ORCA jobs go through the core queue.
"""
import json
import sys
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RX.parents[1]))
from pipeline import orca, thermo                            # noqa: E402
from pipeline.geometry import Geometry, formal_charge, read_xyz  # noqa: E402
from pipeline.steps import Reaction                          # noqa: E402

SPECIES = ["hydroxide", "bromomethane", "methanol", "bromide"]


def gas_keywords(method: str, job: str) -> str:
    return f"{method} {job}".strip()


def main():
    rx = Reaction.load(RX)
    smiles = {s["name"]: s["smiles"] for _, s in rx.species()}
    frag = json.loads((RX / "01_ts/fragments.json").read_text(encoding="utf-8"))
    out = {"rule": "gas phase = same keywords and settings as the SMD jobs, no solvent (J. Picione 2026-11-19)",
           "species": {}}
    for sp in SPECIES:
        q = formal_charge(smiles[sp])
        fd = RX / f"01_ts/fragments/{sp}"     # a one-atom single point writes no .xyz; its start.xyz is the geometry
        g0 = read_xyz(fd / f"{frag[sp]['job']}.xyz" if (fd / f"{frag[sp]['job']}.xyz").exists() else fd / "start.xyz")
        wd = rx.d(f"03_solvation/gas/{sp}")
        one = len(g0.symbols) == 1
        job = "" if one else "Opt Freq"
        inp = orca.make_input(gas_keywords(rx.cfg["method"], job), g0, q, 1, blocks="", nprocs=1,
                              maxcore_mb=rx.cfg["maxcore_mb"])
        name = sp
        r = orca.run(wd, name, inp, rx.log, note=f"solvation check: gas phase {smiles[sp]}, from the SMD geometry")
        if not r.ok or r.opt_converged is False:
            raise orca.OrcaError(f"gas {sp} failed; see {r.out}")
        if r.imaginary:            # same handling as the fragments step: push along the mode, reoptimize once
            h = orca.read_hess(wd / f"{name}.hess")
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            gp = read_xyz(wd / f"{name}.xyz")
            gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
            prev, name = name, f"{sp}_min"
            wd = rx.d(f"03_solvation/gas/{sp}_min")
            r = orca.run(wd, name, orca.make_input(gas_keywords(rx.cfg["method"], job), gp, q, 1, blocks="",
                                                   nprocs=1, maxcore_mb=rx.cfg["maxcore_mb"]), rx.log,
                         note=f"rerun: {prev} converged to a saddle point (imaginary "
                              f"{[round(f, 1) for f in r.imaginary]} cm-1); displaced 0.1 Å along that mode")
            if not r.ok or r.opt_converged is False:
                raise orca.OrcaError(f"gas {name} failed; see {r.out}")
        g = g0 if one else read_xyz(wd / f"{name}.xyz")
        if one:
            a = thermo.atom(g.symbols[0], 1)
            G, gsrc = r.energy_eh + a["G_minus_E_kJ"] / thermo.EH_KJ, a["method"]
        else:
            G, gsrc = r.gibbs_eh, "ORCA, quasi-RRHO, 298.15 K, 1 atm"
        res = {"job": str((wd / name).relative_to(RX)).replace("\\", "/"), "charge": q, "E_eh": r.energy_eh,
               "G_eh": G, "G_source": gsrc, "imaginary_cm": r.imaginary or [], "sp": {}}
        for lv in rx.cfg["diagnostic_levels"]:
            sd = rx.d(f"03_solvation/gas/{sp}/sp/{lv['label']}")
            s = orca.run(sd, "sp", orca.make_input(lv["method"], g, q, 1, blocks="", nprocs=rx.cfg["nprocs"],
                                                   maxcore_mb=rx.cfg["maxcore_mb"]),
                         rx.log, note=f"solvation check: gas-phase single point {lv['label']}")
            if not s.ok:
                raise orca.OrcaError(f"gas single point {sp} {lv['label']} failed; see {s.out}")
            res["sp"][lv["label"]] = s.energy_eh
        out["species"][sp] = res
        (RX / "03_solvation/gas.json").write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
