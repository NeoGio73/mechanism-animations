"""Secondary hydride shift: the "anti" conformer of the secondary cation at the working level (J. Picione, 2026-11-19).

Rule: the energy zero is the lowest-free-energy conformer of the secondary cation (same rule as for separated species,
docs/protocol.md "Reference conformers"), not automatically the form the path reaches (frame 0, `secondary_eclip_irc`).
The anti minimum was optimized only in the gas phase (M06-2X/ma-def2-SVP, `reactions/hshift_3_methyl_2_butyl/
01_check/method/M06-2X/anti`, C2–C3–CH₃ 71.3°). Here it is reoptimized at the working level (Opt Freq, SMD
dichloromethane, DefGrid3), checked to be still the secondary cation (H on C3), and given single points at every path
level with the same inputs as the diagnostics of the path ends. Nothing on the path is rerun.

Output: 01_stationary/secondary_anti/ (optimization), 01_stationary/secondary_anti/sp/<level>/ (single points),
01_stationary/secondary_anti.json (energies, geometry checks).
"""
import json
import sys
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RX.parents[1]))
from pipeline import orca                                   # noqa: E402
from pipeline.geometry import Geometry, read_xyz            # noqa: E402
from pipeline.steps import Reaction                         # noqa: E402

C2, C3, H = 1, 3, 9                                          # 0-based, as in reaction.toml [key_atoms]
SRC = RX.parent / "hshift_3_methyl_2_butyl/01_check/method/M06-2X/anti/anti.xyz"


def d(g, i, j):
    return float(np.linalg.norm(g.coords[i] - g.coords[j]))


def ang(g, i, j, k):
    a, b = g.coords[i] - g.coords[j], g.coords[k] - g.coords[j]
    return float(np.degrees(np.arccos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))


def tors(g, i, j, k, l):
    p = g.coords
    b0, b1, b2 = p[i] - p[j], p[k] - p[j], p[l] - p[k]
    b1 /= np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def describe(g):
    me = [i for i, s in enumerate(g.symbols) if s == "C" and i not in (0, C2, C3)]   # the two methyls on C3
    st = max(me, key=lambda i: d(g, C3, i))
    return {"C3_H_A": round(d(g, C3, H), 4), "C2_H_A": round(d(g, C2, H), 4),
            "C3_CH3_stretched_A": round(d(g, C3, st), 4), "C2_C3_CH3_deg": round(ang(g, C2, C3, st), 2),
            "C1_C2_C3_H_deg": round(tors(g, 0, C2, C3, H), 2)}


def main():
    rx = Reaction.load(RX)
    g0 = read_xyz(SRC)
    wd = rx.d("01_stationary/secondary_anti")
    if any(wd.iterdir()):
        sys.exit(f"{wd} is not empty; one job per clean folder")
    (RX / "00_inputs/secondary_anti_start.xyz").write_text(SRC.read_text())
    g0.write(wd / "start.xyz")
    name = "secondary_anti"
    r = orca.run(wd, name, rx.inp("Opt Freq", g0), rx.log,
                 note="anti conformer of the secondary cation, from the gas-phase M06-2X anti minimum")
    if not r.ok or r.opt_converged is False:
        raise orca.OrcaError(f"{name} failed; see {r.out}")
    if r.imaginary:
        h = orca.read_hess(wd / f"{name}.hess")
        v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
        gp = read_xyz(wd / f"{name}.xyz")
        gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
        prev, name = name, f"{name}_min"
        wd = rx.d("01_stationary/secondary_anti_min")
        r = orca.run(wd, name, rx.inp("Opt Freq", gp), rx.log,
                     note=f"rerun: {prev} converged to a saddle point (imaginary {[round(f, 1) for f in r.imaginary]}"
                          " cm-1); displaced 0.1 Å along that mode, same keywords")
        if not r.ok or r.opt_converged is False:
            raise orca.OrcaError(f"{name} failed; see {r.out}")
    g = read_xyz(wd / f"{name}.xyz")
    res = {"job": str((wd / name).relative_to(RX)).replace("\\", "/"), "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
           "imaginary_cm": r.imaginary, "geometry": describe(g),
           "still_secondary": d(g, C3, H) < d(g, C2, H), "sp": {}}
    for lv in rx.cfg["diagnostic_levels"]:
        sd = rx.d(f"01_stationary/secondary_anti/sp/{lv['label']}")
        inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', g, rx.cfg["charge"], rx.cfg["mult"],
                              blocks=rx.solvent_block, nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
        s = orca.run(sd, "sp", inp, rx.log, note=f"anti conformer, single point {lv['label']}")
        if not s.ok:
            raise orca.OrcaError(f"single point {lv['label']} failed; see {s.out}")
        res["sp"][lv["label"]] = s.energy_eh
    (RX / "01_stationary/secondary_anti.json").write_text(json.dumps(res, indent=2, ensure_ascii=False),
                                                          encoding="utf-8")
    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
