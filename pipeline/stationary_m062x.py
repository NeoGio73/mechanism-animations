"""Transition states and end points with frequencies at the M06-2X working level, VeryTightSCF, for all
six reactions (J. Picione, 2026-11-12: work that does not depend on the path method).

Level: M06-2X-D3(0)/ma-def2-SVP, DefGrid3, RIJCOSX (def2/J), VeryTightSCF, each reaction's SMD solvent.
TS: OptTS Freq with an exact Hessian at the start and every 5 steps; accepted only with exactly one
imaginary frequency (otherwise reported, not fixed). Minima: Opt Freq; a minimum with an imaginary mode
is displaced 0.1 Å along it and re-optimized once with TightOpt (as pipeline/steps.endpoints).
Starting geometries: the M06-2X structures of the level tests where they exist, otherwise the optimized
B3LYP end points. The end points here are optimized from those structures, not from M06-2X IRC ends;
when the M06-2X paths are run, each IRC end is checked against them.
Each job alone in its own folder, 3 cores (two at a time beside the 2-core E2 lane), through the core ledger.
Folders: reactions/<reaction>_m062x/01_stationary/<name>/; summary: docs/stationary_m062x.json
Usage: python -m pipeline.stationary_m062x [reaction ...]
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca
from .geometry import Geometry, read_xyz

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reactions"
LEVEL = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX VeryTightSCF DefGrid3"
TS_BLOCK = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end\n"

H3 = "hshift_2_3_dimethyl_2_butyl/01_ts/level_test/M06-2X"
MB = "hshift_3_methyl_2_butyl/01_check"
JOBS = {  # reaction: (new folder, charge, mult, solvent, [(name, kind, start xyz relative to reactions/)])
    "sn2": ("sn2_hydroxide_bromomethane_m062x", -1, 1, "water", [
        ("ts", "ts", "sn2_hydroxide_bromomethane_masvp/01_ts/level_test/M06-2X/ts/ts.xyz"),
        ("reactant_complex", "min", "sn2_hydroxide_bromomethane/02_irc/endpoints/reactant_complex/reactant_complex.xyz"),
        ("product_complex", "min", "sn2_hydroxide_bromomethane/02_irc/endpoints/product_complex/product_complex_tight.xyz")]),
    "e2": ("e2_ethoxide_2_bromopropane_m062x", -1, 1, "ethanol", [
        ("ts", "ts", "e2_ethoxide_2_bromopropane/01_ts/level_test/M06-2X/ts/ts.xyz"),
        ("reactant_complex", "min", "e2_ethoxide_2_bromopropane/02_irc/endpoints/reactant_complex/reactant_complex_tight.xyz"),
        ("product_complex", "min", "e2_ethoxide_2_bromopropane/02_irc/endpoints/product_complex/product_complex_cont.xyz")]),
    "da": ("da_butadiene_ethylene_m062x", 0, 1, "toluene", [
        ("ts", "ts", "da_butadiene_ethylene/01_ts/level_test/M06-2X/ts/ts.xyz"),
        ("reactant_complex", "min", "da_butadiene_ethylene/02_irc/endpoints/reactant_complex/reactant_complex_tight.xyz"),
        ("product_complex", "min", "da_butadiene_ethylene/02_irc/endpoints/product_complex/product_complex.xyz")]),
    "cyanide": ("nucadd_cyanide_acetone_m062x", -1, 1, "water", [
        ("ts", "ts", "nucadd_cyanide_acetone/01_ts/level_test/M06-2X/ts/ts.xyz"),
        ("reactant_complex", "min", "nucadd_cyanide_acetone/02_irc/endpoints/reactant_complex/reactant_complex_tight.xyz"),
        ("product_complex", "min", "nucadd_cyanide_acetone/02_irc/endpoints/product_complex/product_complex.xyz")]),
    "hshift_tertiary": ("hshift_2_3_dimethyl_2_butyl_m062x", 1, 1, "dichloromethane", [
        ("ts", "ts", f"{H3}/dcm/ts/ts.xyz"),
        ("open_cation", "min", f"{H3}/dcm/min/min.xyz"),
        ("bridged", "min", f"{H3}/bridged/dcm/ts_minus/minus.xyz")]),
    "hshift_methylbutyl": ("hshift_3_methyl_2_butyl_m062x", 1, 1, "dichloromethane", [
        ("ts", "ts", f"{MB}/rotation/M06-2X/eclip/ts.xyz"),               # gas-phase M06-2X rotation TS
        ("secondary_eclip", "min", f"{MB}/method/M06-2X/eclip/eclip_min.xyz"),  # gas-phase M06-2X minimum
        ("tertiary", "min", f"{MB}/dcm/perp/perp.xyz")]),                  # B3LYP DCM, already tertiary
}


def inp(job: str, g: Geometry, q: int, m: int, solvent: str, blocks: str = "") -> str:
    geo = "\n".join(f"{s:2s} {x:16.8f} {y:16.8f} {z:16.8f}" for s, (x, y, z) in zip(g.symbols, g.coords))
    return (f"! {LEVEL} {job} CPCM({solvent})\n%pal nprocs 3 end\n%maxcore 3000\n"
            f'%cpcm smd true SMDsolvent "{solvent.upper()}" end\n{blocks}\n* xyz {q} {m}\n{geo}\n*\n')


def one(key: str, name: str, kind: str, start: str) -> dict:
    folder, q, m, solv, _ = JOBS[key]
    base = R / folder / "01_stationary" / name
    base.mkdir(parents=True, exist_ok=True)
    log = orca.RunLog(R / folder / "run_log.md")
    g = read_xyz(R / start)
    g.write(base / "start.xyz", f"from reactions/{start}")
    note = f"M06-2X VeryTightSCF {'TS' if kind == 'ts' else 'minimum'}, from {start}"
    if kind == "ts":
        r = orca.run(base, "ts", inp("OptTS Freq", g, q, m, solv, TS_BLOCK), log, note=note)
        jn = "ts"
    else:
        r = orca.run(base, name, inp("Opt Freq", g, q, m, solv), log, note=note)
        jn = name
        if r.ok and r.opt_converged is False:
            prev, jn = jn, f"{name}_cont"
            r = orca.run(base, jn, inp("Opt Freq", read_xyz(base / f"{prev}.xyz"), q, m, solv, "%geom MaxIter 200 end\n"),
                         log, note=f"rerun: {prev} reached the cycle limit; continued from its last geometry, MaxIter 200")
        if r.ok and r.opt_converged and r.imaginary:
            h = orca.read_hess(base / f"{jn}.hess")
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            gp = read_xyz(base / f"{jn}.xyz")
            gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
            prev, jn = jn, f"{name}_tight"
            r = orca.run(base, jn, inp("TightOpt Freq", gp, q, m, solv), log,
                         note=f"rerun: {prev} had imaginary mode(s) {[round(f, 1) for f in r.imaginary]} cm-1; "
                              "displaced 0.1 Å along it, TightOpt")
    ok = bool(r.ok and r.opt_converged)
    n_im = len(r.imaginary) if ok else None
    accepted = ok and (n_im == 1 if kind == "ts" else n_im == 0)
    return {"reaction": key, "folder": f"{folder}/01_stationary/{name}", "name": name, "kind": kind,
            "job": jn, "ok": ok, "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
            "imaginary_cm": r.imaginary if ok else None, "accepted": accepted, "start": start}


def main(keys: list[str]) -> list[dict]:
    todo = [(k, n, kind, s) for k in keys for (n, kind, s) in JOBS[k][4]]
    todo.sort(key=lambda t: t[2] != "ts")                    # transition states first
    out_p = ROOT / "docs/stationary_m062x.json"
    res = json.loads(out_p.read_text(encoding="utf-8")) if out_p.exists() else {}

    def run(t):
        try:
            d = one(*t)
        except Exception as e:                                # one failure must not stop the others
            d = {"reaction": t[0], "name": t[1], "kind": t[2], "ok": False, "accepted": False, "error": str(e)}
        res[f"{t[0]}/{t[1]}"] = d
        out_p.write_text(json.dumps(res, indent=1), encoding="utf-8")
        print(json.dumps(d), flush=True)
        return d

    with ThreadPoolExecutor(2) as ex:
        return list(ex.map(run, todo))


if __name__ == "__main__":
    main(sys.argv[1:] or list(JOBS))
