"""Tertiary hydride shift: lowest-free-energy C2–C3 rotamer of the 2,3-dimethyl-2-butyl cation as the zero
(J. Picione, 2026-11-19; same rule as for the secondary shift). No new ORCA job: reads the rotamer jobs run by the
fragments step ([[references]] rotamer_h_perp, rotamer_h_plane in reaction.toml) and the existing structures.

Compared in G at 298.15 K (ORCA's own value, the rule) and at 135 K, the temperature of the Saunders–Kates
measurement (pipeline/thermo.py: the same quasi-RRHO model as ORCA, checked against ORCA at 298.15 K for every
structure). Candidates: the path-end minimum (`open_cation`, frame 0 optimized), the separately optimized cation
(`dimethylbutyl_cation`, the current zero) and the two rotamers. For each temperature: the lowest candidate, how far
frame 0 lies above it in E at every level, and ΔG‡ from it (TS frequencies from 01_stationary/ts).
Output: 01_ts/rotamer_zero.json, 01_ts/rotamer_zero.md.
Run: .venv/Scripts/python.exe -X utf8 reactions/hshift_2_3_dimethyl_2_butyl_m062x/00_inputs/rotamer_zero.py
"""
import csv
import io
import json
import re
import sys
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RX.parents[1]))
from pipeline.geometry import read_xyz                         # noqa: E402
from pipeline.steps import Reaction                            # noqa: E402
from pipeline.thermo import R_KJ, check_against_orca, symmetry_number, thermo   # noqa: E402
from pipeline.units import kj                                  # noqa: E402

TEMPS = (298.15, 135.0)
LEVELS = ["working", "M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]
LNAME = {"working": "M06-2X-D3(0)/ma-def2-SVP", "M06-2X_ma-def2-TZVP": "M06-2X/ma-def2-TZVP",
         "B3LYP_ma-def2-TZVP": "B3LYP/ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP": "DLPNO-CCSD(T)/ma-def2-TZVP"}
EXP = {"dG_kJ": 13.0, "err_kJ": 0.4, "T_K": 135.0, "source": "Saunders and Kates, J. Am. Chem. Soc. 1978, 100, 7082"}


def tors(p, i, j, k, l):
    b0, b1, b2 = p[i] - p[j], p[k] - p[j], p[l] - p[k]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def describe(stem: Path) -> dict:
    """Atoms (reaction.toml): 0 C2 (cation centre), 1 C3, 2 the H on C3, 3–4 methyls on C2, 5–6 methyls on C3.
    Angles of the C3–H and C3–CH₃ bonds to the empty p orbital (normal of the C2 plane): 0° = lined up with it."""
    g = read_xyz(stem.with_suffix(".xyz"))
    p = g.coords
    n = np.cross(p[3] - p[0], p[4] - p[0])
    n /= np.linalg.norm(n)
    to_p = lambda a: float(np.degrees(np.arccos(abs(np.dot((p[a] - p[1]) / np.linalg.norm(p[a] - p[1]), n)))))
    t = stem.with_suffix(".out").read_text(errors="replace")
    pg = re.findall(r"Point Group:\s+(\S+),", t)
    return {"Me_C2_C3_H_deg": round(tors(p, 3, 0, 1, 2), 1), "C3H_to_p_deg": round(to_p(2), 1),
            "C3Me_to_p_deg": sorted([round(to_p(5), 1), round(to_p(6), 1)]),
            "C3_H_A": round(float(np.linalg.norm(p[1] - p[2])), 4), "C2_H_A": round(float(np.linalg.norm(p[0] - p[2])), 4),
            "point_group": pg[-1] if pg else None}


def main():
    rx = Reaction.load(RX)
    frag = json.loads((RX / "01_ts/fragments.json").read_text(encoding="utf-8"))
    diag = json.loads((RX / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["energies_eh"]
    ep = json.loads((RX / "02_irc/endpoints.json").read_text(encoding="utf-8"))["reactant_complex"]
    ts = json.loads((RX / "01_ts/ts.json").read_text(encoding="utf-8"))
    cand = {"open_cation": {"stem": RX / "01_stationary/open_cation/open_cation", "E": ep["E_eh"],
                            "G": ep["G_eh"], "diag": "reactant_complex",
                            "what": "path-end minimum (frame 0 optimized)"},
            "dimethylbutyl_cation": {"what": "separately optimized cation (current zero)"},
            "rotamer_h_perp": {"what": "start: C3–H along the empty p orbital"},
            "rotamer_h_plane": {"what": "start: C3–H in the cation plane"}}
    for k, c in cand.items():
        if k in frag:
            f = frag[k]
            c.update(stem=RX / f"01_ts/fragments/{k}/{f['job']}", E=f["E_eh"], G=f["G_eh"], diag=k,
                     job=f["job"], imaginary=f["imaginary_cm"])
            if f["imaginary_cm"]:
                raise RuntimeError(f"{k} has imaginary frequencies {f['imaginary_cm']}")
    missing = [k for k, c in cand.items() if "stem" not in c]
    if missing and "--partial" not in sys.argv:
        sys.exit(f"not computed yet: {missing} (run the fragments step first)")
    cand = {k: c for k, c in cand.items() if "stem" in c}
    for k, c in cand.items():
        c["geometry"] = describe(c["stem"])
        sig = symmetry_number(c["stem"].with_suffix(".out"))
        c["sigma"] = sig
        c["check_vs_orca_298K_kJ"] = check_against_orca(c["stem"].with_suffix(".hess"), c["stem"].with_suffix(".out"))
        c["GmE_kJ"] = {T: thermo(c["stem"].with_suffix(".hess"), T, sig)["G_minus_E_kJ"] for T in TEMPS}
        c["E_lv"] = {"working": c["E"], **{lv: diag[lv][c["diag"]] for lv in LEVELS[1:]}}
    tstem = RX / "01_stationary/ts/ts"
    tsig = symmetry_number(tstem.with_suffix(".out"))
    ts_gme = {T: thermo(tstem.with_suffix(".hess"), T, 1)["G_minus_E_kJ"] for T in TEMPS}
    ts_E = {"working": ts["E_eh"], **{lv: diag[lv]["ts"] for lv in LEVELS[1:]}}

    # G relative to the lowest, working level: at 298.15 K ORCA's own G (the rule); at 135 K E + G(T) − E
    Grel = {}
    for T in TEMPS:
        g = {k: (kj(c["G"] - c["E"]) if T == 298.15 else c["GmE_kJ"][T]) + kj(c["E"]) for k, c in cand.items()}
        lo = min(g, key=g.get)
        Grel[T] = {"lowest": lo, "G_kJ_vs_lowest": {k: v - g[lo] for k, v in g.items()}}
    # frame 0 at every level
    files = {"working": RX / "02_irc/energies.csv", **{lv: RX / f"02_irc/energies_{lv}.csv" for lv in LEVELS[1:]}}
    f0 = {}
    for lv, p in files.items():
        if p.exists():
            f0[lv] = float(next(csv.DictReader(io.StringIO(p.read_text(encoding="utf-8"))))["E_eh"])
        else:
            o = RX / f"02_irc/path_sp/{lv}/f000/f.out"
            m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", o.read_text(errors="replace")) if o.exists() else []
            f0[lv] = float(m[-1]) if m else None
    res = {"rule": "zero = lowest-free-energy C2–C3 rotamer of the cation, working-level G (J. Picione 2026-11-19)",
           "temperatures_K": list(TEMPS), "candidates": {}, "by_temperature": {}, "experiment": EXP,
           "ts_sigma_orca": tsig, "ts_sigma_used": 1}
    for k, c in cand.items():
        res["candidates"][k] = {x: c.get(x) for x in ("what", "job", "geometry", "sigma", "check_vs_orca_298K_kJ")}
        res["candidates"][k]["E_eh"] = c["E_lv"]
        res["candidates"][k]["G_minus_E_kJ"] = {str(T): v for T, v in c["GmE_kJ"].items()}
    for T in TEMPS:
        lo = cand[Grel[T]["lowest"]]
        d = {"lowest": Grel[T]["lowest"], "G_kJ_vs_lowest": Grel[T]["G_kJ_vs_lowest"], "levels": {}}
        gme_lo = kj(lo["G"] - lo["E"]) if T == 298.15 else lo["GmE_kJ"][T]
        gme_ts = kj(ts["G_eh"] - ts["E_eh"]) if T == 298.15 else ts_gme[T]
        for lv in LEVELS:
            dE = kj(ts_E[lv] - lo["E_lv"][lv])
            d["levels"][lv] = {"frame0_dE_kJ": None if f0[lv] is None else kj(f0[lv] - lo["E_lv"][lv]),
                               "open_cation_dE_kJ": kj(cand["open_cation"]["E_lv"][lv] - lo["E_lv"][lv]),
                               "dE_ts_kJ": dE, "dG_ts_kJ": dE + gme_ts - gme_lo}
        res["by_temperature"][str(T)] = d
    (RX / "01_ts/rotamer_zero.json").write_text(json.dumps(res, indent=2, ensure_ascii=False, default=float),
                                                encoding="utf-8")
    f = lambda x, n=2: f"{x:+.{n}f} ({x / 4.184:+.{n}f})"
    L = ["# Tertiary hydride shift: lowest-free-energy rotamer of the cation (J. Picione, 2026-11-19)", "",
         "kJ/mol (kcal/mol). Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane, quasi-RRHO. 298.15 K: ORCA's "
         "G; 135 K: pipeline/thermo.py (same model, checked against ORCA at 298.15 K).", "",
         "| Structure | Me–C2–C3–H | C3–H to p axis | C3–CH₃ to p axis | point group | G vs lowest, 298.15 K | G vs lowest, 135 K | thermo check vs ORCA (G − E) |",
         "|---|---|---|---|---|---|---|---|"]
    for k, c in cand.items():
        gm = c["geometry"]
        L.append(f"| `{k}` ({c['what']}) | {gm['Me_C2_C3_H_deg']}° | {gm['C3H_to_p_deg']}° | "
                 f"{gm['C3Me_to_p_deg'][0]}°, {gm['C3Me_to_p_deg'][1]}° | {gm['point_group']} | "
                 f"{f(Grel[298.15]['G_kJ_vs_lowest'][k])} | {f(Grel[135.0]['G_kJ_vs_lowest'][k])} | "
                 f"{c['check_vs_orca_298K_kJ']['G_minus_E_kJ']:+.3f} kJ/mol |")
    for T in TEMPS:
        d = res["by_temperature"][str(T)]
        L += ["", f"**{T:g} K: lowest = `{d['lowest']}`.**", "",
              "| Level | Frame 0 above it (E) | Path-end minimum above it (E) | ΔE‡ | ΔG‡ |", "|---|---|---|---|---|"]
        for lv, v in d["levels"].items():
            L.append(f"| {LNAME[lv]} | {'pending' if v['frame0_dE_kJ'] is None else f(v['frame0_dE_kJ'])} | "
                     f"{f(v['open_cation_dE_kJ'])} | {f(v['dE_ts_kJ'], 1)} | {f(v['dG_ts_kJ'], 1)} |")
    L += ["", f"Measured (Saunders and Kates 1978, 135 K, superacid): ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol)."]
    (RX / "01_ts/rotamer_zero.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
