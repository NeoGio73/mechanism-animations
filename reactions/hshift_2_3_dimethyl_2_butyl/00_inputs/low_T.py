"""Comparison with Saunders and Kates 1978 (J. Am. Chem. Soc. 100, 7082; read by J. Picione
2026-11-11): ΔG‡ = 13.0 ± 0.4 kJ/mol and k = 3.1 × 10⁷ s⁻¹ at −138 °C (135 K), superacid medium
(SO₂ClF/SO₂F₂), ¹³C NMR line broadening.

ΔG‡(135 K) from the lowest cation conformer (the optimized IRC end), DCM (SMD) and gas phase, no
new ORCA job: ΔE‡ at each level + G(T) − E from the working-level frequencies (pipeline/thermo.py,
quasi-RRHO, same model as ORCA). Levels other than the working level are added when their single
points exist (01_ts/low_T_levels.json: {"label": {"dcm": {"ts": E, "min": E}, "gas": {...}}}).

Symmetry: the cation minimum is C1 (σ = 1) and chiral; the TS is C2 within 0.0002 Å (ORCA
detected C2 at an early geometry, C1 at the final one). Its C2 axis swaps C2 and C3 and reverses
the reaction coordinate, so it is not counted in σ‡ (σ‡ = 1; equivalently σ‡ = 2 with a path
degeneracy of 2). Each mirror-image minimum crosses its own mirror-image TS, so the factor 2 for
the two enantiomeric conformers appears in reactant and TS and cancels.

Rate: Eyring k = κ (k_B T/h) exp(−ΔG‡/RT), first order (no standard-state term). κ = 1, and the
Wigner factor 1 + (hν‡/k_B T)²/24 as a rough indicator only (outside its reliable range here).
Run: .venv/Scripts/python.exe -X utf8 reactions/hshift_2_3_dimethyl_2_butyl/00_inputs/low_T.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))

from pipeline.thermo import H_PL, KB, R_KJ, check_against_orca, thermo   # noqa: E402
from pipeline.units import kj                                           # noqa: E402

T = 135.0
EXP = {"dG_kJ": 13.0, "dG_err_kJ": 0.4, "k_s": 3.1e7, "T_K": T,
       "k_other": {"145 K (−128 °C)": 4.6e7, "154 K (−119 °C)": 1.4e8}}
C_CM = 2.99792458e10
PH = {"dcm": {"ts": "01_ts/ts/ts", "min": "02_irc/endpoints/reactant_complex/reactant_complex",
              "min_mirror": "02_irc/endpoints/product_complex/product_complex",
              "aligned": "01_ts/fragments/dimethylbutyl_cation/dimethylbutyl_cation_min"},
      "gas": {"ts": "01_ts/gas_phase_opt/ts/ts", "min": "01_ts/gas_phase_opt/cation/cation"}}


def energy(stem):
    import re
    t = (ROOT / f"{stem}.out").read_text(errors="replace")
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])


def imag(stem):
    from pipeline import orca
    f = orca.read_hess(ROOT / f"{stem}.hess")["freqs"]
    return float(f[f < 0][0])


def fmt(x):
    return f"{x:.1f} kJ/mol ({x / 4.184:.1f} kcal/mol)"


def sci(x):
    e = int(np.floor(np.log10(abs(x))))
    return f"{x / 10 ** e:.2f} × 10{str(e).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"


def main():
    out = {"T_K": T, "experiment": EXP, "symmetry": "σ(minimum) = 1, σ‡ = 1 (C2 axis reverses the "
           "reaction coordinate); mirror-image factor 2 cancels", "phases": {}}
    pref = KB * T / H_PL
    levels = {"B3LYP-D3(BJ)/ma-def2-SVP (working)": None}
    extra = ROOT / "01_ts/low_T_levels.json"
    if extra.exists():
        levels.update(json.loads(extra.read_text(encoding="utf-8")))
    rows = []
    for ph, st in PH.items():
        g = {k: thermo(ROOT / f"{v}.hess", T, 1)["G_minus_E_kJ"] for k, v in st.items()}
        g298 = {k: thermo(ROOT / f"{v}.hess", 298.15, 1)["G_minus_E_kJ"] for k, v in st.items()}
        chk = {k: check_against_orca(ROOT / f"{v}.hess", ROOT / f"{v}.out") for k, v in st.items()}
        E = {k: energy(v) for k, v in st.items()}
        nu = abs(imag(st["ts"]))
        u = H_PL * C_CM * nu / (KB * T)
        wig = 1 + u ** 2 / 24
        d = {"G_minus_E_kJ_135K": g, "check_vs_orca_298K_kJ": chk, "imag_cm": -nu,
             "wigner_kappa": wig, "levels": {}}
        if "aligned" in st:
            dGc = kj(E["aligned"] - E["min"]) + g["aligned"] - g["min"]
            d["aligned_conformer_dG_135K_kJ"] = dGc
            d["aligned_conformer_population_135K"] = 1 / (1 + np.exp(dGc / (R_KJ * T)))
        for lab, lv in levels.items():
            e_ts, e_min = (E["ts"], E["min"]) if lv is None else (lv[ph]["ts"], lv[ph]["min"])
            dE = kj(e_ts - e_min)
            dG = dE + g["ts"] - g["min"]
            dG298 = dE + g298["ts"] - g298["min"]
            k = pref * np.exp(-dG / (R_KJ * T))
            d["levels"][lab] = {"dE_kJ": dE, "dG_135K_kJ": dG, "dG_298K_kJ": dG298, "k_135K_s": k,
                                "k_ratio_to_exp": k / EXP["k_s"], "k_wigner_s": k * wig,
                                "dG_minus_exp_kJ": dG - EXP["dG_kJ"]}
            rows.append((lab, ph, dE, dG, k, k / EXP["k_s"], k * wig))
        out["phases"][ph] = d
    (ROOT / "01_ts/low_T_135K.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")

    L = ["Comparison with Saunders and Kates 1978 at 135 K (−138 °C). Barriers from the lowest cation "
         "conformer; thermal corrections from the working-level frequencies (quasi-RRHO); σ‡ = 1, σ = 1 "
         "(see low_T.py); Eyring k = (k_B T/h) exp(−ΔG‡/RT), no tunnelling unless marked.", "",
         "| Level | Phase | ΔE‡ | ΔG‡ (135 K) | k (135 K), s⁻¹ | k / k(exp) | k × Wigner κ, s⁻¹ |",
         "|---|---|---|---|---|---|---|"]
    for lab, ph, dE, dG, k, r, kw in rows:
        L.append(f"| {lab} | {'SMD dichloromethane' if ph == 'dcm' else 'gas phase'} | {fmt(dE)} | {fmt(dG)} "
                 f"| {sci(k)} | {r:.2g} | {sci(kw)} |")
    L += ["", f"Experiment (superacid medium, 50:50 SO₂ClF/SO₂F₂): ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol), "
          f"k = 3.1 × 10⁷ s⁻¹ at 135 K.", ""]
    for ph, d in out["phases"].items():
        L.append(f"- {ph}: imaginary frequency {d['imag_cm']:.0f} cm⁻¹, Wigner κ(135 K) = {d['wigner_kappa']:.2f} "
                 "(rough indicator only; hν‡/k_BT ≈ "
                 f"{H_PL * C_CM * abs(d['imag_cm']) / (KB * T):.1f}, far outside the small-correction range "
                 "where the Wigner formula holds).")
        if "aligned_conformer_dG_135K_kJ" in d:
            L.append(f"- {ph}: the H-aligned conformer lies {fmt(d['aligned_conformer_dG_135K_kJ'])} above the "
                     f"lowest in G(135 K); population {d['aligned_conformer_population_135K']:.1e}.")
    (ROOT / "01_ts/low_T_135K.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    print(json.dumps({ph: d["check_vs_orca_298K_kJ"] for ph, d in out["phases"].items()}, indent=1))


if __name__ == "__main__":
    main()
