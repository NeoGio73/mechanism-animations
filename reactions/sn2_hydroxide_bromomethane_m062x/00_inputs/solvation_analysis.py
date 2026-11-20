"""SN2 per-species solvation check, analysis (J. Picione, 2026-11-19). No ORCA job.

Inputs: 03_solvation/gas.json (gas-phase jobs, 00_inputs/solvation_check.py), 01_ts/fragments.json and
01_ts/diagnostics.json (the SMD water jobs the reaction energies use), 00_inputs/hydration_measured.toml (measured
hydration free energies with their sources).

Standard state: every solvation free energy is 1 M (gas) → 1 M (solution) (Ben-Naim). ORCA's G is computed for an
ideal gas at 1 atm in both phases, so G(SMD, solution geometry) − G(gas, gas geometry) is the 1 M → 1 M value: the
1 atm → 1 M term (RT ln 24.46 = 7.92 kJ/mol at 298.15 K) is the same in both phases and cancels.
For each level: G(phase) = E(level, phase) + [G − E](working-level frequencies, same phase and geometry).

Split of the reaction free energy error (computed − measured, solution):
  ΔrG(solution) = ΔrG(gas) + Σ ν·ΔG_solv  (ν = −1 reactants, +1 products), in both computation and experiment, so
  error = [ΔrG(gas, computed) − ΔrG(gas, measured)] + Σ ν·[ΔG_solv(computed) − ΔG_solv(measured)]
  with ΔrG(gas, measured) = ΔrG(solution, measured) − Σ ν·ΔG_solv(measured). Exact, no remainder.
Output: 03_solvation/solvation.json, 03_solvation/solvation.md.
"""
import json
import sys
import tomllib
from pathlib import Path

RX = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RX.parents[1]))
from pipeline.units import kj                                  # noqa: E402

NU = {"hydroxide": -1, "bromomethane": -1, "methanol": 1, "bromide": 1}
LEVELS = ["working", "M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]
LNAME = {"working": "M06-2X-D3(0)/ma-def2-SVP", "M06-2X_ma-def2-TZVP": "M06-2X/ma-def2-TZVP",
         "B3LYP_ma-def2-TZVP": "B3LYP/ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP": "DLPNO-CCSD(T)/ma-def2-TZVP"}
KCAL = 4.184


def main():
    gas = json.loads((RX / "03_solvation/gas.json").read_text(encoding="utf-8"))["species"]
    frag = json.loads((RX / "01_ts/fragments.json").read_text(encoding="utf-8"))
    diag = json.loads((RX / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["energies_eh"]
    meas = tomllib.loads((RX / "00_inputs/hydration_measured.toml").read_text(encoding="utf-8"))
    dGs_m = {k: meas["species"][k]["dG_solv_kJ"] for k in NU}
    dGr_soln_m = meas["reaction"]["dG_solution_kJ"]
    dGr_gas_m = dGr_soln_m - sum(NU[k] * dGs_m[k] for k in NU)
    res = {"standard_state": "1 M gas -> 1 M solution; 298.15 K", "measured": {"dG_solv_kJ": dGs_m,
           "dG_solution_kJ": dGr_soln_m, "dG_gas_kJ_derived": dGr_gas_m,
           "Br_minus_OH_kJ": dGs_m["bromide"] - dGs_m["hydroxide"]}, "levels": {}}
    import re                                          # the SMD radius ORCA used for Br (bromide SMD job)
    m = re.search(r"Radius for Br used is .*?=\s*([\d.]+) Ang", (RX / f"01_ts/fragments/bromide/{frag['bromide']['job']}.out")
                  .read_text(encoding="utf-8", errors="replace"))
    res["bromide_radius_A"] = float(m.group(1)) if m else None
    m2 = re.search(r"CPCM Dielectric\s*:\s*(-?[\d.]+) Eh", (RX / f"01_ts/fragments/bromide/{frag['bromide']['job']}.out")
                   .read_text(encoding="utf-8", errors="replace"))
    res["bromide_dielectric_kJ"] = kj(float(m2.group(1))) if m2 else None
    for lv in LEVELS:
        Gsol, Ggas, dGs = {}, {}, {}
        for k in NU:
            f, g = frag[k], gas[k]
            e_sol = f["E_eh"] if lv == "working" else diag[lv][k]
            e_gas = g["E_eh"] if lv == "working" else g["sp"][lv]
            Gsol[k] = e_sol + (f["G_eh"] - f["E_eh"])
            Ggas[k] = e_gas + (g["G_eh"] - g["E_eh"])
            dGs[k] = kj(Gsol[k] - Ggas[k])
        dGr_sol = kj(sum(NU[k] * Gsol[k] for k in NU))
        dGr_gas = kj(sum(NU[k] * Ggas[k] for k in NU))
        err = dGr_sol - dGr_soln_m
        e_gas = dGr_gas - dGr_gas_m
        e_sp = {k: NU[k] * (dGs[k] - dGs_m[k]) for k in NU}
        assert abs(err - e_gas - sum(e_sp.values())) < 1e-6
        res["levels"][lv] = {"dG_solv_kJ": dGs, "Br_minus_OH_kJ": dGs["bromide"] - dGs["hydroxide"],
                             "dG_solution_kJ": dGr_sol, "dG_gas_kJ": dGr_gas, "error_total_kJ": err,
                             "error_gas_kJ": e_gas, "error_solvation_kJ": e_sp,
                             "solvation_error_minus_meas_kJ": {k: dGs[k] - dGs_m[k] for k in NU}}
    (RX / "03_solvation/solvation.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    f = lambda x, n=1: f"{x:+.{n}f} ({x / KCAL:+.{n}f})"
    sp = ["hydroxide", "bromomethane", "methanol", "bromide"]
    L = ["kJ/mol (kcal/mol), 298.15 K, 1 M (gas) → 1 M (solution).", "",
         "**Solvation free energies in water** (SMD, solution geometry, minus gas, gas geometry):", "",
         "| Level | Hydroxide | Bromomethane | Methanol | Bromide | Bromide − hydroxide |", "|---|---|---|---|---|---|"]
    for lv in LEVELS:
        d = res["levels"][lv]
        L.append(f"| {LNAME[lv]} | " + " | ".join(f(d["dG_solv_kJ"][k]) for k in sp) + f" | {f(d['Br_minus_OH_kJ'])} |")
    L.append("| Measured | " + " | ".join(f(dGs_m[k]) for k in sp) + f" | {f(res['measured']['Br_minus_OH_kJ'])} |")
    L += ["", f"**Reaction free energy and its error split** (measured: solution {f(dGr_soln_m)}, gas "
          f"{f(dGr_gas_m)} = solution minus the measured hydration terms):", "",
          "| Level | ΔrG, solution | ΔrG, gas | Error, total | (a) gas-phase ΔrG | (b) hydroxide | (b) bromomethane | "
          "(b) methanol | (b) bromide |", "|---|---|---|---|---|---|---|---|---|"]
    for lv in LEVELS:
        d = res["levels"][lv]
        L.append(f"| {LNAME[lv]} | {f(d['dG_solution_kJ'])} | {f(d['dG_gas_kJ'])} | {f(d['error_total_kJ'])} | "
                 f"{f(d['error_gas_kJ'])} | " + " | ".join(f(d["error_solvation_kJ"][k]) for k in sp) + " |")
    L += ["", "(b) is each species' contribution to the reaction error: ν × [computed − measured solvation], ν = −1 for "
          "reactants and +1 for products; (a) + (b) = total."]
    # Findings, from the numbers above (DLPNO-CCSD(T) is the reference level of the summary).
    d = res["levels"]["DLPNO-CCSD(T)_ma-def2-TZVP"]
    es = d["error_solvation_kJ"]
    big = max(es, key=lambda k: abs(es[k]))
    ions = es["hydroxide"] + es["bromide"]
    neut = es["bromomethane"] + es["methanol"]
    br_lo, br_hi = (min(res["levels"][lv]["error_solvation_kJ"]["bromide"] for lv in LEVELS),
                    max(res["levels"][lv]["error_solvation_kJ"]["bromide"] for lv in LEVELS))
    rad = res.get("bromide_radius_A")
    L += ["", "**Findings** (DLPNO-CCSD(T)/ma-def2-TZVP, the reference level, unless stated):", "",
          f"- The solvation of {big} carries most of the error: {f(es[big])}. SMD solvates bromide by "
          f"{f(d['dG_solv_kJ']['bromide'])} against {f(dGs_m['bromide'])} measured; the same at every level "
          f"({br_lo:+.1f} to {br_hi:+.1f} kJ/mol), because the solvation of a single ion barely depends on the level.",
          f"- Hydroxide is also under-solvated, by {f(d['dG_solv_kJ']['hydroxide'] - dGs_m['hydroxide'])}; as a "
          f"reactant it pulls the error the other way ({f(es['hydroxide'])}). The two ions together give "
          f"{f(ions)}: SMD's bromide − hydroxide difference is {f(d['Br_minus_OH_kJ'])} against "
          f"{f(res['measured']['Br_minus_OH_kJ'])} measured, so SMD favours hydroxide over bromide too much.",
          "- That ion total does not depend on the absolute proton value behind the single-ion data (−265.9 "
          "kcal/mol, uncertain by at least 2 kcal/mol = 8.4 kJ/mol): one anion is consumed and one formed, so the "
          "proton value cancels. The split between the two ions does depend on it.",
          f"- The neutrals add {f(neut)} (bromomethane {f(es['bromomethane'])}, methanol {f(es['methanol'])}).",
          f"- The gas-phase reaction free energy is off by {f(d['error_gas_kJ'])}. Without the solvation errors, "
          f"DLPNO-CCSD(T) would give ΔrG in water within {abs(d['error_gas_kJ']):.1f} kJ/mol of experiment.",
          ]
    sv = {lv: sum(res["levels"][lv]["error_solvation_kJ"].values()) for lv in LEVELS}
    L += ["- The total solvation error is " + ", ".join(f"{LNAME[lv]} {f(sv[lv])}" for lv in LEVELS)
          + ". The bromide term hardly changes with the level; the hydroxide term does ("
          + ", ".join(f"{res['levels'][lv]['error_solvation_kJ']['hydroxide']:+.1f}" for lv in LEVELS)
          + " kJ/mol in the same order). The gas-phase part: "
          + ", ".join(f"{LNAME[lv]} {f(res['levels'][lv]['error_gas_kJ'])}" for lv in LEVELS)
          + ". So M06-2X/ma-def2-TZVP's small total error (" + f(res["levels"]["M06-2X_ma-def2-TZVP"]["error_total_kJ"])
          + ") is a cancellation of a gas-phase error against the solvation error, and B3LYP/ma-def2-TZVP has the "
          "gas-phase value closest to experiment but the largest total error."]
    if rad:
        eps, e2 = 78.355, 1389.35                     # water (ORCA's value); e²/(4πε₀) in kJ/mol·Å
        born = -(1 - 1 / eps) * e2 / (2 * rad)
        a_meas = (1 - 1 / eps) * e2 / (2 * -dGs_m["bromide"])
        L += [f"- SMD cavity radius ORCA used for Br: {rad:.2f} Å (ORCA 6.1.1, \"element-dependent radii\"). The "
              f"electrostatic part of the computed bromide solvation, {res['bromide_dielectric_kJ']:.1f} kJ/mol, is "
              f"close to the Born value for a point charge in a sphere of that radius ({born:.1f} kJ/mol). The measured "
              f"value would need a Born sphere of about {a_meas:.1f} Å. So the bromide error follows from the size "
              "of the cavity SMD gives Br; this radius was not checked against the SMD paper (Marenich, Cramer and Truhlar 2009)."]
    L += ["", "**Sources and standard state.** Every solvation free energy is 1 M (gas) → 1 M (solution), 298.15 K. "
          "ORCA's G is for 1 atm in both phases, so the 1 atm → 1 M term (7.92 kJ/mol) cancels in G(solution) − G(gas). "
          "Neutrals: FreeSolv v0.52 (Mobley and Guthrie, J. Comput.-Aided Mol. Des. 2014, 28, 711), entries "
          "mobley_1636752 (methanol, −5.10 kcal/mol) and mobley_8983100 (bromomethane, −0.82 kcal/mol), experimental "
          "reference Rizzo et al., J. Chem. Theory Comput. 2006, 2, 128; FreeSolv's ±0.60 kcal/mol for these two is a "
          "default (\"not presently available\"), not a measured uncertainty, and FreeSolv does not state its standard "
          "state (its methanol value equals the 1 M → 1 M value of Kelly et al., Table 3). Ions: C. P. Kelly, C. J. "
          "Cramer and D. G. Truhlar, J. Phys. Chem. B 2006, 110, 16066, Table 8 ('this work'; p. 58 of the authors' "
          "final manuscript, comp.chem.umn.edu/Truhlar/docs/final744.pdf; standard state in Section 2, pp. 5–6): "
          "OH⁻ −104.7 and Br⁻ −68.3 kcal/mol with ΔG*solv(H⁺) = −265.9 kcal/mol. Cross-check, same table, column "
          "'ref 23' (Pliego and Riveros 2002): OH⁻ −105.0, Br⁻ −68.6 kcal/mol, the same bromide − hydroxide difference.",
          "", "**Caveats.** The measured solution ΔrG (−100.8 kJ/mol) built CH₃Br(aq) from a Henry constant (+4.70 kJ/mol, "
          "1 bar → 1 mol/kg); FreeSolv's value used here for bromomethane corresponds to +4.53 kJ/mol for the same "
          "step. That moves the derived gas-phase value, and the split between (a) and the bromomethane term, by "
          "0.17 kJ/mol; the total error is unchanged. Gas-phase G uses the working-level frequencies of the gas "
          "geometry; solution G those of the SMD geometry. No structure had an imaginary frequency."]
    (RX / "03_solvation/solvation.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
