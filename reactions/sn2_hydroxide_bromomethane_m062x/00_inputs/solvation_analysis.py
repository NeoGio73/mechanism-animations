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
    (RX / "03_solvation/solvation.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
