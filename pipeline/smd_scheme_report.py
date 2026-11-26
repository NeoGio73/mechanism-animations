"""SMD-scheme check, report part (pipeline/smd_scheme.py runs the jobs). J. Picione, 2026-11-23.

Scheme A (current): E = DLPNO-CCSD(T)/ma-def2-TZVP with SMD.
Scheme B: E = DLPNO-CCSD(T)/ma-def2-TZVP gas phase + [E(SMD) − E(gas)] at M06-2X/ma-def2-TZVP, same geometries.
Both get the same working-level thermal corrections, reference species and 1 M conversion as the paper: the energies
go through pipeline/results_summary.main's own quantities() (ΔG‡ = dG_sep, ΔrG = dGr; for one-molecule reactions
the zero is the reactant conformer the summary uses). Scheme A must reproduce the summary's DLPNO numbers exactly
(checked). If any |B − A| > 2 kJ/mol the report says STOP: no table changes before J. Picione has seen it.

Writes docs/results/smd_scheme_check.md and .json.
"""
from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path

from . import results_summary as rs
from .smd_scheme import LEVELS, REACTIONS, jobs

ROOT = Path(__file__).resolve().parents[1]
DL, M6 = LEVELS
EH = rs.EH
LIMIT = 2.0
NAMES = {"sn2_hydroxide_bromomethane_m062x": "SN2, HO⁻ + CH₃Br (water)",
         "e2_ethoxide_2_bromopropane_m062x": "E2, EtO⁻ + 2-bromopropane (ethanol)",
         "nucadd_cyanide_acetone_m062x": "CN⁻ + acetone (water)",
         "da_butadiene_ethylene_m062x": "Diels–Alder, butadiene + ethylene (SMD as in the paper)",
         "hshift_3_methyl_2_butyl_m062x": "1,2-H shift, 3-methyl-2-butyl cation (DCM)",
         "hshift_2_3_dimethyl_2_butyl_m062x": "1,2-H shift, 2,3-dimethyl-2-butyl cation (DCM)"}


def _e(out: Path) -> float:
    tx = out.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in tx:
        raise RuntimeError(f"{out} did not end normally: no result")
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)[-1])


def measured(name: str) -> dict:
    """Measured ΔG‡ and ΔrG (kJ/mol, 298.15 K, 1 M) as the results summary takes them (section 7): card_notes
    [experiment_row], then summary_notes [measured], then card_notes [measured_barrier]."""
    import tomllib
    d = ROOT / "reactions" / name / "00_inputs"
    rd = lambda f: tomllib.loads((d / f).read_text(encoding="utf-8")) if (d / f).exists() else {}
    notes, snotes = rd("card_notes.toml"), rd("summary_notes.toml")
    er = dict(notes.get("experiment_row", {}))
    er.update(snotes.get("measured", {}))
    if "measured_barrier" in notes and "dG_TS_vs_separated_1M_kJ" not in er:
        er["dG_TS_vs_separated_1M_kJ"] = notes["measured_barrier"]["dG_kJ"]
    return {"barrier": er.get("dG_TS_vs_separated_1M_kJ"), "reaction": er.get("dG_reaction_1M_kJ")}


def decision(rows: list) -> list:
    """J. Picione, 2026-11-26: scheme A kept; for SN2, E2 and cyanide, does B move each value toward the measured one?"""
    t = rs.t
    L = ["", "## Decision (J. Picione, 2026-11-26)", "",
         "Scheme A stays in every table, for consistency with the paths, cards and videos, which are all on scheme A. "
         "The B − A columns go into the SI as the uncertainty from the solvation scheme. Scheme B is not uniformly "
         "closer to experiment:", "",
         "| Reaction | Quantity | Measured | A − measured | B − measured | B moves it |", "|---|---|---|---|---|---|"]
    for r in rows:
        if not r["reaction"].startswith(("sn2", "e2", "nucadd")):
            continue
        m = measured(r["reaction"])
        r["measured"] = m
        for q, key, mk in (("ΔG‡", "dG_barrier", "barrier"), ("ΔrG", "dGr", "reaction")):
            if m[mk] is None:
                L.append(f"| {r['label']} | {q} | not measured | – | – | – |")
                r[f"{key}_toward_measured"] = None
                continue
            ea, eb = r[f"{key}_A"] - m[mk], r[f"{key}_B"] - m[mk]
            toward = abs(eb) < abs(ea)
            r[f"{key}_toward_measured"] = toward
            L.append(f"| {r['label']} | {q} | {t(m[mk])} | {t(ea, 1, True)} | {t(eb, 1, True)} | "
                     f"{'toward' if toward else 'away from'} the measured value, by {abs(abs(ea) - abs(eb)):.1f} kJ/mol |")
    return L


def per_species(rows: list) -> list:
    """Per species: B − A = [E(SMD) − E(gas)] at M06-2X/ma-def2-TZVP − [E(SMD) − E(gas)] at DLPNO-CCSD(T)/ma-def2-TZVP.
    A negative value puts the species lower in scheme B. ΔG‡(B) − ΔG‡(A) = TS value − sum of the reactants' values."""
    L = ["", "## Per species", "",
         "B − A for each structure: the M06-2X/ma-def2-TZVP solvation energy, E(SMD) − E(gas), minus the "
         "DLPNO-CCSD(T)/ma-def2-TZVP one, at the same geometry (from the jobs above; no new calculation). Negative: "
         "the species lies lower in scheme B. The barrier difference is the TS value minus the reactants' values; the "
         "reaction difference is the products' values minus the reactants' values. kJ/mol (kcal/mol).", "",
         "| Reaction | Species | Role | Solvation, M06-2X | Solvation, DLPNO-CCSD(T) | B − A |", "|---|---|---|---|---|---|"]
    t = rs.t
    for r in rows:
        sm, sd = r["solvation_M062X_kJ"], r["solvation_DLPNO_kJ"]
        role = {k: "reactant" for k in r["reference"]}
        role.update({k: "product" for k in r["products"]})
        role["ts"] = "transition state"
        order = list(r["reference"]) + ["ts"] + list(r["products"]) + [k for k in sm if k not in role]
        for k in order:
            L.append(f"| {r['label'].split(' (')[0]} | {k} | {role.get(k, 'not in ΔG‡ or ΔrG')} | {t(sm[k])} | "
                     f"{t(sd[k])} | {t(sm[k] - sd[k], 2, True)} |")
    return L


def report() -> dict:
    gas = {}
    for j in jobs():
        gas.setdefault(j["reaction"], {}).setdefault(j["level"], {})[j["species"]] = \
            _e(j["wd"] / f"{j['src'].stem}.out")
    rows, worst = [], 0.0
    for name in REACTIONS:
        with tempfile.TemporaryDirectory() as td:
            rs.main(name, str(Path(td) / "summary.md"))
        L = rs.LAST
        E, q = L["E"], L["quantities"]
        g = gas[name]
        EB = {}
        for k, e_smd_dl in E[DL].items():
            EB[k] = g[DL][k] + (E[M6][k] - g[M6][k])
        E["schemeB"] = EB
        A, B = L["Q"][DL], q("schemeB")
        a_chk = q(DL)
        assert all(abs(a_chk[x] - A[x]) < 1e-9 for x in A), name
        solv = {k: (E[M6][k] - g[M6][k]) * EH for k in EB}
        solv_dl = {k: (E[DL][k] - g[DL][k]) * EH for k in EB}
        # a reactant conformer with its own TS (pipeline/own_barrier.py): the TS shown measured from it is not the
        # barrier of any path, so the barrier here is the path shown, from the form it starts from (J. Picione,
        # 2026-11-26); that conformer's own rotation TS has no scheme-B (gas-phase) single point
        from .own_barrier import own_barriers
        from .steps import Reaction
        bk = "dG_cx" if own_barriers(Reaction.load(ROOT / "reactions" / name)) else "dG_sep"
        r = {"reaction": name, "label": NAMES[name] + (" (path shown, from the eclipsed form)" if bk == "dG_cx" else ""),
             "reference": L["react"], "products": L["prods"], "barrier_from": "path start" if bk == "dG_cx" else "zero",
             "dG_barrier_A": A[bk], "dG_barrier_B": B[bk], "dG_barrier_diff": B[bk] - A[bk],
             "dGr_A": A["dGr"], "dGr_B": B["dGr"], "dGr_diff": B["dGr"] - A["dGr"],
             "dE_barrier_A": A[bk.replace("dG", "dE")], "dE_barrier_B": B[bk.replace("dG", "dE")],
             "solvation_M062X_kJ": solv, "solvation_DLPNO_kJ": solv_dl}
        worst = max(worst, abs(r["dG_barrier_diff"]), abs(r["dGr_diff"]))
        rows.append(r)
    out = {"limit_kJ": LIMIT, "worst_abs_diff_kJ": worst, "stop": worst > LIMIT, "rows": rows}
    t = rs.t
    md = ["# SMD-scheme check (stationary points, all six reactions)", "",
          "Scheme A (as in the paper): DLPNO-CCSD(T)/ma-def2-TZVP with SMD. Scheme B: gas-phase DLPNO-CCSD(T)/ma-def2-TZVP "
          "+ [E(SMD) − E(gas)] at M06-2X/ma-def2-TZVP, at the same geometries. Same working-level thermal corrections, "
          "reference species and 1 M standard state as the paper's barriers (pipeline/results_summary). "
          "kJ/mol (kcal/mol), 298.15 K.", "",
          f"**Largest difference: {worst:.2f} kJ/mol ({worst / 4.184:.2f} kcal/mol). "
          + ("STOP: over the 2 kJ/mol limit; no table changed.**" if out["stop"] else
             "Every difference is within the 2 kJ/mol limit.**"), "",
          "| Reaction | ΔG‡ A | ΔG‡ B | B − A | ΔrG A | ΔrG B | B − A |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['label']} | {t(r['dG_barrier_A'])} | {t(r['dG_barrier_B'])} | {t(r['dG_barrier_diff'], 2, True)} "
                  f"| {t(r['dGr_A'])} | {t(r['dGr_B'])} | {t(r['dGr_diff'], 2, True)} |")
    md += ["", "Reference species (zero) per reaction: " + "; ".join(
        f"{r['label'].split(',')[0]}: {' + '.join(r['reference'])} → {' + '.join(r['products'])}" for r in rows)]
    md += decision(rows) + per_species(rows)
    md += ["", "Code: `pipeline/smd_scheme.py` (jobs), `pipeline/smd_scheme_report.py` (this table). Gas-phase jobs: "
           "`reactions/<reaction>/01_ts/smd_scheme/gas/`."]
    (ROOT / "docs/results/smd_scheme_check.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (ROOT / "docs/results/smd_scheme_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                             encoding="utf-8")
    print("\n".join(md))
    return out
