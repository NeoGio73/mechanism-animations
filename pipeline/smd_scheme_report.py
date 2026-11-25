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
        r = {"reaction": name, "label": NAMES[name], "reference": L["react"], "products": L["prods"],
             "dG_barrier_A": A["dG_sep"], "dG_barrier_B": B["dG_sep"], "dG_barrier_diff": B["dG_sep"] - A["dG_sep"],
             "dGr_A": A["dGr"], "dGr_B": B["dGr"], "dGr_diff": B["dGr"] - A["dGr"],
             "dE_barrier_A": A["dE_sep"], "dE_barrier_B": B["dE_sep"],
             "solvation_M062X_kJ": solv, "solvation_DLPNO_kJ": solv_dl}
        worst = max(worst, abs(r["dG_barrier_diff"]), abs(r["dGr_diff"]))
        rows.append(r)
    out = {"limit_kJ": LIMIT, "worst_abs_diff_kJ": worst, "stop": worst > LIMIT, "rows": rows}
    (ROOT / "docs/results/smd_scheme_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                             encoding="utf-8")
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
        f"{r['label'].split(',')[0]}: {' + '.join(r['reference'])} → {' + '.join(r['products'])}" for r in rows),
        "", "Code: `pipeline/smd_scheme.py` (jobs), `pipeline/smd_scheme_report.py` (this table). Gas-phase jobs: "
        "`reactions/<reaction>/01_ts/smd_scheme/gas/`."]
    (ROOT / "docs/results/smd_scheme_check.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return out
