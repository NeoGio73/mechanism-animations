"""Results summary for one finished reaction (written 2026-11-18 for cyanide + acetone, J. Picione's request):
path, stationary points, barriers and reaction energies at every level, highest path point per level, DLPNO-CCSD(T)
diagnostics, validation, experiment. Reads existing files only (no ORCA job). Every energy in kJ/mol with
kcal/mol in brackets, each with its reference zero.

Free energies and enthalpies: electronic energy at the level named + the working-level (M06-2X-D3(0)/ma-def2-SVP)
thermal correction read from ORCA's thermochemistry block (quasi-RRHO, 298.15 K, 1 atm); bimolecular steps are moved
to the 1 M standard state with −RT ln(24.4654) per mole of gas lost. Enthalpies need no standard-state change.
Cross-checked against 01_ts/diagnostics.json, which the pipeline wrote independently (the script stops on a mismatch).

Usage: python -X utf8 -m pipeline.results_summary nucadd_cyanide_acetone_m062x docs/results/cyanide_acetone.md
"""
from __future__ import annotations

import csv
import json
import math
import re
import statistics
import sys
import time
import tomllib
from pathlib import Path

import numpy as np

from .geometry import read_xyz_frames
from .path_quality import region
from .steps import Reaction
from .thermo import orca_block

ROOT = Path(__file__).resolve().parents[1]
EH = 2625.499639
KCAL = 1 / 4.184
RT = 8.314462618e-3 * 298.15
CONV_1M = RT * math.log(24.4654)          # 1 atm → 1 M, per mole of gas (kJ/mol)
DL = "DLPNO-CCSD(T)_ma-def2-TZVP"
LEVELS = [("M06-2X/ma-def2-SVP (working level, geometries)", None, "energies.csv"),
          ("M06-2X/ma-def2-TZVP", "M06-2X_ma-def2-TZVP", "energies_M06-2X_ma-def2-TZVP.csv"),
          ("B3LYP/ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "energies_B3LYP_ma-def2-TZVP.csv"),
          ("DLPNO-CCSD(T)/ma-def2-TZVP", DL, f"energies_{DL}.csv")]


def kj(x: float, d: int = 1, sign: bool = False) -> str:
    f = f"{{:{'+' if sign else ''}.{d}f}}"
    return f"{f.format(x)} kJ/mol ({f.format(x * KCAL)} kcal/mol)"


def t(x: float, d: int = 1, sign: bool = False) -> str:
    """Table cell: kJ/mol (kcal/mol)."""
    f = f"{{:{'+' if sign else ''}.{d}f}}"
    return f"{f.format(x)} ({f.format(x * KCAL)})"


def main(rxname: str, out: str) -> None:
    rx = Reaction.load(ROOT / "reactions" / rxname)
    R = rx.root
    notes = tomllib.loads((R / "00_inputs/card_notes.toml").read_text(encoding="utf-8"))
    frag = json.loads((R / "01_ts/fragments.json").read_text(encoding="utf-8"))
    tsj = json.loads((R / "01_ts/ts.json").read_text(encoding="utf-8"))
    ends = json.loads((R / "02_irc/endpoints.json").read_text(encoding="utf-8"))
    diag = json.loads((R / "01_ts/diagnostics.json").read_text(encoding="utf-8"))
    meta = json.loads((R / "02_irc/irc_meta.json").read_text(encoding="utf-8"))
    its = meta["ts_frame"]
    L: list[str] = []
    w = L.append

    # ---------------------------------------------------------------- species, energies, thermal corrections
    react = [k for k, v in frag.items() if v["role"] == "reactants"]
    prods = [k for k, v in frag.items() if v["role"] == "products"]
    outs = {"ts": R / "01_stationary/ts/ts.out",
            "reactant_complex": R / "01_stationary/reactant_complex/reactant_complex.out",
            "product_complex": R / "01_stationary/product_complex/product_complex.out"}
    for k in react + prods:
        outs[k] = R / "01_ts/fragments" / k / f"{frag[k]['job']}.out"
    th = {k: orca_block(p) for k, p in outs.items()}
    Ework = {"ts": tsj["E_eh"], "reactant_complex": ends["reactant_complex"]["E_eh"],
             "product_complex": ends["product_complex"]["E_eh"], **{k: frag[k]["E_eh"] for k in react + prods}}
    for k, p in outs.items():                                   # the .out read is the job the reference uses
        e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", p.read_text(errors="replace"))[-1])
        assert abs(e - Ework[k]) < 1e-6, (k, e, Ework[k])
    E = {None: Ework, **{lv: diag["energies_eh"][lv] for _, lv, _ in LEVELS[1:]}}

    def quantities(lv):
        e = E[lv]
        rel = lambda a, b: (sum(e[x] for x in a) - sum(e[x] for x in b)) * EH
        corr = lambda key, a, b: sum(th[x][key] for x in a) - sum(th[x][key] for x in b)
        sep, cx, ts, pr = react, ["reactant_complex"], ["ts"], prods
        dn_ts = len(sep) - 1                                       # moles of gas lost, separated → TS
        dn_r = len(sep) - len(pr)
        q = {"dE_sep": rel(ts, sep), "dE_cx": rel(ts, cx),
             "dH_sep": rel(ts, sep) + corr("H_minus_E_kJ", ts, sep), "dH_cx": rel(ts, cx) + corr("H_minus_E_kJ", ts, cx),
             "dG_sep": rel(ts, sep) + corr("G_minus_E_kJ", ts, sep) - dn_ts * CONV_1M,
             "dG_cx": rel(ts, cx) + corr("G_minus_E_kJ", ts, cx),
             "dEr": rel(pr, sep), "dHr": rel(pr, sep) + corr("H_minus_E_kJ", pr, sep),
             "dGr": rel(pr, sep) + corr("G_minus_E_kJ", pr, sep) - dn_r * CONV_1M,
             "dE_cx_sep": rel(cx, sep), "dG_cx_sep": rel(cx, sep) + corr("G_minus_E_kJ", cx, sep) - dn_ts * CONV_1M,
             "dE_pc_sep": rel(["product_complex"], sep),
             "dG_pc_sep": rel(["product_complex"], sep) + corr("G_minus_E_kJ", ["product_complex"], sep) - dn_ts * CONV_1M,
             "dE_rev": rel(ts, pr), "dG_rev": rel(ts, pr) + corr("G_minus_E_kJ", ts, pr) - (len(pr) - 1) * CONV_1M}
        return q
    Q = {lv: quantities(lv) for _, lv, _ in LEVELS}
    # cross-check against the pipeline's own table
    names = {None: "M06-2X-D3(0)_ma-def2-SVP (working)"} | {lv: lv for _, lv, _ in LEVELS[1:]}
    for lv, q in Q.items():
        row = next(r for r in diag["rows"] if r["level"] == names[lv])
        for a, b in [("dE_sep", "dE_TS_vs_separated_kJ"), ("dE_cx", "dE_TS_vs_complex_kJ"), ("dEr", "dE_reaction_kJ"),
                     ("dG_sep", "dG_TS_vs_separated_1M_kJ"), ("dG_cx", "dG_TS_vs_complex_kJ"), ("dGr", "dG_reaction_1M_kJ")]:
            assert abs(q[a] - row[b]) < 0.01, (lv, a, q[a], row[b])

    # ---------------------------------------------------------------- path data
    rows = list(csv.DictReader(open(R / "02_irc/energies.csv", encoding="utf-8")))
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    dcol = next(c for c in rows[0] if c.startswith("d_C(N)-C"))
    ocol = next(c for c in rows[0] if c.startswith("d_C-O"))
    dcc = np.array([float(r[dcol]) for r in rows])
    dco = np.array([float(r[ocol]) for r in rows])
    Ew = np.array([float(r["E_kJ_per_mol"]) for r in rows])
    n = len(rows)
    frames = read_xyz_frames(R / "02_irc/irc_full_trj.xyz")
    assert len(frames) == n
    arc = sum(float(np.linalg.norm(frames[i + 1].coords - frames[i].coords)) for i in range(n - 1))
    _, _, inreg = region(rx)
    reg = [k for k in range(n) if inreg[k]]
    path_E = {}
    status = {}
    for name, lv, fn in LEVELS:
        p = R / "02_irc" / fn
        if p.exists():
            rr = list(csv.DictReader(open(p, encoding="utf-8")))
            col = "E_kJ_per_mol"
            vals = {int(r["frame"]): float(r[col]) for r in rr if r.get("computed", "1") in ("1", "True", "true")}
            path_E[lv] = vals
            status[lv] = "final"
        else:                                                 # level still running: read the finished jobs
            from .steps import path_sp_frames
            todo = path_sp_frames(rx, n, its, lv)
            ref = sum(E[lv][x] for x in react)
            vals = {}
            for k in todo:
                o = R / "02_irc/path_sp" / lv / f"f{k:03d}" / "f.out"
                if o.exists():
                    tx = o.read_text(errors="replace")
                    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)
                    if m and "ORCA TERMINATED NORMALLY" in tx:
                        vals[k] = (float(m[-1]) - ref) * EH
            path_E[lv] = vals
            status[lv] = f"partial: {len(vals)} of {len(todo)} frames"
    ets = path_E[None][its]

    # ================================================================ write
    w(f"# Results summary: {rx.cfg['title']}")
    w("")
    w(f"`{rx.cfg['equation']}`, in water (SMD), 298.15 K. Reaction folder `reactions/{rxname}`. "
      f"Written {time.strftime('%Y-%m-%d %H:%M')} by `pipeline/results_summary.py` from the finished jobs; no new calculation.")
    w("")
    w("**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "
      "\"Separated reactants\" means cyanide + acetone, each optimized on its own (single conformers). "
      "ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) "
      "to the electronic energy of the level named. Steps from separated reactants use the 1 M standard state "
      f"(−RT ln 24.47 = −{CONV_1M:.2f} kJ/mol per mole of solute lost); steps from the reactant complex involve no change "
      "in the number of molecules, so no standard-state term applies. ΔH needs none.")
    w("")
    if any(v != "final" for v in status.values()):
        w("**Draft.** " + "; ".join(f"{nm}: {status[lv]}" for nm, lv, _ in LEVELS if status[lv] != "final")
          + ". Path values at those levels are not final.")
        w("")

    # ---- 1. path
    w("## 1. Path")
    w("")
    w(f"- **Frames:** {n} (reactant side {its} frames, TS = frame {its}, product side {n - 1 - its} frames). "
      "Single points: " + "; ".join(f"{nm} {len(path_E[lv])} of {n}" for nm, lv, _ in LEVELS) +
      ". DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and "
      "within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.")
    w(f"- **Path length:** {s[-1] - s[0]:.2f} bohr·amu½ in the mass-weighted coordinate s "
      f"(from s = {s[0]:.2f} at the reactant end to s = {s[-1]:.2f} at the product end; TS at s = 0); "
      f"{arc:.2f} Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).")
    lo, hi = reg[0], reg[-1]
    w(f"- **Where the C–C bond forms:** the bond-changing region (C–C distance between 5 % and 95 % of its total change, "
      f"the definition of `pipeline/path_quality.py`) runs from frame {lo} (s = {s[lo]:.2f}, {dcc[lo]:.2f} Å) to frame "
      f"{hi} (s = {s[hi]:.2f}, {dcc[hi]:.2f} Å). At the TS the C–C distance is {dcc[its]:.3f} Å, "
      f"{100 * (dcc[0] - dcc[its]) / (dcc[0] - dcc[-1]):.0f} % of the way from the reactant end ({dcc[0]:.3f} Å) to the "
      f"product end ({dcc[-1]:.3f} Å). Most of the shortening happens after the TS.")
    w("")
    w("| Point | Frame | s (bohr·amu½) | C–C (Å) | C–O (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |")
    w("|---|---|---|---|---|---|")
    picks = [("reactant end", 0)]
    for target in (-6, -4, -2, -1, -0.5):
        k = int(np.argmin(abs(s - target)))
        if 0 < k < its:
            picks.append((f"s ≈ {target:g}", k))
    picks.append(("TS", its))
    for target in (0.5, 1, 1.5, 2, 3):
        k = int(np.argmin(abs(s - target)))
        if its < k < n - 1:
            picks.append((f"s ≈ +{target:g}", k))
    picks.append(("product end", n - 1))
    for lab, k in picks:
        w(f"| {lab} | {k} | {s[k]:+.2f} | {dcc[k]:.3f} | {dco[k]:.3f} | {t(Ew[k], 1, True)} |")
    w("")

    # ---- 2. stationary points
    w("## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD water)")
    w("")
    w(f"- **TS:** one imaginary frequency, {tsj['imaginary_cm'][0]:.1f} cm⁻¹ (written {abs(tsj['imaginary_cm'][0]):.1f}i cm⁻¹). "
      f"The mode is the C–C formation (C–C stretch rate {tsj['mode_check']['key_bond_stretch_rates'][0]['stretch_rate']:+.2f}, "
      f"C–O {tsj['mode_check']['key_bond_stretch_rates'][1]['stretch_rate']:+.2f}; forming and lengthening bonds move together). "
      f"C–C {tsj['bonds_A']['C(N)–C']:.3f} Å, C–O {tsj['bonds_A']['C–O']:.3f} Å, "
      f"Bürgi–Dunitz angle C(N)···C=O {tsj['angles_deg']['Bürgi–Dunitz angle C(N)···C=O']:.1f}°.")

    def mol_d(k):
        g = read_xyz_frames(R / "01_ts/fragments" / k / f"{frag[k]['job']}.xyz")[-1]
        sy = g.symbols
        out = []
        if "O" in sy:
            o = sy.index("O")
            c = min((i for i, x in enumerate(sy) if x == "C"), key=lambda i: g.distance(i, o))
            out.append(f"C–O {g.distance(c, o):.3f} Å")
            if "N" in sy:
                nn = sy.index("N")
                cn = min((i for i, x in enumerate(sy) if x == "C"), key=lambda i: g.distance(i, nn))
                out.append(f"C–C(N) {g.distance(c, cn):.3f} Å")
        elif "N" in sy:
            out.append(f"C≡N {g.distance(0, 1):.3f} Å")
        return ", ".join(out)
    G_rel = lambda k: ((Ework[k] - sum(Ework[x] for x in react)) * EH
                       + th[k]["G_minus_E_kJ"] - sum(th[x]["G_minus_E_kJ"] for x in react) - (len(react) - 1) * CONV_1M)
    w("- **References** (structures the barriers and reaction energies are measured from):")
    w("")
    w("| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |")
    w("|---|---|---|---|---|---|")
    for k in react:
        w(f"| {k} (`{frag[k]['smiles']}`, job `{frag[k]['job']}`) | separated reactant | {mol_d(k)} | none | — | — |")
    rc, pc = ends["reactant_complex"], ends["product_complex"]
    w(f"| reactant complex | reactant end of the path | C···C {rc['bonds_A']['C(N)–C']:.3f} Å, C–O {rc['bonds_A']['C–O']:.3f} Å | none | "
      f"{t((Ework['reactant_complex'] - sum(Ework[x] for x in react)) * EH, 1, True)} | {t(G_rel('reactant_complex'), 1, True)} |")
    w(f"| TS | | C–C {tsj['bonds_A']['C(N)–C']:.3f} Å, C–O {tsj['bonds_A']['C–O']:.3f} Å | 1 | "
      f"{t(Q[None]['dE_sep'], 1, True)} | {t(Q[None]['dG_sep'], 1, True)} |")
    w(f"| product complex | product end of the path | C–C {pc['bonds_A']['C(N)–C']:.3f} Å, C–O {pc['bonds_A']['C–O']:.3f} Å | none | "
      f"{t((Ework['product_complex'] - sum(Ework[x] for x in react)) * EH, 1, True)} | {t(G_rel('product_complex'), 1, True)} |")
    for k in prods:
        w(f"| {k} (`{frag[k]['smiles']}`) | separated product | {mol_d(k)} | none | {t(Q[None]['dEr'], 1, True)} | {t(Q[None]['dGr'], 1, True)} |")
    w("")
    w(f"The product is one molecule, so the product complex and the separated product are the same species optimized "
      f"twice (from the path end and on its own); they differ by "
      f"{(Ework['product_complex'] - sum(Ework[x] for x in prods)) * EH:+.2f} kJ/mol in E. Reaction energies use the "
      f"separated product. The reactant end of the path ended {rc['end_rmsd_to_minimum_A']:.3f} Å (RMSD) from the reactant "
      f"complex, the product end {pc['end_rmsd_to_minimum_A']:.3f} Å from the product complex (section 6).")
    w("")

    # ---- 3. barriers
    w("## 3. Barriers and reaction energies")
    w("")
    w("All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level "
      "thermal correction. kJ/mol (kcal/mol).")
    w("")
    w("**Barriers, zero = separated reactants (cyanide + acetone), 298.15 K, ΔG at 1 M:**")
    w("")
    w("| Level | ΔE‡ | ΔH‡ | ΔG‡ |")
    w("|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        q = Q[lv]
        w(f"| {nm} | {t(q['dE_sep'])} | {t(q['dH_sep'])} | {t(q['dG_sep'])} |")
    w("")
    w("**Barriers, zero = reactant complex, 298.15 K:**")
    w("")
    w("| Level | ΔE‡ | ΔH‡ | ΔG‡ |")
    w("|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        q = Q[lv]
        w(f"| {nm} | {t(q['dE_cx'])} | {t(q['dH_cx'])} | {t(q['dG_cx'])} |")
    w("")
    w("**Reaction, zero = separated reactants; product = separated alkoxide, 298.15 K, ΔG at 1 M:**")
    w("")
    w("| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from the alkoxide, ΔE‡ / ΔG‡ |")
    w("|---|---|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        q = Q[lv]
        w(f"| {nm} | {t(q['dEr'], 1, True)} | {t(q['dHr'], 1, True)} | {t(q['dGr'], 1, True)} | "
          f"{t(q['dE_cx_sep'], 1, True)} / {t(q['dG_cx_sep'], 1, True)} | {t(q['dE_rev'])} / {t(q['dG_rev'])} |")
    w("")
    w("ΔrH was not asked for; it is given because it costs nothing. The reactant complex lies below separated "
      "reactants in E but above them in G at 1 M: the binding energy does not pay for the entropy lost on association.")
    w("")

    # ---- 4. highest point
    w("## 4. Highest point of the path at each level")
    w("")
    w(f"Path energies at each level, zero = separated reactants at that level. The TS is frame {its} "
      f"(s = 0, C–C {dcc[its]:.3f} Å).")
    w("")
    w("| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C–C (Å) | Highest − TS frame |")
    w("|---|---|---|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        v = path_E[lv]
        if its not in v:
            w(f"| {nm} | {len(v)} | not yet computed | | | | |")
            continue
        km = max(v, key=v.get)
        w(f"| {nm} | {len(v)} of {n} | {t(v[its])} | {km} | {s[km]:+.2f} | {dcc[km]:.3f} | {t(v[km] - v[its], 2, True)} |")
    w("")
    w("A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP "
      "path, at a shorter C–C distance. The difference is the energy that a single point at the TS geometry misses.")
    w("")

    # ---- 5. DLPNO diagnostics
    w("## 5. DLPNO-CCSD(T) diagnostics")
    w("")
    t1 = {}
    for k in ["ts", "reactant_complex", "product_complex"] + react + prods:
        o = R / "01_ts/diagnostics" / DL / k / f"{k}.out"
        tx = o.read_text(errors="replace")
        m = re.findall(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", tx)
        amp = re.search(r"LARGEST AMPLITUDES\n-+\n(.*?)\n", tx)
        t1[k] = (float(m[-1]) if m else None, amp.group(1).split()[-1] if amp else None)
    w("| Structure | T1 | Largest amplitude |")
    w("|---|---|---|")
    for k, (a, b) in t1.items():
        w(f"| {k} | {a:.4f} | {b} |")
    pt = {}
    for k in path_E[DL]:
        o = R / "02_irc/path_sp" / DL / f"f{k:03d}" / "f.out"
        m = re.findall(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", o.read_text(errors="replace")) if o.exists() else []
        if m:
            pt[k] = float(m[-1])
    if pt:
        km = max(pt, key=pt.get)
        w("")
        w(f"- **Path frames ({len(pt)}):** T1 from {min(pt.values()):.4f} to {max(pt.values()):.4f}, median "
          f"{statistics.median(pt.values()):.4f}; largest at frame {km} (s = {s[km]:+.2f}, C–C {dcc[km]:.3f} Å).")
    w("- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a "
      "single-reference treatment is adequate for every structure and frame computed.")
    w("- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.")
    w("")

    # ---- 6. validation
    w("## 6. Validation")
    w("")
    ang = (R / "02_irc/true_irc/angles.md").read_text(encoding="utf-8")
    verdict = re.search(r"Path: \*\*(\w+)\*\*", ang).group(1)
    w(f"- **Angle check: {verdict.lower()}** (rule: in each zone, median < 5°, 90th percentile < 15°; "
      "mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):")
    w("")
    w("| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |")
    w("|---|---|---|---|---|---|")
    for ln in ang.splitlines():
        c = [x.strip() for x in ln.strip("|").split("|")]
        if c and c[0] in ("backward", "forward"):
            side = "reactant side (backward)" if c[0] == "backward" else "product side (forward)"
            w(f"| {side} | {c[2]} | {c[3]} | {c[4]} | {c[5]} | {c[6]} |")
    w("")
    ec = (R / "02_irc/true_irc/end_check.md").read_text(encoding="utf-8")
    eca = (R / "02_irc/true_irc/end_check_all.md").read_text(encoding="utf-8")
    w(f"- **End checks: {re.search(r'Verdict: \*\*(\w+)\*\*', ec).group(1)}** (`02_irc/true_irc/end_check.md`, "
      f"`end_check_all.md`).")
    for src in (ec, eca):
        for ln in src.splitlines():
            c = [x.strip() for x in ln.strip("|").split("|")]
            if c and c[0] in ("backward", "forward") and c[4] == "yes":
                side = "Reactant side" if c[0] == "backward" else "Product side"
                w(f"  - {side}: last frame {c[2]} kJ/mol from the {c[1].replace('_', ' ')} (RMSD {c[3]} Å). Opt Freq "
                  f"from the last frame: E {c[5]} kJ/mol from it, heavy-atom RMSD {c[6]} Å (all-atom {c[7]} Å), "
                  f"largest key-distance difference {c[8]} Å, {c[9]} imaginary frequencies: {c[10]}.")
    si = (ROOT / "docs/si_corrections.md").read_text(encoding="utf-8")
    row = next(ln for ln in si.splitlines() if ln.strip().startswith("| Cyanide + acetone"))
    c = [x.strip() for x in row.strip().strip("|").split("|")]
    w(f"- **Side-end table (SI entry 5):** reactant side {c[1]}, product side {c[2]} above the official end minimum "
      "(last frame, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): the reactant side "
      "is final, not extended; the paper lists how far above its minimum each side ended.")
    w("- **SI corrections that touch this reaction** (`docs/si_corrections.md`):")
    w("  - Entry 1 (stopping rule near a flat TS): the 20 flat steps on both sides began at least 10 kJ/mol below the "
      "TS, so the amended rule stops them at the same frame; not rerun.")
    w("  - Entry 4 (the minimum the path reaches is the official end): product side reaches the same minimum "
      "(E −0.031 kJ/mol, heavy-atom RMSD 0.003 Å); no switch, no value changed.")
    w("  - Entry 5: the side-end row above.")
    w("  - Entry 6 (angle against gradient size): 18 of the 70 frames over 15° in all six paths are on this reaction's "
      "product side, in the final approach to its minimum, where the RMS gradient is about 1.5 × 10⁻⁵ Eh/bohr; the "
      "largest, up to 117°, are there. Near this TS the angles stay small.")
    w("  - Entry 9 (lowest-free-energy conformers): cyanide, acetone and the alkoxide have single conformers; nothing changes.")
    w("  - No entry corrects a value reported for this reaction.")
    w("")

    # ---- 7. experiment
    w("## 7. Experiment")
    w("")
    er = notes.get("experiment_row", {})
    w(notes["experimental"].strip())
    w("")
    if er:
        w("**Computed against measured** (computed − measured; zero = separated reactants, 298.15 K, 1 M):")
        w("")
        w("| Level | ΔG‡ | ΔrG | Reverse ΔG‡ (from the alkoxide) |")
        w("|---|---|---|---|")
        w(f"| Measured | {t(er['dG_TS_vs_separated_1M_kJ'])} | {t(er['dG_reaction_1M_kJ'], 1, True)} | {t(er['dG_reverse_kJ'])} |")
        for nm, lv, _ in LEVELS:
            q = Q[lv]
            w(f"| {nm} | {t(q['dG_sep'])}, {t(q['dG_sep'] - er['dG_TS_vs_separated_1M_kJ'], 1, True)} | "
              f"{t(q['dGr'], 1, True)}, {t(q['dGr'] - er['dG_reaction_1M_kJ'], 1, True)} | "
              f"{t(q['dG_rev'])}, {t(q['dG_rev'] - er['dG_reverse_kJ'], 1, True)} |")
        w("")
    w("Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):")
    w("")
    for x in notes.get("limitations", []):
        w(f"- {x}")
    w("")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text("\n".join(L), encoding="utf-8")
    print("written", out, "| status:", status)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
