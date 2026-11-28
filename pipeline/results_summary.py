"""Results summary for one finished reaction (written 2026-11-18 for cyanide + acetone; made general for every reaction
2026-11-19, J. Picione: one summary per reaction when its single points finish): path, stationary points, barriers and
reaction energies at every level, highest path point per level, DLPNO-CCSD(T) diagnostics, validation (angle check,
end checks, SI corrections), experiment, job times. Reads existing files only (no ORCA job). Every energy in kJ/mol
with kcal/mol in brackets, each with its reference zero.

Free energies and enthalpies: electronic energy at the level named + the working-level (M06-2X-D3(0)/ma-def2-SVP)
thermal correction (ORCA's quasi-RRHO block; a single atom from pipeline/thermo.atom); steps that change the number of
molecules are moved to the 1 M standard state with RT ln(24.4654) per mole of solute. Enthalpies need no change.
Reference species follow pipeline.steps.reference_names (lowest free energy conformer, SI entry 9). A reaction with
no separated species (one molecule throughout, fragments.json empty) uses its reactant-end minimum as the zero.
Cross-checked against 01_ts/diagnostics.json, which the pipeline wrote independently (the script stops on a mismatch).

Hand-written parts per reaction: 00_inputs/card_notes.toml (experiment, limitations) and
00_inputs/summary_notes.toml (si_row = the row label in the side-end table of docs/si_corrections.md entry 5;
si_notes = the SI corrections that touch this reaction).

Usage: python -X utf8 -m pipeline.results_summary nucadd_cyanide_acetone_m062x docs/results/cyanide_acetone.md
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import statistics
import sys
import time
import tomllib
from pathlib import Path

import numpy as np

from .geometry import read_xyz_frames
from .path_quality import region
from .steps import Reaction, reference_names
from .thermo import atom, orca_block

ROOT = Path(__file__).resolve().parents[1]
EH = 2625.499639
KCAL = 1 / 4.184
RT = 8.314462618e-3 * 298.15
CONV_1M = RT * math.log(24.4654)          # 1 atm → 1 M, per mole of solute (kJ/mol)
DL = "DLPNO-CCSD(T)_ma-def2-TZVP"
LEVELS = [("M06-2X/ma-def2-SVP (working level, geometries)", None, "energies.csv"),
          ("M06-2X/ma-def2-TZVP", "M06-2X_ma-def2-TZVP", "energies_M06-2X_ma-def2-TZVP.csv"),
          ("B3LYP/ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "energies_B3LYP_ma-def2-TZVP.csv"),
          ("DLPNO-CCSD(T)/ma-def2-TZVP", DL, f"energies_{DL}.csv")]
RCOV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Cl": 1.02, "Br": 1.20, "I": 1.39}
SIDE = {"backward": "reactant side (backward)", "forward": "product side (forward)"}
LAST: dict = {}          # the last main() call's energies and quantities function (see main)


def kj(x: float, d: int = 1, sign: bool = False) -> str:
    f = f"{{:{'+' if sign else ''}.{d}f}}"
    return f"{f.format(x)} kJ/mol ({f.format(x * KCAL)} kcal/mol)"


def t(x: float, d: int = 1, sign: bool = False) -> str:
    """Table cell: kJ/mol (kcal/mol)."""
    x = 0.0 if abs(x) < 1e-6 else x                    # a degenerate reaction: no "-0.0"
    f = f"{{:{'+' if sign else ''}.{d}f}}"
    return f"{f.format(x)} ({f.format(x * KCAL)})"


def nice(k: str) -> str:
    return {"ts": "TS", "reactant_complex": "reactant complex", "product_complex": "product complex"}.get(k, k)


def run_seconds(out: Path) -> float | None:
    m = re.findall(r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec",
                   out.read_text(errors="replace"))
    if not m:
        return None
    d, h, mi, s, ms = map(int, m[-1])
    return d * 86400 + h * 3600 + mi * 60 + s + ms / 1000


def bond_types(g) -> str:
    """Distinct bonded atom-pair types of a molecule with their lengths (heavy atoms; all atoms if none)."""
    sy, n = g.symbols, len(g.symbols)
    if n == 1:
        return "single atom"
    pairs: dict[str, list[float]] = {}
    for heavy in (True, False):
        for i in range(n):
            for j in range(i + 1, n):
                if heavy and "H" in (sy[i], sy[j]):
                    continue
                d = g.distance(i, j)
                if d < 1.2 * (RCOV.get(sy[i], 0.8) + RCOV.get(sy[j], 0.8)):
                    a, b = sorted((sy[i], sy[j]), key=lambda x: (x == "H", x != "C", x))
                    pairs.setdefault(f"{a}–{b}", []).append(d)
        if pairs:
            break
    out = []
    for k, v in sorted(pairs.items()):
        out.append(f"{k} {min(v):.3f} Å" if max(v) - min(v) < 0.0015 else f"{k} {min(v):.3f}–{max(v):.3f} Å")
    return ", ".join(out)


def geom(p: Path):
    """Last geometry of a job: its .xyz, else the start.xyz beside it (a single atom is not optimized)."""
    x = p.with_suffix(".xyz")
    return read_xyz_frames(x if x.exists() else p.parent / "start.xyz")[-1]


def main(rxname: str, out: str) -> None:
    rx = Reaction.load(ROOT / "reactions" / rxname)
    R = rx.root
    notes = tomllib.loads((R / "00_inputs/card_notes.toml").read_text(encoding="utf-8"))
    sn_p = R / "00_inputs/summary_notes.toml"
    snotes = tomllib.loads(sn_p.read_text(encoding="utf-8")) if sn_p.exists() else {}
    frag = json.loads((R / "01_ts/fragments.json").read_text(encoding="utf-8"))
    tsj = json.loads((R / "01_ts/ts.json").read_text(encoding="utf-8"))
    ends = json.loads((R / "02_irc/endpoints.json").read_text(encoding="utf-8"))
    diag = json.loads((R / "01_ts/diagnostics.json").read_text(encoding="utf-8"))
    meta = json.loads((R / "02_irc/irc_meta.json").read_text(encoding="utf-8"))
    its = meta["ts_frame"]
    solvent = rx.cfg["solvent"]
    L: list[str] = []
    w = L.append

    # ---------------------------------------------------------------- species, energies, thermal corrections
    uni = not frag
    zc = {}                                    # one molecule: other conformers of the reactant (pipeline/unimolecular_zero.py)
    if uni:                                    # one molecule throughout: zero = lowest-G conformer of the reactant
        from .unimolecular_zero import zero as uzero
        uz = uzero(rx)
        zc = {k: v for k, v in uz["candidates"].items() if k != "reactant_complex"}
        react, prods = [uz["name"]], ["product_complex"]
    else:
        react, prods = reference_names(rx, frag)
    refs_other = [k for k in frag if k not in react + prods]      # conformers kept on record, not the zero
    outs = {"ts": Path(f"{rx.ts_stem()}.out")}
    for k in ("reactant_complex", "product_complex"):
        outs[k] = (R / ends[k]["xyz"]).with_suffix(".out")
    for k in frag:
        outs[k] = R / "01_ts/fragments" / k / f"{frag[k]['job']}.out"
    for k, v in zc.items():
        outs[k] = R / f"{v['job']}.out"
    Ework = {"ts": tsj["E_eh"], **{k: ends[k]["E_eh"] for k in ("reactant_complex", "product_complex")},
             **{k: v["E_eh"] for k, v in frag.items()}, **{k: v["E_eh"] for k, v in zc.items()}}
    Gwork = {"ts": tsj["G_eh"], **{k: ends[k]["G_eh"] for k in ("reactant_complex", "product_complex")},
             **{k: v["G_eh"] for k, v in frag.items()}, **{k: v["G_eh"] for k, v in zc.items()}}
    th = {}
    for k, p in outs.items():
        tx = p.read_text(errors="replace")
        e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)[-1])
        assert abs(e - Ework[k]) < 1e-6, (k, e, Ework[k])          # the .out read is the job the reference uses
        gk = (Gwork[k] - Ework[k]) * EH
        if "Final Gibbs free energy" in tx:
            b = orca_block(p)
            assert abs(b["G_minus_E_kJ"] - gk) < 0.01, (k, b["G_minus_E_kJ"], gk)
            th[k] = {"H_minus_E_kJ": b["H_minus_E_kJ"], "G_minus_E_kJ": gk}
        else:                                                      # a single atom (bromide): standard formulas
            sym = geom(p).symbols
            assert len(sym) == 1, (k, "no thermochemistry block")
            a = atom(sym[0], 1)
            assert abs(a["G_minus_E_kJ"] - gk) < 0.01, (k, a["G_minus_E_kJ"], gk)
            th[k] = {"H_minus_E_kJ": a["H_minus_E_kJ"], "G_minus_E_kJ": gk}
    E = {None: Ework, **{lv: {**diag["energies_eh"][lv], **{k: v["sp"][lv] for k, v in zc.items()}}
                         for _, lv, _ in LEVELS[1:]}}
    zoff = uni and react[0] != "reactant_complex"      # the zero is not the form the path reaches
    from .own_barrier import own_barriers               # that form's own TS, when it has one (J. Picione, 2026-11-26)
    own = own_barriers(rx) if zoff else None

    def quantities(lv):
        e = E[lv]
        rel = lambda a, b: (sum(e[x] for x in a) - sum(e[x] for x in b)) * EH
        corr = lambda key, a, b: sum(th[x][key] for x in a) - sum(th[x][key] for x in b)
        sep, cx, ts, pr, pc = react, ["reactant_complex"], ["ts"], prods, ["product_complex"]
        dn_ts = len(sep) - 1                                       # moles of solute lost, separated → TS
        dn_r = len(sep) - len(pr)
        rev = pr if len(pr) == 1 else pc                           # reverse barrier from one molecule
        q = {"dE_sep": rel(ts, sep), "dE_cx": rel(ts, cx),
             "dH_sep": rel(ts, sep) + corr("H_minus_E_kJ", ts, sep), "dH_cx": rel(ts, cx) + corr("H_minus_E_kJ", ts, cx),
             "dG_sep": rel(ts, sep) + corr("G_minus_E_kJ", ts, sep) - dn_ts * CONV_1M,
             "dG_cx": rel(ts, cx) + corr("G_minus_E_kJ", ts, cx),
             "dEr": rel(pr, sep), "dHr": rel(pr, sep) + corr("H_minus_E_kJ", pr, sep),
             "dGr": rel(pr, sep) + corr("G_minus_E_kJ", pr, sep) - dn_r * CONV_1M,
             "dE_cx_sep": rel(cx, sep), "dG_cx_sep": rel(cx, sep) + corr("G_minus_E_kJ", cx, sep) - dn_ts * CONV_1M,
             "dE_pc_sep": rel(pc, sep),
             "dG_pc_sep": rel(pc, sep) + corr("G_minus_E_kJ", pc, sep) - dn_ts * CONV_1M,
             "dE_rev": rel(ts, rev), "dH_rev": rel(ts, rev) + corr("H_minus_E_kJ", ts, rev),
             "dG_rev": rel(ts, rev) + corr("G_minus_E_kJ", ts, rev)}
        return q
    Q = {lv: quantities(lv) for _, lv, _ in LEVELS}
    # for other checks on the same reference states and thermal corrections (pipeline/smd_scheme_report.py): add an
    # energy set to E under a new key, then call quantities(key)
    LAST.clear()
    LAST.update({"E": E, "quantities": quantities, "Q": Q, "react": react, "prods": prods})
    # cross-check against the pipeline's own table (working level = first row)
    for (nm, lv, _), first in zip(LEVELS, [True, False, False, False]):
        row = diag["rows"][0] if first else next(r for r in diag["rows"] if r["level"] == lv)
        pairs = [("dE_cx", "dE_TS_vs_complex_kJ"), ("dG_cx", "dG_TS_vs_complex_kJ")]
        if not uni:                                # (one molecule: only the values from the path end are in the table)
            pairs += [("dE_sep", "dE_TS_vs_separated_kJ"), ("dE_sep", "dE_TS_vs_separated_kJ"), ("dEr", "dE_reaction_kJ"),
                      ("dG_sep", "dG_TS_vs_separated_1M_kJ"), ("dGr", "dG_reaction_1M_kJ")]
        for a, b in pairs:
            if row.get(b) is not None:
                assert abs(Q[lv][a] - row[b]) < 0.01, (lv, a, Q[lv][a], row[b])

    # ---------------------------------------------------------------- path data
    rows = list(csv.DictReader(open(R / "02_irc/energies.csv", encoding="utf-8")))
    n = len(rows)
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    kb = [(b["bond"], b["role"], tuple(b["atoms"]), b["stretch_rate"]) for b in tsj["mode_check"]["key_bond_stretch_rates"]]
    col = {}
    for lab, _, (a, b), _ in kb:
        col[lab] = next(c for c in rows[0] if c.startswith("d_") and
                        (c.endswith(f"_{a}_{b}_A") or c.endswith(f"_{b}_{a}_A")))
    D = {lab: np.array([float(r[col[lab]]) for r in rows]) for lab in col}
    prim, prim_role = next(((lab, ro) for lab, ro, _, _ in kb if ro == "forms"), (kb[0][0], kb[0][1]))
    dprim = D[prim]
    frames = read_xyz_frames(R / "02_irc/irc_full_trj.xyz")
    assert len(frames) == n
    arc = sum(float(np.linalg.norm(frames[i + 1].coords - frames[i].coords)) for i in range(n - 1))
    _, _, inreg = region(rx)
    reg = [k for k in range(n) if inreg[k]]
    zero = {lv: sum(E[lv][x] for x in react) for _, lv, _ in LEVELS}
    path_E, status, copied, jobt = {}, {}, {}, {}
    from .steps import path_sp_frames
    for name, lv, fn in LEVELS:
        p = R / "02_irc" / fn
        if p.exists():
            rr = list(csv.DictReader(open(p, encoding="utf-8")))
            path_E[lv] = {int(r["frame"]): (float(r["E_eh"]) - zero[lv]) * EH for r in rr if r.get("E_eh")}
            copied[lv] = sum(1 for r in rr if r.get("copied_from_frame"))
            status[lv] = "final"
        else:                                                 # level still running: read the finished jobs
            vals = {}
            for k in path_sp_frames(rx, n, its, lv):
                o = R / "02_irc/path_sp" / lv / f"f{k:03d}" / "f.out"
                if o.exists():
                    tx = o.read_text(errors="replace")
                    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)
                    if m and "ORCA TERMINATED NORMALLY" in tx:
                        vals[k] = (float(m[-1]) - zero[lv]) * EH
            path_E[lv] = vals
            copied[lv] = 0
            status[lv] = f"partial: {len(vals)} of {len(path_sp_frames(rx, n, its, lv))} frames"
        if lv is not None:                                     # job times of the path single points
            ts_ = []
            for o in sorted((R / "02_irc/path_sp" / lv).glob("f*/f.out")):
                sec = run_seconds(o)
                if sec is not None:
                    ts_.append((int(o.parent.name[1:]), sec, os.path.getmtime(o)))
            jobt[lv] = ts_
    Ew = (np.array([float(r["E_eh"]) for r in rows]) - zero[None]) * EH
    ets = path_E[None][its]
    rnm = snotes.get("reactant_name", "reactant-end minimum")      # e.g. "secondary cation" for a one-molecule reaction
    pnm = snotes.get("product_name", "product-end minimum")
    zero_name = (f"the {rnm} (reactant-end minimum, `{ends['reactant_complex']['minimum']}`)" if uni and not zoff
                 else f"the {rnm}, lowest-free-energy conformer (`{react[0]}`)" if uni
                 else "separated reactants (" + " + ".join(react) + ")")
    path_form = f"the {rnm} as the path reaches it (`{ends['reactant_complex']['minimum']}`, frame 0)"

    # ================================================================ write
    w(f"# Results summary: {rx.cfg['title']}")
    w("")
    w(f"`{rx.cfg['equation']}`, in {solvent} (SMD), 298.15 K. Reaction folder `reactions/{rxname}`. "
      f"Written {time.strftime('%Y-%m-%d %H:%M')} by `pipeline/results_summary.py` from the finished jobs; no new calculation.")
    w("")
    if uni and zoff:
        zl = {None: "M06-2X/ma-def2-SVP", **{lv: nm.split(" (")[0] for nm, lv, _ in LEVELS[1:]}}
        offs = "; ".join(f"{zl[lv]} {t(Q[lv]['dE_cx_sep'], 2, True)} kJ/mol" for _, lv, _ in LEVELS)
        zt = (f"The reaction has one molecule throughout, so there are no separated reactants, and every step keeps the "
              f"number of molecules, so no standard-state term applies. The zero is the lowest-free-energy conformer of "
              f"the {rnm} (`{react[0]}`, working-level G; rule of J. Picione, 2026-11-19, the same as for separated "
              f"species), not the form the path reaches. The path starts from {path_form}, which lies above the zero by "
              f"{t(Q[None]['dG_cx_sep'], 2, True)} kJ/mol in G (working level) and, in E, by {offs}. Barriers from {path_form} "
              "are given as well." + (
              f" The `{own['name']}` form does not reach the TS of the path shown; it has its own TS (`{own['ts_key']}`), "
              f"so its barrier is given through that TS, and the TS shown measured from `{own['name']}` is an energy "
              "difference, not the barrier of any path (J. Picione, 2026-11-26; section 3)." if own else ""))
    elif uni:
        zt = (f"The zero is the {rnm} at the reactant end of the path (`{ends['reactant_complex']['minimum']}`): the reaction has one molecule "
              "throughout, so there are no separated reactants, and every step keeps the number of molecules, so no "
              "standard-state term applies.")
    else:
        conf = "single conformers" if not refs_other else "each the conformer lowest in free energy (SI entry 9)"
        zt = ((f"\"Separated reactants\" means {' + '.join(react)}, each optimized on its own ({conf}). " if len(react) > 1
               else f"\"Separated reactant\" means {react[0]}, optimized on its own ({conf.replace('single conformers', 'single conformer').replace('each the', 'the')}). ") +
              "Steps that change the number of molecules use the 1 M standard state "
              f"(RT ln 24.47 = {CONV_1M:.2f} kJ/mol per mole of solute: subtracted for each molecule lost, added for each "
              "gained); steps from the reactant complex involve no change in the number of molecules, so no "
              "standard-state term applies.")
    w("**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. " + zt +
      " ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) "
      "to the electronic energy of the level named. ΔH needs no standard-state term.")
    w("")
    if any(v != "final" for v in status.values()):
        w("**Draft.** " + "; ".join(f"{nm}: {status[lv]}" for nm, lv, _ in LEVELS if status[lv] != "final")
          + ". Path values at those levels are not final.")
        w("")

    # ---- 1. path
    w("## 1. Path")
    w("")
    sparse = DL in rx.cfg.get("sparse_path_levels", [])
    cp = "; ".join(f"{nm} {len(path_E[lv])} of {n}" + (f" ({len(path_E[lv]) - copied[lv]} computed, {copied[lv]} "
                   "copied from their mirror partner, SI entry 7)" if copied[lv] else "") for nm, lv, _ in LEVELS)
    w(f"- **Frames:** {n} (reactant side {its} frames, TS = frame {its}, product side {n - 1 - its} frames). "
      f"Single points: {cp}." + (" DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the "
      "bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends." if sparse else ""))
    w(f"- **Path length:** {s[-1] - s[0]:.2f} bohr·amu½ in the mass-weighted coordinate s "
      f"(from s = {s[0]:.2f} at the reactant end to s = {s[-1]:.2f} at the product end; TS at s = 0); "
      f"{arc:.2f} Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).")
    lo, hi = reg[0], reg[-1]
    fb = [lab for lab, ro, _, _ in kb if ro in ("forms", "breaks")]
    frac = 100 * (dprim[0] - dprim[its]) / (dprim[0] - dprim[-1])
    verb = "shortening" if dprim[-1] < dprim[0] else "lengthening"
    w(f"- **Where the bonds change:** the bond-changing region ({' or '.join(fb)} between 5 % and 95 % of its total "
      f"change, the definition of `pipeline/path_quality.py`) runs from frame {lo} (s = {s[lo]:.2f}, {prim} "
      f"{dprim[lo]:.2f} Å) to frame {hi} (s = {s[hi]:.2f}, {dprim[hi]:.2f} Å). At the TS the {prim} distance is "
      f"{dprim[its]:.3f} Å, {frac:.0f} % of the way from the reactant end ({dprim[0]:.3f} Å) to the product end "
      f"({dprim[-1]:.3f} Å): {'most' if frac > 50 else 'less than half'} of the {verb} has happened by the TS, and the "
      f"remaining {100 - frac:.0f} % ({abs(dprim[its] - dprim[-1]):.2f} Å) comes after it, over the product side "
      f"(s = 0 to {s[-1]:.2f}).")
    w("")
    w("| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |")
    w("|---|---|---|---|---|---|")
    for lab, ro, _, _ in kb:
        d = D[lab]
        tot = d[-1] - d[0]
        sh = f"{100 * (d[its] - d[0]) / tot:.0f} %" if abs(tot) > 0.02 else f"— (changes by {tot:+.3f} Å only)"
        w(f"| {lab} | {ro} | {d[0]:.3f} | {d[its]:.3f} | {d[-1]:.3f} | {sh} |")
    w("")
    w("| Point | Frame | s (bohr·amu½) | " + " | ".join(f"{lab} (Å)" for lab, *_ in kb) +
      f" | E, M06-2X/ma-def2-SVP, zero = {((rnm + f' (`{react[0]}`)') if zoff else rnm) if uni else 'separated reactants'} |")
    w("|---|---|---|" + "---|" * len(kb) + "---|")
    picks = [("reactant end", 0)]
    for target in (-10, -6, -4, -2, -1, -0.5):
        k = int(np.argmin(abs(s - target)))
        if 0 < k < its and abs(s[k] - target) < 0.3:
            picks.append((f"s ≈ {target:g}", k))
    picks.append(("TS", its))
    for target in (0.5, 1, 1.5, 2, 3, 6, 10):
        k = int(np.argmin(abs(s - target)))
        if its < k < n - 1 and abs(s[k] - target) < 0.3:
            picks.append((f"s ≈ +{target:g}", k))
    picks.append(("product end", n - 1))
    for lab, k in picks:
        w(f"| {lab} | {k} | {s[k]:+.2f} | " + " | ".join(f"{D[b][k]:.3f}" for b, *_ in kb) + f" | {t(Ew[k], 1, True)} |")
    w("")

    # ---- 2. stationary points
    w(f"## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD {solvent})")
    w("")
    rates = ", ".join(f"{lab} {r:+.2f}" for lab, _, _, r in kb)
    tog = tsj["mode_check"].get("forming_and_breaking_bonds_move_together")
    geo = [f"{lab} {v:.3f} Å" for lab, v in tsj["bonds_A"].items()]
    geo += [f"{lab} {v:.1f}°" for lab, v in (tsj.get("angles_deg") or {}).items()]
    geo += [f"{lab} {v:.1f}°" for lab, v in (tsj.get("dihedrals_deg") or {}).items()]
    w(f"- **TS:** {tsj['n_imaginary']} imaginary frequency, {tsj['imaginary_cm'][0]:.1f} cm⁻¹ (written "
      f"{abs(tsj['imaginary_cm'][0]):.1f}i cm⁻¹). Stretch rates of the key bonds in the mode: {rates}"
      + ("; forming and breaking bonds move together" if tog else ("; forming and breaking bonds do not move together"
         if tog is not None else "")) + f". {', '.join(geo)}.")
    w("- **References** (structures the barriers and reaction energies are measured from):")
    w("")
    ztab = (f"{rnm} (`{react[0]}`)" if zoff else rnm) if uni else "separated reactants"
    w(f"| Structure | Role | Key distances | Imaginary frequencies | E, zero = {ztab} | G (298.15 K, 1 M), zero = {ztab} |")
    w("|---|---|---|---|---|---|")
    e0 = sum(Ework[x] for x in react)
    g0 = sum(Gwork[x] for x in react)

    def EGrel(k, nmol):
        de = (Ework[k] - e0) * EH
        dg = (Gwork[k] - g0) * EH - (len(react) - nmol) * CONV_1M
        return de, dg
    kd = lambda b: ", ".join(f"{lab} {v:.3f} Å" for lab, v in b.items())
    if not uni:
        for k in react:
            g = geom(outs[k])
            w(f"| {k} (`{frag[k]['smiles']}`, job `{frag[k]['job']}`) | separated reactant | {bond_types(g)} | none | — | — |")
    elif zoff:
        g = geom(outs[react[0]])
        bz = {lab: float(np.linalg.norm(g.coords[a] - g.coords[b])) for lab, _, (a, b), _ in kb}
        w(f"| {rnm}, `{react[0]}` (job `{zc[react[0]]['job']}`) | lowest-free-energy conformer; the zero | {kd(bz)} | none | 0 | 0 |")
    for k, role in (("reactant_complex", "reactant end of the path"), ("ts", ""), ("product_complex", "product end of the path")):
        b = tsj["bonds_A"] if k == "ts" else ends[k]["bonds_A"]
        nim = len(tsj["imaginary_cm"]) if k == "ts" else (len(ends[k]["imaginary_cm"]) or "none")
        rname = ({"reactant_complex": rnm, "product_complex": pnm}.get(k, nice(k)) if uni else nice(k)) +             ("" if k == "ts" or ends[k]["minimum"] == k else f" (`{ends[k]['minimum']}`)")
        if uni and k == "reactant_complex" and not zoff:
            w(f"| {rname} | {role}; the zero | {kd(b)} | {nim} | 0 | 0 |")
            continue
        de, dg = EGrel(k, 1)
        w(f"| {rname} | {role} | {kd(b)} | {nim} | {t(de, 1, True)} | {t(dg, 1, True)} |")
    if not uni:
        for k in prods:
            g = geom(outs[k])
            if len(prods) == 1:
                w(f"| {k} (`{frag[k]['smiles']}`) | separated product | {bond_types(g)} | none | {t(Q[None]['dEr'], 1, True)} | "
                  f"{t(Q[None]['dGr'], 1, True)} |")
            else:
                w(f"| {k} (`{frag[k]['smiles']}`, job `{frag[k]['job']}`) | separated product | {bond_types(g)} | none | — | — |")
        if len(prods) > 1:
            w(f"| {' + '.join(prods)} | separated products | | | {t(Q[None]['dEr'], 1, True)} | {t(Q[None]['dGr'], 1, True)} |")
    w("")
    if uni and zoff:
        w(f"The reaction has one molecule throughout. The zero is the {rnm} conformer lowest in free energy (`{react[0]}`); "
          f"the path starts from another conformer, {path_form}, {t(Q[None]['dE_cx_sep'], 2, True)} kJ/mol above it in E and "
          f"{t(Q[None]['dG_cx_sep'], 2, True)} kJ/mol in G. The product is the {pnm} at the product end.")
    elif uni:
        w(f"The reaction has one molecule throughout: the zero is the {rnm} at the reactant end of the path, and the "
          f"product is the {pnm} at the product end.")
    elif len(prods) == 1:
        dpe = (Ework['product_complex'] - sum(Ework[x] for x in prods)) * EH
        if abs(dpe) < 1:
            w(f"The product is one molecule, so the product complex and the separated product are the same species "
              f"optimized twice (from the path end and on its own); they differ by {dpe:+.2f} kJ/mol in E. Reaction "
              "energies use the separated product.")
        else:
            w(f"The product is one molecule. The product complex is the form the path reaches, the separated product the "
              f"same molecule optimized on its own; here they are different conformers, {dpe:+.2f} kJ/mol apart in E "
              "(product complex above), so the reaction energies (separated product) and the energy of the path's "
              "product end differ by that much (card notes, section 7).")
    else:
        w(f"The products are {len(prods)} separate molecules ({' + '.join(prods)}); the product complex holds them together "
          "at the end of the path. Reaction energies use the separated products.")
    if len(react) == 1 and not uni:
        dre = Q[None]['dE_cx_sep']
        w(f"The reactant is one molecule, so the reactant complex and the separated reactant are the same species "
          f"optimized twice (from the path end and on its own); they differ by {dre:+.2f} kJ/mol in E."
          if abs(dre) < 1 else
          f"The reactant is one molecule. The reactant complex is the form the path reaches, the separated reactant the "
          f"same molecule optimized on its own; here they are different conformers, {dre:+.2f} kJ/mol apart in E.")
    for k in refs_other:
        of = frag[k].get("conformer_of")
        of = of[0] if isinstance(of, list) else of
        partner = of or next((x for x in frag if k in ([frag[x].get("conformer_of")] if isinstance(
            frag[x].get("conformer_of"), str) else frag[x].get("conformer_of") or [])), None)
        if partner in frag:
            w(f"Conformer kept on record, not the zero: {k} (`{frag[k]['smiles']}`) lies "
              f"{(Ework[k] - Ework[partner]) * EH:+.2f} kJ/mol in E and {(Gwork[k] - Gwork[partner]) * EH:+.2f} kJ/mol in G "
              f"from {partner} (SI entry 9).")
    w(f"The reactant end of the path ended {ends['reactant_complex']['end_rmsd_to_minimum_A']:.3f} Å (RMSD) from the "
      f"reactant-end minimum, the product end {ends['product_complex']['end_rmsd_to_minimum_A']:.3f} Å from the product-end "
      "minimum (section 6).")
    w("")

    # ---- 3. barriers
    w("## 3. Barriers and reaction energies")
    w("")
    w("All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level "
      "thermal correction. kJ/mol (kcal/mol).")
    w("")
    if uni and own:
        w(f"**The TS of the path shown, measured from {zero_name}, 298.15 K. Not the barrier of any path:** "
          f"`{own['name']}` does not reach this TS (its own barrier is in the table after next).")
    elif uni:
        w(f"**Barriers, zero = {zero_name}, 298.15 K:**")
    else:
        w(f"**Barriers, zero = {zero_name}, 298.15 K, ΔG at 1 M:**")
    w("")
    w("| Level | ΔE‡ | ΔH‡ | ΔG‡ |")
    w("|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        q = Q[lv]
        w(f"| {nm} | {t(q['dE_sep'])} | {t(q['dH_sep'])} | {t(q['dG_sep'])} |")
    w("")
    if not uni or zoff:
        w(f"**Barriers, zero = {path_form if uni else 'reactant complex'}, 298.15 K:**")
        w("")
        w("| Level | ΔE‡ | ΔH‡ | ΔG‡ |")
        w("|---|---|---|---|")
        for nm, lv, _ in LEVELS:
            q = Q[lv]
            w(f"| {nm} | {t(q['dE_cx'])} | {t(q['dH_cx'])} | {t(q['dG_cx'])} |")
        w("")
        if uni and own:
            w(f"**Barrier of `{own['name']}` through its own TS (`{own['ts_key']}`, {abs(own['imaginary_cm'][0]):.1f}i "
              f"cm⁻¹, `{own['ts_folder']}`; "
              + (f"connects `{own['name']}` with the {pnm}: {own['irc']}; the optimizations displaced 0.10 Å both "
                 f"ways along the imaginary mode agree (`01_stationary/rotation_ts.md`)" if own.get("irc")
                 else f"confirmed to connect `{own['name']}` with the {pnm}") + "), zero = "
              f"`{own['name']}`, 298.15 K:**")
            w("")
            w("| Level | ΔE‡ | ΔH‡ | ΔG‡ |")
            w("|---|---|---|---|")
            for nm, lv, _ in LEVELS:
                q = own["rows"][lv]
                w(f"| {nm} | {t(q['dE'])} | {t(q['dH'])} | {t(q['dG'])} |")
            w("")
        for i, (nm1, lv1, _) in enumerate(LEVELS):
            for nm2, lv2, _ in LEVELS[i + 1:]:
                a, b = Q[lv1], Q[lv2]
                if all(round(a[k], 1) == round(b[k], 1) for k in ("dE_cx", "dH_cx", "dG_cx")):
                    w(f"{nm1} and {nm2} give the same rounded barriers from the complex by coincidence, not by copying: "
                      f"ΔE‡ {a['dE_cx']:.2f} and {b['dE_cx']:.2f} kJ/mol, ΔG‡ {a['dG_cx']:.2f} and {b['dG_cx']:.2f} kJ/mol "
                      f"(their barriers from separated reactants differ by {abs(a['dE_sep'] - b['dE_sep']):.1f} kJ/mol, and so "
                      f"do their complex binding energies).")
                    w("")
    if uni:
        pn = f"{pnm} (product-end minimum, `{ends['product_complex']['minimum']}`)"
        w(f"**Reaction, zero = {zero_name}; product = the {pn}, 298.15 K:**")
    else:
        pn = "separated " + " + ".join(prods)
        w(f"**Reaction, zero = separated reactants; product{'s' if len(prods) > 1 else ''} = {pn}, 298.15 K, ΔG at 1 M:**")
    w("")
    revn = (prods[0] if len(prods) == 1 else "the product complex") if not uni else f"the {pnm}"
    cxc = not uni or zoff
    cxh = " Reactant complex, ΔE / ΔG |" if not uni else f" {rnm[0].upper() + rnm[1:]} as the path reaches it, ΔE / ΔG |"
    head = "| Level | ΔrE | ΔrH | ΔrG |" + (cxh if cxc else "") + f" Reverse barrier from {revn}, ΔE‡ / ΔG‡ |"
    w(head)
    w("|---|---|---|---|" + ("---|" if cxc else "") + "---|")
    for nm, lv, _ in LEVELS:
        q = Q[lv]
        w(f"| {nm} | {t(q['dEr'], 1, True)} | {t(q['dHr'], 1, True)} | {t(q['dGr'], 1, True)} | " +
          (f"{t(q['dE_cx_sep'], 2 if uni else 1, True)} / {t(q['dG_cx_sep'], 2 if uni else 1, True)} | " if cxc else "") +
          f"{t(q['dE_rev'])} / {t(q['dG_rev'])} |")
    w("")
    w("ΔrH was not asked for; it is given because it costs nothing.")
    if not uni and len(react) > 1:
        eb = [Q[lv]["dE_cx_sep"] < 0 for _, lv, _ in LEVELS]
        gb = [Q[lv]["dG_cx_sep"] < 0 for _, lv, _ in LEVELS]
        if all(eb) and not any(gb):
            L[-1] += (" The reactant complex lies below separated reactants in E but above them in G at 1 M at every level: "
                      "the binding energy does not pay for the entropy lost on association.")
        elif all(eb) and all(gb):
            L[-1] += " The reactant complex lies below separated reactants in both E and G (1 M) at every level."
        else:
            L[-1] += (" The reactant complex lies " + "; ".join(
                f"{'below' if Q[lv]['dE_cx_sep'] < 0 else 'above'} in E and {'below' if Q[lv]['dG_cx_sep'] < 0 else 'above'} "
                f"in G at {nm.split(' (')[0]}" for nm, lv, _ in LEVELS) + " (separated reactants as zero, 1 M).")
    w("")

    # ---- 4. highest point
    w("## 4. Highest point of the path at each level")
    w("")
    w(f"Path energies at each level, zero = {zero_name} at that level. The TS is frame {its} "
      f"(s = 0, {prim} {dprim[its]:.3f} Å).")
    w("")
    w(f"| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its {prim} (Å) | Highest − TS frame |")
    w("|---|---|---|---|---|---|---|")
    for nm, lv, _ in LEVELS:
        v = path_E[lv]
        if its not in v:
            w(f"| {nm} | {len(v)} | not yet computed | | | | |")
            continue
        km = max(v, key=v.get)
        w(f"| {nm} | {len(v)} of {n} | {t(v[its])} | {km} | {s[km]:+.2f} | {dprim[km]:.3f} | {t(v[km] - v[its], 2, True)} |")
    w("")
    later = "shorter" if dprim[-1] < dprim[0] else "longer"
    w(f"A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP "
      f"path, at a {later} {prim} distance. The difference is the energy that a single point at the TS geometry misses.")
    w("")

    # ---- 5. DLPNO diagnostics
    w("## 5. DLPNO-CCSD(T) diagnostics")
    w("")
    t1 = {}
    for k in dict.fromkeys(["ts", "reactant_complex", "product_complex"] + react + prods + refs_other):
        o = R / "01_ts/diagnostics" / DL / k / f"{k}.out"
        if not o.exists():
            continue
        tx = o.read_text(errors="replace")
        m = re.findall(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", tx)
        amp = re.search(r"LARGEST AMPLITUDES\n-+\n(.*?)\n", tx)
        t1[k] = (float(m[-1]) if m else None, amp.group(1).split()[-1] if amp else None)
    w("| Structure | T1 | Largest amplitude |")
    w("|---|---|---|")
    for k, (a, b) in t1.items():
        w(f"| {k} | {'—' if a is None else f'{a:.4f}'} | {b or '—'} |")
    pt = {}
    for k in path_E[DL]:
        o = R / "02_irc/path_sp" / DL / f"f{k:03d}" / "f.out"
        m = re.findall(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", o.read_text(errors="replace")) if o.exists() else []
        if m:
            pt[k] = float(m[-1])
    w("")
    if pt:
        km = max(pt, key=pt.get)
        w(f"- **Path frames ({len(pt)} computed):** T1 from {min(pt.values()):.4f} to {max(pt.values()):.4f}, median "
          f"{statistics.median(pt.values()):.4f}; largest at frame {km} (s = {s[km]:+.2f}, {prim} {dprim[km]:.3f} Å).")
    allt = [a for a, _ in t1.values() if a is not None] + list(pt.values())
    hi_t1 = [f"{k} ({a:.4f})" for k, (a, _) in t1.items() if a is not None and a >= 0.02] + \
            [f"frame {k} ({v:.4f})" for k, v in pt.items() if v >= 0.02]
    if allt and not hi_t1:
        w("- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a "
          "single-reference treatment is adequate for every structure and frame computed.")
    elif hi_t1:
        w(f"- **Reading:** T1 reaches 0.02 or more, the usual limit for closed-shell molecules (Lee and Taylor 1989), for "
          f"{', '.join(hi_t1[:8])}{' and others' if len(hi_t1) > 8 else ''}: the coupled-cluster energies there need care.")
    w("- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.")
    w("")

    # ---- 6. validation
    w("## 6. Validation")
    w("")
    ang = (R / "02_irc/true_irc/angles.md").read_text(encoding="utf-8")
    verdict = re.search(r"Path: \*\*(.+?)\*\*", ang).group(1)
    w(f"- **Angle check: {verdict.lower()}** (rule: in each zone, median < 5°, 90th percentile < 15°; "
      "mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):")
    w("")
    w("| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |")
    w("|---|---|---|---|---|---|")
    for ln in ang.splitlines():
        c = [x.strip() for x in ln.strip().strip("|").split("|")]
        if c and c[0] in SIDE:
            w(f"| {SIDE[c[0]]} | {c[2]} | {c[3]} | {c[4]} | {c[5]} | {c[6].strip('*').lower()} |")
    w("")
    ecps = [p for p in (R / "02_irc/true_irc/end_check_all.md", R / "02_irc/true_irc/end_check.md") if p.exists()]
    ecs = [p.read_text(encoding="utf-8") for p in ecps]
    vm = re.search(r"Verdict: \*\*(\w+)\*\*", ecs[0])
    w(f"- **End checks: {vm.group(1) if vm else 'see file'}** ({', '.join(f'`02_irc/true_irc/{p.name}`' for p in ecps)}; "
      "from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end "
      "structure, SI entry 4):")
    best = {}
    for ec in ecs:                                   # per side: the latest file's optimized check, else the older one
        for ln in ec.splitlines():
            c = [x.strip() for x in ln.strip().strip("|").split("|")]
            if c and c[0] in SIDE and len(c) >= 11:
                if c[0] not in best or (best[c[0]][4] != "yes" and c[4] == "yes"):
                    best[c[0]] = c
    for sd in SIDE:
        if sd in best:
            c = best[sd]
            side = SIDE[c[0]][0].upper() + SIDE[c[0]][1:]
            if c[4] == "yes":
                w(f"  - {side}: last frame {c[2]} kJ/mol from `{c[1]}` (RMSD {c[3]} Å). Opt Freq from the last frame: "
                  f"E {c[5]} kJ/mol from it, heavy-atom RMSD {c[6]} Å (all-atom {c[7]} Å), largest key-distance "
                  f"difference {c[8]} Å, {c[9]} imaginary frequencies: {c[10]}.")
            else:
                w(f"  - {side}: {c[10]} (last frame {c[2]} kJ/mol from `{c[1]}`).")
    w(f"  - Official end structures now: reactant end `{ends['reactant_complex']['minimum']}`, product end "
      f"`{ends['product_complex']['minimum']}` (`02_irc/endpoints.json`).")
    si = (ROOT / "docs/si_corrections.md").read_text(encoding="utf-8")
    if snotes.get("si_row"):
        row = next(ln for ln in si.splitlines() if ln.strip().startswith(f"| {snotes['si_row']} |"))
        c = [x.strip() for x in row.strip().strip("|").split("|")]
        w(f"- **Side-end table (SI entry 5):** reactant side: {c[1]}; product side: {c[2]} (last frame above the official "
          "end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not "
          "extended to it; the paper lists how far above its minimum each side ended.")
    if snotes.get("si_notes"):
        w("- **SI corrections that touch this reaction** (`docs/si_corrections.md`):")
        for x in snotes["si_notes"]:
            w(f"  - {x}")
    w("")

    # ---- 7. experiment
    w("## 7. Experiment")
    w("")
    er = dict(notes.get("experiment_row", {}))
    er.update(snotes.get("measured", {}))          # measured values given only in the card text (summary_notes.toml)
    if "measured_barrier" in notes and "dG_TS_vs_separated_1M_kJ" not in er:
        er["dG_TS_vs_separated_1M_kJ"] = notes["measured_barrier"]["dG_kJ"]
    if notes.get("experimental"):
        w(notes["experimental"].strip())
        if not notes.get("experimental_confirmed", False):
            w("")
            w("These measured values are not yet confirmed by J. Picione (card notes).")
    else:
        w("No measured value is recorded for this reaction (`00_inputs/card_notes.toml` has none).")
    w("")
    cols = [(k, h) for k, h in (("dG_TS_vs_separated_1M_kJ", "ΔG‡"), ("dG_reaction_1M_kJ", "ΔrG"),
                                ("dG_reverse_kJ", f"Reverse ΔG‡ (from {revn})")) if k in er]
    if cols:
        qk = {"dG_TS_vs_separated_1M_kJ": "dG_sep", "dG_reaction_1M_kJ": "dGr", "dG_reverse_kJ": "dG_rev"}
        w(f"**Computed against measured** (computed − measured; zero = {zero_name}, 298.15 K, 1 M):")
        w("")
        w("| Level | " + " | ".join(h for _, h in cols) + " |")
        w("|---|" + "---|" * len(cols))
        w("| Measured | " + " | ".join(t(er[k], 1, k == "dG_reaction_1M_kJ") for k, _ in cols) + " |")
        for nm, lv, _ in LEVELS:
            q = Q[lv]
            w(f"| {nm} | " + " | ".join(f"{t(q[qk[k]], 1, k == 'dG_reaction_1M_kJ')}, {t(q[qk[k]] - er[k], 1, True)}"
                                        for k, _ in cols) + " |")
        w("")
    if notes.get("limitations"):
        w("Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):")
        w("")
        for x in notes["limitations"]:
            w(f"- {x}")
        w("")

    # ---- 8. job times
    w("## 8. Job times")
    w("")
    w(f"ORCA job time (TOTAL RUN TIME in each output), {rx.cfg['nprocs']} cores per job; several jobs ran side by side, "
      "so the wall-clock span of a level is longer than its job time divided by the number of lanes.")
    w("")
    w("| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |")
    w("|---|---|---|---|---|---|---|---|")
    for nm, lv, _ in LEVELS[1:]:
        jt = jobt.get(lv) or []
        if not jt:
            w(f"| {nm} | 0 | | | | | | |")
            continue
        secs = [x[1] for x in jt]
        tsf = next((x[1] for x in jt if x[0] == its), None)
        slow = max(jt, key=lambda x: x[1])
        st = min(x[2] - x[1] for x in jt)
        en = max(x[2] for x in jt)
        w(f"| {nm} | {len(jt)} | {sum(secs) / 3600:.1f} | {statistics.median(secs):.0f} | "
          f"{'—' if tsf is None else f'{tsf:.0f}'} | {slow[1]:.0f} (frame {slow[0]}) | "
          f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(st))} | {time.strftime('%Y-%m-%d %H:%M', time.localtime(en))} |")
    w("")
    sides = []
    for side in ("backward", "forward"):
        secs = []
        dirs = sorted(p.name for p in (R / "02_irc/true_irc").glob(f"{side}*") if p.is_dir())
        for d in dirs:
            for o in (R / "02_irc/true_irc" / d).glob("g*/j.out"):
                sec = run_seconds(o)
                if sec is not None:
                    secs.append(sec)
        if secs:
            sides.append(f"{SIDE[side]}: {len(secs)} gradient jobs, {sum(secs) / 3600:.1f} h job time "
                         f"(folders {', '.join(f'`{d}`' for d in dirs)})")
    if sides:
        w("- **Path (IRC) gradient jobs, working level:** " + "; ".join(sides) + ". The Hessian jobs and the TS and end-"
          "minimum optimizations are not counted.")
        w("")
    irc_rel = next((c["own_ts"].get("irc") for c in rx.cfg.get("zero_conformers", []) if c.get("own_ts")), None)
    if irc_rel:                                  # a separate check, never added to the path (J. Picione, 2026-11-27)
        ic = json.loads((R / irc_rel).read_text(encoding="utf-8"))
        idir = (R / irc_rel).parent
        n_g = {s: len(list((idir / s).glob("g*/j.out"))) for s in ("backward", "forward")}
        w(f"- **Separate check, not part of the path or its total:** IRC from `{own['ts_key'] if own else '?'}` "
          f"(`{Path(irc_rel).parent.as_posix()}`), {ic['job_time_h']:.2f} h of job time, {ic['nprocs']} cores: "
          f"{n_g['backward']} + {n_g['forward']} gradient jobs (backward + forward) and one Opt Freq at each end.")
        w("")
    # ---- 9+. extra sections (summary_notes.toml: extra_sections = [{title, file}], file relative to the reaction
    # folder, written by that check's own script; copied in verbatim)
    for n, xs in enumerate(snotes.get("extra_sections", []), start=9):
        w(f"## {n}. {xs['title']}")
        w("")
        w((R / xs["file"]).read_text(encoding="utf-8").strip())
        w("")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text("\n".join(L), encoding="utf-8")
    print("written", out, "| status:", status)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
