"""Diels–Alder gas-phase comparison, butadiene conformer ensemble (J. Picione, 2026-11-21).

Usage:  python -X utf8 -m pipeline.gas_conf <reaction folder> run|md

Follows pipeline/gas_opt.py (03_gas_phase). At 760–921 K the gauche s-cis form of butadiene is not a small
minority, so the reactant is the equilibrium conformer mixture, as in the experiment:
  1. gas-phase single points on the already optimized gauche butadiene at the three single-point levels
     (run only after his machine block is released: the job waits for docs/level_test_queue/window_release);
  2. conformer-weighted butadiene partition function: s-trans + both mirror-image gauche forms,
       G(butadiene, mixture) = −RT ln[ exp(−G_trans/RT) + 2 exp(−G_gauche/RT) ],
     each G = E(level) + [G(T) − E] from the working-level gas-phase frequencies (quasi-RRHO, 1 atm), with
     σ(s-trans, C2h) = 2 and σ(gauche, C2) = 2 for one enantiomer; the factor 2 counts the two enantiomers;
  3. k(TST) = (kB T/h)(RT/p°) exp(−ΔG‡/RT) at 760, 800, 840, 880 and 921 K, then a least-squares fit of
     ln k against 1/T (as Rowley and Steiner did) gives Ea = −R·slope and A = exp(intercept), for the mixture
     and, with the same fit, for s-trans butadiene alone.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from . import gas_opt, orca, thermo
from .gas_opt import SUB, _sci, levels
from .geometry import read_xyz
from .steps import Reaction, _read_json, _write_json
from .units import fmt_kj, kj

TEMPS = [760.0, 800.0, 840.0, 880.0, 921.0]
FRACTION_T = [760.0, 840.0, 921.0]
R_KJ = thermo.R_KJ
Q = Path(__file__).resolve().parents[1] / "docs/level_test_queue"
WORK = "M06-2X_ma-def2-SVP (working)"
ORDER = ["M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]


def say(msg: str) -> None:
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def c2_check(xyz_path: Path) -> dict:
    """Best C2 axis among the principal axes: RMSD (Å) between the structure and its 180° rotation."""
    g = read_xyz(xyz_path)
    X = np.array(g.coords, float)
    m = g.masses()
    X = X - (m[:, None] * X).sum(0) / m.sum()
    I = sum(mi * (np.dot(r, r) * np.eye(3) - np.outer(r, r)) for mi, r in zip(m, X))
    _, ax = np.linalg.eigh(I)
    best = None
    for k in range(3):
        u = ax[:, k]
        Rm = 2 * np.outer(u, u) - np.eye(3)            # 180° rotation about u
        Y = X @ Rm.T
        cost = np.linalg.norm(X[:, None] - Y[None], axis=2)
        cost[np.array(g.symbols)[:, None] != np.array(g.symbols)[None]] = 1e6
        r, c = linear_sum_assignment(cost)
        rmsd = float(np.sqrt((cost[r, c] ** 2).mean()))
        if best is None or rmsd < best:
            best = rmsd
    dih = _dihedral(X[0], X[1], X[2], X[3])
    return {"c2_rmsd_A": best, "CCCC_dihedral_deg": dih}


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def run(rx: Reaction) -> None:
    root = rx.root / SUB
    rel = Q / "window_release"
    if not rel.exists():
        say("waiting for the machine block to be released (docs/level_test_queue/window_release)")
        (root / "HOLD_conformers.txt").write_text(
            f"{time.strftime('%Y-%m-%d %H:%M')}: gauche butadiene single points held until J. Picione's machine "
            "block is released (his request: run after the block ends).\n", encoding="utf-8")
        while not rel.exists():
            time.sleep(60)
        time.sleep(120)
        (root / "HOLD_conformers.txt").unlink(missing_ok=True)
        say("block released")
    opt = _read_json(root / "opt.json")
    row = opt["butadiene_gauche"]
    g = read_xyz(rx.root / row["dir"] / f"{row['job']}.xyz")
    out_p = root / "sp_gauche.json"
    E = _read_json(out_p) if out_p.exists() else {}
    for lv in sorted(levels(rx), key=lambda l: not l["label"].startswith("DLPNO")):
        if lv["label"] in E:
            continue
        wd = rx.d(f"{SUB}/sp/{lv['label']}/butadiene_gauche")
        inp = orca.make_input(lv["method"], g, rx.cfg["charge"], rx.cfg["mult"], nprocs=rx.cfg["nprocs"],
                              maxcore_mb=rx.cfg["maxcore_mb"])
        r = orca.run(wd, "butadiene_gauche", inp, rx.log, front=True,
                     note=f"gas-phase conformer ensemble: {lv['label']} single point on the gas-phase gauche "
                          "butadiene, no solvent model")
        if not r.ok:
            (root / "STOPPED_conformers.md").write_text(
                f"Stopped: gauche {lv['label']} single point failed; see {r.out}. Not a result.\n", encoding="utf-8")
            say(f"STOPPED: gauche {lv['label']} failed")
            return
        E[lv["label"]] = r.energy_eh
        _write_json(out_p, E)
        say(f"done  {lv['label']} gauche ({r.wall_s:.0f} s)")
    analyse(rx)
    say("finished")


def _th(rx: Reaction, row: dict, T: float, sigma: int) -> dict:
    hp = rx.root / row["dir"] / f"{row['job']}.hess"
    return thermo.thermo(hp, T, sigma)


def _fit(Ts: list[float], ks: list[float]) -> tuple[float, float, float]:
    x = 1.0 / np.array(Ts)
    y = np.log(np.array(ks))
    slope, icpt = np.polyfit(x, y, 1)
    resid = y - (slope * x + icpt)
    return -slope * R_KJ, math.exp(icpt), float(np.abs(resid).max())


def analyse(rx: Reaction) -> dict:
    root = rx.root / SUB
    opt = _read_json(root / "opt.json")
    E = _read_json(root / "sp.json")
    Eg = _read_json(root / "sp_gauche.json") if (root / "sp_gauche.json").exists() else {}
    gp = _read_json(root / "gas_phase.json")
    x = gp["experiment"]
    A_exp, Ea_exp = x["A_L_mol_s"], x["Ea_kJ"]
    sig = {"ts": 1, "ethylene": 4, "butadiene_s-trans": 2, "butadiene_gauche": 2}
    orca_sig = {sp: thermo.symmetry_number(rx.root / opt[sp]["dir"] / f"{opt[sp]['job']}.out") for sp in sig}
    c2 = c2_check(rx.root / opt["butadiene_gauche"]["dir"] / "butadiene_gauche.xyz")
    th = {T: {sp: _th(rx, opt[sp], T, sig[sp]) for sp in sig} for T in TEMPS + [298.15]}
    tabs = [(WORK, {sp: opt[sp]["E_eh"] for sp in opt})]
    for lab in ORDER:
        if lab in E and lab in Eg:
            tabs.append((lab, {**E[lab], "butadiene_gauche": Eg[lab]}))
    rows = []
    for lab, En in tabs:
        G = lambda sp, T: kj(En[sp]) + th[T][sp]["G_minus_E_kJ"]
        H = lambda sp, T: kj(En[sp]) + th[T][sp]["H_minus_E_kJ"]
        per_T = []
        for T in TEMPS + [298.15]:
            RT = R_KJ * T
            gt, gg = G("butadiene_s-trans", T), G("butadiene_gauche", T)
            z = 1.0 + 2.0 * math.exp(-(gg - gt) / RT)
            G_mix = gt - RT * math.log(z)
            frac_g = 1.0 - 1.0 / z
            kT_h = thermo.KB * T / thermo.H_PL
            v0 = R_KJ * 1000.0 * T / thermo.P_ATM * 1000.0
            dG_tr = G("ts", T) - G("ethylene", T) - gt
            dG_mix = G("ts", T) - G("ethylene", T) - G_mix
            per_T.append({"T_K": T, "gauche_fraction": frac_g,
                          "dG_gauche_minus_trans_one_enantiomer_kJ": gg - gt,
                          "dH_gauche_minus_trans_kJ": H("butadiene_gauche", T) - H("butadiene_s-trans", T),
                          "dG_act_s_trans_kJ": dG_tr, "dG_act_mixture_kJ": dG_mix,
                          "k_s_trans_L_mol_s": kT_h * v0 * math.exp(-dG_tr / RT),
                          "k_mixture_L_mol_s": kT_h * v0 * math.exp(-dG_mix / RT),
                          "k_exp_L_mol_s": A_exp * math.exp(-Ea_exp / RT)})
        fitT = [p for p in per_T if p["T_K"] in TEMPS]
        Ea_t, A_t, r_t = _fit(TEMPS, [p["k_s_trans_L_mol_s"] for p in fitT])
        Ea_m, A_m, r_m = _fit(TEMPS, [p["k_mixture_L_mol_s"] for p in fitT])
        p840 = next(p for p in per_T if p["T_K"] == 840.0)
        rows.append({"level": lab, "per_T": per_T,
                     "s_trans": {"Ea_kJ": Ea_t, "A_L_mol_s": A_t, "A_ratio": A_t / A_exp, "Ea_minus_exp_kJ": Ea_t - Ea_exp,
                                 "max_fit_residual_ln_k": r_t, "k840_ratio": p840["k_s_trans_L_mol_s"] / p840["k_exp_L_mol_s"]},
                     "mixture": {"Ea_kJ": Ea_m, "A_L_mol_s": A_m, "A_ratio": A_m / A_exp, "Ea_minus_exp_kJ": Ea_m - Ea_exp,
                                 "max_fit_residual_ln_k": r_m, "k840_ratio": p840["k_mixture_L_mol_s"] / p840["k_exp_L_mol_s"]}})
    out = {"temperatures_K": TEMPS, "symmetry_numbers_used": sig, "orca_symmetry_numbers": orca_sig,
           "gauche_c2_check": c2, "rows": rows,
           "experiment": {"Ea_kJ": Ea_exp, "Ea_err_kJ": x["Ea_err_kJ"], "A_L_mol_s": A_exp},
           "levels_missing_gauche_sp": [l for l in ORDER if l not in Eg]}
    _write_json(root / "conformers.json", out)
    write_md(out)
    gas_opt.write_md(rx, gp)                            # gas_phase.md now ends with the conformer section
    return out


def _pm(v: float) -> str:
    s = fmt_kj(v)
    return ("+" + s.replace("(", "(+")) if v > 0 else s.replace("-", "−")


def write_md(out: dict) -> None:
    x = out["experiment"]
    c2 = out["gauche_c2_check"]
    L = ["**Butadiene conformer ensemble (J. Picione, 2026-11-21).** At 760–921 K the reactant is the equilibrium "
         "mixture of s-trans butadiene and the two mirror-image gauche (s-cis) forms, as in the experiment. "
         "k(TST) was computed at 760, 800, 840, 880 and 921 K and ln k fitted against 1/T, as Rowley and Steiner did; "
         "Ea = −R × slope, A = exp(intercept), in L mol⁻¹ s⁻¹. The s-trans-only values below use the same fit, so "
         "they differ slightly from the 840 K point values (Ea = ΔH‡ + 2RT) in the table above.", "",
         "| Level | Reactant butadiene | Ea (fit, 760–921 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | k(840 K) / measured k |",
         "|---|---|---|---|---|---|---|"]
    for r in out["rows"]:
        for key, name in (("s_trans", "s-trans only"), ("mixture", "conformer mixture")):
            f = r[key]
            L.append(f"| {r['level'] if key == 's_trans' else ''} | {name} | {fmt_kj(f['Ea_kJ'])} | {_pm(f['Ea_minus_exp_kJ'])} | "
                     f"{_sci(f['A_L_mol_s'])} | {f['A_ratio']:.2f} | {f['k840_ratio']:.2f} |")
    L.append(f"| Measured (Rowley and Steiner 1951) | equilibrium mixture | {x['Ea_kJ']:.1f} ± {x['Ea_err_kJ']} kJ/mol "
             f"({x['Ea_kJ'] / 4.184:.1f} ± {x['Ea_err_kJ'] / 4.184:.1f} kcal/mol) | | {_sci(x['A_L_mol_s'])} | 1 | 1 |")
    L += ["", "**Gauche fraction of butadiene** (both enantiomers together; equilibrium in the gas phase):", "",
          "| Level | 298 K | 760 K | 840 K | 921 K | G(gauche, one enantiomer) − G(s-trans), 840 K |", "|---|---|---|---|---|---|"]
    for r in out["rows"]:
        p = {q["T_K"]: q for q in r["per_T"]}
        L.append(f"| {r['level']} | {p[298.15]['gauche_fraction'] * 100:.1f} % | {p[760.0]['gauche_fraction'] * 100:.1f} % | "
                 f"{p[840.0]['gauche_fraction'] * 100:.1f} % | {p[921.0]['gauche_fraction'] * 100:.1f} % | "
                 f"{fmt_kj(p[840.0]['dG_gauche_minus_trans_one_enantiomer_kJ'])} |")
    os_ = out["orca_symmetry_numbers"]
    L += ["", "**How the symmetry numbers and the two enantiomers are counted.** G(butadiene, mixture) = "
          "−RT ln[exp(−G(s-trans)/RT) + 2 exp(−G(gauche)/RT)]. s-trans butadiene is C2h, σ = 2. Gauche butadiene is "
          f"C2 (CCCC dihedral {abs(c2['CCCC_dihedral_deg']):.1f}°; the optimized structure matches its own 180° rotation "
          f"within {c2['c2_rmsd_A']:.4f} Å RMSD), so one enantiomer has σ = 2; the factor 2 in front of its term counts "
          "the two mirror-image forms, which are separate minima with equal energy. Together the pair therefore "
          "counts like one structure with σ = 1 (2/σ = 1). "
          f"ORCA's point-group detection labelled the optimized gauche form C1 (σ = {os_['butadiene_gauche']}), "
          "so its own G (and the 9.2 kJ/mol gauche–s-trans ΔG quoted above) already equals the pair value; here σ = 2 "
          "is set explicitly and the factor 2 added. Ethylene σ = 4 (D2h), TS σ = 1 (Cs, achiral, one TS)."]
    if out["levels_missing_gauche_sp"]:
        L += ["", "Not yet computed: gauche single points at " + ", ".join(out["levels_missing_gauche_sp"]) + "."]
    root_md = Path(_ROOT[0]) / "conformers.md"
    fnd = root_md.with_name("findings_conformers.md")
    if fnd.exists():
        L += ["", fnd.read_text(encoding="utf-8").strip()]
    root_md.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


_ROOT: list[str] = [""]


def main() -> None:
    rx = Reaction.load(Path(sys.argv[1]))
    _ROOT[0] = str(rx.root / SUB)
    if sys.argv[2] == "run":
        run(rx)
    else:
        analyse(rx)


if __name__ == "__main__":
    main()
