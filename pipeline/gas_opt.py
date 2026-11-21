"""Gas-phase Diels–Alder comparison with Rowley and Steiner, on gas-phase geometries (J. Picione, 2026-11-20).

Usage:  python -X utf8 -m pipeline.gas_opt <reaction folder> run|md

Gas phase, no solvent model. Working level (the folder's `method`): Opt Freq for ethylene and both butadiene
conformers (s-trans, and the gauche s-cis form the path reaches, to confirm which is lowest), OptTS Freq for the
TS from the toluene TS (one imaginary frequency, same synchronous mode as the toluene TS). Then single points at
M06-2X/ma-def2-TZVP, B3LYP/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP on the gas-phase geometries.

Same method as the earlier folder's gas_phase step (pipeline/gas_phase.py), except that geometries and thermal
corrections now come from gas-phase optimizations and frequencies:
  ΔH‡(T)  = ΔE‡(level) + Δ[H(T) − E] (working-level gas-phase frequencies, harmonic, ideal gas)
  Ea      = ΔH‡(T) + 2RT                     (bimolecular, gas phase, Δn‡ = −1)
  ΔG‡(T)  = ΔE‡(level) + Δ[G(T) − E]          (quasi-RRHO entropy, 1 atm standard state)
  k(TST)  = (kB T/h)(RT/p°) exp(−ΔG‡/RT)      (concentration units, κ = 1, symmetry numbers from ORCA)
  A       = k(TST) exp(Ea/RT)                 (= e² (kB T/h)(RT/p°) exp(ΔS‡/R))
All jobs go to the head of the core queue (J. Picione: ahead of the path single points). The run must finish
before the 03:00 freeze for his machine block; a job that could not finish before then waits until the block is
released.
"""
from __future__ import annotations

import json
import math
import sys
import threading
import time
from pathlib import Path

import numpy as np

from . import orca, thermo
from .geometry import Geometry, read_xyz
from .steps import Reaction, _read_json, _write_json
from .units import fmt_kj, kj

SUB = "03_gas_phase"
T = 840.0
R_KJ = thermo.R_KJ
Q = Path(__file__).resolve().parents[1] / "docs/level_test_queue"
FREEZE = "2026-11-21 02:55"           # five minutes before the reservation (no new job from 03:00)
EST_S = {"opt_ts": 1800, "opt": 600, "DLPNO-CCSD(T)_ma-def2-TZVP": 1500, "dft": 300}
MINIMA = ["ethylene", "butadiene_s-trans", "butadiene_gauche"]
START = {"ethylene": "ethylene", "butadiene_s-trans": "butadiene_s-trans", "butadiene_gauche": "butadiene_s-cis"}
SIGMA_NOTE = "rotational symmetry numbers from ORCA's point-group detection"
_print_lock = threading.Lock()


def say(msg: str) -> None:
    with _print_lock:
        print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def before_freeze(est_s: float, what: str) -> None:
    """Start a job only if it can finish before the freeze; otherwise wait until his block is released."""
    tf = time.mktime(time.strptime(FREEZE, "%Y-%m-%d %H:%M"))
    if time.time() + est_s <= tf:
        return
    rel = Q / "window_release"
    say(f"HOLD {what}: it may not finish before {FREEZE} (estimate {est_s / 60:.0f} min); waiting for the block release")
    (Path(SUB_ROOT[0]) / "HOLD.txt").write_text(
        f"{time.strftime('%Y-%m-%d %H:%M')}: {what} held until the machine block is released "
        f"(estimate {est_s / 60:.0f} min would pass {FREEZE}).\n", encoding="utf-8")
    while not rel.exists():
        time.sleep(60)
    time.sleep(120)
    say(f"block released; continuing with {what}")


SUB_ROOT: list[str] = [""]


def _ts_stem(rx: Reaction) -> Path:
    return Path(f"{rx.ts_stem()}")


def _start_geom(rx: Reaction, sp: str) -> Geometry:
    if sp == "ts":
        return read_xyz(Path(f"{_ts_stem(rx)}.xyz"))
    frag = _read_json(rx.root / "01_ts/fragments.json")
    n = START[sp]
    return read_xyz(rx.root / f"01_ts/fragments/{n}/{frag[n].get('job', n)}.xyz")


def _opt_input(rx: Reaction, sp: str, g: Geometry) -> str:
    if sp == "ts":   # as the toluene TS, without the solvent model
        kw, blocks = f'{rx.cfg["method"]} OptTS Freq', "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end"
    else:            # as the toluene fragments, without the solvent model
        kw, blocks = f'{rx.cfg["method"]} Opt Freq', ""
    return orca.make_input(kw, g, rx.cfg["charge"], rx.cfg["mult"], blocks=blocks, nprocs=rx.cfg["nprocs"],
                           maxcore_mb=rx.cfg["maxcore_mb"])


def _kabsch(P: np.ndarray, Q_: np.ndarray) -> np.ndarray:
    Pc, Qc = P - P.mean(0), Q_ - Q_.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    dd = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, dd]) @ Vt


def ts_check(rx: Reaction, wd: Path, name: str) -> dict:
    """One imaginary frequency; same mode as the toluene TS (overlap ≥ 0.8); synchronous by the folder's rules."""
    g = read_xyz(wd / f"{name}.xyz")
    h = orca.read_hess(wd / f"{name}.hess")
    k = int(np.argmin(h["freqs"]))
    v = h["modes"][:, k].reshape(-1, 3)
    g0 = read_xyz(Path(f"{_ts_stem(rx)}.xyz"))
    h0 = orca.read_hess(Path(f"{_ts_stem(rx)}.hess"))
    v0 = h0["modes"][:, int(np.argmin(h0["freqs"]))].reshape(-1, 3)
    R = _kabsch(g.coords, g0.coords)
    vs = v @ R
    ov = float(abs(np.sum(vs * v0)) / (np.linalg.norm(vs) * np.linalg.norm(v0)))
    bonds = {b["label"]: b["atoms"] for b in rx.cfg["key_bonds"]}
    dist, rate = {}, {}
    for lab, (a, b) in bonds.items():
        u = g.coords[a] - g.coords[b]
        dist[lab] = float(np.linalg.norm(u))
        rate[lab] = float(np.dot(v[a] - v[b], u / np.linalg.norm(u)))
    f1, f2 = rx.cfg["ts_mode"]["synchronous"]
    spread = abs(dist[f1] - dist[f2])
    ratio = min(abs(rate[f1]), abs(rate[f2])) / max(abs(rate[f1]), abs(rate[f2]))
    same_sign = np.sign(rate[f1]) == np.sign(rate[f2])
    in_range = all(lo <= dist[b] <= hi for b, (lo, hi) in rx.cfg["ts_required"].items())
    n_imag = int(sum(1 for f in h["freqs"] if f < 0))
    ok = (n_imag == 1 and ov >= 0.8 and spread <= rx.cfg["ts_mode"]["synchronous_max_spread_A"]
          and ratio >= rx.cfg["ts_mode"]["synchronous_min_rate_ratio"] and bool(same_sign) and in_range)
    return {"n_imaginary": n_imag, "freq_cm": float(h["freqs"][k]), "overlap_with_toluene_mode": round(ov, 3),
            "distances_A": {b: round(d, 4) for b, d in dist.items()}, "forming_spread_A": round(spread, 4),
            "forming_rate_ratio": round(ratio, 3), "forming_same_sign": bool(same_sign),
            "forming_in_required_range": bool(in_range), "accepted": bool(ok)}


def opt_one(rx: Reaction, sp: str) -> dict:
    wd = rx.d(f"{SUB}/opt/{sp}")
    before_freeze(EST_S["opt_ts" if sp == "ts" else "opt"], f"{sp} optimization")
    say(f"start {sp} optimization")
    r = orca.run(wd, sp, _opt_input(rx, sp, _start_geom(rx, sp)), rx.log, front=True,
                 note="gas-phase comparison (J. Picione 2026-11-20): " + ("OptTS Freq" if sp == "ts" else "Opt Freq")
                      + ", no solvent model")
    row = {"species": sp, "job": sp, "wall_s": r.wall_s, "ok": bool(r.ok and r.opt_converged is not False)}
    if row["ok"] and sp != "ts" and r.imaginary:
        h = orca.read_hess(wd / f"{sp}.hess")
        v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
        gp = read_xyz(wd / f"{sp}.xyz")
        gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
        prev, name2 = r.imaginary, f"{sp}_min"
        before_freeze(EST_S["opt"], f"{sp} reoptimization")
        r = orca.run(wd, name2, _opt_input(rx, sp, gp), rx.log, front=True,
                     note=f"rerun: {sp} stopped on a saddle point ({[round(f, 1) for f in prev]} cm-1); "
                          "displaced 0.1 Å along that mode, same keywords")
        row.update(job=name2, first_imaginary_cm=prev, wall_s=row["wall_s"] + r.wall_s,
                   ok=bool(r.ok and r.opt_converged is not False))
    if row["ok"]:
        row.update(E_eh=r.energy_eh, G298_eh=r.gibbs_eh, imaginary_cm=r.imaginary)
        if sp == "ts":
            row["mode"] = ts_check(rx, wd, row["job"])
            row["usable"] = row["mode"]["accepted"]
        else:
            row["usable"] = not r.imaginary
    else:
        row["usable"] = False
    row["dir"] = str(wd.relative_to(rx.root)).replace("\\", "/")
    say(f"done  {sp} optimization: ok={row['ok']} usable={row['usable']} ({r.wall_s:.0f} s)")
    return row


def levels(rx: Reaction) -> list[dict]:
    return [{"label": d["label"], "method": d["method"]} for d in rx.cfg["diagnostic_levels"]]


def sp_one(rx: Reaction, lv: dict, row: dict) -> tuple[float, float]:
    g = read_xyz(rx.root / row["dir"] / f"{row['job']}.xyz")
    wd = rx.d(f"{SUB}/sp/{lv['label']}/{row['species']}")
    est = EST_S.get(lv["label"], EST_S["dft"]) * (1.0 if row["species"] == "ts" else 0.4)
    before_freeze(est, f"{lv['label']} single point on {row['species']}")
    inp = orca.make_input(lv["method"], g, rx.cfg["charge"], rx.cfg["mult"], nprocs=rx.cfg["nprocs"],
                          maxcore_mb=rx.cfg["maxcore_mb"])
    r = orca.run(wd, row["species"], inp, rx.log, front=True,
                 note=f"gas-phase comparison: {lv['label']} single point on the gas-phase {row['job']}, no solvent model")
    if not r.ok:
        raise orca.OrcaError(f"{row['species']} {lv['label']} failed; see {r.out}")
    say(f"done  {lv['label']} {row['species']} ({r.wall_s:.0f} s)")
    return r.energy_eh, r.wall_s


def run(rx: Reaction) -> None:
    root = rx.d(SUB)
    SUB_ROOT[0] = str(root)
    opt: dict[str, dict] = {}
    E: dict[str, dict[str, float]] = {lv["label"]: {} for lv in levels(rx)}
    err: list[str] = []

    def lane(species: list[str]):
        try:
            for sp in species:
                opt[sp] = opt_one(rx, sp)
                _write_json(root / "opt.json", opt)
                if not opt[sp]["usable"]:
                    err.append(f"{sp}: optimization not usable; see opt.json")
                    return
            for sp in [s for s in species if s != "butadiene_gauche"]:
                for lv in sorted(levels(rx), key=lambda l: not l["label"].startswith("DLPNO")):
                    E[lv["label"]][sp] = sp_one(rx, lv, opt[sp])[0]
                    _write_json(root / "sp.json", E)
        except Exception as e:                      # recorded, never a result
            err.append(f"{species}: {type(e).__name__}: {e}")

    th = [threading.Thread(target=lane, args=(["ts"],)), threading.Thread(target=lane, args=(MINIMA,))]
    for t in th:
        t.start()
    for t in th:
        t.join()
    if err:
        (root / "STOPPED.md").write_text("Stopped:\n" + "\n".join(f"- {e}" for e in err) + "\n", encoding="utf-8")
        say("STOPPED: " + "; ".join(err))
        return
    analyse(rx)
    say("finished")


def _thermo(rx: Reaction, row: dict, T_: float) -> dict:
    hp = rx.root / row["dir"] / f"{row['job']}.hess"
    sigma = thermo.symmetry_number(hp.with_suffix(".out"))
    t = thermo.thermo(hp, T_, sigma)
    t["check_298_vs_orca_kJ"] = thermo.check_against_orca(hp, hp.with_suffix(".out"))
    return t


def analyse(rx: Reaction) -> dict:
    root = rx.root / SUB
    opt = _read_json(root / "opt.json")
    E = _read_json(root / "sp.json")
    g = rx.cfg["gas_phase"]
    A_exp_L = g["experimental_A_cm3_molecule_s"] * thermo.NA / 1000.0
    Ea_exp = g["experimental_Ea_kJ"]
    k_exp = A_exp_L * math.exp(-Ea_exp / (R_KJ * T))
    ref = ["butadiene_s-trans", "ethylene"]
    th = {sp: _thermo(rx, opt[sp], T) for sp in ["ts", *ref]}
    th298 = {sp: _thermo(rx, opt[sp], 298.15) for sp in ["ts", *ref]}
    d = lambda tab, key: tab["ts"][key] - sum(tab[s][key] for s in ref)
    dH_corr, dG_corr = d(th, "H_minus_E_kJ"), d(th, "G_minus_E_kJ")
    dG298_corr = d(th298, "G_minus_E_kJ")
    E_work = {sp: opt[sp]["E_eh"] for sp in opt}
    order = ["M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]
    tab = [{"label": "M06-2X_ma-def2-SVP (working)", "E": E_work}] + [
        {"label": k, "E": E[k]} for k in sorted(E, key=lambda k: order.index(k) if k in order else 99)]
    kT_h = thermo.KB * T / thermo.H_PL
    v0 = R_KJ * 1000.0 * T / thermo.P_ATM * 1000.0            # L/mol at 1 atm
    rows = []
    for t in tab:
        dE = kj(t["E"]["ts"] - sum(t["E"][s] for s in ref))
        dH = dE + dH_corr
        Ea = dH + 2 * R_KJ * T
        dG = dE + dG_corr
        k = kT_h * v0 * math.exp(-dG / (R_KJ * T))
        A = k * math.exp(Ea / (R_KJ * T))
        rows.append({"level": t["label"], "dE_kJ": dE, "dH_840_kJ": dH, "Ea_840_kJ": Ea,
                     "Ea_minus_exp_kJ": Ea - Ea_exp, "dG_840_1atm_kJ": dG, "dS_840_1atm_J_mol_K": (dH - dG) / T * 1000,
                     "dG_298_1atm_kJ": dE + dG298_corr,
                     "k_TST_L_mol_s": k, "k_ratio_to_exp": k / k_exp,
                     "A_L_mol_s": A, "A_cm3_mol_s": A * 1000.0, "A_cm3_molecule_s": A * 1000.0 / thermo.NA,
                     "A_ratio_to_exp": A / A_exp_L})
    gauche = opt.get("butadiene_gauche", {})
    conf = None
    if gauche.get("usable"):
        conf = {"dE_kJ": kj(gauche["E_eh"] - opt["butadiene_s-trans"]["E_eh"]),
                "dG298_kJ": kj(gauche["G298_eh"] - opt["butadiene_s-trans"]["G298_eh"]),
                "lowest": "butadiene_s-trans" if gauche["G298_eh"] > opt["butadiene_s-trans"]["G298_eh"] else "butadiene_gauche"}
    out = {"T_K": T, "experiment": {"source": g["experimental_source"], "Ea_kJ": Ea_exp,
                                    "Ea_err_kJ": g.get("experimental_Ea_err_kJ"), "A_L_mol_s": A_exp_L,
                                    "A_cm3_mol_s": A_exp_L * 1000.0,
                                    "A_cm3_molecule_s": g["experimental_A_cm3_molecule_s"], "k_840_L_mol_s": k_exp,
                                    "T_range_K": g.get("T_exp_range_K")},
           "reference": ref, "thermal": {sp: {k_: v for k_, v in th[sp].items()} for sp in th},
           "sigma_note": SIGMA_NOTE, "conformers": conf,
           "ts": opt["ts"]["mode"], "ts_imaginary_cm": opt["ts"]["imaginary_cm"],
           "delta_thermal_840_kJ": {"H_minus_E": dH_corr, "G_minus_E": dG_corr}, "rows": rows,
           "A_units": "A and k(TST) are second-order rate constants in concentration units: L mol⁻¹ s⁻¹ "
                      "(= 10³ cm³ mol⁻¹ s⁻¹, Rowley and Steiner's unit), also given per molecule "
                      "(cm³ molecule⁻¹ s⁻¹). ΔG‡ and ΔS‡ use the 1 atm standard state; the factor RT/p° turns that "
                      "into concentration units, so A and k do not depend on the choice of standard state and "
                      "compare directly with the measured A."}
    _write_json(root / "gas_phase.json", out)
    write_md(rx, out)
    return out


def _sci(x: float, digits: int = 2) -> str:
    m, e = f"{x:.{digits}e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def write_md(rx: Reaction, out: dict) -> None:
    x = out["experiment"]
    ts = out["ts"]
    L = ["Gas-phase comparison with Rowley and Steiner (J. Picione, 2026-11-20). Gas phase, no solvent model, gas-phase geometries and frequencies at M06-2X-D3(0)/ma-def2-SVP; "
         f"barrier from s-trans butadiene + ethylene; T = {T:.0f} K (midpoint of the measured "
         f"{x['T_range_K'][0]}–{x['T_range_K'][1]} K).", "",
         "| Level | ΔE‡ | ΔH‡ (840 K) | Ea = ΔH‡ + 2RT (840 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | "
         "k(TST, 840 K) (L mol⁻¹ s⁻¹) | k / measured k |", "|---|---|---|---|---|---|---|---|---|"]
    for r in out["rows"]:
        L.append(f"| {r['level']} | {fmt_kj(r['dE_kJ'])} | {fmt_kj(r['dH_840_kJ'])} | {fmt_kj(r['Ea_840_kJ'])} | "
                 f"{(('+' + fmt_kj(r['Ea_minus_exp_kJ']).replace('(', '(+')) if r['Ea_minus_exp_kJ'] > 0 else fmt_kj(r['Ea_minus_exp_kJ']).replace('-', '−'))} | {_sci(r['A_L_mol_s'])} | {r['A_ratio_to_exp']:.2f} | "
                 f"{r['k_TST_L_mol_s']:.3g} | {r['k_ratio_to_exp']:.2f} |")
    L.append(f"| Measured (Rowley and Steiner 1951) | | | {x['Ea_kJ']:.1f} ± {x['Ea_err_kJ']} kJ/mol ({x['Ea_kJ'] / 4.184:.1f} ± {x['Ea_err_kJ'] / 4.184:.1f} kcal/mol) | | "
             f"{_sci(x['A_L_mol_s'])} | 1 | {x['k_840_L_mol_s']:.3g} | 1 |")
    r0 = out["rows"][0]
    L += ["", "**Units and standard state of A.** " + out["A_units"] +
          f" Measured A = {_sci(x['A_cm3_mol_s'])} cm³ mol⁻¹ s⁻¹ = {_sci(x['A_L_mol_s'])} L mol⁻¹ s⁻¹ = "
          f"{_sci(x['A_cm3_molecule_s'])} cm³ molecule⁻¹ s⁻¹; computed A = {_sci(r0['A_cm3_mol_s'])} cm³ mol⁻¹ s⁻¹ = "
          f"{_sci(r0['A_cm3_molecule_s'])} cm³ molecule⁻¹ s⁻¹.",
          "",
          f"A is the same at every level: all levels share the working-level gas-phase frequencies, so ΔS‡ "
          f"(840 K, 1 atm) = {r0['dS_840_1atm_J_mol_K']:.1f} J mol⁻¹ K⁻¹".replace("= -", "= −") + f" at every level; the levels differ only in Ea.",
          "",
          f"Method: ΔH‡(840 K) = ΔE‡ + Δ[H(840 K) − E]; Ea = ΔH‡ + 2RT; ΔG‡(840 K, 1 atm) = ΔE‡ + Δ[G(840 K) − E] "
          f"(quasi-RRHO entropy); k(TST) = (kB T/h)(RT/p°) exp(−ΔG‡/RT), κ = 1; A = k(TST) exp(Ea/RT). "
          f"Thermal terms from the gas-phase working-level frequencies ({out['sigma_note']}: "
          + ", ".join(f"{sp} σ = {out['thermal'][sp]['sigma']}" for sp in out["thermal"]) + "). "
          "Same method as the earlier B3LYP-geometry folder, except that geometries and frequencies are now gas phase.",
          "",
          f"**Gas-phase TS:** one imaginary frequency, " + f"{ts['freq_cm']:.0f}".replace("-", "−") + f" cm⁻¹; overlap with the toluene TS mode "
          f"{ts['overlap_with_toluene_mode']:.2f}; forming bonds C1–C6 {ts['distances_A']['C1–C6']:.3f} Å and "
          f"C4–C5 {ts['distances_A']['C4–C5']:.3f} Å (spread {ts['forming_spread_A']:.3f} Å, stretch-rate ratio "
          f"{ts['forming_rate_ratio']:.2f}): " + ("synchronous, accepted by the folder's TS rules." if ts['accepted'] else "NOT accepted by the folder's TS rules.") + " Toluene TS: 2.263 Å."]
    c = out["conformers"]
    if c:
        L += ["", f"**Butadiene conformers (gas phase, working level):** the gauche s-cis form lies "
                  f"{fmt_kj(c['dE_kJ'])} in E and {fmt_kj(c['dG298_kJ'])} in G (298.15 K) above s-trans; "
                  f"the lowest conformer is {c['lowest'].replace('butadiene_', '').replace('gauche', 'gauche s-cis')}."]
    fnd = rx.root / SUB / "findings.md"
    if fnd.exists():
        L += ["", fnd.read_text(encoding="utf-8").strip()]
    (rx.root / SUB / "gas_phase.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


def main() -> None:
    rx = Reaction.load(Path(sys.argv[1]))
    if sys.argv[2] == "run":
        run(rx)
    elif sys.argv[2] == "md":
        analyse(rx)


if __name__ == "__main__":
    main()
