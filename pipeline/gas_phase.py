"""Gas-phase stationary-point energies and an Arrhenius activation energy for comparison with a
gas-phase experiment (step `gas_phase`, run with the barrier single points).

reaction.toml:

    [gas_phase]
    species = ["ts", "butadiene_s-cis", "butadiene_s-trans", "ethylene", "cyclohexene"]
    reference = ["butadiene_s-trans", "ethylene"]   # reactants the measurement refers to
    T_exp_K = 800                                   # mean temperature of the measurements
    T_exp_provisional = true                        # until the range is taken from the paper
    experimental_Ea_kJ = 115
    experimental_source = "Rowley and Steiner 1951"

Every species is recomputed without the solvent model at the working level, the protocol level
(`sp_method`) and each `[[diagnostic_levels]]` entry, on the solution (working-level) geometries.
Thermal corrections come from the working-level frequencies (SMD solvent), as for all barriers:
H(T) = E + ZPE + vibrational thermal energy (harmonic) + 4RT (3RT/2 translation, 3RT/2 rotation,
RT for pV; non-linear molecules), checked against ORCA's own enthalpy at 298.15 K. ΔG‡(298.15 K)
uses ORCA's free energies (quasi-RRHO entropy), 1 atm standard state in the gas phase and 1 M in
solution. For a gas-phase bimolecular reaction, Ea = ΔH‡(T) + 2RT (transition-state theory,
Δn‡ = −1).
"""
from __future__ import annotations

import json
import math
import re

import numpy as np

from . import orca
from .geometry import formal_charge, read_xyz
from .steps import Reaction, _read_json, _stationary_points, _write_json
from .units import fmt_kj, kj

R_KJ = 8.314462618e-3            # kJ/(mol K)
H_CM_TO_KJ = 0.01196266          # kJ/mol per cm⁻¹ (N_A h c)
STD_1M_KJ = lambda T: R_KJ * T * math.log(0.0820573661 * T)   # 1 atm → 1 M, per species


def _hess_path(rx: Reaction, name: str):
    if name == "ts":
        return rx.root / "01_ts/ts/ts.hess"
    job = _read_json(rx.root / "01_ts/fragments.json").get(name, {}).get("job", name)
    return rx.root / f"01_ts/fragments/{name}/{job}.hess"


def thermal_H_kj(rx: Reaction, name: str, T: float) -> float:
    """H(T) − E in kJ/mol from the working-level harmonic frequencies (imaginary mode left out)."""
    h = orca.read_hess(_hess_path(rx, name))
    f = h["freqs"]
    nat = len(h["symbols"])
    xyz = read_xyz(_hess_path(rx, name).with_suffix(".xyz")).coords
    c = xyz - xyz.mean(0)
    linear = nat == 2 or np.linalg.matrix_rank(c, tol=1e-3) <= 1
    nzero = 5 if linear else 6
    real = np.sort(f[f > 0])[-(3 * nat - nzero - int((f < 0).sum())):] if nat > 1 else np.array([])
    e = real * H_CM_TO_KJ                         # quantum energies, kJ/mol
    zpe = 0.5 * e.sum()
    x = e / (R_KJ * T)
    evib = (e / np.expm1(x)).sum()
    trans_rot_pv = (1.5 + (1.0 if linear else 1.5) + 1.0) * R_KJ * T
    return zpe + evib + trans_rot_pv


def _orca_thermo(rx: Reaction, name: str) -> dict:
    """ORCA's last thermochemistry block (298.15 K): E, H, G in Eh."""
    out = (_hess_path(rx, name).with_suffix(".out")).read_text(errors="replace")
    E = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1])
    H = float(re.findall(r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)", out)[-1])
    G = float(re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", out)[-1])
    return {"E": E, "H": H, "G": G}


def _levels(rx: Reaction) -> list[dict]:
    lv = [{"label": rx.level().replace("/", "_") + " (working)", "method": rx.cfg["method"], "key": "working"},
          {"label": rx.level(sp=True).replace("/", "_") + " (protocol)", "method": rx.cfg["sp_method"], "key": "protocol"}]
    lv += [{"label": d["label"], "method": d["method"], "key": d["label"]} for d in rx.cfg.get("diagnostic_levels", [])]
    return lv


def _solution_energies(rx: Reaction, key: str) -> dict:
    """Solution-phase single-point energies (Eh) per species at one level, from earlier steps."""
    if key == "working":
        frag = _read_json(rx.root / "01_ts/fragments.json")
        E = {n: v["E_eh"] for n, v in frag.items()}
        E["ts"] = _read_json(rx.root / "01_ts/ts.json")["E_eh"]
        return E
    if key == "protocol":
        return _read_json(rx.root / "01_ts/sp_tzvp.json")
    return _read_json(rx.root / "01_ts/diagnostics.json")["energies_eh"][key]


def gas_phase(rx: Reaction) -> None:
    cfg = rx.cfg["gas_phase"]
    names = cfg["species"]
    pts = _stationary_points(rx)
    T = float(cfg.get("T_exp_K", 800.0))
    ref = cfg["reference"]
    reac = [s["name"] for role, s in rx.species() if role == "reactants"]

    # thermal corrections (kJ/mol) from the working-level frequencies
    Hcorr_T = {n: thermal_H_kj(rx, n, T) for n in names}
    Hcorr_298 = {n: thermal_H_kj(rx, n, 298.15) for n in names}
    check = {}
    for n in names:
        o = _orca_thermo(rx, n)
        check[n] = round(Hcorr_298[n] - kj(o["H"] - o["E"]), 3)
        Hcorr_298[n] = kj(o["H"] - o["E"])          # use ORCA's own value at 298.15 K
    Gcorr_298 = {n: kj(_orca_thermo(rx, n)["G"] - _orca_thermo(rx, n)["E"]) for n in names}
    print("H(298.15 K) from frequencies minus ORCA's, kJ/mol:", check)

    rows = []
    for lv in _levels(rx):
        E_sol = _solution_energies(rx, lv["key"])
        E_gas = {}
        for n in names:
            xyz, q, m = pts[n]
            wd = rx.d(f"01_ts/gas_phase/{lv['key']}/{n}")
            inp = orca.make_input(lv["method"], read_xyz(xyz), q, m, nprocs=rx.cfg["nprocs"],
                                  maxcore_mb=rx.cfg["maxcore_mb"])
            r = orca.run(wd, n, inp, rx.log, note=f"gas-phase single point, {lv['label']}")
            if not r.ok:
                raise orca.OrcaError(f"gas-phase single point {lv['key']} {n} failed; see {r.out}")
            E_gas[n] = r.energy_eh
        for phase, E in (("SMD " + rx.cfg["solvent"], E_sol), ("gas phase", E_gas)):
            def barrier(reactants, corr=None):
                d = kj(E["ts"] - sum(E[n] for n in reactants))
                if corr:
                    d += corr["ts"] - sum(corr[n] for n in reactants)
                return d
            dG = barrier(ref, Gcorr_298) - (STD_1M_KJ(298.15) if phase != "gas phase" else 0.0)
            dH_T = barrier(ref, Hcorr_T)
            rows.append({
                "level": lv["label"], "phase": phase,
                "dE_from_reactants_kJ": barrier(reac),
                "dE_from_reference_kJ": barrier(ref),
                "dH_298_kJ": barrier(ref, Hcorr_298),
                "dG_298_kJ": dG,
                "standard_state": "1 atm" if phase == "gas phase" else "1 M",
                "dH_T_kJ": dH_T,
                "Ea_T_kJ": dH_T + 2 * R_KJ * T,
                "s_trans_to_s_cis_dE_kJ": kj(E[reac[0]] - E[ref[0]]) if reac[0] in E and ref[0] in E else None,
            })

    out = {"T_K": T, "T_provisional": bool(cfg.get("T_exp_provisional", True)),
           "reference": ref, "reactants": reac,
           "experimental_Ea_kJ": cfg.get("experimental_Ea_kJ"),
           "thermal_note": "working-level (SMD) harmonic frequencies; H(298.15 K) and G(298.15 K) from ORCA "
                           "(G with quasi-RRHO entropy); H(T) from the same frequencies, harmonic",
           "H298_check_kJ": check, "rows": rows}
    _write_json(rx.root / "01_ts/gas_phase.json", out)

    tag = f"{T:.0f} K" + (" (provisional)" if out["T_provisional"] else "")
    ref_txt = " + ".join(ref)
    md = [f"Barriers relative to {ref_txt} unless marked; solution geometries (SMD {rx.cfg['solvent']}, "
          f"{rx.basis()}); thermal corrections from the working-level frequencies.", "",
          f"| Level | Phase | ΔE‡ from {' + '.join(reac)} | ΔE‡ | ΔH‡ (298.15 K) | ΔG‡ (298.15 K) | "
          f"ΔH‡ ({tag}) | Ea = ΔH‡ + 2RT ({tag}) |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['level']} | {r['phase']} | {fmt_kj(r['dE_from_reactants_kJ'])} | "
                  f"{fmt_kj(r['dE_from_reference_kJ'])} | {fmt_kj(r['dH_298_kJ'])} | "
                  f"{fmt_kj(r['dG_298_kJ'])} ({r['standard_state']}) | {fmt_kj(r['dH_T_kJ'])} | "
                  f"{fmt_kj(r['Ea_T_kJ'])} |")
    if out["experimental_Ea_kJ"]:
        md += ["", f"Experiment ({cfg.get('experimental_source', '')}, gas phase): "
                   f"Ea ≈ {fmt_kj(out['experimental_Ea_kJ'])}."]
    (rx.root / "01_ts/gas_phase.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
