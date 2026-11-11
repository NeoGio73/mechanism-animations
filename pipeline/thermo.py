"""Ideal-gas thermochemistry at any temperature from an ORCA .hess file (no new ORCA job).

Same model as ORCA's own thermochemistry block, so the 298.15 K values can be checked against
ORCA's output (`check_against_orca`):
  - translation: Sackur–Tetrode at 1 atm;
  - rotation: rigid rotor with the point-group symmetry number ORCA reports;
  - vibration: harmonic energies (ZPE + thermal); entropy by Grimme's quasi-RRHO
    (Chem. Eur. J. 2012, 18, 9955): each mode's entropy is interpolated between the harmonic
    value and a free-rotor value, weight 1/(1 + (ω0/ν)⁴), ω0 = 100 cm⁻¹, Bav = 1e-44 kg m²;
  - the imaginary mode of a transition state is left out; electronic entropy zero (singlets).
All results in kJ/mol, per mole of molecules.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from . import orca

H_PL = 6.62607015e-34       # J s
KB = 1.380649e-23           # J/K
C_CM = 2.99792458e10        # cm/s
NA = 6.02214076e23
R_KJ = KB * NA / 1000.0     # kJ/(mol K)
AMU = 1.66053906660e-27     # kg
BOHR = 0.529177210903e-10   # m
P_ATM = 101325.0            # Pa
EH_KJ = 2625.499639
OMEGA0_CM = 100.0
BAV = 1.0e-44               # kg m²


def _atoms_bohr(path: Path):
    lines = Path(path).read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.strip() == "$atoms")
    n = int(lines[i + 1])
    rows = [lines[i + 2 + k].split() for k in range(n)]
    return np.array([float(r[1]) for r in rows]), np.array([[float(v) for v in r[2:5]] for r in rows]) * BOHR


def symmetry_number(out_path: Path) -> int:
    m = re.findall(r"Symmetry Number:\s+(\d+)", Path(out_path).read_text(errors="replace"))
    if not m:
        raise ValueError(f"no symmetry number in {out_path}")
    return int(m[-1])


def thermo(hess_path: Path, T: float, sigma: int, p_pa: float = P_ATM) -> dict:
    """H − E, T·S terms and G − E (kJ/mol) at temperature T and pressure p (default 1 atm)."""
    h = orca.read_hess(hess_path)
    masses, xyz = _atoms_bohr(hess_path)
    nat = len(masses)
    m_kg = masses * AMU
    com = (m_kg[:, None] * xyz).sum(0) / m_kg.sum()
    r = xyz - com
    inertia = np.zeros((3, 3))
    for mi, ri in zip(m_kg, r):
        inertia += mi * (np.dot(ri, ri) * np.eye(3) - np.outer(ri, ri))
    I = np.sort(np.linalg.eigvalsh(inertia))
    linear = nat == 2 or I[0] < 1e-3 * I[2]

    f = h["freqs"]
    nzero = 5 if linear else 6
    nvib = 3 * nat - nzero - int((f < 0).sum())
    nu = np.sort(f[f > 0])[-nvib:]
    RT = R_KJ * T

    # energies
    e = nu * H_PL * C_CM * NA / 1000.0                     # kJ/mol per quantum
    x = e / RT
    zpe = 0.5 * e.sum()
    evib = (e / np.expm1(x)).sum()
    H_minus_E = zpe + evib + (1.5 + (1.0 if linear else 1.5) + 1.0) * RT

    # entropies (J/(mol K) then kJ)
    M = m_kg.sum()
    s_trans = R_KJ * (np.log((2 * np.pi * M * KB * T / H_PL ** 2) ** 1.5 * KB * T / p_pa) + 2.5)
    theta = H_PL ** 2 / (8 * np.pi ** 2 * KB * I[(1 if linear else 0):])
    if linear:
        s_rot = R_KJ * (np.log(T / (sigma * theta[-1])) + 1.0)
    else:
        s_rot = R_KJ * (np.log(np.sqrt(np.pi) / sigma * np.sqrt(T ** 3 / np.prod(theta))) + 1.5)
    s_harm = R_KJ * (x / np.expm1(x) - np.log(-np.expm1(-x)))
    mu = H_PL / (8 * np.pi ** 2 * nu * C_CM)
    mu_eff = mu * BAV / (mu + BAV)
    s_free = R_KJ * (0.5 + np.log(np.sqrt(8 * np.pi ** 3 * mu_eff * KB * T / H_PL ** 2)))
    w = 1.0 / (1.0 + (OMEGA0_CM / nu) ** 4)
    s_vib = (w * s_harm + (1 - w) * s_free).sum()
    S = s_trans + s_rot + s_vib
    return {"T_K": T, "H_minus_E_kJ": H_minus_E, "TS_trans_kJ": T * s_trans, "TS_rot_kJ": T * s_rot,
            "TS_vib_kJ": T * s_vib, "G_minus_E_kJ": H_minus_E - T * S, "sigma": sigma, "linear": bool(linear)}


def orca_block(out_path: Path) -> dict:
    """ORCA's last thermochemistry block (kJ/mol): H − E, T·S terms, G − E."""
    t = Path(out_path).read_text(errors="replace")
    g = lambda pat: float(re.findall(pat + r"\s+\.\.\.\s+(-?\d+\.\d+)", t)[-1])
    E = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
    return {"H_minus_E_kJ": (g("Total Enthalpy") - E) * EH_KJ,
            "TS_trans_kJ": g("Translational entropy") * EH_KJ,
            "TS_rot_kJ": g("Rotational entropy") * EH_KJ,
            "TS_vib_kJ": g("Vibrational entropy") * EH_KJ,
            "G_minus_E_kJ": (g("Final Gibbs free energy") - E) * EH_KJ}


def check_against_orca(hess_path: Path, out_path: Path) -> dict:
    """Differences (this module − ORCA) at 298.15 K, kJ/mol."""
    mine = thermo(hess_path, 298.15, symmetry_number(out_path))
    ref = orca_block(out_path)
    return {k: round(mine[k] - ref[k], 3) for k in ref}
