"""Measured rate parameters converted to activation quantities (J. Picione, 2026-11-27; first use: the reverse
direction of the secondary hydride shift, Saunders and Hagen 1968).

From Arrhenius parameters (Ea, log10 A, A in s⁻¹ for a first-order step) at temperature T:
  k = A exp(−Ea/RT);  ΔG‡ = −RT ln(k h / kB T)  (Eyring, transmission coefficient 1);
  ΔH‡ = Ea − RT  (first order, condensed phase);  ΔS‡ = (ΔH‡ − ΔG‡)/T = R (ln A − ln(kB T/h) − 1).
The uncertainty of ΔS‡ is R ln10 × σ(log10 A); that of ΔH‡ is σ(Ea). Ea and log A from one fit are correlated, so
no uncertainty is given for k or ΔG‡ here.
"""
from __future__ import annotations

import math

R = 8.314462618          # J/(mol K)
KB = 1.380649e-23        # J/K
H = 6.62607015e-34       # J s
KCAL = 4.184


def arrhenius_to_activation(Ea_kcal: float, log10_A: float, T: float = 298.15, Ea_err_kcal: float | None = None,
                            log10_A_err: float | None = None) -> dict:
    Ea = Ea_kcal * KCAL * 1000.0
    A = 10.0 ** log10_A
    k = A * math.exp(-Ea / (R * T))
    dG = -R * T * math.log(k * H / (KB * T))
    dH = Ea - R * T
    dS = (dH - dG) / T
    return {"T_K": T, "k_s": k, "dG_kJ": dG / 1000.0, "dH_kJ": dH / 1000.0, "dS_J_molK": dS,
            "dH_err_kJ": Ea_err_kcal * KCAL if Ea_err_kcal is not None else None,
            "dS_err_J_molK": R * math.log(10.0) * log10_A_err if log10_A_err is not None else None,
            "Ea_kJ": Ea / 1000.0}


def eyring_dG(k: float, T: float = 298.15) -> float:
    """ΔG‡ in kJ/mol from a rate constant (s⁻¹, or L mol⁻¹ s⁻¹ for a second-order step: 1 M standard state)."""
    return -R * T * math.log(k * H / (KB * T)) / 1000.0


if __name__ == "__main__":
    print(arrhenius_to_activation(15.3, 13.2, 298.15, 0.2, 0.2))
