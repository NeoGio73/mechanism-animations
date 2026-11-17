"""Eyring fit of measured second-order rate constants (J. Picione, 2026-11-17).

SN2 HO⁻ + CH₃Br in water: E. A. Moelwyn-Hughes, Proc. R. Soc. Lond. A 1949, 196, 540–553,
doi:10.1098/rspa.1949.0044, Table 5 (NaOH in water, rates corrected for the competing reaction with water).
Least squares of ln(k/T) against 1/T: ln(k/T) = ln(kB/h) + ΔS‡/R − ΔH‡/(RT); k in L mol⁻¹ s⁻¹, so
ΔG‡ refers to the 1 M standard state. Usage: python -X utf8 -m pipeline.eyring_fit
"""
from __future__ import annotations

import numpy as np

KB, H, R = 1.380649e-23, 6.62607015e-34, 8.314462618
SN2_TABLE5 = [(292.69, 0.708e-4), (297.62, 1.37e-4), (303.41, 2.87e-4),
              (308.46, 5.33e-4), (313.29, 9.31e-4), (322.91, 29.1e-4)]   # (T / K, k2 / L mol⁻¹ s⁻¹)


def eyring(points: list[tuple[float, float]], T0: float = 298.15) -> dict:
    T, k = np.array(points).T
    X = np.vstack([1 / T, np.ones_like(T)]).T
    y = np.log(k / T)
    (m, b), res, *_ = np.linalg.lstsq(X, y, rcond=None)
    cov = res[0] / (len(T) - 2) * np.linalg.inv(X.T @ X)
    g = np.array([-R, -T0 * R])                       # dΔG‡/d(m, b)
    dH, dS = -m * R, (b - np.log(KB / H)) * R
    return {"dH_kJ": dH / 1e3, "dH_sd_kJ": np.sqrt(cov[0, 0]) * R / 1e3,
            "dS_J_K": dS, "dS_sd_J_K": np.sqrt(cov[1, 1]) * R,
            "dG_kJ": (dH - T0 * dS) / 1e3, "dG_sd_kJ": np.sqrt(g @ cov @ g) / 1e3, "T0": T0,
            "T_range_K": (T.min(), T.max())}


if __name__ == "__main__":
    f = eyring(SN2_TABLE5)
    c = 4.184
    print(f"ΔH‡ = {f['dH_kJ']:.1f} ± {f['dH_sd_kJ']:.1f} kJ/mol ({f['dH_kJ'] / c:.1f} kcal/mol)")
    print(f"ΔS‡ = {f['dS_J_K']:.1f} ± {f['dS_sd_J_K']:.1f} J mol⁻¹ K⁻¹ ({f['dS_J_K'] / c:.1f} cal mol⁻¹ K⁻¹)")
    print(f"ΔG‡({f['T0']} K, 1 M) = {f['dG_kJ']:.2f} ± {f['dG_sd_kJ']:.2f} kJ/mol ({f['dG_kJ'] / c:.1f} kcal/mol); "
          f"measured range {f['T_range_K'][0]}–{f['T_range_K'][1]} K")
