# Tertiary hydride shift: lowest-free-energy rotamer of the cation (J. Picione, 2026-11-19)

kJ/mol (kcal/mol). Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane, quasi-RRHO. 298.15 K: ORCA's G; 135 K: pipeline/thermo.py (same model, checked against ORCA at 298.15 K).

| Structure | Me–C2–C3–H | C3–H to p axis | C3–CH₃ to p axis | point group | G vs lowest, 298.15 K | G vs lowest, 135 K | thermo check vs ORCA (G − E) |
|---|---|---|---|---|---|---|---|
| `open_cation` (path-end minimum (frame 0 optimized)) | 150.2° | 61.6° | 14.3°, 69.9° | C1 | +1.19 (+0.28) | +0.94 (+0.22) | -0.001 kJ/mol |
| `dimethylbutyl_cation` (separately optimized cation (current zero)) | 150.5° | 61.9° | 13.9°, 69.9° | C1 | +0.00 (+0.00) | +0.00 (+0.00) | -0.003 kJ/mol |
| `rotamer_h_perp` (start: C3–H along the empty p orbital) | 150.4° | 61.9° | 14.5°, 69.2° | C1 | +0.64 (+0.15) | +0.50 (+0.12) | -0.002 kJ/mol |
| `rotamer_h_plane` (start: C3–H in the cation plane) | 149.9° | 61.3° | 14.7°, 70.0° | C1 | +0.40 (+0.09) | +0.29 (+0.07) | -0.002 kJ/mol |

**298.15 K: lowest = `dimethylbutyl_cation`.**

| Level | Frame 0 above it (E) | Path-end minimum above it (E) | ΔE‡ | ΔG‡ |
|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP | +1.27 (+0.30) | -0.20 (-0.05) | +15.3 (+3.7) | +14.1 (+3.4) |
| M06-2X/ma-def2-TZVP | +1.86 (+0.44) | -0.19 (-0.04) | +19.2 (+4.6) | +18.0 (+4.3) |
| B3LYP/ma-def2-TZVP | +1.67 (+0.40) | -0.23 (-0.05) | +28.5 (+6.8) | +27.3 (+6.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +2.22 (+0.53) | -0.18 (-0.04) | +18.2 (+4.3) | +17.0 (+4.1) |

**135 K: lowest = `dimethylbutyl_cation`.**

| Level | Frame 0 above it (E) | Path-end minimum above it (E) | ΔE‡ | ΔG‡ |
|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP | +1.27 (+0.30) | -0.20 (-0.05) | +15.3 (+3.7) | +13.2 (+3.2) |
| M06-2X/ma-def2-TZVP | +1.86 (+0.44) | -0.19 (-0.04) | +19.2 (+4.6) | +17.1 (+4.1) |
| B3LYP/ma-def2-TZVP | +1.67 (+0.40) | -0.23 (-0.05) | +28.5 (+6.8) | +26.4 (+6.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +2.22 (+0.53) | -0.18 (-0.04) | +18.2 (+4.3) | +16.1 (+3.8) |

Measured (Saunders and Kates 1978, 135 K, superacid): ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol).
