# SMD-scheme check (stationary points, all six reactions)

Scheme A (as in the paper): DLPNO-CCSD(T)/ma-def2-TZVP with SMD. Scheme B: gas-phase DLPNO-CCSD(T)/ma-def2-TZVP + [E(SMD) − E(gas)] at M06-2X/ma-def2-TZVP, at the same geometries. Same working-level thermal corrections, reference species and 1 M standard state as the paper's barriers (pipeline/results_summary). kJ/mol (kcal/mol), 298.15 K.

**Largest difference: 6.58 kJ/mol (1.57 kcal/mol). STOP: over the 2 kJ/mol limit; no table changed.**

| Reaction | ΔG‡ A | ΔG‡ B | B − A | ΔrG A | ΔrG B | B − A |
|---|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br (water) | 100.4 (24.0) | 103.4 (24.7) | +3.01 (+0.72) | -85.3 (-20.4) | -78.8 (-18.8) | +6.58 (+1.57) |
| E2, EtO⁻ + 2-bromopropane (ethanol) | 99.6 (23.8) | 102.8 (24.6) | +3.16 (+0.75) | -108.3 (-25.9) | -104.3 (-24.9) | +3.99 (+0.95) |
| CN⁻ + acetone (water) | 79.3 (19.0) | 78.2 (18.7) | -1.16 (-0.28) | 36.4 (8.7) | 36.3 (8.7) | -0.12 (-0.03) |
| Diels–Alder, butadiene + ethylene (SMD as in the paper) | 144.4 (34.5) | 144.8 (34.6) | +0.48 (+0.11) | -123.3 (-29.5) | -122.6 (-29.3) | +0.63 (+0.15) |
| 1,2-H shift, 3-methyl-2-butyl cation (DCM) (path shown, from the eclipsed form) | 11.1 (2.7) | 11.5 (2.8) | +0.40 (+0.10) | -45.4 (-10.9) | -45.2 (-10.8) | +0.29 (+0.07) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation (DCM) | 17.0 (4.1) | 17.4 (4.1) | +0.33 (+0.08) | 0.0 (0.0) | 0.0 (0.0) | +0.00 (+0.00) |

Reference species (zero) per reaction: SN2: hydroxide + bromomethane → methanol + bromide; E2: ethoxide + 2-bromopropane → propene + ethanol + bromide; CN⁻ + acetone (water): cyanide + acetone → acetone_cyanohydrin_alkoxide; Diels–Alder: butadiene_s-trans + ethylene → cyclohexene; 1: secondary_anti → product_complex; 1: dimethylbutyl_cation → dimethylbutyl_cation_shifted

## Decision (J. Picione, 2026-11-26)

Scheme A stays in every table, for consistency with the paths, cards and videos, which are all on scheme A. The B − A columns go into the SI as the uncertainty from the solvation scheme. Scheme B is not uniformly closer to experiment:

| Reaction | Quantity | Measured | A − measured | B − measured | B moves it |
|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br (water) | ΔG‡ | 94.9 (22.7) | +5.5 (+1.3) | +8.5 (+2.0) | away from the measured value, by 3.0 kJ/mol |
| SN2, HO⁻ + CH₃Br (water) | ΔrG | -100.8 (-24.1) | +15.5 (+3.7) | +22.0 (+5.3) | away from the measured value, by 6.6 kJ/mol |
| E2, EtO⁻ + 2-bromopropane (ethanol) | ΔG‡ | 105.1 (25.1) | -5.5 (-1.3) | -2.3 (-0.6) | toward the measured value, by 3.2 kJ/mol |
| E2, EtO⁻ + 2-bromopropane (ethanol) | ΔrG | not measured | – | – | – |
| CN⁻ + acetone (water) | ΔG‡ | 69.6 (16.6) | +9.7 (+2.3) | +8.6 (+2.1) | toward the measured value, by 1.2 kJ/mol |
| CN⁻ + acetone (water) | ΔrG | 10.5 (2.5) | +25.9 (+6.2) | +25.8 (+6.2) | toward the measured value, by 0.1 kJ/mol |

## Per species

B − A for each structure: the M06-2X/ma-def2-TZVP solvation energy, E(SMD) − E(gas), minus the DLPNO-CCSD(T)/ma-def2-TZVP one, at the same geometry (from the jobs above; no new calculation). Negative: the species lies lower in scheme B. The barrier difference is the TS value minus the reactants' values; the reaction difference is the products' values minus the reactants' values. kJ/mol (kcal/mol).

| Reaction | Species | Role | Solvation, M06-2X | Solvation, DLPNO-CCSD(T) | B − A |
|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br | hydroxide | reactant | -400.5 (-95.7) | -394.0 (-94.2) | -6.55 (-1.57) |
| SN2, HO⁻ + CH₃Br | bromomethane | reactant | -4.9 (-1.2) | -4.2 (-1.0) | -0.78 (-0.19) |
| SN2, HO⁻ + CH₃Br | ts | transition state | -254.3 (-60.8) | -250.0 (-59.8) | -4.33 (-1.03) |
| SN2, HO⁻ + CH₃Br | methanol | product | -19.8 (-4.7) | -19.1 (-4.6) | -0.76 (-0.18) |
| SN2, HO⁻ + CH₃Br | bromide | product | -222.9 (-53.3) | -222.9 (-53.3) | +0.00 (+0.00) |
| SN2, HO⁻ + CH₃Br | reactant_complex | not in ΔG‡ or ΔrG | -365.2 (-87.3) | -360.7 (-86.2) | -4.56 (-1.09) |
| SN2, HO⁻ + CH₃Br | product_complex | not in ΔG‡ or ΔrG | -225.2 (-53.8) | -224.7 (-53.7) | -0.49 (-0.12) |
| E2, EtO⁻ + 2-bromopropane | ethoxide | reactant | -318.0 (-76.0) | -313.1 (-74.8) | -4.87 (-1.16) |
| E2, EtO⁻ + 2-bromopropane | 2-bromopropane | reactant | -22.3 (-5.3) | -21.3 (-5.1) | -0.97 (-0.23) |
| E2, EtO⁻ + 2-bromopropane | ts | transition state | -216.7 (-51.8) | -214.0 (-51.1) | -2.69 (-0.64) |
| E2, EtO⁻ + 2-bromopropane | propene | product | -8.3 (-2.0) | -7.2 (-1.7) | -1.12 (-0.27) |
| E2, EtO⁻ + 2-bromopropane | ethanol | product | -27.5 (-6.6) | -26.8 (-6.4) | -0.74 (-0.18) |
| E2, EtO⁻ + 2-bromopropane | bromide | product | -224.4 (-53.6) | -224.4 (-53.6) | +0.01 (+0.00) |
| E2, EtO⁻ + 2-bromopropane | reactant_complex | not in ΔG‡ or ΔrG | -290.5 (-69.4) | -286.0 (-68.4) | -4.43 (-1.06) |
| E2, EtO⁻ + 2-bromopropane | product_complex | not in ΔG‡ or ΔrG | -230.3 (-55.0) | -229.8 (-54.9) | -0.48 (-0.12) |
| E2, EtO⁻ + 2-bromopropane | ethanol_trans | not in ΔG‡ or ΔrG | -27.4 (-6.5) | -26.4 (-6.3) | -0.97 (-0.23) |
| CN⁻ + acetone | cyanide | reactant | -266.9 (-63.8) | -265.7 (-63.5) | -1.22 (-0.29) |
| CN⁻ + acetone | acetone | reactant | -22.0 (-5.3) | -17.8 (-4.3) | -4.19 (-1.00) |
| CN⁻ + acetone | ts | transition state | -251.3 (-60.1) | -244.8 (-58.5) | -6.57 (-1.57) |
| CN⁻ + acetone | acetone_cyanohydrin_alkoxide | product | -285.8 (-68.3) | -280.3 (-67.0) | -5.54 (-1.32) |
| CN⁻ + acetone | reactant_complex | not in ΔG‡ or ΔrG | -280.2 (-67.0) | -272.9 (-65.2) | -7.30 (-1.74) |
| CN⁻ + acetone | product_complex | not in ΔG‡ or ΔrG | -285.4 (-68.2) | -279.9 (-66.9) | -5.53 (-1.32) |
| Diels–Alder, butadiene + ethylene | butadiene_s-trans | reactant | -12.4 (-3.0) | -11.6 (-2.8) | -0.83 (-0.20) |
| Diels–Alder, butadiene + ethylene | ethylene | reactant | -4.9 (-1.2) | -4.4 (-1.1) | -0.51 (-0.12) |
| Diels–Alder, butadiene + ethylene | ts | transition state | -18.6 (-4.4) | -17.7 (-4.2) | -0.86 (-0.21) |
| Diels–Alder, butadiene + ethylene | cyclohexene | product | -21.1 (-5.1) | -20.4 (-4.9) | -0.71 (-0.17) |
| Diels–Alder, butadiene + ethylene | reactant_complex | not in ΔG‡ or ΔrG | -17.8 (-4.2) | -16.3 (-3.9) | -1.50 (-0.36) |
| Diels–Alder, butadiene + ethylene | product_complex | not in ΔG‡ or ΔrG | -20.9 (-5.0) | -20.2 (-4.8) | -0.69 (-0.16) |
| Diels–Alder, butadiene + ethylene | butadiene_s-cis | not in ΔG‡ or ΔrG | -12.4 (-3.0) | -11.5 (-2.8) | -0.84 (-0.20) |
| 1,2-H shift, 3-methyl-2-butyl cation | secondary_anti | reactant | -246.0 (-58.8) | -246.6 (-58.9) | +0.56 (+0.13) |
| 1,2-H shift, 3-methyl-2-butyl cation | ts | transition state | -249.0 (-59.5) | -249.9 (-59.7) | +0.86 (+0.21) |
| 1,2-H shift, 3-methyl-2-butyl cation | product_complex | product | -245.0 (-58.6) | -245.8 (-58.8) | +0.85 (+0.20) |
| 1,2-H shift, 3-methyl-2-butyl cation | reactant_complex | not in ΔG‡ or ΔrG | -245.7 (-58.7) | -246.2 (-58.8) | +0.46 (+0.11) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | dimethylbutyl_cation | reactant | -236.6 (-56.6) | -237.6 (-56.8) | +0.93 (+0.22) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | ts | transition state | -237.2 (-56.7) | -238.5 (-57.0) | +1.26 (+0.30) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | dimethylbutyl_cation_shifted | product | -236.6 (-56.6) | -237.6 (-56.8) | +0.93 (+0.22) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | reactant_complex | not in ΔG‡ or ΔrG | -236.8 (-56.6) | -237.8 (-56.8) | +0.94 (+0.22) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | product_complex | not in ΔG‡ or ΔrG | -236.8 (-56.6) | -237.8 (-56.8) | +0.94 (+0.22) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | rotamer_h_perp | not in ΔG‡ or ΔrG | -236.6 (-56.5) | -237.5 (-56.8) | +0.93 (+0.22) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation | rotamer_h_plane | not in ΔG‡ or ΔrG | -236.7 (-56.6) | -237.6 (-56.8) | +0.95 (+0.23) |

Code: `pipeline/smd_scheme.py` (jobs), `pipeline/smd_scheme_report.py` (this table). Gas-phase jobs: `reactions/<reaction>/01_ts/smd_scheme/gas/`.
