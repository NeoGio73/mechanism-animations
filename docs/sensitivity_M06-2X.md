# Functional sensitivity test: B3LYP vs M06-2X geometries (J. Picione, 2026-11-11)

Geometries: B3LYP-D3(BJ)/ma-def2-SVP (working level) and M06-2X-D3(0)/ma-def2-SVP (DefGrid3), same SMD solvent, Opt Freq; TS from the B3LYP TS. Energies: DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO) single points on each geometry set, each with its own level's thermal corrections. Solution reactions: ΔG‡(298.15 K) from separated reactants, 1 M. Diels–Alder: gas-phase Ea = ΔH‡(840 K) + 2RT from s-trans-butadiene + ethylene (gas single points on the SMD geometries, SMD frequencies; σ ts 1, s-trans 2, ethylene 4), as on the card. Every M06-2X TS passes the pipeline's TS acceptance test (one imaginary mode, key bonds moving the expected way); all M06-2X minima have only real frequencies.

## Transition states

| Reaction | Bond (role) | B3LYP | M06-2X | Imaginary freq. B3LYP / M06-2X | Mode overlap |
|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br (water) | C–O (forms) | 2.170 Å | 2.112 Å | 487i / 619i cm⁻¹ | 0.986 |
| | C–Br (breaks) | 2.372 Å | 2.323 Å | | |
| E2, EtO⁻ + 2-bromopropane (ethanol) | O–H (forms) | 1.392 Å | 1.294 Å | 888i / 1390i cm⁻¹ | 0.975 |
| | Cβ–H (breaks) | 1.279 Å | 1.348 Å | | |
| | Cα–Br (breaks) | 2.461 Å | 2.258 Å | | |
| Cyanide + acetone (water) | C(N)–C (forms) | 2.034 Å | 2.080 Å | 319i / 367i cm⁻¹ | 0.999 |
| Diels–Alder, butadiene + ethylene (gas Ea, 840 K) | C1–C6 (forms) | 2.265 Å | 2.263 Å | 516i / 512i cm⁻¹ | 0.998 |
| | C4–C5 (forms) | 2.265 Å | 2.263 Å | | |

## Coupled-cluster barrier on each geometry set vs measured

| Reaction | Measured | CC on B3LYP geometry | − measured | CC on M06-2X geometry | − measured |
|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br (water) | ΔG‡ 95.0 kJ/mol (22.7 kcal/mol) | 100.6 kJ/mol (24.0 kcal/mol) | +5.6 kJ/mol | 100.3 kJ/mol (24.0 kcal/mol) | +5.3 kJ/mol |
| E2, EtO⁻ + 2-bromopropane (ethanol) | ΔG‡ 105.1 kJ/mol (25.1 kcal/mol) | 107.2 kJ/mol (25.6 kcal/mol) | +2.1 kJ/mol | 99.8 kJ/mol (23.8 kcal/mol) | -5.3 kJ/mol |
| Cyanide + acetone (water) | ΔG‡ 69.6 kJ/mol (16.6 kcal/mol) | 79.5 kJ/mol (19.0 kcal/mol) | +9.9 kJ/mol | 79.3 kJ/mol (18.9 kcal/mol) | +9.7 kJ/mol |
| Diels–Alder, butadiene + ethylene (gas Ea, 840 K) | Ea (840 K, gas) 115.1 kJ/mol (27.5 kcal/mol) | 119.8 kJ/mol (28.6 kcal/mol) | +4.7 kJ/mol | 119.3 kJ/mol (28.5 kcal/mol) | +4.2 kJ/mol |

Diels–Alder rate at 840 K, k(TST)/k(exp): 0.29 on B3LYP geometries, 0.25 on M06-2X geometries.

Notes:
- SN2: the B3LYP-geometry value is for the ma-def2-SVP rerun TS (TS frame, separated reactants), computed here because the rerun lane is held before its barrier steps; the def2-SVP pilot card gives 95.6 kJ/mol.
- E2: the warning 'Cα–Br 2.26 Å outside the expected 2.3–2.6 Å' (reaction.toml [ts_expected]) is a range warning, not a rejection.
- Sources of measured values: SN2: water, 25 °C; Mabey and Mill 1978; E2: ethanol, 25 °C, E2 share; Dhar, Hughes and Ingold 1948; Cyanide + acetone (water): water, 25 °C; Schlesinger and Miller 1973; Diels–Alder: gas phase, 760–921 K.

Data: docs/sensitivity_M06-2X.json (pipeline/level_test.py sens_report).
