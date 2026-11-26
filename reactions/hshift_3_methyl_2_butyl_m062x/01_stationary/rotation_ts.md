# Rotation TSs from the anti minimum, secondary hydride shift (J. Picione, 2026-11-22)

Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane; G at 298.15 K (ORCA quasi-RRHO, working-level thermal correction); zero = secondary_anti. kJ/mol (kcal/mol). Written by `00_inputs/rotation_ts_m062x.py` from finished jobs.

## ts_anti_rotation
- Start: gas-phase M06-2X anti TS
- TS job: {'ok': True, 'opt_converged': True, 'imaginary_cm': [-162.44023604481774]}
- TS geometry: {'C3_H_A': 1.113, 'C2_H_A': 2.02, 'C1_C2_C3_H_deg': -55.4, 'H_on': 'C3 (secondary)'}; core RMSD to the production TS 0.509 Å
- plus side: reached secondary_anti; {'C3_H_A': 1.094, 'C2_H_A': 2.124, 'C1_C2_C3_H_deg': -1.3, 'H_on': 'C3 (secondary)'}; secondary_anti RMSD 0.001 Å, ΔE +0.00 kJ/mol; secondary_eclip_irc RMSD 1.084 Å, ΔE -3.05 kJ/mol; tertiary_irc RMSD 0.784 Å, ΔE +39.39 kJ/mol
- minus side: reached tertiary_irc; {'C3_H_A': 1.995, 'C2_H_A': 1.116, 'C1_C2_C3_H_deg': -122.4, 'H_on': 'C2 (tertiary)'}; secondary_anti RMSD 0.783 Å, ΔE -39.50 kJ/mol; secondary_eclip_irc RMSD 0.948 Å, ΔE -42.55 kJ/mol; tertiary_irc RMSD 0.007 Å, ΔE -0.11 kJ/mol
- Connects ['secondary_anti', 'tertiary_irc']; wanted ['secondary_anti', 'tertiary_irc']: **confirmed**

## ts_anti_eclip_pos
- Start: highest point 15 of scan arm pos
- TS job: {'ok': True, 'opt_converged': True, 'imaginary_cm': [-156.4631423491728]}
- TS geometry: {'C3_H_A': 1.095, 'C2_H_A': 2.131, 'C1_C2_C3_H_deg': 91.4, 'H_on': 'C3 (secondary)'}; core RMSD to the production TS 1.225 Å
- plus side: reached None; {'C3_H_A': 1.096, 'C2_H_A': 2.14, 'C1_C2_C3_H_deg': 94.2, 'H_on': 'C3 (secondary)'}; secondary_anti RMSD 0.727 Å, ΔE +0.06 kJ/mol; secondary_eclip_irc RMSD 1.077 Å, ΔE -3.00 kJ/mol; tertiary_irc RMSD 1.232 Å, ΔE +39.45 kJ/mol
- minus side: reached None; {'C3_H_A': 1.095, 'C2_H_A': 2.13, 'C1_C2_C3_H_deg': 92.0, 'H_on': 'C3 (secondary)'}; secondary_anti RMSD 0.738 Å, ΔE +1.48 kJ/mol; secondary_eclip_irc RMSD 1.102 Å, ΔE -1.58 kJ/mol; tertiary_irc RMSD 1.243 Å, ΔE +40.87 kJ/mol
- Connects ['None', 'None']; wanted ['secondary_anti', 'secondary_eclip_irc']: **NOT confirmed**

## ts_methyl_shift
- Start: scan_pos point 2 (bridging methyl between C3 and C2)
- TS job: {'ok': True, 'opt_converged': True, 'imaginary_cm': [-165.53667151450983]}
- TS geometry: {'C3_H_A': 1.096, 'C2_H_A': 2.132, 'C1_C2_C3_H_deg': 8.4, 'H_on': 'C3 (secondary)', 'C2_Me4_A': 1.733, 'C3_Me4_A': 1.769, 'C2_Me5_A': 2.544, 'C3_Me5_A': 1.493, 'Me_moved_to_C2': [4]}; core RMSD to the production TS 0.817 Å
- plus side: reached None; {'C3_H_A': 1.096, 'C2_H_A': 2.133, 'C1_C2_C3_H_deg': 10.2, 'H_on': 'C3 (secondary)', 'C2_Me4_A': 1.718, 'C3_Me4_A': 1.781, 'C2_Me5_A': 2.545, 'C3_Me5_A': 1.492, 'Me_moved_to_C2': [4]}; secondary_anti RMSD 0.053 Å, ΔE +1.55 kJ/mol; secondary_eclip_irc RMSD 1.092 Å, ΔE -1.51 kJ/mol; tertiary_irc RMSD 0.808 Å, ΔE +40.94 kJ/mol; secondary_anti_swapped RMSD 0.037 Å, ΔE +1.55 kJ/mol
- minus side: reached secondary_anti; {'C3_H_A': 1.094, 'C2_H_A': 2.124, 'C1_C2_C3_H_deg': -1.2, 'H_on': 'C3 (secondary)', 'C2_Me4_A': 1.804, 'C3_Me4_A': 1.674, 'C2_Me5_A': 2.55, 'C3_Me5_A': 1.503, 'Me_moved_to_C2': []}; secondary_anti RMSD 0.002 Å, ΔE +0.01 kJ/mol; secondary_eclip_irc RMSD 1.084 Å, ΔE -3.05 kJ/mol; tertiary_irc RMSD 0.784 Å, ΔE +39.40 kJ/mol; secondary_anti_swapped RMSD 0.069 Å, ΔE +0.01 kJ/mol
- Connects ['None', 'secondary_anti']; wanted ['secondary_anti', 'secondary_anti_swapped']: **NOT confirmed**

## Relaxed scan, arm neg (-1.2° → -161.9°)
| Point | C1–C2–C3–H (°) | E rel. secondary_anti (kJ/mol) | C3–H (Å) | C2–H (Å) | H on |
|---|---|---|---|---|---|
| 1 | -1.2 | +0.0 | 1.094 | 2.124 | C3 (secondary) |
| 2 | -10.7 | +1.6 | 1.095 | 2.123 | C3 (secondary) |
| 3 | -20.1 | +6.4 | 1.095 | 2.125 | C3 (secondary) |
| 4 | -29.6 | +14.0 | 1.096 | 2.128 | C3 (secondary) |
| 5 | -39.0 | +20.5 | 1.1 | 2.106 | C3 (secondary) |
| 6 | -48.5 | +24.2 | 1.106 | 2.066 | C3 (secondary) |
| 7 | -57.9 | +24.9 | 1.116 | 1.998 | C3 (secondary) |
| 8 | -67.4 | +13.2 | 1.28 | 1.333 | C3 (secondary) |
| 9 | -76.8 | +0.1 | 1.336 | 1.284 | C2 (tertiary) |
| 10 | -86.3 | -10.0 | 1.433 | 1.231 | C2 (tertiary) |
| 11 | -95.7 | -20.2 | 1.705 | 1.156 | C2 (tertiary) |
| 12 | -105.2 | -31.7 | 1.862 | 1.132 | C2 (tertiary) |
| 13 | -114.6 | -38.2 | 1.939 | 1.123 | C2 (tertiary) |
| 14 | -124.1 | -39.4 | 2.008 | 1.115 | C2 (tertiary) |
| 15 | -133.5 | -39.4 | 2.09 | 1.104 | C2 (tertiary) |
| 16 | -143.0 | -37.1 | 2.115 | 1.101 | C2 (tertiary) |
| 17 | -152.4 | -31.8 | 2.129 | 1.099 | C2 (tertiary) |
| 18 | -161.9 | -24.9 | 2.135 | 1.098 | C2 (tertiary) |

- Hydride moved to C2 at point: 9; highest point 7 (interior: True); TS guess: no

## Relaxed scan, arm pos (-1.2° → 161.9°)
| Point | C1–C2–C3–H (°) | E rel. secondary_anti (kJ/mol) | C3–H (Å) | C2–H (Å) | H on |
|---|---|---|---|---|---|
| 1 | -1.2 | +0.0 | 1.094 | 2.124 | C3 (secondary) |
| 2 | 8.4 | +1.6 | 1.094 | 2.127 | C3 (secondary) |
| 3 | 18.0 | +0.9 | 1.096 | 2.144 | C3 (secondary) |
| 4 | 27.6 | +4.5 | 1.097 | 2.15 | C3 (secondary) |
| 5 | 37.2 | +11.2 | 1.097 | 2.158 | C3 (secondary) |
| 6 | 46.8 | +20.0 | 1.098 | 2.164 | C3 (secondary) |
| 7 | 56.4 | +21.9 | 1.1 | 2.175 | C3 (secondary) |
| 8 | 66.0 | +22.8 | 1.1 | 2.17 | C3 (secondary) |
| 9 | 75.6 | +9.1 | 1.097 | 2.152 | C3 (secondary) |
| 10 | 85.1 | +2.5 | 1.096 | 2.145 | C3 (secondary) |
| 11 | 94.7 | +0.1 | 1.095 | 2.138 | C3 (secondary) |
| 12 | 104.3 | +1.1 | 1.094 | 2.122 | C3 (secondary) |
| 13 | 113.9 | +8.4 | 1.095 | 2.123 | C3 (secondary) |
| 14 | 123.5 | +23.8 | 1.096 | 2.133 | C3 (secondary) |
| 15 | 133.1 | +47.3 | 1.097 | 2.153 | C3 (secondary) |
| 16 | 142.7 | +2.4 | 1.095 | 2.116 | C3 (secondary) |
| 17 | 152.3 | +0.1 | 1.094 | 2.122 | C3 (secondary) |
| 18 | 161.9 | +0.9 | 1.094 | 2.127 | C3 (secondary) |

- Hydride moved to C2 at point: None; highest point 15 (interior: True); TS guess: yes

## ΔG‡ and ΔE‡ from secondary_anti

*Note added 2026-11-26 (J. Picione's decision): the four "production TS" rows are energy differences only. `secondary_anti` does not reach the production TS, so they are not the barrier of any path. The anti form's barrier is the `ts_anti_rotation` rows; the barrier of the path shown is measured from the eclipsed form (ΔG‡ 11.1 kJ/mol, 2.7 kcal/mol, at DLPNO-CCSD(T)/ma-def2-TZVP, 298 K).*

| TS | Level | ΔE‡ | ΔG‡ |
|---|---|---|---|
| production TS (eclip side, 01_ts) | M06-2X_ma-def2-SVP (working) | 19.0 (4.5) | 15.8 (3.8) |
| production TS (eclip side, 01_ts) | M06-2X_ma-def2-TZVP | 20.5 (4.9) | 17.2 (4.1) |
| production TS (eclip side, 01_ts) | DLPNO-CCSD(T)_ma-def2-TZVP | 17.3 (4.1) | 14.0 (3.3) |
| production TS (eclip side, 01_ts) | B3LYP_ma-def2-TZVP | 3.7 (0.9) | 0.5 (0.1) |
| ts_anti_rotation | M06-2X_ma-def2-SVP (working) | 25.0 (6.0) | 20.9 (5.0) |
| ts_anti_rotation | M06-2X_ma-def2-TZVP | 26.2 (6.3) | 22.2 (5.3) |
| ts_anti_rotation | DLPNO-CCSD(T)_ma-def2-TZVP | 22.6 (5.4) | 18.6 (4.4) |
| ts_anti_rotation | B3LYP_ma-def2-TZVP | 9.3 (2.2) | 5.2 (1.2) |
| ts_anti_eclip_pos | M06-2X_ma-def2-SVP (working) | 1.6 (0.4) | 3.2 (0.8) |
| ts_anti_eclip_pos | M06-2X_ma-def2-TZVP | 0.6 (0.1) | 2.2 (0.5) |
| ts_anti_eclip_pos | DLPNO-CCSD(T)_ma-def2-TZVP | 0.5 (0.1) | 2.1 (0.5) |
| ts_anti_eclip_pos | B3LYP_ma-def2-TZVP | 0.7 (0.2) | 2.3 (0.5) |
| ts_methyl_shift | M06-2X_ma-def2-SVP (working) | 1.7 (0.4) | 3.1 (0.7) |
| ts_methyl_shift | M06-2X_ma-def2-TZVP | 0.7 (0.2) | 2.1 (0.5) |
| ts_methyl_shift | DLPNO-CCSD(T)_ma-def2-TZVP | 0.7 (0.2) | 2.0 (0.5) |
| ts_methyl_shift | B3LYP_ma-def2-TZVP | 0.8 (0.2) | 2.2 (0.5) |
