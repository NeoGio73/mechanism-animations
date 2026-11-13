# E2 arbitration (EtO⁻ + 2-bromopropane, SMD ethanol): parts a, b and c

DLPNO-CCSD(T) TightPNO/ma-def2-TZVP SMD ethanol single points throughout. ΔG‡ at 298.15 K, 1 M, from separated reactants.

## a. The 7.4 kJ/mol between the CC barriers on B3LYP and M06-2X geometries

| Geometry set | CC ΔE‡ (electronic) | thermal part (own frequencies) | 1 atm → 1 M | CC ΔG‡ |
|---|---|---|---|---|
| B3LYP | 78.6 kJ/mol (18.8 kcal/mol) | 36.5 kJ/mol (8.7 kcal/mol) | −7.92 kJ/mol | 107.2 kJ/mol (25.6 kcal/mol) |
| M06-2X | 72.2 kJ/mol (17.2 kcal/mol) | 35.5 kJ/mol (8.5 kcal/mol) | −7.92 kJ/mol | 99.8 kJ/mol (23.8 kcal/mol) |
| B3LYP − M06-2X | 6.4 kJ/mol (1.5 kcal/mol) | 1.0 kJ/mol (0.2 kcal/mol) | 0 | 7.4 kJ/mol (1.8 kcal/mol) |

## b. ωB97X-D3(BJ) and RI-MP2 TS (ma-def2-SVP, SMD ethanol), CC barrier on each

| Level | TS from | O–H (Å) | Cβ–H (Å) | Cα–Br (Å) | imaginary (cm⁻¹) | accepted | CC ΔE‡ | thermal | CC ΔG‡ |
|---|---|---|---|---|---|---|---|---|---|
| B3LYP-D3(BJ) | (production) | 1.392 | 1.279 | 2.461 | 888i | yes | 78.6 kJ/mol (18.8 kcal/mol) | 36.5 kJ/mol (8.7 kcal/mol) | 107.2 kJ/mol (25.6 kcal/mol) |
| M06-2X | (level test) | 1.294 | 1.348 | 2.258 | 1390i | yes | 72.2 kJ/mol (17.2 kcal/mol) | 35.5 kJ/mol (8.5 kcal/mol) | 99.8 kJ/mol (23.8 kcal/mol) |
| wB97X-D3(BJ) | B3LYP | 1.263 | 1.378 | 2.308 | 1610i | yes | 73.5 kJ/mol (17.6 kcal/mol) | 35.2 kJ/mol (8.4 kcal/mol) | 100.8 kJ/mol (24.1 kcal/mol) |
| wB97X-D3(BJ) | M06-2X | 1.263 | 1.378 | 2.308 | 1610i | yes | 73.5 kJ/mol (17.6 kcal/mol) | 35.2 kJ/mol (8.4 kcal/mol) | 100.8 kJ/mol (24.1 kcal/mol) |
| RI-MP2 | B3LYP | 1.250 | 1.378 | 2.250 | 1481i | yes | 73.0 kJ/mol (17.5 kcal/mol) | 32.8 kJ/mol (7.8 kcal/mol) | 97.9 kJ/mol (23.4 kcal/mol) |
| RI-MP2 | M06-2X | 1.250 | 1.378 | 2.250 | 1490i, 34i | NO: 2 imaginary frequencies, need exactly 1 | 73.0 kJ/mol (17.5 kcal/mol) | 37.6 kJ/mol (9.0 kcal/mol) | 102.7 kJ/mol (24.5 kcal/mol) |

Reactants at each level: wB97X-D3(BJ): ethoxide imaginary 0, 2-bromopropane imaginary 0; RI-MP2: ethoxide imaginary 0, 2-bromopropane imaginary 0.

## c. Where CC peaks near the TS (energies relative to the TS geometry at the same level)

s in bohr·amu½ from the TS, + toward the reactants on both paths. Peak from the parabola through the highest computed frame and its two neighbours.

| Path | points, abs(s) ≤ 1.2 | neighbours used (s) | peak s | peak above TS geometry | O–H at peak (Å) | Cβ–H (Å) | Cα–Br (Å) | at the TS: O–H / Cβ–H / Cα–Br (Å) |
|---|---|---|---|---|---|---|---|---|
| B3LYP production path | 9 | -0.366, +0.366 | +0.116 | +0.87 kJ/mol (+0.21 kcal/mol) | 1.446 | 1.232 | 2.454 | 1.392 / 1.280 / 2.461 |
| M06-2X short IRC, 0.10 step | 21 | -0.233, +0.000 | -0.147 | +0.98 kJ/mol (+0.23 kcal/mol) | 1.237 | 1.403 | 2.266 | 1.294 / 1.348 / 2.258 |

M06-2X short IRC frames (CC relative to the TS geometry):

| s | E (kJ/mol) | O–H (Å) | Cβ–H (Å) | Cα–Br (Å) |
|---|---|---|---|---|
| -1.621 | -23.269 | 1.084 | 1.563 | 2.402 |
| -1.505 | -20.397 | 1.094 | 1.552 | 2.391 |
| -1.382 | -17.665 | 1.103 | 1.540 | 2.380 |
| -1.259 | -14.860 | 1.113 | 1.528 | 2.369 |
| -1.137 | -12.193 | 1.124 | 1.515 | 2.357 |
| -1.026 | -9.808 | 1.135 | 1.503 | 2.347 |
| -0.907 | -7.502 | 1.147 | 1.490 | 2.336 |
| -0.799 | -5.473 | 1.159 | 1.478 | 2.326 |
| -0.699 | -3.784 | 1.171 | 1.466 | 2.317 |
| -0.560 | -2.100 | 1.186 | 1.450 | 2.304 |
| -0.461 | -0.922 | 1.198 | 1.439 | 2.295 |
| -0.331 | +0.177 | 1.214 | 1.424 | 2.284 |
| -0.233 | +0.645 | 1.226 | 1.413 | 2.274 |
| -0.100 | +0.880 | 1.243 | 1.398 | 2.262 |
| +0.000 | +0.000 | 1.294 | 1.348 | 2.258 |
| +0.100 | -3.180 | 1.346 | 1.297 | 2.255 |
| +0.223 | -5.049 | 1.362 | 1.285 | 2.243 |
| +0.324 | -6.897 | 1.376 | 1.276 | 2.234 |
| +0.452 | -9.125 | 1.393 | 1.263 | 2.223 |
| +0.548 | -11.240 | 1.406 | 1.254 | 2.214 |
| +0.677 | -13.602 | 1.423 | 1.242 | 2.203 |
| +0.776 | -16.034 | 1.436 | 1.233 | 2.195 |
| +0.887 | -18.335 | 1.451 | 1.224 | 2.185 |
| +0.993 | -20.675 | 1.464 | 1.215 | 2.176 |
| +1.101 | -23.007 | 1.478 | 1.207 | 2.167 |
| +1.210 | -25.312 | 1.492 | 1.199 | 2.157 |
| +1.319 | -27.578 | 1.505 | 1.192 | 2.148 |
| +1.428 | -29.755 | 1.518 | 1.185 | 2.139 |
| +1.539 | -31.958 | 1.532 | 1.178 | 2.130 |
