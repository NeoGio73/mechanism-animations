Sensitivity test only: every other number in this summary stays on SMD water without explicit waters. Cluster-continuum model: n = 1, 2, 3 explicit waters, the same n on both sides (CN⁻·(H₂O)n + acetone → TS·(H₂O)n → alkoxide·(H₂O)n), all in SMD water; lowest-G arrangement for each species and n (working level); single points on those at the other levels with working-level G − E; 298.15 K, 1 M. kJ/mol (kcal/mol).

**ΔG‡ from separated reactants:**

| Level | n = 0 (production) | n = 1 | n = 2 | n = 3 | n = 3, earlier reference (CCC) |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 67.4 (16.1) | 73.5 (17.6) | 71.0 (17.0) | 70.2 (16.8) | 62.4 (14.9) |
| M06-2X/ma-def2-TZVP | 72.1 (17.2) | 77.0 (18.4) | 74.2 (17.7) | 74.5 (17.8) | 67.1 (16.0) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 79.3 (19.0) | 85.4 (20.4) | 85.8 (20.5) | 85.0 (20.3) | 73.5 (17.6) |
| B3LYP/ma-def2-TZVP | 76.6 (18.3) | 84.2 (20.1) | 84.5 (20.2) | 82.9 (19.8) | 72.4 (17.3) |
| Measured | 69.6 (16.6) | | | | |

**ΔrG:**

| Level | n = 0 (production) | n = 1 | n = 2 | n = 3 | n = 3, earlier reference (CCC) |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | +26.3 (+6.3) | +20.4 (+4.9) | +10.9 (+2.6) | +0.2 (+0.0) | -7.7 (-1.8) |
| M06-2X/ma-def2-TZVP | +33.0 (+7.9) | +27.6 (+6.6) | +19.3 (+4.6) | +10.3 (+2.5) | +2.9 (+0.7) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +36.4 (+8.7) | +33.2 (+7.9) | +24.7 (+5.9) | +14.7 (+3.5) | +3.2 (+0.8) |
| B3LYP/ma-def2-TZVP | +45.8 (+11.0) | +43.8 (+10.5) | +35.6 (+8.5) | +27.0 (+6.4) | +16.4 (+3.9) |
| Measured | +10.5 (+2.5) | | | | |

**Retry of two cyanide·(H₂O)3 arrangements** (J. Picione 2026-11-27, fixed before the run: one attempt each, started 0.1 Å along the remaining imaginary mode, same level and keywords (Opt Freq); no second attempt, no change of settings. A true minimum joins the usable set and the n = 3 reactant is again the lowest working-level G of all usable arrangements; single points on any new kept structure and the n = 3 rows rebuilt. If neither is a minimum, CCC stays and the n = 3 values are lower bounds; no energy from a structure with an imaginary mode is used.)

- cyanide_w3_NNC: was 85i cm⁻¹; after the retry a true minimum, usable.
- cyanide_w3_NCC: was 108i cm⁻¹; after the retry still not a minimum (33i cm⁻¹), not used.

The n = 3 reactant is now cyanide·(H₂O)3 NNC, the lowest working-level free energy of all usable arrangements; it replaces CCC, which was kept only because the two arrangements with water on both ends of cyanide still had an imaginary frequency. The last column keeps the earlier n = 3 values.

**Lowest structures and their hydrogen bonds** (H···X and O(water)···X, Å):

| Structure | Kept arrangement | Water → acceptor: H···X / O···X |
|---|---|---|
| cyanide·(H₂O)1 | C | 1 → C (cyanide): 1.91 / 2.90 |
| ts·(H₂O)1 | O (-308i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.75 / 2.73; 1 → N (nitrile): 2.48 / 3.34 |
| alkoxide·(H₂O)1 | O | 1 → O (alkoxide/carbonyl): 1.48 / 2.52 |
| cyanide·(H₂O)2 | NC | 1 → C (cyanide): 1.92 / 2.91; 2 → N (nitrile): 1.80 / 2.78 |
| ts·(H₂O)2 | OO (-294i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.79 / 2.76; 1 → N (nitrile): 2.48 / 3.35; 2 → O (alkoxide/carbonyl): 1.80 / 2.76; 2 → N (nitrile): 2.47 / 3.34 |
| alkoxide·(H₂O)2 | OO | 1 → O (alkoxide/carbonyl): 1.56 / 2.58; 2 → O (alkoxide/carbonyl): 1.55 / 2.57 |
| cyanide·(H₂O)3 | NNC | 1 → C (cyanide): 1.94 / 2.92; 2 → N (nitrile): 1.93 / 2.91; 3 → N (nitrile): 1.91 / 2.89 |
| ts·(H₂O)3 | OON (-315i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → N (nitrile): 1.84 / 2.82; 2 → O (alkoxide/carbonyl): 1.79 / 2.75; 2 → N (nitrile): 2.51 / 3.37; 3 → O (alkoxide/carbonyl): 1.79 / 2.76; 3 → N (nitrile): 2.53 / 3.38 |
| alkoxide·(H₂O)3 | OOO | 1 → O (alkoxide/carbonyl): 1.61 / 2.61; 2 → O (alkoxide/carbonyl): 1.61 / 2.61; 3 → O (alkoxide/carbonyl): 1.61 / 2.61 |

**Lowest frequency of every structure** (read by script from each final .hess; negative = imaginary), cm⁻¹:

| Structure | Job | Lowest frequency | Imaginary |
|---|---|---|---|
| alkoxide_w1_N | alkoxide_w1_N_min | 26.7 | 0 |
| alkoxide_w1_O | alkoxide_w1_O | 50.9 | 0 |
| alkoxide_w2_NN | alkoxide_w2_NN | 24.7 | 0 |
| alkoxide_w2_ON | alkoxide_w2_ON | 25.1 | 0 |
| alkoxide_w2_OO | alkoxide_w2_OO_min | 11.8 | 0 |
| alkoxide_w3_ONN | alkoxide_w3_ONN | 28.5 | 0 |
| alkoxide_w3_OON | alkoxide_w3_OON | 21.1 | 0 |
| alkoxide_w3_OOO | alkoxide_w3_OOO_min | 24.7 | 0 |
| cyanide_w1_C | cyanide_w1_C | 80.7 | 0 |
| cyanide_w1_N | cyanide_w1_N | 57.9 | 0 |
| cyanide_w2_CC | cyanide_w2_CC_min | -82.0 | 1 |
| cyanide_w2_NC | cyanide_w2_NC | 19.0 | 0 |
| cyanide_w2_NN | cyanide_w2_NN_min | 23.5 | 0 |
| cyanide_w3_CCC | cyanide_w3_CCC | 53.9 | 0 |
| cyanide_w3_NCC | cyanide_w3_NCC_retry | -33.2 | 1 |
| cyanide_w3_NNC | cyanide_w3_NNC_retry | 17.4 | 0 |
| cyanide_w3_NNN | cyanide_w3_NNN | 65.7 | 0 |
| ts_w1_N | ts_w1_N | -13.7 | 1 |
| ts_w1_O | ts_w1_O | -308.0 | 1 |
| ts_w2_NN | ts_w2_NN | -342.8 | 1 |
| ts_w2_ON | ts_w2_ON | -352.8 | 1 |
| ts_w2_OO | ts_w2_OO | -294.1 | 1 |
| ts_w3_ONN | ts_w3_ONN | -357.5 | 1 |
| ts_w3_OON | ts_w3_OON | -315.3 | 1 |
| ts_w3_OOO | ts_w3_OOO | -270.7 | 1 |

Arrangements tried: 25; usable: 22. Script `pipeline/explicit_water.py`; jobs in `04_explicit_water/`.
