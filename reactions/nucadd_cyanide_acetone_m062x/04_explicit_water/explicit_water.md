Sensitivity test only: every other number in this summary stays on SMD water without explicit waters. Cluster-continuum model: n = 1, 2, 3 explicit waters, the same n on both sides (CN⁻·(H₂O)n + acetone → TS·(H₂O)n → alkoxide·(H₂O)n), all in SMD water; lowest-G arrangement for each species and n (working level); single points on those at the other levels with working-level G − E; 298.15 K, 1 M. kJ/mol (kcal/mol).

**ΔG‡ from separated reactants:**

| Level | n = 0 (production) | n = 1 | n = 2 | n = 3 |
|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 67.4 (16.1) | 73.5 (17.6) | 71.0 (17.0) | 62.4 (14.9) |
| M06-2X/ma-def2-TZVP | 72.1 (17.2) | 77.0 (18.4) | 74.2 (17.7) | 67.1 (16.0) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 79.3 (19.0) | 85.4 (20.4) | 85.8 (20.5) | 73.5 (17.6) |
| B3LYP/ma-def2-TZVP | 76.6 (18.3) | 84.2 (20.1) | 84.5 (20.2) | 72.4 (17.3) |
| Measured | 69.6 (16.6) | | | |

**ΔrG:**

| Level | n = 0 (production) | n = 1 | n = 2 | n = 3 |
|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | +26.3 (+6.3) | +20.4 (+4.9) | +10.9 (+2.6) | -7.7 (-1.8) |
| M06-2X/ma-def2-TZVP | +33.0 (+7.9) | +27.6 (+6.6) | +19.3 (+4.6) | +2.9 (+0.7) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +36.4 (+8.7) | +33.2 (+7.9) | +24.7 (+5.9) | +3.2 (+0.8) |
| B3LYP/ma-def2-TZVP | +45.8 (+11.0) | +43.8 (+10.5) | +35.6 (+8.5) | +16.4 (+3.9) |
| Measured | +10.5 (+2.5) | | | |

**Lowest structures and their hydrogen bonds** (H···X and O(water)···X, Å):

| Structure | Kept arrangement | Water → acceptor: H···X / O···X |
|---|---|---|
| cyanide·(H₂O)1 | C | 1 → C (cyanide): 1.91 / 2.90 |
| ts·(H₂O)1 | O (-308i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.75 / 2.73; 1 → N (nitrile): 2.48 / 3.34 |
| alkoxide·(H₂O)1 | O | 1 → O (alkoxide/carbonyl): 1.48 / 2.52 |
| cyanide·(H₂O)2 | NC | 1 → C (cyanide): 1.92 / 2.91; 2 → N (nitrile): 1.80 / 2.78 |
| ts·(H₂O)2 | OO (-294i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.79 / 2.76; 1 → N (nitrile): 2.48 / 3.35; 2 → O (alkoxide/carbonyl): 1.80 / 2.76; 2 → N (nitrile): 2.47 / 3.34 |
| alkoxide·(H₂O)2 | OO | 1 → O (alkoxide/carbonyl): 1.56 / 2.58; 2 → O (alkoxide/carbonyl): 1.55 / 2.57 |
| cyanide·(H₂O)3 | CCC | 1 → C (cyanide): 2.06 / 3.03; 1 → O of water 2: 2.34 / 3.10; 2 → O of water 3: 1.91 / 2.85; 3 → C (cyanide): 2.00 / 2.97 |
| ts·(H₂O)3 | OON (-315i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → N (nitrile): 1.84 / 2.82; 2 → O (alkoxide/carbonyl): 1.79 / 2.75; 2 → N (nitrile): 2.51 / 3.37; 3 → O (alkoxide/carbonyl): 1.79 / 2.76; 3 → N (nitrile): 2.53 / 3.38 |
| alkoxide·(H₂O)3 | OOO | 1 → O (alkoxide/carbonyl): 1.61 / 2.61; 2 → O (alkoxide/carbonyl): 1.61 / 2.61; 3 → O (alkoxide/carbonyl): 1.61 / 2.61 |

Arrangements tried: 25; usable: 21. Script `pipeline/explicit_water.py`; jobs in `04_explicit_water/`.
