A sensitivity test only: every other number on this card is without added waters. One to three water molecules are placed around the ions, the same number on both sides (CN⁻·(H₂O)n + acetone → TS·(H₂O)n → alkoxide·(H₂O)n), all in SMD water. For each species and n, the true minimum (or, for the transition state, the structure with the reaction's imaginary mode) lowest in free energy at the working level is kept; the other levels are single points on it, with working-level thermal corrections. 298.15 K, 1 M; kJ/mol (kcal/mol).

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

**Why the n = 3 values changed.** The n = 3 cyanide·water cluster with two waters on N and one on C kept a small imaginary frequency after its first optimization (85i cm⁻¹, the free O–H hydrogen of one water), so at first it could not be used and the cluster with waters on C was the reference. One more optimization, started along that mode with the same settings, reached a true minimum (lowest frequency 17.4 cm⁻¹). It lies 7.8 kJ/mol (1.9 kcal/mol) below the earlier reference in working-level free energy (0.7 kJ/mol (0.2 kcal/mol) above it in electronic energy) and is now the n = 3 reactant. The last column keeps the earlier n = 3 values.

**Lowest structures and their hydrogen bonds** (H···X and O(water)···X, Å):

| Structure | Kept arrangement | Water → acceptor: H···X / O···X |
|---|---|---|
| cyanide·(H₂O)1 | C | 1 → C (cyanide): 1.91 / 2.90 |
| ts·(H₂O)1 | O (308i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.75 / 2.73; 1 → N (nitrile): 2.48 / 3.34 |
| alkoxide·(H₂O)1 | O | 1 → O (alkoxide/carbonyl): 1.48 / 2.52 |
| cyanide·(H₂O)2 | NC | 1 → C (cyanide): 1.92 / 2.91; 2 → N (nitrile): 1.80 / 2.78 |
| ts·(H₂O)2 | OO (294i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → O (alkoxide/carbonyl): 1.79 / 2.76; 1 → N (nitrile): 2.48 / 3.35; 2 → O (alkoxide/carbonyl): 1.80 / 2.76; 2 → N (nitrile): 2.47 / 3.34 |
| alkoxide·(H₂O)2 | OO | 1 → O (alkoxide/carbonyl): 1.56 / 2.58; 2 → O (alkoxide/carbonyl): 1.55 / 2.57 |
| cyanide·(H₂O)3 | NNC | 1 → C (cyanide): 1.94 / 2.92; 2 → N (nitrile): 1.93 / 2.91; 3 → N (nitrile): 1.91 / 2.89 |
| ts·(H₂O)3 | OON (315i cm⁻¹, overlap with the bare TS mode 0.99) | 1 → N (nitrile): 1.84 / 2.82; 2 → O (alkoxide/carbonyl): 1.79 / 2.75; 2 → N (nitrile): 2.51 / 3.37; 3 → O (alkoxide/carbonyl): 1.79 / 2.76; 3 → N (nitrile): 2.53 / 3.38 |
| alkoxide·(H₂O)3 | OOO | 1 → O (alkoxide/carbonyl): 1.61 / 2.61; 2 → O (alkoxide/carbonyl): 1.61 / 2.61; 3 → O (alkoxide/carbonyl): 1.61 / 2.61 |

**Arrangements not used** (25 tried, 22 usable):

- ts·(H₂O)1, water on N: not a transition state. It moved back to a loose complex of acetone and cyanide·water (C···C 3.30 Å) with a residual 14i cm⁻¹ mode, overlap 0.06 with the reaction mode.
- cyanide·(H₂O)2, both waters on C: one imaginary frequency left (82i cm⁻¹, the free O–H hydrogen of one water). It lies 4.2 kJ/mol (1.0 kcal/mol) above the kept arrangement in electronic energy and 4.7 kJ/mol (1.1 kcal/mol) above it in free energy, so it would not have been kept.
- cyanide·(H₂O)3, one water on N and two on C: an imaginary frequency left after two optimizations (33i cm⁻¹). Not used; it lies 0.1 kJ/mol (0.0 kcal/mol) below the kept arrangement in electronic energy but 2.4 kJ/mol (0.6 kcal/mol) above it in free energy, so it would not have been kept.

Files: `04_explicit_water/explicit_water.md` and `.json` (with the lowest frequency of every structure); jobs in `04_explicit_water/`.
