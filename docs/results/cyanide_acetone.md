# Results summary: Nucleophilic addition, cyanide + acetone

`N≡C⁻ + (CH₃)₂C=O → (CH₃)₂C(O⁻)C≡N`, in water (SMD), 298.15 K. Reaction folder `reactions/nucadd_cyanide_acetone_m062x`. Written 2026-11-19 08:43 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "Separated reactants" means cyanide + acetone, each optimized on its own (single conformers). Steps that change the number of molecules use the 1 M standard state (RT ln 24.47 = 7.93 kJ/mol per mole of solute: subtracted for each molecule lost, added for each gained); steps from the reactant complex involve no change in the number of molecules, so no standard-state term applies. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 253 (reactant side 158 frames, TS = frame 158, product side 94 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 253 of 253; M06-2X/ma-def2-TZVP 253 of 253; B3LYP/ma-def2-TZVP 253 of 253; DLPNO-CCSD(T)/ma-def2-TZVP 227 of 253. DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 12.30 bohr·amu½ in the mass-weighted coordinate s (from s = -7.89 at the reactant end to s = 4.41 at the product end; TS at s = 0); 2.30 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (C(N)–C between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 17 (s = -7.05, C(N)–C 3.11 Å) to frame 215 (s = 2.85, 1.60 Å). At the TS the C(N)–C distance is 2.080 Å, 66 % of the way from the reactant end (3.196 Å) to the product end (1.513 Å): most of the shortening has happened by the TS, and the remaining 34 % (0.57 Å) comes after it, over the product side (s = 0 to 4.41).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| C(N)–C | forms | 3.196 | 2.080 | 1.513 | 66 % |
| C–O | lengthens | 1.222 | 1.261 | 1.376 | 25 % |

| Point | Frame | s (bohr·amu½) | C(N)–C (Å) | C–O (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |
|---|---|---|---|---|---|
| reactant end | 0 | -7.89 | 3.196 | 1.222 | -10.8 (-2.6) |
| s ≈ -6 | 38 | -6.00 | 2.983 | 1.222 | -8.8 (-2.1) |
| s ≈ -4 | 78 | -4.00 | 2.716 | 1.224 | -1.5 (-0.4) |
| s ≈ -2 | 118 | -2.00 | 2.416 | 1.231 | +14.4 (+3.4) |
| s ≈ -1 | 138 | -1.00 | 2.250 | 1.242 | +23.1 (+5.5) |
| s ≈ -0.5 | 148 | -0.50 | 2.166 | 1.250 | +25.9 (+6.2) |
| TS | 158 | +0.00 | 2.080 | 1.261 | +27.0 (+6.5) |
| s ≈ +0.5 | 168 | +0.50 | 1.994 | 1.274 | +25.4 (+6.1) |
| s ≈ +1 | 178 | +1.00 | 1.908 | 1.290 | +20.8 (+5.0) |
| s ≈ +1.5 | 188 | +1.50 | 1.822 | 1.308 | +13.2 (+3.1) |
| s ≈ +2 | 198 | +2.00 | 1.738 | 1.327 | +3.1 (+0.7) |
| s ≈ +3 | 218 | +3.00 | 1.582 | 1.364 | -17.3 (-4.1) |
| product end | 252 | +4.41 | 1.513 | 1.376 | -22.5 (-5.4) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD water)

- **TS:** 1 imaginary frequency, -366.8 cm⁻¹ (written 366.8i cm⁻¹). Stretch rates of the key bonds in the mode: C(N)–C +1.04, C–O -0.14; forming and breaking bonds move together. C(N)–C 2.080 Å, C–O 1.261 Å, Bürgi–Dunitz angle C(N)···C=O 106.2°.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |
|---|---|---|---|---|---|
| cyanide (`[C-]#N`, job `cyanide`) | separated reactant | C–N 1.174 Å | none | — | — |
| acetone (`CC(C)=O`, job `acetone_min`) | separated reactant | C–C 1.500 Å, C–O 1.221 Å | none | — | — |
| reactant complex | reactant end of the path | C(N)–C 3.208 Å, C–O 1.221 Å | none | -11.8 (-2.8) | +18.3 (+4.4) |
| TS |  | C(N)–C 2.080 Å, C–O 1.261 Å | 1 | +27.1 (+6.5) | +67.4 (+16.1) |
| product complex | product end of the path | C(N)–C 1.514 Å, C–O 1.376 Å | none | -22.5 (-5.4) | +26.4 (+6.3) |
| acetone_cyanohydrin_alkoxide (`CC(C)(C#N)[O-]`) | separated product | C–C 1.513–1.540 Å, C–N 1.159 Å, C–O 1.376 Å | none | -22.8 (-5.5) | +26.3 (+6.3) |

The product is one molecule, so the product complex and the separated product are the same species optimized twice (from the path end and on its own); they differ by +0.36 kJ/mol in E. Reaction energies use the separated product.
The reactant end of the path ended 0.245 Å (RMSD) from the reactant-end minimum, the product end 0.003 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = separated reactants (cyanide + acetone), 298.15 K, ΔG at 1 M:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 27.1 (6.5) | 28.8 (6.9) | 67.4 (16.1) |
| M06-2X/ma-def2-TZVP | 31.8 (7.6) | 33.5 (8.0) | 72.1 (17.2) |
| B3LYP/ma-def2-TZVP | 36.3 (8.7) | 38.0 (9.1) | 76.6 (18.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 39.0 (9.3) | 40.7 (9.7) | 79.3 (19.0) |

**Barriers, zero = reactant complex, 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 38.8 (9.3) | 36.8 (8.8) | 49.1 (11.7) |
| M06-2X/ma-def2-TZVP | 43.0 (10.3) | 41.0 (9.8) | 53.3 (12.7) |
| B3LYP/ma-def2-TZVP | 43.0 (10.3) | 41.0 (9.8) | 53.3 (12.7) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 47.1 (11.3) | 45.1 (10.8) | 57.4 (13.7) |

M06-2X/ma-def2-TZVP and B3LYP/ma-def2-TZVP give the same rounded barriers from the complex by coincidence, not by copying: ΔE‡ 42.98 and 42.99 kJ/mol, ΔG‡ 53.26 and 53.27 kJ/mol (their barriers from separated reactants differ by 4.5 kJ/mol, and so do their complex binding energies).

**Reaction, zero = separated reactants; product = separated acetone_cyanohydrin_alkoxide, 298.15 K, ΔG at 1 M:**

| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from acetone_cyanohydrin_alkoxide, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | -22.8 (-5.5) | -16.1 (-3.8) | +26.3 (+6.3) | -11.8 (-2.8) / +18.3 (+4.4) | 49.9 (11.9) / 41.1 (9.8) |
| M06-2X/ma-def2-TZVP | -16.1 (-3.8) | -9.3 (-2.2) | +33.0 (+7.9) | -11.2 (-2.7) / +18.9 (+4.5) | 47.9 (11.4) / 39.1 (9.4) |
| B3LYP/ma-def2-TZVP | -3.2 (-0.8) | +3.5 (+0.8) | +45.8 (+11.0) | -6.7 (-1.6) / +23.4 (+5.6) | 39.5 (9.4) / 30.8 (7.4) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -12.6 (-3.0) | -5.9 (-1.4) | +36.4 (+8.7) | -8.1 (-1.9) / +21.9 (+5.2) | 51.6 (12.3) / 42.9 (10.3) |

ΔrH was not asked for; it is given because it costs nothing. The reactant complex lies below separated reactants in E but above them in G at 1 M at every level: the binding energy does not pay for the entropy lost on association.

## 4. Highest point of the path at each level

Path energies at each level, zero = separated reactants (cyanide + acetone) at that level. The TS is frame 158 (s = 0, C(N)–C 2.080 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C(N)–C (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 253 of 253 | 27.0 (6.5) | 158 | +0.00 | 2.080 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 253 of 253 | 31.8 (7.6) | 161 | +0.15 | 2.054 | +0.15 (+0.04) |
| B3LYP/ma-def2-TZVP | 253 of 253 | 36.3 (8.7) | 168 | +0.50 | 1.994 | +1.42 (+0.34) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 227 of 253 | 39.0 (9.3) | 160 | +0.10 | 2.063 | +0.17 (+0.04) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter C(N)–C distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0159 | 0.052950 |
| reactant_complex | 0.0148 | 0.048768 |
| product_complex | 0.0134 | 0.053322 |
| cyanide | 0.0166 | 0.047862 |
| acetone | 0.0140 | 0.048100 |
| acetone_cyanohydrin_alkoxide | 0.0134 | 0.053687 |

- **Path frames (227 computed):** T1 from 0.0134 to 0.0159, median 0.0151; largest at frame 152 (s = -0.30, C(N)–C 2.132 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: passes** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 158 | 0 | 0.22° / 1.64° (20) | 0.06° / 5.18° (137) | yes |
| product side (forward) | 94 | 17 | 0.12° / 1.47° (20) | 0.02° / 0.25° (56) | yes |

- **End checks: ok** (`02_irc/true_irc/end_check_all.md`, `02_irc/true_irc/end_check.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame +0.941 kJ/mol from `reactant_complex` (RMSD 0.245 Å). Opt Freq from the last frame: E -0.072 kJ/mol from it, heavy-atom RMSD 0.110 Å (all-atom 0.393 Å), largest key-distance difference +0.007 (C(N)–C) Å, 0 imaginary frequencies: reaches the same minimum.
  - Product side (forward): last frame -0.028 kJ/mol from `product_complex` (RMSD 0.003 Å). Opt Freq from the last frame: E -0.031 kJ/mol from it, heavy-atom RMSD 0.003 Å (all-atom 0.003 Å), largest key-distance difference -0.000 (C(N)–C) Å, 0 imaginary frequencies: reaches the same minimum.
  - Official end structures now: reactant end `reactant_complex`, product end `product_complex` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: 0.94 kJ/mol; product side: −0.03 kJ/mol (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 1 (stopping rule near a flat TS): the 20 flat steps on both sides began at least 10 kJ/mol below the TS, so the amended rule stops them at the same frame; not rerun.
  - Entry 4 (the minimum the path reaches is the official end): product side reaches the same minimum (E −0.031 kJ/mol, heavy-atom RMSD 0.003 Å); no switch, no value changed.
  - Entry 5: the side-end row above.
  - Entry 6 (angle against gradient size): 18 of the 70 frames over 15° in all six paths are on this reaction's product side, in the final approach to its minimum, where the RMS gradient is about 1.5 × 10⁻⁵ Eh/bohr; the largest, up to 117°, are there. Near this TS the angles stay small.
  - Entry 9 (lowest-free-energy conformers): cyanide, acetone and the alkoxide have single conformers; nothing changes.
  - No entry corrects a value reported for this reaction.

## 7. Experiment

Experiment (water, 25 °C): cyanide addition to acetone, alkoxide formation. Equilibrium K = [R₂C(O⁻)CN]/([R₂C=O][CN⁻]), log K = −1.84, K = 0.0145 M⁻¹, ΔG° = +10.5 kJ/mol (+2.5 kcal/mol); rate = k[CN⁻][R₂C=O], log k = 0.60, k = 4.0 M⁻¹ s⁻¹, ΔG‡ = 69.6 kJ/mol (16.6 kcal/mol) (Eyring, 298.15 K, 1 M standard state, no transmission correction). All three measured values are in the tables above, next to the computed ones. Source: Schlesinger and Miller 1973 (G. Schlesinger and S. L. Miller, J. Am. Chem. Soc. 1973, 95, 3729–3735, doi:10.1021/ja00792a043), Table IV, acetone row: k₂ = 3.97 M⁻¹ s⁻¹ (CN⁻ + acetone → alkoxide), K₂ = 68.5 M (alkoxide → acetone + CN⁻; 1/K₂ = 0.0146 M⁻¹), k₋₂ = 272 s⁻¹ (alkoxide → acetone + CN⁻), giving ΔG‡ = 69.6 kJ/mol (16.6 kcal/mol), ΔG° = +10.5 kJ/mol (+2.5 kcal/mol) and ΔG‡(reverse) = 59.1 kJ/mol (14.1 kcal/mol). The same values, as log k = 0.60 and log K = −1.84, are compiled in Guthrie 1998 (J. P. Guthrie, J. Am. Chem. Soc. 1998, 120, 1688–1694, doi:10.1021/ja970506q, Table 1). Values confirmed by J. Picione from Guthrie's Table 1 and from the original paper (2026-11-11).

**Computed against measured** (computed − measured; zero = separated reactants (cyanide + acetone), 298.15 K, 1 M):

| Level | ΔG‡ | ΔrG | Reverse ΔG‡ (from acetone_cyanohydrin_alkoxide) |
|---|---|---|---|
| Measured | 69.6 (16.6) | +10.5 (+2.5) | 59.1 (14.1) |
| M06-2X/ma-def2-SVP (working level, geometries) | 67.4 (16.1), -2.2 (-0.5) | +26.3 (+6.3), +15.8 (+3.8) | 41.1 (9.8), -18.0 (-4.3) |
| M06-2X/ma-def2-TZVP | 72.1 (17.2), +2.5 (+0.6) | +33.0 (+7.9), +22.5 (+5.4) | 39.1 (9.4), -20.0 (-4.8) |
| B3LYP/ma-def2-TZVP | 76.6 (18.3), +7.0 (+1.7) | +45.8 (+11.0), +35.3 (+8.4) | 30.8 (7.4), -28.3 (-6.8) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 79.3 (19.0), +9.7 (+2.3) | +36.4 (+8.7), +25.9 (+6.2) | 42.9 (10.3), -16.2 (-3.9) |

Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):

- The path ends at the alkoxide, (CH₃)₂C(CN)O⁻, which lies above separated cyanide + acetone in free energy (ΔG = +26.3 kJ/mol (+6.3 kcal/mol) at the working level and +36.4 kJ/mol (+8.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP, 298 K, 1 M; electronically it is 22.8 kJ/mol (5.5 kcal/mol) and 12.6 kJ/mol (3.0 kcal/mol) below): the addition step on its own is uphill. The neutral cyanohydrin forms only after a later proton transfer (from HCN or water) that this path does not model.
- Implicit solvent under-stabilizes the alkoxide: SMD has no explicit water molecules, so the hydrogen bonds that water donates to the alkoxide oxygen are missing. Against experiment, DLPNO-CCSD(T)/ma-def2-TZVP puts the barrier 9.7 kJ/mol (2.3 kcal/mol) too high (79.3 against 69.6 kJ/mol) and the reaction free energy 25.9 kJ/mol (6.2 kcal/mol) too high (+36.4 against +10.5 kJ/mol). The test with added water molecules (`04_explicit_water/explicit_water.md`) supports the missing hydrogen bonds for the reaction free energy. It does not explain the error in the barrier: the added waters raise the barrier by about 6 kJ/mol, away from the measured value. The barrier error is left unassigned.
  _2026-11-27 (J. Picione): the argument that the barrier error is about 38 % of the product error, in proportion to the transition state's pyramidalization, was replaced after the explicit-water test._
- The working-level barrier, ΔG‡ = 67.4 kJ/mol (16.1 kcal/mol), is closer to experiment (69.6 kJ/mol) than coupled cluster only through cancellation of errors: with the same functional, the larger ma-def2-TZVP basis raises it to 72.1 kJ/mol (17.2 kcal/mol).
- The measured equilibrium constant rests on an estimated acidity of the cyanohydrin O–H (K₃ = 6.3 × 10⁻¹³, pKa 12.2, from a substituent-constant correlation, Ballinger and Long 1960, not a direct measurement), so ΔG° = +10.5 kJ/mol carries an uncertainty of a few kJ/mol. The forward rate constant does not depend on it. The reverse rate constant does: Schlesinger and Miller derive it rather than measure it, k₋₂ = k₂K₂ = k₂K′K_HCN/K₃ (their equation 24, p. 3734; 3.97 × 68.5 = 272 s⁻¹), so the reverse barrier from it (59.1 kJ/mol) carries the same uncertainty.
- Table IV of Schlesinger and Miller does not state where the acetone rate constant comes from. They cite Svirbely and Roth 1953 as the earlier kinetic study of acetone, so that is its probable origin; this is not confirmed.
- The measured rate constant refers to free cyanide (rate = k[CN⁻][R₂C=O]); observed rates depend on pH, because only CN⁻, not HCN, adds to the carbonyl.
- Reactant end of the path: the stopping rule ended the reactant side 0.94 kJ/mol (0.22 kcal/mol) above the optimized reactant complex, at an RMSD of 0.245 Å, while the energy was still falling slowly in this loosely bound complex. An optimization from that last frame reaches the same complex: −0.07 kJ/mol (−0.02 kcal/mol) from it, no imaginary frequency, all-atom RMSD 0.393 Å (heavy atoms 0.110 Å), the C(N)···C distance within 0.007 Å; it differs only in the rotation of the methyl groups (accepted by J. Picione, 2026-11-14).

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 253 | 2.4 | 34 | 38 | 40 (frame 247) | 2026-11-16 22:19 | 2026-11-17 07:11 |
| B3LYP/ma-def2-TZVP | 253 | 1.2 | 17 | 17 | 20 (frame 43) | 2026-11-17 07:15 | 2026-11-17 14:30 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 227 | 26.3 | 423 | 465 | 491 (frame 232) | 2026-11-17 14:33 | 2026-11-18 21:13 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 361 gradient jobs, 5.0 h job time (folders `backward`); product side (forward): 295 gradient jobs, 4.2 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.
