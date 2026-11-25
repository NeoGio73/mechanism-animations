# Results summary: E2, ethoxide + 2-bromopropane

`EtO⁻ + (CH₃)₂CHBr → CH₂=CHCH₃ + EtOH + Br⁻`, in ethanol (SMD), 298.15 K. Reaction folder `reactions/e2_ethoxide_2_bromopropane_m062x`. Written 2026-11-25 06:27 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "Separated reactants" means ethoxide + 2-bromopropane, each optimized on its own (each the conformer lowest in free energy (SI entry 9)). Steps that change the number of molecules use the 1 M standard state (RT ln 24.47 = 7.93 kJ/mol per mole of solute: subtracted for each molecule lost, added for each gained); steps from the reactant complex involve no change in the number of molecules, so no standard-state term applies. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 532 (reactant side 132 frames, TS = frame 132, product side 399 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 532 of 532; M06-2X/ma-def2-TZVP 532 of 532; B3LYP/ma-def2-TZVP 532 of 532; DLPNO-CCSD(T)/ma-def2-TZVP 507 of 532. DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 26.54 bohr·amu½ in the mass-weighted coordinate s (from s = -6.59 at the reactant end to s = 19.95 at the product end; TS at s = 0); 5.82 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (O–H or Cβ–H or Cα–Br between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 18 (s = -5.70, O–H 2.33 Å) to frame 497 (s = 18.25, 0.97 Å). At the TS the O–H distance is 1.295 Å, 77 % of the way from the reactant end (2.409 Å) to the product end (0.970 Å): most of the shortening has happened by the TS, and the remaining 23 % (0.32 Å) comes after it, over the product side (s = 0 to 19.95).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| O–H | forms | 2.409 | 1.295 | 0.970 | 77 % |
| Cβ–H | breaks | 1.103 | 1.348 | 2.446 | 18 % |
| Cα–Br | breaks | 1.997 | 2.259 | 3.857 | 14 % |
| Cα–Cβ | shortens | 1.512 | 1.441 | 1.337 | 41 % |

| Point | Frame | s (bohr·amu½) | O–H (Å) | Cβ–H (Å) | Cα–Br (Å) | Cα–Cβ (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |
|---|---|---|---|---|---|---|---|
| reactant end | 0 | -6.59 | 2.409 | 1.103 | 1.997 | 1.512 | -8.5 (-2.0) |
| s ≈ -6 | 12 | -5.99 | 2.365 | 1.103 | 1.997 | 1.512 | -8.2 (-2.0) |
| s ≈ -4 | 52 | -4.00 | 2.100 | 1.104 | 1.997 | 1.510 | -5.0 (-1.2) |
| s ≈ -2 | 92 | -2.00 | 1.766 | 1.115 | 2.088 | 1.503 | +9.8 (+2.3) |
| s ≈ -1 | 112 | -1.00 | 1.604 | 1.135 | 2.181 | 1.485 | +27.8 (+6.6) |
| s ≈ -0.5 | 122 | -0.50 | 1.513 | 1.159 | 2.229 | 1.467 | +38.8 (+9.3) |
| TS | 132 | +0.00 | 1.295 | 1.348 | 2.259 | 1.441 | +50.0 (+12.0) |
| s ≈ +0.5 | 142 | +0.50 | 1.056 | 1.593 | 2.283 | 1.414 | +34.0 (+8.1) |
| s ≈ +1 | 152 | +1.00 | 1.015 | 1.680 | 2.344 | 1.380 | +21.2 (+5.1) |
| s ≈ +1.5 | 162 | +1.50 | 1.002 | 1.731 | 2.412 | 1.364 | +9.7 (+2.3) |
| s ≈ +2 | 172 | +2.00 | 0.994 | 1.770 | 2.481 | 1.354 | -1.2 (-0.3) |
| s ≈ +3 | 192 | +3.00 | 0.985 | 1.832 | 2.615 | 1.343 | -21.3 (-5.1) |
| s ≈ +6 | 252 | +6.00 | 0.977 | 1.999 | 2.943 | 1.335 | -61.7 (-14.7) |
| s ≈ +10 | 332 | +10.00 | 0.972 | 2.230 | 3.272 | 1.335 | -86.7 (-20.7) |
| product end | 531 | +19.95 | 0.970 | 2.446 | 3.857 | 1.337 | -104.9 (-25.1) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD ethanol)

- **TS:** 1 imaginary frequency, -1389.4 cm⁻¹ (written 1389.4i cm⁻¹). Stretch rates of the key bonds in the mode: O–H +1.04, Cβ–H -1.01, Cα–Br -0.07, Cα–Cβ +0.09; forming and breaking bonds move together. O–H 1.295 Å, Cβ–H 1.348 Å, Cα–Br 2.258 Å, Cα–Cβ 1.441 Å, H–Cβ–Cα–Br dihedral 176.9°.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |
|---|---|---|---|---|---|
| ethoxide (`CC[O-]`, job `ethoxide`) | separated reactant | C–C 1.532 Å, C–O 1.377 Å | none | — | — |
| 2-bromopropane (`CC(C)Br`, job `2-bromopropane`) | separated reactant | C–Br 1.994 Å, C–C 1.513 Å | none | — | — |
| reactant complex (`reactant_complex_irc`) | reactant end of the path | O–H 2.623 Å, Cβ–H 1.102 Å, Cα–Br 1.997 Å, Cα–Cβ 1.513 Å | none | -9.7 (-2.3) | +29.2 (+7.0) |
| TS |  | O–H 1.295 Å, Cβ–H 1.348 Å, Cα–Br 2.258 Å, Cα–Cβ 1.441 Å | 1 | +50.0 (+12.0) | +77.6 (+18.5) |
| product complex (`product_complex_irc`) | product end of the path | O–H 0.969 Å, Cβ–H 2.544 Å, Cα–Br 4.593 Å, Cα–Cβ 1.337 Å | none | -109.9 (-26.3) | -79.4 (-19.0) |
| propene (`C=CC`, job `propene`) | separated product | C–C 1.335–1.497 Å | none | — | — |
| ethanol (`CCO`, job `ethanol`) | separated product | C–C 1.516 Å, C–O 1.421 Å | none | — | — |
| bromide (`[Br-]`, job `bromide`) | separated product | single atom | none | — | — |
| propene + ethanol + bromide | separated products | | | -92.6 (-22.1) | -119.3 (-28.5) |

The products are 3 separate molecules (propene + ethanol + bromide); the product complex holds them together at the end of the path. Reaction energies use the separated products.
Conformer kept on record, not the zero: ethanol_trans (`CCO`) lies -0.08 kJ/mol in E and +0.90 kJ/mol in G from ethanol (SI entry 9).
The reactant end of the path ended 0.203 Å (RMSD) from the reactant-end minimum, the product end 0.402 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = separated reactants (ethoxide + 2-bromopropane), 298.15 K, ΔG at 1 M:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 50.0 (12.0) | 39.9 (9.5) | 77.6 (18.5) |
| M06-2X/ma-def2-TZVP | 56.8 (13.6) | 46.7 (11.2) | 84.4 (20.2) |
| B3LYP/ma-def2-TZVP | 48.6 (11.6) | 38.5 (9.2) | 76.2 (18.2) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 72.1 (17.2) | 62.0 (14.8) | 99.6 (23.8) |

**Barriers, zero = reactant complex, 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 59.7 (14.3) | 44.7 (10.7) | 48.4 (11.6) |
| M06-2X/ma-def2-TZVP | 66.4 (15.9) | 51.4 (12.3) | 55.1 (13.2) |
| B3LYP/ma-def2-TZVP | 54.3 (13.0) | 39.3 (9.4) | 43.0 (10.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 79.0 (18.9) | 64.0 (15.3) | 67.7 (16.2) |

**Reaction, zero = separated reactants; products = separated propene + ethanol + bromide, 298.15 K, ΔG at 1 M:**

| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from the product complex, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | -92.6 (-22.1) | -90.5 (-21.6) | -119.3 (-28.5) | -9.7 (-2.3) / +29.2 (+7.0) | 159.9 (38.2) / 157.0 (37.5) |
| M06-2X/ma-def2-TZVP | -88.4 (-21.1) | -86.3 (-20.6) | -115.1 (-27.5) | -9.6 (-2.3) / +29.3 (+7.0) | 161.2 (38.5) / 158.3 (37.8) |
| B3LYP/ma-def2-TZVP | -88.4 (-21.1) | -86.3 (-20.6) | -115.1 (-27.5) | -5.6 (-1.3) / +33.2 (+7.9) | 151.8 (36.3) / 148.9 (35.6) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -81.5 (-19.5) | -79.5 (-19.0) | -108.3 (-25.9) | -6.9 (-1.6) / +31.9 (+7.6) | 166.9 (39.9) / 164.0 (39.2) |

ΔrH was not asked for; it is given because it costs nothing. The reactant complex lies below separated reactants in E but above them in G at 1 M at every level: the binding energy does not pay for the entropy lost on association.

## 4. Highest point of the path at each level

Path energies at each level, zero = separated reactants (ethoxide + 2-bromopropane) at that level. The TS is frame 132 (s = 0, O–H 1.295 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its O–H (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 532 of 532 | 50.0 (12.0) | 132 | +0.00 | 1.295 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 532 of 532 | 56.8 (13.6) | 133 | +0.05 | 1.269 | +0.05 (+0.01) |
| B3LYP/ma-def2-TZVP | 532 of 532 | 48.6 (11.6) | 134 | +0.10 | 1.243 | +0.90 (+0.22) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 507 of 532 | 72.1 (17.2) | 134 | +0.10 | 1.243 | +1.07 (+0.26) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter O–H distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0122 | 0.035977 |
| reactant_complex | 0.0100 | 0.060172 |
| product_complex | 0.0084 | 0.052191 |
| ethoxide | 0.0123 | 0.033339 |
| 2-bromopropane | 0.0082 | 0.060280 |
| propene | 0.0104 | 0.056331 |
| ethanol | 0.0094 | 0.033050 |
| bromide | 0.0041 | 0.031531 |
| ethanol_trans | 0.0093 | 0.032436 |

- **Path frames (507 computed):** T1 from 0.0084 to 0.0122, median 0.0085; largest at frame 132 (s = +0.00, O–H 1.295 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: passes** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 132 | 0 | 0.28° / 0.65° (20) | 0.69° / 7.29° (111) | yes |
| product side (forward) | 399 | 0 | 0.24° / 0.95° (20) | 1.06° / 1.53° (378) | yes |

- **End checks: decision** (`02_irc/true_irc/end_check.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame +0.899 kJ/mol from `reactant_complex` (RMSD 0.315 Å). Opt Freq from the last frame: E -0.267 kJ/mol from it, heavy-atom RMSD 0.226 Å (all-atom 0.277 Å), largest key-distance difference -0.065 (O–H) Å, 0 imaginary frequencies: fails: energy, heavy-atom RMSD, key distances.
  - Product side (forward): last frame +4.730 kJ/mol from `product_complex` (RMSD 0.526 Å). Opt Freq from the last frame: E -0.279 kJ/mol from it, heavy-atom RMSD 0.150 Å (all-atom 0.173 Å), largest key-distance difference +0.101 (Cβ–H) Å, 0 imaginary frequencies: fails: energy, key distances.
  - Official end structures now: reactant end `reactant_complex_irc`, product end `product_complex_irc` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: 1.17 kJ/mol; product side: 5.01 kJ/mol (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 3 (decision of J. Picione, 2026-11-16, `02_irc/true_irc/end_check_decision.md`): on each side Opt Freq from the last frame reached a slightly different form of the hand-built complex (covalent key distances within 0.001 Å; reactant side ethoxide O···Hβ 2.623 vs 2.688 Å, product side Cβ···H(ethanol) 2.544 vs 2.443 Å; E lower by 0.267 and 0.279 kJ/mol). The reached minima are the official end structures (`reactant_complex_irc`, `product_complex_irc`); the hand-built complexes stay on record. Barriers from separated reactants are unaffected.
  - Entry 4: the same rule for every path.
  - Entry 5: the side-end row above.
  - Entry 9 (lowest-free-energy conformers): trans ethanol is 0.90 kJ/mol (0.22 kcal/mol) above gauche in G, so gauche ethanol stays the reference; nothing changes.

## 7. Experiment

Experiment: 2-bromopropane + 1 N NaOEt in dry ethanol, 25.0 °C (M. L. Dhar, E. D. Hughes and C. K. Ingold, J. Chem. Soc. 1948, 2058–2065, doi:10.1039/JR9480002058; Mechanism of elimination reactions, Part X). Table I: overall second-order 10⁵k₂ = 0.295 L mol⁻¹ s⁻¹ (k₂ = 2.95 × 10⁻⁶), corrected for solvent expansion, with 80.3 % olefin (two determinations, 82.7 and 78.0 %), the rest SN2 ether. Table VII's 10⁶k₂ = 2.88 is the uncorrected mean of a single run and is not used. Table II: 10⁵k(E2) = 0.237 L mol⁻¹ s⁻¹, i.e. k(E2) = 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹. Caveats stated in the paper: the first-order (solvolysis) processes contribute up to 2–3 % of the rate, and the olefin proportions are good to within 3 %, with a correction for loss of the volatile propene. Values checked against the original paper by J. Picione (2026-11-13; again 2026-11-20, with the table details above). See also E. D. Hughes, C. K. Ingold, S. Masterman and B. J. McNulty, J. Chem. Soc. 1940, 899–912, doi:10.1039/JR9400000899 (Part V), and M. L. Dhar, E. D. Hughes, C. K. Ingold, A. M. M. Mandour, G. A. Maw and L. I. Woolf, J. Chem. Soc. 1948, 2093, doi:10.1039/JR9480002093 (Part XVI). Eyring at 298.15 K (1 M standard state): k(E2) = 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹ gives ΔG‡(E2) = 105.1 kJ/mol (25.1 kcal/mol); the SN2 share (0.58 × 10⁻⁶ L mol⁻¹ s⁻¹) gives ΔG‡(SN2) = 108.6 kJ/mol (26.0 kcal/mol). The two pathways differ by only about 3.5 kJ/mol (0.8 kcal/mol), which is why a secondary bromide with ethoxide gives mostly, not only, alkene.

**Computed against measured** (computed − measured; zero = separated reactants (ethoxide + 2-bromopropane), 298.15 K, 1 M):

| Level | ΔG‡ |
|---|---|
| Measured | 105.1 (25.1) |
| M06-2X/ma-def2-SVP (working level, geometries) | 77.6 (18.5), -27.5 (-6.6) |
| M06-2X/ma-def2-TZVP | 84.4 (20.2), -20.7 (-5.0) |
| B3LYP/ma-def2-TZVP | 76.2 (18.2), -28.9 (-6.9) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 99.6 (23.8), -5.5 (-1.3) |

Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):

- The single-point curves peak just after the working-level transition state because the single points sit on a path optimized at another level: M06-2X/ma-def2-TZVP one frame after it (frame 133, s = +0.05 bohr·amu½, 0.05 kJ/mol (0.01 kcal/mol) above the TS frame) and DLPNO-CCSD(T)/ma-def2-TZVP two frames after it (frame 134, s = +0.10 bohr·amu½, O–H 1.243 Å, 1.1 kJ/mol (0.3 kcal/mol) above the TS frame). (Corrected 2026-11-25 from the final path; the earlier text, frame 28 and 0.6 kJ/mol with the coupled-cluster maximum on the TS frame, came from an earlier folder.)
- Only the anti-periplanar E2 path is computed; the competing SN2 path (about 20 % of the product in experiment) is not shown.

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 532 | 9.8 | 65 | 72 | 129 (frame 45) | 2026-11-16 20:13 | 2026-11-17 15:20 |
| B3LYP/ma-def2-TZVP | 532 | 6.0 | 42 | 47 | 81 (frame 358) | 2026-11-17 15:28 | 2026-11-21 01:59 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 507 | 149.6 | 1061 | 1295 | 1406 (frame 137) | 2026-11-21 02:10 | 2026-11-25 06:11 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 337 gradient jobs, 5.1 h job time (folders `backward`); product side (forward): 900 gradient jobs, 12.8 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.

## 9. Sensitivity check: SMD18 bromine radius (J. Picione, 2026-11-20)

Sensitivity check only: every other number in this summary stays on standard SMD (Br radius 3.06 Å). SMD18 (Engelage, Mohamed and Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652) changed the SMD Coulomb radius of Br from 3.06 to 2.60 Å. The 5 bromine-containing structures (TS, Reactant complex, Product complex, 2-bromopropane, Bromide) were recomputed at all four levels with only that radius changed: single points on the existing M06-2X/ma-def2-SVP geometries, SMD ethanol, working-level thermal corrections. Species without bromine do not change. kJ/mol (kcal/mol).

**Radius check.** In all 20 new outputs ORCA 6.1.1 printed "Radius for Br used is 4.9133 Bohr (= 2.6000 Ang.)", against 3.0600 Å in the standard outputs. The SMD-CDS term did not change (largest change 0.000 kJ/mol): only the electrostatic cavity uses this radius. Bromide's cavity area went from 420.2 to 303.4 Å². At the working level the standard single points reproduce the stored optimization energies within 0.017 kJ/mol.

**Energy change with the 2.60 Å radius, E(2.60 Å) − E(3.06 Å):**

| Structure | M06-2X/ma-def2-SVP (working level) | M06-2X/ma-def2-TZVP | B3LYP/ma-def2-TZVP | DLPNO-CCSD(T)/ma-def2-TZVP |
|---|---|---|---|---|
| TS | -9.4 (-2.2) | -8.8 (-2.1) | -8.6 (-2.0) | -8.2 (-2.0) |
| Reactant complex | -2.9 (-0.7) | -2.6 (-0.6) | -2.5 (-0.6) | -2.4 (-0.6) |
| Product complex | -29.8 (-7.1) | -28.8 (-6.9) | -27.8 (-6.7) | -29.0 (-6.9) |
| 2-bromopropane | -2.8 (-0.7) | -2.6 (-0.6) | -2.5 (-0.6) | -2.4 (-0.6) |
| Bromide | -32.3 (-7.7) | -31.2 (-7.5) | -30.1 (-7.2) | -31.3 (-7.5) |

**ΔG‡ from separated reactants (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 77.6 (18.5) | 71.0 (17.0) | -6.6 (-1.6) |
| M06-2X/ma-def2-TZVP | 84.4 (20.2) | 78.1 (18.7) | -6.3 (-1.5) |
| B3LYP/ma-def2-TZVP | 76.2 (18.2) | 70.1 (16.8) | -6.1 (-1.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 99.6 (23.8) | 93.8 (22.4) | -5.8 (-1.4) |
| Measured | 105.1 (25.1) | | |

**ΔG‡ from the reactant complex (298.15 K):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 48.4 (11.6) | 41.9 (10.0) | -6.5 (-1.6) |
| M06-2X/ma-def2-TZVP | 55.1 (13.2) | 48.9 (11.7) | -6.2 (-1.5) |
| B3LYP/ma-def2-TZVP | 43.0 (10.3) | 37.0 (8.8) | -6.0 (-1.4) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 67.7 (16.2) | 61.9 (14.8) | -5.8 (-1.4) |

**Reaction free energy ΔrG (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | -119.3 (-28.5) | -148.8 (-35.6) | -29.5 (-7.1) |
| M06-2X/ma-def2-TZVP | -115.1 (-27.5) | -143.7 (-34.4) | -28.7 (-6.8) |
| B3LYP/ma-def2-TZVP | -115.1 (-27.5) | -142.7 (-34.1) | -27.6 (-6.6) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -108.3 (-25.9) | -137.2 (-32.8) | -28.9 (-6.9) |

**Findings.**

- **Barrier.** With the SMD18 radius every barrier drops by 5.8–6.6 kJ/mol (1.4–1.6 kcal/mol), a little less than in the SN2 check (7.3–8.0 kJ/mol). The TS gains about 26 % of what free bromide gains (−8.2 against −31.3 kJ/mol at DLPNO-CCSD(T)), while 2-bromopropane and the reactant complex gain only 2.4 kJ/mol. DLPNO-CCSD(T)/ma-def2-TZVP ΔG‡ from separated reactants goes from 99.6 to 93.8 kJ/mol, against 105.1 measured: from 5.5 kJ/mol below to 11.3 kJ/mol below. For this reaction the new radius moves the barrier away from experiment.
- **Reaction free energy.** ΔrG drops by 27.6–29.5 kJ/mol (6.6–7.1 kcal/mol), almost all of it from bromide. There is no measured ΔrG for this reaction to compare with.
- **SN2 and E2 together.** The radius moves the two barriers by similar amounts in the same direction. It brings SN2 closer to experiment (+5.5 → −1.8 kJ/mol) and E2 further away (−5.5 → −11.3 kJ/mol), so it does not improve both at once.
- **Geometries** are the standard-SMD ones; no structure was reoptimized with the new radius.
