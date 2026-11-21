# Results summary: Diels–Alder, butadiene + ethylene

`CH₂=CH–CH=CH₂ + CH₂=CH₂ → cyclohexene`, in toluene (SMD), 298.15 K. Reaction folder `reactions/da_butadiene_ethylene_m062x`. Written 2026-11-21 05:11 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "Separated reactants" means butadiene_s-trans + ethylene, each optimized on its own (each the conformer lowest in free energy (SI entry 9)). Steps that change the number of molecules use the 1 M standard state (RT ln 24.47 = 7.93 kJ/mol per mole of solute: subtracted for each molecule lost, added for each gained); steps from the reactant complex involve no change in the number of molecules, so no standard-state term applies. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 302 (reactant side 179 frames, TS = frame 179, product side 122 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 302 of 302; M06-2X/ma-def2-TZVP 302 of 302; B3LYP/ma-def2-TZVP 302 of 302; DLPNO-CCSD(T)/ma-def2-TZVP 275 of 302. DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 14.91 bohr·amu½ in the mass-weighted coordinate s (from s = -8.94 at the reactant end to s = 5.97 at the product end; TS at s = 0); 3.49 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (C1–C6 or C4–C5 between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 20 (s = -7.95, C1–C6 3.31 Å) to frame 264 (s = 4.25, 1.64 Å). At the TS the C1–C6 distance is 2.263 Å, 61 % of the way from the reactant end (3.401 Å) to the product end (1.542 Å): most of the shortening has happened by the TS, and the remaining 39 % (0.72 Å) comes after it, over the product side (s = 0 to 5.97).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| C1–C6 | forms | 3.401 | 2.263 | 1.542 | 61 % |
| C4–C5 | forms | 3.401 | 2.263 | 1.542 | 61 % |
| C2–C3 | shortens | 1.473 | 1.413 | 1.339 | 45 % |
| C1–C2 | lengthens | 1.339 | 1.378 | 1.505 | 23 % |
| C3–C4 | lengthens | 1.339 | 1.378 | 1.505 | 23 % |
| C5–C6 | lengthens | 1.332 | 1.381 | 1.547 | 23 % |

| Point | Frame | s (bohr·amu½) | C1–C6 (Å) | C4–C5 (Å) | C2–C3 (Å) | C1–C2 (Å) | C3–C4 (Å) | C5–C6 (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |
|---|---|---|---|---|---|---|---|---|---|
| reactant end | 0 | -8.94 | 3.401 | 3.401 | 1.473 | 1.339 | 1.339 | 1.332 | +3.3 (+0.8) |
| s ≈ -6 | 59 | -6.00 | 3.084 | 3.084 | 1.470 | 1.340 | 1.340 | 1.332 | +7.9 (+1.9) |
| s ≈ -4 | 99 | -4.00 | 2.829 | 2.829 | 1.466 | 1.341 | 1.341 | 1.333 | +19.9 (+4.8) |
| s ≈ -2 | 139 | -2.00 | 2.556 | 2.556 | 1.453 | 1.346 | 1.346 | 1.339 | +48.5 (+11.6) |
| s ≈ -1 | 159 | -1.00 | 2.411 | 2.411 | 1.437 | 1.357 | 1.357 | 1.352 | +68.2 (+16.3) |
| s ≈ -0.5 | 169 | -0.50 | 2.337 | 2.337 | 1.425 | 1.366 | 1.366 | 1.365 | +75.5 (+18.0) |
| TS | 179 | +0.00 | 2.263 | 2.263 | 1.413 | 1.378 | 1.378 | 1.381 | +78.4 (+18.7) |
| s ≈ +0.5 | 189 | +0.50 | 2.189 | 2.189 | 1.399 | 1.393 | 1.393 | 1.401 | +74.6 (+17.8) |
| s ≈ +1 | 199 | +1.00 | 2.116 | 2.116 | 1.385 | 1.410 | 1.410 | 1.424 | +61.2 (+14.6) |
| s ≈ +1.5 | 209 | +1.50 | 2.042 | 2.042 | 1.373 | 1.428 | 1.428 | 1.450 | +36.1 (+8.6) |
| s ≈ +2 | 219 | +2.00 | 1.968 | 1.968 | 1.362 | 1.445 | 1.445 | 1.475 | +1.2 (+0.3) |
| s ≈ +3 | 239 | +3.00 | 1.820 | 1.820 | 1.349 | 1.476 | 1.476 | 1.519 | -81.0 (-19.4) |
| product end | 301 | +5.97 | 1.542 | 1.542 | 1.339 | 1.505 | 1.505 | 1.547 | -195.1 (-46.6) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD toluene)

- **TS:** 1 imaginary frequency, -511.8 cm⁻¹ (written 511.8i cm⁻¹). Stretch rates of the key bonds in the mode: C1–C6 +0.82, C4–C5 +0.82, C2–C3 +0.15, C1–C2 -0.15, C3–C4 -0.15, C5–C6 -0.20; forming and breaking bonds move together. C1–C6 2.263 Å, C4–C5 2.263 Å, C2–C3 1.413 Å, C1–C2 1.378 Å, C3–C4 1.378 Å, C5–C6 1.381 Å, C1=C2–C3=C4 dihedral 0.0°.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |
|---|---|---|---|---|---|
| butadiene_s-trans (`C=CC=C`, job `butadiene_s-trans`) | separated reactant | C–C 1.339–1.463 Å | none | — | — |
| ethylene (`C=C`, job `ethylene`) | separated reactant | C–C 1.331 Å | none | — | — |
| reactant complex (`reactant_complex_irc`) | reactant end of the path | C1–C6 3.473 Å, C4–C5 3.664 Å, C2–C3 1.473 Å, C1–C2 1.339 Å, C3–C4 1.339 Å, C5–C6 1.332 Å | none | +2.7 (+0.7) | +29.9 (+7.1) |
| TS |  | C1–C6 2.263 Å, C4–C5 2.263 Å, C2–C3 1.413 Å, C1–C2 1.378 Å, C3–C4 1.378 Å, C5–C6 1.381 Å | 1 | +78.4 (+18.7) | +126.0 (+30.1) |
| product complex | product end of the path | C1–C6 1.541 Å, C4–C5 1.542 Å, C2–C3 1.338 Å, C1–C2 1.505 Å, C3–C4 1.505 Å, C5–C6 1.547 Å | none | -195.1 (-46.6) | -130.4 (-31.2) |
| cyclohexene (`C1CCC=CC1`) | separated product | C–C 1.339–1.529 Å | none | -218.1 (-52.1) | -151.7 (-36.2) |

The product is one molecule. The product complex is the form the path reaches, the separated product the same molecule optimized on its own; here they are different conformers, +22.96 kJ/mol apart in E (product complex above), so the reaction energies (separated product) and the energy of the path's product end differ by that much (card notes, section 7).
Conformer kept on record, not the zero: butadiene_s-cis (`C=CC=C`) lies +11.71 kJ/mol in E and +9.48 kJ/mol in G from butadiene_s-trans (SI entry 9).
The reactant end of the path ended 0.596 Å (RMSD) from the reactant-end minimum, the product end 0.002 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = separated reactants (butadiene_s-trans + ethylene), 298.15 K, ΔG at 1 M:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 78.4 (18.7) | 82.8 (19.8) | 126.0 (30.1) |
| M06-2X/ma-def2-TZVP | 84.2 (20.1) | 88.7 (21.2) | 131.9 (31.5) |
| B3LYP/ma-def2-TZVP | 79.8 (19.1) | 84.2 (20.1) | 127.4 (30.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 96.7 (23.1) | 101.1 (24.2) | 144.4 (34.5) |

**Barriers, zero = reactant complex, 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 75.6 (18.1) | 75.2 (18.0) | 96.2 (23.0) |
| M06-2X/ma-def2-TZVP | 80.0 (19.1) | 79.5 (19.0) | 100.5 (24.0) |
| B3LYP/ma-def2-TZVP | 72.5 (17.3) | 72.0 (17.2) | 93.0 (22.2) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 90.2 (21.6) | 89.7 (21.4) | 110.7 (26.5) |

**Reaction, zero = separated reactants; product = separated cyclohexene, 298.15 K, ΔG at 1 M:**

| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from cyclohexene, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | -218.1 (-52.1) | -199.6 (-47.7) | -151.7 (-36.2) | +2.7 (+0.7) / +29.9 (+7.1) | 296.4 (70.8) / 277.7 (66.4) |
| M06-2X/ma-def2-TZVP | -198.2 (-47.4) | -179.7 (-42.9) | -131.8 (-31.5) | +4.2 (+1.0) / +31.4 (+7.5) | 282.4 (67.5) / 263.7 (63.0) |
| B3LYP/ma-def2-TZVP | -171.3 (-41.0) | -152.8 (-36.5) | -104.9 (-25.1) | +7.2 (+1.7) / +34.4 (+8.2) | 251.1 (60.0) / 232.4 (55.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -189.7 (-45.3) | -171.2 (-40.9) | -123.3 (-29.5) | +6.5 (+1.6) / +33.7 (+8.0) | 286.3 (68.4) / 267.6 (64.0) |

ΔrH was not asked for; it is given because it costs nothing. The reactant complex lies above in E and above in G at M06-2X/ma-def2-SVP; above in E and above in G at M06-2X/ma-def2-TZVP; above in E and above in G at B3LYP/ma-def2-TZVP; above in E and above in G at DLPNO-CCSD(T)/ma-def2-TZVP (separated reactants as zero, 1 M).

## 4. Highest point of the path at each level

Path energies at each level, zero = separated reactants (butadiene_s-trans + ethylene) at that level. The TS is frame 179 (s = 0, C1–C6 2.263 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C1–C6 (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 302 of 302 | 78.4 (18.7) | 179 | +0.00 | 2.263 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 302 of 302 | 84.2 (20.1) | 182 | +0.15 | 2.241 | +0.40 (+0.10) |
| B3LYP/ma-def2-TZVP | 302 of 302 | 79.8 (19.1) | 183 | +0.20 | 2.234 | +0.44 (+0.11) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 275 of 302 | 96.7 (23.1) | 183 | +0.20 | 2.234 | +0.82 (+0.19) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter C1–C6 distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0104 | 0.060839 |
| reactant_complex | 0.0109 | 0.061350 |
| product_complex | 0.0097 | 0.058802 |
| butadiene_s-trans | 0.0113 | 0.042107 |
| ethylene | 0.0111 | 0.060250 |
| cyclohexene | 0.0097 | 0.042666 |
| butadiene_s-cis | 0.0111 | 0.046809 |

- **Path frames (275 computed):** T1 from 0.0097 to 0.0108, median 0.0107; largest at frame 0 (s = -8.94, C1–C6 3.401 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: passes** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 179 | 0 | 0.03° / 1.65° (20) | 0.04° / 2.53° (158) | yes |
| product side (forward) | 122 | 2 | 0.03° / 1.56° (20) | 0.01° / 4.57° (99) | yes |

- **End checks: ok** (`02_irc/true_irc/end_check_all.md`, `02_irc/true_irc/end_check.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame +0.327 kJ/mol from `reactant_complex` (RMSD 0.178 Å). Opt Freq from the last frame: E -0.297 kJ/mol from it, heavy-atom RMSD 0.250 Å (all-atom 0.518 Å), largest key-distance difference +0.185 (C4–C5) Å, 0 imaginary frequencies: a different form of reactant_complex with the same bonds (fails: energy, heavy-atom RMSD, key distances); now the official end structure reactant_complex_irc.
  - Product side (forward): last frame +0.011 kJ/mol from `product_complex` (RMSD 0.002 Å). Opt Freq from the last frame: E +0.004 kJ/mol from it, heavy-atom RMSD 0.001 Å (all-atom 0.002 Å), largest key-distance difference -0.000 (C1–C2) Å, 0 imaginary frequencies: reaches the same minimum.
  - Official end structures now: reactant end `reactant_complex_irc`, product end `product_complex` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: 0.33 kJ/mol; product side: 0.01 kJ/mol (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 1 (stopping rule near a flat TS): the 20 flat steps began at least 10 kJ/mol below the TS, so the amended rule stops the reactant side at the same frame; not rerun. The product side ended by the energy-rise stop and was accepted by J. Picione (2026-11-14); the rule change does not apply to it.
  - Entry 4 (the minimum the path reaches is the official end): on the reactant side the reached form (Opt, then TightOpt after an imaginary mode) has the C4···C5 contact 0.185 Å longer, heavy-atom RMSD 0.250 Å, same bonds, E −0.297 kJ/mol and G +2.245 kJ/mol from the earlier reference; it is now the official reactant end (`reactant_complex_irc`). The single points were run on it from the start. Product side: same minimum.
  - Entry 5: the side-end row above.
  - Entry 6 (angle against gradient size): 5 of the 70 frames over 15° in all six paths are on this path, in the final approach to a minimum, where the gradient is nearly zero. Near this TS the angles stay small.
  - Entry 9 (lowest-free-energy conformers): s-trans butadiene is 9.48 kJ/mol (2.27 kcal/mol) lower in G than the gauche form the path reaches, so it is the zero; barriers from separated reactants include that gap. Barriers from the reactant complex are unaffected.

## 7. Experiment

Experiment (gas phase): C₂H₄ + 1,3-butadiene → cyclohexene, second order, 760–921 K: k = 3.0 × 10¹⁰ exp(−27,500 cal mol⁻¹/RT) cm³ mol⁻¹ s⁻¹ (eq. 5), i.e. A = 3.0 × 10⁷ L mol⁻¹ s⁻¹ (4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹); a refit of the 13 runs in Table VI gives Ea = 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol). Source: D. Rowley and H. Steiner, Discuss. Faraday Soc. 1951, 10, 198–213, doi:10.1039/DF9511000198; confirmed from the original paper by J. Picione (2026-11-13). Cross-check: NIST Chemical Kinetics Database, squib 1951ROW/STE198-213 (same A and Ea). The measured activation energy refers to equilibrium butadiene in the gas phase at 760–921 K. That is almost all s-trans at room temperature (about 2 % gauche s-cis) but not at the temperatures of the measurement: about 27 % gauche at 760 K and 34 % at 921 K (DLPNO-CCSD(T)/ma-def2-TZVP, gas phase; corrected 2026-11-21, J. Picione). This calculation is in SMD toluene and starts from s-cis butadiene, so the s-trans → s-cis energy is reported separately; the comparison with this measurement uses separate gas-phase calculations, from s-trans butadiene alone and from the equilibrium conformer mixture, each with ethylene. Concerted, synchronous transition state: K. N. Houk, Y.-T. Lin and F. K. Brown, J. Am. Chem. Soc. 1986, 108, 554–556, doi:10.1021/ja00263a059; E. Goldstein, B. Beno and K. N. Houk, J. Am. Chem. Soc. 1996, 118, 6036–6043, doi:10.1021/ja9601494.

Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):

- The path ends at the boat form of cyclohexene (mirror-symmetric, the new CH₂–CH₂ bond eclipsed): a synchronous path keeps the mirror plane of the TS. Optimized from the IRC end it stays in the boat form (lowest vibration 47 cm⁻¹, a very flat ring motion), 23.0 kJ/mol (5.5 kcal/mol) above the half-chair in E (M06-2X/ma-def2-SVP; numbers rechecked for this folder 2026-11-20), the ground-state conformer. Reaction energies use half-chair cyclohexene; the boat is the end of the path shown in the viewer and video.
- Gas-phase comparison with Rowley and Steiner (760–921 K), on gas-phase geometries and frequencies (gas-phase section of the card and of the results summary). With the equilibrium butadiene conformer mixture as the reactant (s-trans plus both gauche forms, about 30 % gauche at 840 K; computed 2026-11-21) and an ln k against 1/T fit over 760–921 K as in the experiment, DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 114.9 kJ/mol (27.5 kcal/mol), 0.2 kJ/mol (0.05 kcal/mol) below the measured 115.1 ± 3.4 kJ/mol, and M06-2X/ma-def2-TZVP 13.0 kJ/mol (3.1 kcal/mol) below. The computed pre-exponential factor, 5.1 × 10⁶ L mol⁻¹ s⁻¹, is 0.17 times the measured one (0.17–0.19 at every level), so DLPNO-CCSD(T) gives 0.17 times the measured rate at 840 K; M06-2X/ma-def2-TZVP gives 1.08 times only because its low Ea cancels the low A. From s-trans butadiene alone (2026-11-20, 840 K point values) the numbers were +4.2 kJ/mol (DLPNO-CCSD(T)), −8.6 kJ/mol (M06-2X/ma-def2-TZVP) and A 0.45 times the measured one; the earlier B3LYP-geometry folder gave −8.8 and +4.7 kJ/mol and an A ratio of 0.57. SI corrections entry 11.
- The product side of the path ended at the minimum (0.011 kJ/mol above the optimized cyclohexene, RMS gradient 1.0 × 10⁻⁶ Eh/bohr) on the old energy-rise stop, before the 20 flat steps of the stopping rule; not rerun (accepted by J. Picione, 2026-11-14).

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 302 | 3.8 | 45 | 46 | 51 (frame 254) | 2026-11-16 22:39 | 2026-11-17 09:31 |
| B3LYP/ma-def2-TZVP | 302 | 1.8 | 21 | 23 | 28 (frame 270) | 2026-11-17 09:33 | 2026-11-17 23:11 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 275 | 61.2 | 698 | 1192 | 1759 (frame 210) | 2026-11-17 23:20 | 2026-11-20 17:16 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 432 gradient jobs, 7.0 h job time (folders `backward`); product side (forward): 324 gradient jobs, 5.4 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.

## 9. Gas-phase comparison with Rowley and Steiner (840 K and the 760–921 K fit)

Gas-phase comparison with Rowley and Steiner (J. Picione, 2026-11-20). Gas phase, no solvent model, gas-phase geometries and frequencies at M06-2X-D3(0)/ma-def2-SVP; barrier from s-trans butadiene + ethylene; T = 840 K (midpoint of the measured 760–921 K).

| Level | ΔE‡ | ΔH‡ (840 K) | Ea = ΔH‡ + 2RT (840 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | k(TST, 840 K) (L mol⁻¹ s⁻¹) | k / measured k |
|---|---|---|---|---|---|---|---|---|
| M06-2X_ma-def2-SVP (working) | 80.0 kJ/mol (19.1 kcal/mol) | 86.9 kJ/mol (20.8 kcal/mol) | 100.9 kJ/mol (24.1 kcal/mol) | −14.2 kJ/mol (−3.4 kcal/mol) | 1.36 × 10⁷ | 0.45 | 7.26 | 3.48 |
| M06-2X_ma-def2-TZVP | 85.6 kJ/mol (20.5 kcal/mol) | 92.5 kJ/mol (22.1 kcal/mol) | 106.5 kJ/mol (25.4 kcal/mol) | −8.6 kJ/mol (−2.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 3.27 | 1.57 |
| B3LYP_ma-def2-TZVP | 81.3 kJ/mol (19.4 kcal/mol) | 88.2 kJ/mol (21.1 kcal/mol) | 102.2 kJ/mol (24.4 kcal/mol) | −12.9 kJ/mol (−3.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 6.01 | 2.88 |
| DLPNO-CCSD(T)_ma-def2-TZVP | 98.4 kJ/mol (23.5 kcal/mol) | 105.3 kJ/mol (25.2 kcal/mol) | 119.3 kJ/mol (28.5 kcal/mol) | +4.2 kJ/mol (+1.0 kcal/mol) | 1.36 × 10⁷ | 0.45 | 0.524 | 0.25 |
| Measured (Rowley and Steiner 1951) | | | 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol) | | 3.00 × 10⁷ | 1 | 2.09 | 1 |

**Units and standard state of A.** A and k(TST) are second-order rate constants in concentration units: L mol⁻¹ s⁻¹ (= 10³ cm³ mol⁻¹ s⁻¹, Rowley and Steiner's unit), also given per molecule (cm³ molecule⁻¹ s⁻¹). ΔG‡ and ΔS‡ use the 1 atm standard state; the factor RT/p° turns that into concentration units, so A and k do not depend on the choice of standard state and compare directly with the measured A. Measured A = 3.00 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 3.00 × 10⁷ L mol⁻¹ s⁻¹ = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹; computed A = 1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 2.27 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹.

A is the same at every level: all levels share the working-level gas-phase frequencies, so ΔS‡ (840 K, 1 atm) = −168.8 J mol⁻¹ K⁻¹ at every level; the levels differ only in Ea.

Method: ΔH‡(840 K) = ΔE‡ + Δ[H(840 K) − E]; Ea = ΔH‡ + 2RT; ΔG‡(840 K, 1 atm) = ΔE‡ + Δ[G(840 K) − E] (quasi-RRHO entropy); k(TST) = (kB T/h)(RT/p°) exp(−ΔG‡/RT), κ = 1; A = k(TST) exp(Ea/RT). Thermal terms from the gas-phase working-level frequencies (rotational symmetry numbers from ORCA's point-group detection: ts σ = 1, butadiene_s-trans σ = 2, ethylene σ = 4). Same method as the earlier B3LYP-geometry folder, except that geometries and frequencies are now gas phase.

**Gas-phase TS:** one imaginary frequency, −514 cm⁻¹; overlap with the toluene TS mode 1.00; forming bonds C1–C6 2.265 Å and C4–C5 2.265 Å (spread 0.000 Å, stretch-rate ratio 1.00): synchronous, accepted by the folder's TS rules. Toluene TS: 2.263 Å.

**Butadiene conformers (gas phase, working level):** the gauche s-cis form lies 11.6 kJ/mol (2.8 kcal/mol) in E and 9.2 kJ/mol (2.2 kcal/mol) in G (298.15 K) above s-trans; the lowest conformer is s-trans.

**Findings (2026-11-20 20:52).**

- **Activation energy.** DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 119.3 kJ/mol (28.5 kcal/mol), 4.2 kJ/mol (1.0 kcal/mol) above the measured 115.1 ± 3.4 kJ/mol, just outside the fit error. M06-2X/ma-def2-TZVP is 8.6 kJ/mol (2.1 kcal/mol) below, B3LYP/ma-def2-TZVP 12.9 kJ/mol (3.1 kcal/mol) below, and the working level 14.2 kJ/mol (3.4 kcal/mol) below.
- **Pre-exponential factor.** A = 1.36 × 10⁷ L mol⁻¹ s⁻¹ (1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹), 0.45 times the measured 3.0 × 10⁷ L mol⁻¹ s⁻¹, at every level (ΔS‡ = −168.8 J mol⁻¹ K⁻¹ at 840 K, 1 atm). A and Ea from one Arrhenius fit over 760–921 K are correlated: the ±3.4 kJ/mol in Ea corresponds to about a factor of 1.6 in A at 840 K, so the computed A is low by somewhat more than the fit uncertainty.
- **Rate at 840 K.** DLPNO-CCSD(T) gives 0.25 times the measured k: 0.45 from A times 0.55 from the higher Ea. M06-2X/ma-def2-TZVP gives 1.57 times: 0.45 from A times 3.4 from its low Ea. As in the earlier folder, M06-2X matches the rate better only because its low Ea cancels the low A.
- **Against the withdrawn numbers** from the earlier B3LYP-geometry folder (solution geometries and SMD frequencies, gas-phase single points): Ea difference −8.8 → −8.6 kJ/mol (M06-2X/ma-def2-TZVP), +4.7 → +4.2 kJ/mol (DLPNO-CCSD(T)); A ratio 0.57 → 0.45.
- **Gas phase against toluene** (same levels, s-trans butadiene + ethylene): the DLPNO-CCSD(T) ΔE‡ is 98.4 kJ/mol in the gas phase and 96.7 kJ/mol in SMD toluene; the forming C–C bonds are 2.265 Å and 2.263 Å.
- **Checks.** Every job ended normally. The minima have no imaginary frequency. The TS has one (−514 cm⁻¹), synchronous, with overlap 1.00 with the toluene TS mode. The thermal terms reproduce ORCA's 298.15 K H and G within 0.001 kJ/mol. The largest DLPNO-CCSD(T) T1 is 0.0113.
- **Job time:** 2092 s in total (4 optimizations 820 s, 9 single points 1272 s, 4 cores each), 20:24–20:52, ahead of the path single points.

**Butadiene conformer ensemble (J. Picione, 2026-11-21).** At 760–921 K the reactant is the equilibrium mixture of s-trans butadiene and the two mirror-image gauche (s-cis) forms, as in the experiment. k(TST) was computed at 760, 800, 840, 880 and 921 K and ln k fitted against 1/T, as Rowley and Steiner did; Ea = −R × slope, A = exp(intercept), in L mol⁻¹ s⁻¹. The s-trans-only values below use the same fit, so they differ slightly from the 840 K point values (Ea = ΔH‡ + 2RT) in the table above.

| Level | Reactant butadiene | Ea (fit, 760–921 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | k(840 K) / measured k |
|---|---|---|---|---|---|---|
| M06-2X_ma-def2-SVP (working) | s-trans only | 100.2 kJ/mol (24.0 kcal/mol) | −14.9 kJ/mol (−3.6 kcal/mol) | 1.24 × 10⁷ | 0.41 | 3.48 |
|  | conformer mixture | 96.5 kJ/mol (23.1 kcal/mol) | −18.6 kJ/mol (−4.4 kcal/mol) | 4.94 × 10⁶ | 0.16 | 2.34 |
| M06-2X_ma-def2-TZVP | s-trans only | 105.8 kJ/mol (25.3 kcal/mol) | −9.3 kJ/mol (−2.2 kcal/mol) | 1.24 × 10⁷ | 0.41 | 1.57 |
|  | conformer mixture | 102.1 kJ/mol (24.4 kcal/mol) | −13.0 kJ/mol (−3.1 kcal/mol) | 5.08 × 10⁶ | 0.17 | 1.08 |
| B3LYP_ma-def2-TZVP | s-trans only | 101.5 kJ/mol (24.3 kcal/mol) | −13.6 kJ/mol (−3.2 kcal/mol) | 1.24 × 10⁷ | 0.41 | 2.88 |
|  | conformer mixture | 98.1 kJ/mol (23.4 kcal/mol) | −17.0 kJ/mol (−4.1 kcal/mol) | 5.68 × 10⁶ | 0.19 | 2.16 |
| DLPNO-CCSD(T)_ma-def2-TZVP | s-trans only | 118.6 kJ/mol (28.3 kcal/mol) | +3.5 kJ/mol (+0.8 kcal/mol) | 1.24 × 10⁷ | 0.41 | 0.25 |
|  | conformer mixture | 114.9 kJ/mol (27.5 kcal/mol) | −0.2 kJ/mol (−0.05 kcal/mol) | 5.13 × 10⁶ | 0.17 | 0.17 |
| Measured (Rowley and Steiner 1951) | equilibrium mixture | 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol) | | 3.00 × 10⁷ | 1 | 1 |

**Gauche fraction of butadiene** (both enantiomers together; equilibrium in the gas phase):

| Level | 298 K | 760 K | 840 K | 921 K | G(gauche, one enantiomer) − G(s-trans), 840 K |
|---|---|---|---|---|---|
| M06-2X_ma-def2-SVP (working) | 2.4 % | 28.9 % | 32.6 % | 35.8 % | 9.9 kJ/mol (2.4 kcal/mol) |
| M06-2X_ma-def2-TZVP | 1.9 % | 27.2 % | 30.9 % | 34.2 % | 10.5 kJ/mol (2.5 kcal/mol) |
| B3LYP_ma-def2-TZVP | 0.9 % | 21.2 % | 24.9 % | 28.4 % | 12.5 kJ/mol (3.0 kcal/mol) |
| DLPNO-CCSD(T)_ma-def2-TZVP | 1.8 % | 26.6 % | 30.4 % | 33.7 % | 10.6 kJ/mol (2.5 kcal/mol) |

**How the symmetry numbers and the two enantiomers are counted.** G(butadiene, mixture) = −RT ln[exp(−G(s-trans)/RT) + 2 exp(−G(gauche)/RT)]. s-trans butadiene is C2h, σ = 2. Gauche butadiene is C2 (CCCC dihedral 33.5°; the optimized structure matches its own 180° rotation within 0.0002 Å RMSD), so one enantiomer has σ = 2; the factor 2 in front of its term counts the two mirror-image forms, which are separate minima with equal energy. Together the pair therefore counts like one structure with σ = 1 (2/σ = 1). ORCA's point-group detection labelled the optimized gauche form C1 (σ = 1), so its own G (and the 9.2 kJ/mol gauche–s-trans ΔG quoted above) already equals the pair value; here σ = 2 is set explicitly and the factor 2 added. Ethylene σ = 4 (D2h), TS σ = 1 (Cs, achiral, one TS).

**Findings, conformer ensemble (2026-11-21 05:09).**

- **Activation energy.** With the equilibrium conformer mixture as the reactant, DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 114.9 kJ/mol (27.5 kcal/mol) from the 760–921 K fit, 0.2 kJ/mol (0.05 kcal/mol) below the measured 115.1 ± 3.4 kJ/mol. From s-trans butadiene alone the same fit gives 118.6 kJ/mol (28.3 kcal/mol), 3.5 kJ/mol above. Every level drops by 3.4–3.7 kJ/mol (0.8–0.9 kcal/mol): at these temperatures about 30 % of the butadiene is gauche, which lies 12.2 kJ/mol (2.9 kcal/mol) above s-trans in H at 840 K (DLPNO-CCSD(T)), so the average reactant enthalpy is higher. M06-2X/ma-def2-TZVP is now 13.0 kJ/mol (3.1 kcal/mol) below the measured value, B3LYP/ma-def2-TZVP 17.0 kJ/mol (4.1 kcal/mol) below.
- **Pre-exponential factor.** A falls from 1.24 × 10⁷ (s-trans only, same fit) to 5.1 × 10⁶ L mol⁻¹ s⁻¹ at DLPNO-CCSD(T), 0.17 times the measured 3.0 × 10⁷ L mol⁻¹ s⁻¹ (0.17–0.19 at every level). Only the s-trans part reacts in this model, and its share falls with temperature, which lowers both the fitted Ea and A. The ±3.4 kJ/mol fit error in Ea corresponds to about a factor of 1.6 in A, so the computed A is low by clearly more than the measurement's uncertainty.
- **Rate at 840 K.** The mixture lowers k(TST) to 0.70 of the s-trans-only value (the s-trans fraction). DLPNO-CCSD(T) gives 0.17 times the measured k, a free-energy gap of 12.4 kJ/mol (3.0 kcal/mol) at 840 K that now comes almost entirely from the entropy of activation (computed ΔS‡ more negative than measured by about R ln 6 ≈ 15 J mol⁻¹ K⁻¹). Likely contributors, not tested here: the harmonic (quasi-RRHO) treatment of the low-frequency modes of the TS and of the butadiene torsion at 840 K.
- **What changes.** The earlier conclusion (s-trans butadiene only, 2026-11-20) was that DLPNO-CCSD(T) puts Ea 4.2 kJ/mol above experiment, just outside the error, with A 0.45 times the measured value. With the conformer mixture the coupled-cluster Ea agrees with experiment, and the remaining disagreement is in A (0.17 times). M06-2X/ma-def2-TZVP comes within 1.08 times the measured k only because its Ea, 13.0 kJ/mol too low, cancels the low A.
- **Gauche fraction** (DLPNO-CCSD(T)): 1.8 % at 298 K, 26.6 % at 760 K, 30.4 % at 840 K, 33.7 % at 921 K. B3LYP gives fewer gauche (24.9 % at 840 K), because it puts gauche 2.0 kJ/mol higher.
- **Checks.** The three gauche single points ended normally (largest DLPNO-CCSD(T) T1 0.0111). The ln k against 1/T fits are straight within 0.007 in ln k at every level. Job time: 3 single points, 05:06–05:09, right after J. Picione's machine block.
