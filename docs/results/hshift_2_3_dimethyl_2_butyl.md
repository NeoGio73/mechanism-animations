# Results summary: 1,2-Hydride shift, 2,3-dimethyl-2-butyl cation (degenerate)

`(CH₃)₂C⁺–CH(CH₃)₂ → (CH₃)₂CH–C⁺(CH₃)₂`, in dichloromethane (SMD), 298.15 K. Reaction folder `reactions/hshift_2_3_dimethyl_2_butyl_m062x`. Written 2026-11-22 17:20 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "Separated reactant" means dimethylbutyl_cation, optimized on its own (the conformer lowest in free energy (SI entry 9)). Steps that change the number of molecules use the 1 M standard state (RT ln 24.47 = 7.93 kJ/mol per mole of solute: subtracted for each molecule lost, added for each gained); steps from the reactant complex involve no change in the number of molecules, so no standard-state term applies. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 543 (reactant side 264 frames, TS = frame 264, product side 278 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 543 of 543; M06-2X/ma-def2-TZVP 543 of 543; B3LYP/ma-def2-TZVP 543 of 543 (272 computed, 271 copied from their mirror partner, SI entry 7); DLPNO-CCSD(T)/ma-def2-TZVP 489 of 543 (245 computed, 244 copied from their mirror partner, SI entry 7). DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 27.25 bohr·amu½ in the mass-weighted coordinate s (from s = -13.29 at the reactant end to s = 13.97 at the product end; TS at s = 0); 10.28 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (C2–H or C3–H between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 55 (s = -10.55, C2–H 2.02 Å) to frame 487 (s = 11.22, 1.11 Å). At the TS the C2–H distance is 1.373 Å, 72 % of the way from the reactant end (2.065 Å) to the product end (1.105 Å): most of the shortening has happened by the TS, and the remaining 28 % (0.27 Å) comes after it, over the product side (s = 0 to 13.97).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| C2–H | forms | 2.065 | 1.373 | 1.105 | 72 % |
| C3–H | breaks | 1.105 | 1.277 | 2.065 | 18 % |

| Point | Frame | s (bohr·amu½) | C2–H (Å) | C3–H (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |
|---|---|---|---|---|---|
| reactant end | 0 | -13.29 | 2.065 | 1.105 | +1.3 (+0.3) |
| s ≈ -10 | 66 | -9.99 | 2.005 | 1.114 | +2.8 (+0.7) |
| s ≈ -6 | 143 | -5.99 | 1.931 | 1.125 | +3.7 (+0.9) |
| s ≈ -4 | 183 | -3.99 | 1.906 | 1.128 | +5.5 (+1.3) |
| s ≈ -2 | 223 | -1.99 | 1.842 | 1.137 | +8.5 (+2.0) |
| s ≈ -1 | 243 | -0.99 | 1.653 | 1.169 | +13.0 (+3.1) |
| s ≈ -0.5 | 253 | -0.49 | 1.507 | 1.211 | +14.9 (+3.6) |
| TS | 264 | +0.00 | 1.373 | 1.277 | +15.3 (+3.7) |
| s ≈ +0.5 | 273 | +0.49 | 1.302 | 1.340 | +15.3 (+3.7) |
| s ≈ +1 | 286 | +1.02 | 1.229 | 1.463 | +15.1 (+3.6) |
| s ≈ +1.5 | 296 | +1.51 | 1.179 | 1.611 | +13.7 (+3.3) |
| s ≈ +2 | 306 | +2.01 | 1.153 | 1.735 | +11.3 (+2.7) |
| s ≈ +3 | 326 | +3.01 | 1.132 | 1.873 | +7.5 (+1.8) |
| s ≈ +6 | 386 | +6.01 | 1.126 | 1.925 | +4.0 (+1.0) |
| s ≈ +10 | 463 | +10.01 | 1.116 | 1.991 | +2.9 (+0.7) |
| product end | 542 | +13.97 | 1.105 | 2.065 | +1.3 (+0.3) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane)

- **TS:** 1 imaginary frequency, -115.2 cm⁻¹ (written 115.2i cm⁻¹). Stretch rates of the key bonds in the mode: C2–H -0.58, C3–H +0.39; forming and breaking bonds move together. C2–H 1.373 Å, C3–H 1.277 Å, C2–H–C3 angle 64.1°.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |
|---|---|---|---|---|---|
| dimethylbutyl_cation (`CC(C)[C+](C)C`, job `dimethylbutyl_cation`) | separated reactant | C–C 1.457–1.566 Å | none | — | — |
| reactant complex (`open_cation`) | reactant end of the path | C2–H 2.086 Å, C3–H 1.102 Å | none | -0.2 (-0.0) | +1.2 (+0.3) |
| TS |  | C2–H 1.373 Å, C3–H 1.277 Å | 1 | +15.3 (+3.7) | +14.1 (+3.4) |
| product complex (`open_cation`) | product end of the path | C2–H 2.086 Å, C3–H 1.102 Å | none | -0.2 (-0.0) | +1.2 (+0.3) |
| dimethylbutyl_cation_shifted (`CC(C)[C+](C)C`) | separated product | C–C 1.457–1.566 Å | none | +0.0 (+0.0) | +0.0 (+0.0) |

The product is one molecule, so the product complex and the separated product are the same species optimized twice (from the path end and on its own); they differ by -0.20 kJ/mol in E. Reaction energies use the separated product.
The reactant is one molecule, so the reactant complex and the separated reactant are the same species optimized twice (from the path end and on its own); they differ by -0.20 kJ/mol in E.
Conformer kept on record, not the zero: rotamer_h_perp (`CC(C)[C+](C)C`) lies +0.07 kJ/mol in E and +0.64 kJ/mol in G from dimethylbutyl_cation (SI entry 9).
Conformer kept on record, not the zero: rotamer_h_plane (`CC(C)[C+](C)C`) lies +0.01 kJ/mol in E and +0.40 kJ/mol in G from dimethylbutyl_cation (SI entry 9).
The reactant end of the path ended 0.260 Å (RMSD) from the reactant-end minimum, the product end 0.712 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = separated reactants (dimethylbutyl_cation), 298.15 K, ΔG at 1 M:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 15.3 (3.7) | 10.9 (2.6) | 14.1 (3.4) |
| M06-2X/ma-def2-TZVP | 19.2 (4.6) | 14.8 (3.5) | 18.0 (4.3) |
| B3LYP/ma-def2-TZVP | 28.5 (6.8) | 24.1 (5.8) | 27.3 (6.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 18.2 (4.3) | 13.8 (3.3) | 17.0 (4.1) |

**Barriers, zero = reactant complex, 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 15.5 (3.7) | 10.7 (2.6) | 13.0 (3.1) |
| M06-2X/ma-def2-TZVP | 19.4 (4.6) | 14.6 (3.5) | 16.8 (4.0) |
| B3LYP/ma-def2-TZVP | 28.7 (6.9) | 24.0 (5.7) | 26.2 (6.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 18.4 (4.4) | 13.6 (3.2) | 15.8 (3.8) |

**Reaction, zero = separated reactants; product = separated dimethylbutyl_cation_shifted, 298.15 K, ΔG at 1 M:**

| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from dimethylbutyl_cation_shifted, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | +0.0 (+0.0) | +0.0 (+0.0) | +0.0 (+0.0) | -0.2 (-0.0) / +1.2 (+0.3) | 15.3 (3.7) / 14.1 (3.4) |
| M06-2X/ma-def2-TZVP | +0.0 (+0.0) | +0.0 (+0.0) | +0.0 (+0.0) | -0.2 (-0.0) / +1.2 (+0.3) | 19.2 (4.6) / 18.0 (4.3) |
| B3LYP/ma-def2-TZVP | +0.0 (+0.0) | +0.0 (+0.0) | +0.0 (+0.0) | -0.2 (-0.1) / +1.2 (+0.3) | 28.5 (6.8) / 27.3 (6.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +0.0 (+0.0) | +0.0 (+0.0) | +0.0 (+0.0) | -0.2 (-0.0) / +1.2 (+0.3) | 18.2 (4.3) / 17.0 (4.1) |

ΔrH was not asked for; it is given because it costs nothing.

## 4. Highest point of the path at each level

Path energies at each level, zero = separated reactants (dimethylbutyl_cation) at that level. The TS is frame 264 (s = 0, C2–H 1.373 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C2–H (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 543 of 543 | 15.3 (3.7) | 264 | +0.00 | 1.373 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 543 of 543 | 19.2 (4.6) | 263 | -0.05 | 1.385 | +0.00 (+0.00) |
| B3LYP/ma-def2-TZVP | 543 of 543 | 28.5 (6.8) | 271 | +0.28 | 1.325 | +0.52 (+0.12) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 489 of 543 | 18.2 (4.3) | 268 | +0.14 | 1.347 | +0.01 (+0.00) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter C2–H distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0105 | 0.050418 |
| reactant_complex | 0.0127 | 0.033675 |
| product_complex | 0.0127 | 0.033675 |
| dimethylbutyl_cation | 0.0127 | 0.033701 |
| dimethylbutyl_cation_shifted | 0.0127 | 0.033701 |
| rotamer_h_perp | 0.0127 | 0.033743 |
| rotamer_h_plane | 0.0127 | 0.033965 |

- **Path frames (245 computed):** T1 from 0.0104 to 0.0127, median 0.0127; largest at frame 22 (s = -12.19, C2–H 2.048 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: fails (flagged)** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 264 | 0 | 9.61° / 39.28° (21) | 5.74° / 10.18° (242) | no |
| product side (forward) | 129 | 0 | 25.16° / 44.82° (20) | 4.27° / 8.38° (108) | no |

- **End checks: ok** (`02_irc/true_irc/end_check.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame +1.468 kJ/mol from `open_cation` (RMSD 0.260 Å). Opt Freq from the last frame: E -0.059 kJ/mol from it, heavy-atom RMSD 0.009 Å (all-atom 0.029 Å), largest key-distance difference +0.000 (C2–H) Å, 0 imaginary frequencies: reaches the same minimum.
  - Product side (forward): last frame -0.003 kJ/mol from `bridged_c2` (RMSD 0.004 Å). Opt Freq from the last frame: E -0.002 kJ/mol from it, heavy-atom RMSD 0.001 Å (all-atom 0.004 Å), largest key-distance difference -0.005 (C3–H) Å, 0 imaginary frequencies: last kept frame matches the bridged reference; optimization lands nearest bridged_c2 (RMSD 0.004 Å, E − bridged reference -0.002 kJ/mol).
  - Official end structures now: reactant end `open_cation`, product end `open_cation` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: open side: 1.47 kJ/mol (end check ok); product side: bridged side: at the bridged structure (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 2 (join that picked the wrong bridged side): the first join (2026-11-15 12:13) built the mirror path from the wrong half; the join was redone at 12:18 with the bridged side chosen as the side that comes closest to `bridged_c2` at any frame. No IRC calculation was affected.
  - Entry 5 (stopping rule near the end minimum): the open side is the earlier side's 113 frames plus 151 frames continued from its last frame (continuation test: identical to a rerun from the TS within 0.000013 Å and 0.00007 kJ/mol); it ends 1.47 kJ/mol (0.35 kcal/mol) above the open cation. The angle check beyond |s| = 1 on that side has median 5.74° (rule < 5°) and 90th percentile 10.18°; kept as computed, both parts stated (J. Picione, 2026-11-16).
  - Entry 6 (angle against gradient size): 21 of the 70 frames over 15° in all six paths lie near this very flat TS and 3 in the slow final approach to the open cation, all where the gradient is nearly zero.
  - Entry 7 (mirror frames): the second half of the path is the mirror image of the first; at M06-2X/ma-def2-TZVP every frame was computed and all 271 mirror frames agree with their partners (largest |ΔE| 1.40 × 10⁻³ kJ/mol), so B3LYP and DLPNO-CCSD(T) mirror frames take their partner's energy.
  - Entry 9 (lowest-free-energy conformers; tertiary rotamer check, J. Picione 2026-11-19): the two other C2–C3 rotamer starts (C3–H along the empty p orbital; C3–H in the cation plane) both relaxed to the rotamer already computed (RMSD ≤ 0.018 Å), so the cation has one rotamer at this level; `dimethylbutyl_cation` is lowest in G at 298.15 K and at 135 K and stays the zero. Frame 0 lies above it in E by 1.27 / 1.86 / 1.67 / 2.22 kJ/mol (working level / M06-2X/TZVP / B3LYP/TZVP / DLPNO-CCSD(T)/TZVP). ΔG‡ at 135 K: 13.2 / 17.1 / 26.4 / 16.1 kJ/mol against 13.0 ± 0.4 measured. The working-level profile now takes its zero by G like the other levels: −0.195 kJ/mol (TS frame 15.5 → 15.3 kJ/mol).

## 7. Experiment

Experiment: Saunders and Kates 1978 (M. Saunders and M. R. Kates, J. Am. Chem. Soc. 1978, 100, 7082–7083, doi:10.1021/ja00490a060), compound I, the degenerate 1,2-hydride shift in the 2,3-dimethyl-2-butyl cation: ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol) at −138 °C (135 K); k = 3.1 × 10⁷ s⁻¹ at 135 K, 4.6 × 10⁷ s⁻¹ at −128 °C (145 K) and 1.4 × 10⁸ s⁻¹ at −119 °C (154 K). Solvent 50:50 SO₂ClF/SO₂F₂ (superacid medium); ¹³C NMR line broadening at 67.9 MHz in the fast-exchange limit; temperatures ±2 °C. The authors did not fit an Arrhenius equation, so there is no measured activation energy. Checked by J. Picione against the original paper (2026-11-11).

The comparison is made at 135 K, not 298 K, and in a different medium: the experiment is in a superacid medium, the calculation in SMD dichloromethane (and in the gas phase). Computed ΔG‡ is measured from the lowest cation conformer. At 135 K, RT is only 1.12 kJ/mol, so a small error in the barrier is a large error in rate. The working-level ΔG‡ at 135 K in dichloromethane, 13.2 kJ/mol (3.2 kcal/mol), is 0.2 kJ/mol above the reported 13.0 kJ/mol, a factor of about 1.2 in rate; the DLPNO-CCSD(T)/ma-def2-TZVP value, 16.1 kJ/mol (3.8 kcal/mol), is 3.1 kJ/mol (0.7 kcal/mol) above it, a factor of about 16. The rate ratio therefore looks much worse than the agreement in energy (numbers corrected 2026-11-22 from the finished path folder; the earlier text had 18.4 kJ/mol and a factor of 120 from an earlier folder).

The authors infer from their data that hyperconjugation is weaker at the transition state than in the cation.

At M06-2X/ma-def2-SVP (DefGrid3; the working level for the two cation reactions) the top of the path is flat within 0.06 kJ/mol (0.014 kcal/mol): the transition state between the open cation and the hydrogen-bridged structure lies 0.010 kJ/mol above the bridged structure in dichloromethane and 0.056 kJ/mol in the gas phase. The bridged well is shallower than the zero-point energy of the transfer mode (0.38 kJ/mol in dichloromethane, 0.54 kJ/mol in the gas phase). The production transition state (115.2i cm⁻¹) gives ΔG‡ = 13.2 kJ/mol (3.2 kcal/mol) at 135 K in dichloromethane at this level, against the measured 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol). (The level test of 2026-11-12, with an earlier optimization of the same flat TS, 108.9i cm⁻¹, gave 12.4 kJ/mol to the transition state and 12.7 kJ/mol to the bridged structure.)

Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):

- Solvent: the experiments on these cations are in superacid media, which the SMD solvent model has no parameters for. Dichloromethane (SMD) is used as a polar, non-nucleophilic stand-in; the gas phase is computed alongside for comparison.
- Tunnelling by the migrating hydrogen is not included. It may matter at the low temperatures of the NMR experiments (135–154 K). The Wigner factor at 135 K is about 1.06 (imaginary frequency 115.2i cm⁻¹ at the working level in DCM, hν‡/k_BT = 1.23). It is small because the top of this barrier is very flat; a flat-topped, wide barrier gives little tunnelling, so the true correction is probably also small, but the Wigner formula is only a rough indicator (corrected 2026-11-22: the earlier text used the B3LYP 563i cm⁻¹, factor about 2.5).
- Only half of the path is computed: the IRC runs from the open cation (H on C3) through the TS to the C₂-symmetric bridged structure. The second half, bridged → TS′ → open cation (H on C2), is its mirror image under the C2 ⇄ C3 swap; those frames are labelled as mirror images in the viewer and are not IRC points. The bridged side is limited by SCF noise: the bridged structure lies only about 0.01 kJ/mol below the TS in DCM, below the SCF noise floor of 0.034 kJ/mol (J. Picione, 2026-11-14).
- At this level (M06-2X/ma-def2-SVP, SMD dichloromethane) the bridged structure and the transition state are not distinguishable within numerical precision: the bridged structure lies 0.012 kJ/mol below the TS, below the SCF noise floor of 0.034 kJ/mol and the optimizer's energy tolerance. The top of the barrier is therefore a flat plateau, not a sharp maximum (J. Picione, 2026-11-14).

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 543 | 8.8 | 59 | 62 | 71 (frame 466) | 2026-11-16 22:41 | 2026-11-17 20:52 |
| B3LYP/ma-def2-TZVP | 272 | 2.3 | 31 | 31 | 34 (frame 131) | 2026-11-17 21:00 | 2026-11-19 14:26 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 245 | 46.5 | 670 | 633 | 870 (frame 6) | 2026-11-19 14:38 | 2026-11-22 16:34 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 598 gradient jobs, 6.5 h job time (folders `backward`); product side (forward): 445 gradient jobs, 4.9 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.
