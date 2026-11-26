# Card details: Diels–Alder, butadiene + ethylene

## Energies at every level (same M06-2X-D3(0)/ma-def2-SVP geometries, SMD toluene)

| Level | ΔE‡ from separated reactants | ΔE‡ from the reactant complex | Reaction energy ΔrE | ΔG‡ from separated reactants (298 K, 1 M) | ΔG‡ from the reactant complex (298 K) | Reaction free energy ΔrG (298 K, 1 M) |
|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (geometries) | 78.4 kJ/mol (18.7 kcal/mol) | 75.6 kJ/mol (18.1 kcal/mol) | −218.1 kJ/mol (−52.1 kcal/mol) | 126.0 kJ/mol (30.1 kcal/mol) | 96.2 kJ/mol (23.0 kcal/mol) | −151.7 kJ/mol (−36.2 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP | 84.2 kJ/mol (20.1 kcal/mol) | 80.0 kJ/mol (19.1 kcal/mol) | −198.2 kJ/mol (−47.4 kcal/mol) | 131.9 kJ/mol (31.5 kcal/mol) | 100.5 kJ/mol (24.0 kcal/mol) | −131.8 kJ/mol (−31.5 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP (default curve) | 96.7 kJ/mol (23.1 kcal/mol) | 90.2 kJ/mol (21.6 kcal/mol) | −189.7 kJ/mol (−45.3 kcal/mol) | 144.4 kJ/mol (34.5 kcal/mol) | 110.7 kJ/mol (26.5 kcal/mol) | −123.3 kJ/mol (−29.5 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP | 79.8 kJ/mol (19.1 kcal/mol) | 72.5 kJ/mol (17.3 kcal/mol) | −171.3 kJ/mol (−41.0 kcal/mol) | 127.4 kJ/mol (30.5 kcal/mol) | 93.0 kJ/mol (22.2 kcal/mol) | −104.9 kJ/mol (−25.1 kcal/mol) |

ΔG adds the M06-2X-D3(0)/ma-def2-SVP thermal corrections (298.15 K, quasi-RRHO) to each level's electronic energies, and the 1 atm → 1 M standard-state correction for steps that change the number of molecules; the reactant complex → TS step needs none.

## Comparison with the gas-phase experiment

Gas phase, no solvent model, gas-phase geometries and frequencies at M06-2X-D3(0)/ma-def2-SVP; barrier from s-trans butadiene + ethylene; T = 840 K (midpoint of the measured 760–921 K).

| Level | ΔE‡ | ΔH‡ (840 K) | Ea = ΔH‡ + 2RT (840 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | k(TST, 840 K) (L mol⁻¹ s⁻¹) | k / measured k |
|---|---|---|---|---|---|---|---|---|
| M06-2X_ma-def2-SVP (working) | 80.0 kJ/mol (19.1 kcal/mol) | 86.9 kJ/mol (20.8 kcal/mol) | 100.9 kJ/mol (24.1 kcal/mol) | −14.2 kJ/mol (−3.4 kcal/mol) | 1.36 × 10⁷ | 0.45 | 7.26 | 3.48 |
| M06-2X_ma-def2-TZVP | 85.6 kJ/mol (20.5 kcal/mol) | 92.5 kJ/mol (22.1 kcal/mol) | 106.5 kJ/mol (25.4 kcal/mol) | −8.6 kJ/mol (−2.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 3.27 | 1.57 |
| B3LYP_ma-def2-TZVP | 81.3 kJ/mol (19.4 kcal/mol) | 88.2 kJ/mol (21.1 kcal/mol) | 102.2 kJ/mol (24.4 kcal/mol) | −12.9 kJ/mol (−3.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 6.01 | 2.88 |
| DLPNO-CCSD(T)_ma-def2-TZVP | 98.4 kJ/mol (23.5 kcal/mol) | 105.3 kJ/mol (25.2 kcal/mol) | 119.3 kJ/mol (28.5 kcal/mol) | +4.2 kJ/mol (+1.0 kcal/mol) | 1.36 × 10⁷ | 0.45 | 0.524 | 0.25 |
| Measured (Rowley and Steiner 1951) | | | 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol) | | 3.00 × 10⁷ | 1 | 2.09 | 1 |

**Units and standard state of A.** A and k(TST) are second-order rate constants in concentration units: L mol⁻¹ s⁻¹ (= 10³ cm³ mol⁻¹ s⁻¹, Rowley and Steiner's unit), also given per molecule (cm³ molecule⁻¹ s⁻¹). ΔG‡ and ΔS‡ use the 1 atm standard state; the factor RT/p° turns that into concentration units, so A and k do not depend on the choice of standard state and compare directly with the measured A. Measured A = 3.00 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 3.00 × 10⁷ L mol⁻¹ s⁻¹ = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹; computed A = 1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 2.27 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹.

A is the same at every level: all levels share the working-level gas-phase frequencies, so ΔS‡ (840 K, 1 atm) = −168.8 J mol⁻¹ K⁻¹ at every level; the levels differ only in Ea.

Method: ΔH‡(840 K) = ΔE‡ + Δ[H(840 K) − E]; Ea = ΔH‡ + 2RT; ΔG‡(840 K, 1 atm) = ΔE‡ + Δ[G(840 K) − E] (quasi-RRHO entropy); k(TST) = (kB T/h)(RT/p°) exp(−ΔG‡/RT), κ = 1; A = k(TST) exp(Ea/RT). Thermal terms from the gas-phase working-level frequencies (rotational symmetry numbers from ORCA's point-group detection: ts σ = 1, butadiene_s-trans σ = 2, ethylene σ = 4).

**Gas-phase TS:** one imaginary frequency, −514 cm⁻¹; overlap with the toluene TS mode 1.00; forming bonds C1–C6 2.265 Å and C4–C5 2.265 Å (spread 0.000 Å, stretch-rate ratio 1.00): synchronous, accepted by the folder's TS rules. Toluene TS: 2.263 Å.

**Butadiene conformers (gas phase, working level):** the gauche s-cis form lies 11.6 kJ/mol (2.8 kcal/mol) in E and 9.2 kJ/mol (2.2 kcal/mol) in G (298.15 K) above s-trans; the lowest conformer is s-trans.

**Findings, s-trans butadiene only.**

- **Activation energy.** DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 119.3 kJ/mol (28.5 kcal/mol), 4.2 kJ/mol (1.0 kcal/mol) above the measured 115.1 ± 3.4 kJ/mol, just outside the fit error. M06-2X/ma-def2-TZVP is 8.6 kJ/mol (2.1 kcal/mol) below, B3LYP/ma-def2-TZVP 12.9 kJ/mol (3.1 kcal/mol) below, and the working level 14.2 kJ/mol (3.4 kcal/mol) below.
- **Pre-exponential factor.** A = 1.36 × 10⁷ L mol⁻¹ s⁻¹ (1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹), 0.45 times the measured 3.0 × 10⁷ L mol⁻¹ s⁻¹, at every level (ΔS‡ = −168.8 J mol⁻¹ K⁻¹ at 840 K, 1 atm). A and Ea from one Arrhenius fit over 760–921 K are correlated: the ±3.4 kJ/mol in Ea corresponds to about a factor of 1.6 in A at 840 K, so the computed A is low by somewhat more than the fit uncertainty.
- **Rate at 840 K.** DLPNO-CCSD(T) gives 0.25 times the measured k: 0.45 from A times 0.55 from the higher Ea. M06-2X/ma-def2-TZVP gives 1.57 times: 0.45 from A times 3.4 from its low Ea. M06-2X matches the rate better only because its low Ea cancels the low A.
- **Gas phase against toluene** (same levels, s-trans butadiene + ethylene): the DLPNO-CCSD(T) ΔE‡ is 98.4 kJ/mol in the gas phase and 96.7 kJ/mol in SMD toluene; the forming C–C bonds are 2.265 Å and 2.263 Å.
- **Checks.** Every job ended normally. The minima have no imaginary frequency. The TS has one (−514 cm⁻¹), synchronous, with overlap 1.00 with the toluene TS mode. The thermal terms reproduce ORCA's 298.15 K H and G within 0.001 kJ/mol. The largest DLPNO-CCSD(T) T1 is 0.0113.

**Butadiene conformer ensemble.** At 760–921 K the reactant is the equilibrium mixture of s-trans butadiene and the two mirror-image gauche (s-cis) forms, as in the experiment. k(TST) was computed at 760, 800, 840, 880 and 921 K and ln k fitted against 1/T, as Rowley and Steiner did; Ea = −R × slope, A = exp(intercept), in L mol⁻¹ s⁻¹. The s-trans-only values below use the same fit, so they differ slightly from the 840 K point values (Ea = ΔH‡ + 2RT) in the table above.

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

**Findings, conformer ensemble.**

- **Activation energy.** With the equilibrium conformer mixture as the reactant, DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 114.9 kJ/mol (27.5 kcal/mol) from the 760–921 K fit, 0.2 kJ/mol (0.05 kcal/mol) below the measured 115.1 ± 3.4 kJ/mol. From s-trans butadiene alone the same fit gives 118.6 kJ/mol (28.3 kcal/mol), 3.5 kJ/mol above. Every level drops by 3.4–3.7 kJ/mol (0.8–0.9 kcal/mol): at these temperatures about 30 % of the butadiene is gauche, which lies 12.2 kJ/mol (2.9 kcal/mol) above s-trans in H at 840 K (DLPNO-CCSD(T)), so the average reactant enthalpy is higher. M06-2X/ma-def2-TZVP is now 13.0 kJ/mol (3.1 kcal/mol) below the measured value, B3LYP/ma-def2-TZVP 17.0 kJ/mol (4.1 kcal/mol) below.
- **Pre-exponential factor.** A falls from 1.24 × 10⁷ (s-trans only, same fit) to 5.1 × 10⁶ L mol⁻¹ s⁻¹ at DLPNO-CCSD(T), 0.17 times the measured 3.0 × 10⁷ L mol⁻¹ s⁻¹ (0.17–0.19 at every level). Only the s-trans part reacts in this model, and its share falls with temperature, which lowers both the fitted Ea and A. The ±3.4 kJ/mol fit error in Ea corresponds to about a factor of 1.6 in A, so the computed A is low by clearly more than the measurement's uncertainty.
- **Rate at 840 K.** The mixture lowers k(TST) to 0.70 of the s-trans-only value (the s-trans fraction). DLPNO-CCSD(T) gives 0.17 times the measured k, a free-energy gap of 12.4 kJ/mol (3.0 kcal/mol) at 840 K that now comes almost entirely from the entropy of activation (computed ΔS‡ more negative than measured by about R ln 6 ≈ 15 J mol⁻¹ K⁻¹). Likely contributors, not tested here: the harmonic (quasi-RRHO) treatment of the low-frequency modes of the TS and of the butadiene torsion at 840 K.
- **What the conformer mixture changes.** From s-trans butadiene alone (840 K point values) DLPNO-CCSD(T) puts Ea 4.2 kJ/mol above experiment, just outside the error, with A 0.45 times the measured value. With the conformer mixture the coupled-cluster Ea agrees with experiment, and the remaining disagreement is in A (0.17 times). M06-2X/ma-def2-TZVP comes within 1.08 times the measured k only because its Ea, 13.0 kJ/mol too low, cancels the low A.
- **Gauche fraction** (DLPNO-CCSD(T)): 1.8 % at 298 K, 26.6 % at 760 K, 30.4 % at 840 K, 33.7 % at 921 K. B3LYP gives fewer gauche (24.9 % at 840 K), because it puts gauche 2.0 kJ/mol higher.
- **Checks.** The three gauche single points ended normally (largest DLPNO-CCSD(T) T1 0.0111). The ln k against 1/T fits are straight within 0.007 in ln k at every level.

## The measured value: source and derivation

Experiment (gas phase): C₂H₄ + 1,3-butadiene → cyclohexene, second order, 760–921 K: k = 3.0 × 10¹⁰ exp(−27,500 cal mol⁻¹/RT) cm³ mol⁻¹ s⁻¹ (eq. 5), i.e. A = 3.0 × 10⁷ L mol⁻¹ s⁻¹ (4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹); a refit of the 13 runs in Table VI gives Ea = 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol). Source: D. Rowley and H. Steiner, Discuss. Faraday Soc. 1951, 10, 198–213, doi:10.1039/DF9511000198. Cross-check: NIST Chemical Kinetics Database, squib 1951ROW/STE198-213 (same A and Ea). The measured activation energy refers to equilibrium butadiene in the gas phase at 760–921 K. That is almost all s-trans at room temperature (about 2 % gauche s-cis) but not at the temperatures of the measurement: about 27 % gauche at 760 K and 34 % at 921 K (DLPNO-CCSD(T)/ma-def2-TZVP, gas phase). This calculation is in SMD toluene and starts from s-cis butadiene, so the s-trans → s-cis energy is reported separately; the comparison with this measurement uses separate gas-phase calculations, from s-trans butadiene alone and from the equilibrium conformer mixture, each with ethylene. Concerted, synchronous transition state: K. N. Houk, Y.-T. Lin and F. K. Brown, J. Am. Chem. Soc. 1986, 108, 554–556, doi:10.1021/ja00263a059; E. Goldstein, B. Beno and K. N. Houk, J. Am. Chem. Soc. 1996, 118, 6036–6043, doi:10.1021/ja9601494.

## Path and orbitals

- IRC frames: 302; reaction coordinate s from −8.94 to 5.97 bohr·amu½ (TS at s = 0).
- Path quality: bond-changing region IRC frames 20–264, median step there 0.050 bohr·amu½; steps longer than 2 × that median: none.
- Orbitals: IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^−3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital, all at the step from IRC frame 188 to 189 (s = +0.45 to +0.50 bohr·amu^1/2; TS at frame 179): π(C3=C4) → σ(C4–C5) 0.714, π(C5=C6) → σ(C1–C6) 0.725, π(C1=C2) → π(C2=C3) 0.716; the neighbouring steps are 1.00–1.00; for comparison, lowest overlap of the shown orbitals in E2 1.000, Tertiary hydride shift 1.000, Secondary hydride shift 1.000, Cyanide + acetone 1.000, SN2 1.000); canonical HOMO/LUMO optional in the viewer.
- Direction of the arrows: In this pericyclic reaction the clockwise and counter-clockwise arrow drawings are mirror images of the same mechanism, and the mirror-symmetric transition state favours neither. The orbitals shown are the computed set, which localized clockwise on this path. The two are equally valid: the direction is the localization's choice, not a physical preference. At 6 frames between the transition state and s = +1.0 bohr·amu½ (frames 181, 182, 185, 187, 188, 198) the localization landed on the counter-clockwise set; there the mirror image of the computed orbitals is shown, which continues the clockwise set. This is valid because the geometry at those frames is mirror-symmetric within 0.001 Å. The viewer's "mirror direction" switch shows the mirror image of the displayed set at every frame (the counter-clockwise drawing). The path is mirror-symmetric within 0.0008 Å except at frames 0–2, 4, 6 and 299 (up to 0.0040 Å), and from the transition state to s = +1.0 bohr·amu½ (frames 179–199) the reflected orbitals lie in the computed occupied space (smallest singular value of the 3×17 overlap matrix 0.99982 or higher).
- Display rotation: Near the reactant end the IBO localization gave the isolated C5=C6 double bond as two bent bonds, not equivalent (s-character 0.43 and −0.00 at IRC frame 0); for display only, the pair was rotated within itself into the σ + π pair that spans the same two orbitals, and the π orbital is shown, in IRC frames 0–139 (where the computed orbital is more than 0.5° from that split; 87° at frame 0), joining the unrotated frames with overlap 1.0000; a rotation inside an occupied pair changes no energy or density, and the unrotated orbitals are kept in the data files (`03_orbitals/display_rotation.json`).

## Technical notes

- M06-2X underestimates this barrier: ΔG‡ = 131.9 kJ/mol (31.5 kcal/mol) at M06-2X-D3(0)/ma-def2-TZVP versus 144.4 kJ/mol (34.5 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (12 kJ/mol lower). The comparison with experiment is the gas-phase Ea table above.
- Reactant complex (a weak, dispersion-bound pair of the neutral reactants) lies 2.7 kJ/mol (0.7 kcal/mol) above separated reactants on the working-level curve and 4.2 kJ/mol (1.0 kcal/mol) above at ma-def2-TZVP, and 6.5 kJ/mol (1.6 kcal/mol) above at DLPNO-CCSD(T)/ma-def2-TZVP; the barriers above are also given from separated reactants.
- Implicit solvent only: no explicit toluene molecules.
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = 5.97 bohr·amu½ (IRC step 122), after the bonds had finished changing; the product-end minimum was optimized from that last point.
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = −8.94 bohr·amu½ (IRC step 179), after the bonds had finished changing; beyond it only the separating fragments drift apart. The reactant complex was optimized from that last point.
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- The canonical HOMO changes character abruptly along the path (near frame 178); canonical orbitals are delocalized and re-ordered by energy, which is why the localized IBOs are the default view.
- The path ends at the boat form of cyclohexene (mirror-symmetric, the new CH₂–CH₂ bond eclipsed): a synchronous path keeps the mirror plane of the TS. Optimized from the IRC end it stays in the boat form (lowest vibration 47 cm⁻¹, a very flat ring motion), 23.0 kJ/mol (5.5 kcal/mol) above the half-chair in E (M06-2X/ma-def2-SVP), the ground-state conformer. Reaction energies use half-chair cyclohexene; the boat is the end of the path shown in the viewer and video.
- Gas-phase comparison with Rowley and Steiner (760–921 K), on gas-phase geometries and frequencies (gas-phase section of the card and of the results summary). With the equilibrium butadiene conformer mixture as the reactant (s-trans plus both gauche forms, about 30 % gauche at 840 K) and an ln k against 1/T fit over 760–921 K as in the experiment, DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 114.9 kJ/mol (27.5 kcal/mol), 0.2 kJ/mol (0.05 kcal/mol) below the measured 115.1 ± 3.4 kJ/mol, and M06-2X/ma-def2-TZVP 13.0 kJ/mol (3.1 kcal/mol) below. The computed pre-exponential factor, 5.1 × 10⁶ L mol⁻¹ s⁻¹, is 0.17 times the measured one (0.17–0.19 at every level), so DLPNO-CCSD(T) gives 0.17 times the measured rate at 840 K; M06-2X/ma-def2-TZVP gives 1.08 times only because its low Ea cancels the low A. From s-trans butadiene alone (840 K point values) the numbers were +4.2 kJ/mol (DLPNO-CCSD(T)), −8.6 kJ/mol (M06-2X/ma-def2-TZVP) and A 0.45 times the measured one.
- The product side of the path ended at the minimum (0.011 kJ/mol above the optimized cyclohexene, RMS gradient 1.0 × 10⁻⁶ Eh/bohr), stopped when the energy began to rise, before the 20 flat steps of the stopping rule.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu½), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/da_butadiene_ethylene.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`da_butadiene_ethylene_<level>.mp4`: the same video with the other energy curves)
- `web/viewer/index.html?r=da_butadiene_ethylene`: interactive viewer
- `run_log.md`: every ORCA job and rerun

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
