# Reaction card: Diels–Alder, butadiene + ethylene

| | |
|---|---|
| Reaction | Diels–Alder, butadiene + ethylene |
| Model system | CH₂=CH–CH=CH₂ + CH₂=CH₂ → cyclohexene |
| Level of theory | B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX; barrier at B3LYP-D3(BJ)/ma-def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, toluene |
| Barrier (ma-def2-TZVP, electronic) | **72.5 kJ/mol (17.3 kcal/mol)** relative to the reactant complex; 65.2 kJ/mol (15.6 kcal/mol) relative to separated reactants |
| Barrier (ma-def2-SVP, electronic, the path shown) | 66.3 kJ/mol (15.9 kcal/mol) relative to the reactant complex; 58.0 kJ/mol (13.9 kcal/mol) relative to separated reactants |
| Reaction energy (ma-def2-TZVP, electronic) | -186.0 kJ/mol (-44.5 kcal/mol) |
| Imaginary frequency at TS | 516i cm⁻¹ (exactly one) |
| TS geometry | C1–C6 2.27 Å (forming), C4–C5 2.26 Å (forming), C2–C3 1.41 Å (becoming double), C1–C2 1.39 Å (becoming single), C3–C4 1.39 Å (becoming single), C5–C6 1.39 Å (becoming single), C1=C2–C3=C4 dihedral 0.0° |
| IRC frames | 110; reaction coordinate s from -5.72 to 5.78 bohr·amu^1/2 (TS at s = 0) |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital, all at the step from IRC frame 59 to 60 (s = +0.63 to +0.86 bohr·amu^1/2; TS at frame 58): π(C1=C2) → σ(C1–C6) 0.835, π(C5=C6) → σ(C4–C5) 0.838, π(C3=C4) → π(C2=C3) 0.827; the neighbouring steps are 0.94–0.96; for comparison, lowest overlap of the shown orbitals in E2 0.990, Nucleophilic addition 0.992, SN2 0.993; tested with intermediate IBO single points at ¼ and ½ of step 59→60 (Cartesian interpolation, `03_orbitals/step_test_59_60/`): the localized orbitals switch abruptly within its first quarter, overlap 0.79 from frame 59 to the ¼ point, then 0.996–0.997 on to frame 60, so more IRC frames would not spread the change; E2, Nucleophilic addition and SN2 show no such switch: every frame-to-frame overlap of their shown orbitals is at least 0.99); canonical HOMO/LUMO optional in the viewer |
| Display rotation | Near the reactant end the IBO localization gave the isolated C5=C6 double bond as two bent bonds, not equivalent (s-character 0.18 and 0.24 at IRC frame 0); for display only, the pair was rotated within itself into the σ + π pair that spans the same two orbitals, and the π orbital is shown, in IRC frames 0–47 (where the computed orbital is more than 0.5° from that split; 41° at frame 0), joining the unrotated frames with overlap 0.9998; a rotation inside an occupied pair changes no energy or density, and the unrotated orbitals are kept in the data files (`03_orbitals/display_rotation.json`) |
| Orbitals that change | π(C1=C2) → σ(C1–C6); π(C5=C6) → σ(C4–C5); π(C3=C4) → π(C2=C3). In this pericyclic reaction the clockwise and counter-clockwise arrow drawings are mirror images of the same mechanism, and the mirror-symmetric transition state favours neither. The direction shown is the one the orbital localization settled on, not a physical preference (the viewer's "mirror direction" switch shows the other one). |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in by J. Picione) |

## Barrier at other levels (same geometries)

| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |
|---|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP (working) | 58.0 kJ/mol (13.9 kcal/mol) | 66.3 kJ/mol (15.9 kcal/mol) | -209.6 kJ/mol (-50.1 kcal/mol) | 108.1 kJ/mol (25.8 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | 65.2 kJ/mol (15.6 kcal/mol) | 72.5 kJ/mol (17.3 kcal/mol) | -186.0 kJ/mol (-44.5 kcal/mol) | 115.3 kJ/mol (27.6 kcal/mol) |
| M06-2X/ma-def2-TZVP | 71.4 kJ/mol (17.1 kcal/mol) | 78.7 kJ/mol (18.8 kcal/mol) | -211.5 kJ/mol (-50.5 kcal/mol) | 121.5 kJ/mol (29.0 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 84.2 kJ/mol (20.1 kcal/mol) | 89.2 kJ/mol (21.3 kcal/mol) | -202.2 kJ/mol (-48.3 kcal/mol) | 134.3 kJ/mol (32.1 kcal/mol) |

ΔG adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic energies.

## Comparison with the gas-phase experiment

Stationary points recomputed without the solvent model at every level (same geometries). Barriers from butadiene s-trans + ethylene (the measured reference) unless marked; Ea = ΔH‡(T) + 2RT for a gas-phase bimolecular reaction, with ΔH‡ recomputed at T = 840 K, the midpoint of the measured range 760–921 K.

| Level | Phase | ΔE‡ from butadiene s-cis + ethylene | ΔE‡ | ΔH‡ (298.15 K) | ΔG‡ (298.15 K) | Ea (840 K) |
|---|---|---|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP (working) | SMD toluene | 58.0 kJ/mol (13.9 kcal/mol) | 71.8 kJ/mol (17.2 kcal/mol) | 76.3 kJ/mol (18.2 kcal/mol) | 119.1 kJ/mol (28.5 kcal/mol) (1 M) | 92.9 kJ/mol (22.2 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-SVP (working) | gas phase | 59.8 kJ/mol (14.3 kcal/mol) | 73.7 kJ/mol (17.6 kcal/mol) | 78.2 kJ/mol (18.7 kcal/mol) | 128.9 kJ/mol (30.8 kcal/mol) (1 atm) | **94.8 kJ/mol (22.7 kcal/mol)** |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | SMD toluene | 65.2 kJ/mol (15.6 kcal/mol) | 79.5 kJ/mol (19.0 kcal/mol) | 84.0 kJ/mol (20.1 kcal/mol) | 126.8 kJ/mol (30.3 kcal/mol) (1 M) | 100.6 kJ/mol (24.0 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | gas phase | 66.6 kJ/mol (15.9 kcal/mol) | 80.9 kJ/mol (19.3 kcal/mol) | 85.4 kJ/mol (20.4 kcal/mol) | 136.1 kJ/mol (32.5 kcal/mol) (1 atm) | **102.0 kJ/mol (24.4 kcal/mol)** |
| M06-2X/ma-def2-TZVP | SMD toluene | 71.4 kJ/mol (17.1 kcal/mol) | 84.0 kJ/mol (20.1 kcal/mol) | 88.5 kJ/mol (21.1 kcal/mol) | 131.2 kJ/mol (31.4 kcal/mol) (1 M) | 105.1 kJ/mol (25.1 kcal/mol) |
| M06-2X/ma-def2-TZVP | gas phase | 72.5 kJ/mol (17.3 kcal/mol) | 85.2 kJ/mol (20.4 kcal/mol) | 89.7 kJ/mol (21.4 kcal/mol) | 140.4 kJ/mol (33.5 kcal/mol) (1 atm) | **106.3 kJ/mol (25.4 kcal/mol)** |
| DLPNO-CCSD(T)/ma-def2-TZVP | SMD toluene | 84.2 kJ/mol (20.1 kcal/mol) | 97.0 kJ/mol (23.2 kcal/mol) | 101.5 kJ/mol (24.3 kcal/mol) | 144.3 kJ/mol (34.5 kcal/mol) (1 M) | 118.1 kJ/mol (28.2 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | gas phase | 85.9 kJ/mol (20.5 kcal/mol) | 98.7 kJ/mol (23.6 kcal/mol) | 103.2 kJ/mol (24.7 kcal/mol) | 153.9 kJ/mol (36.8 kcal/mol) (1 atm) | **119.8 kJ/mol (28.6 kcal/mol)** |

Experiment (gas phase): Ea ≈ 115.0 kJ/mol (27.5 kcal/mol).
Comparison with the gas-phase measurement (gas-phase computed Ea, bold above): DLPNO-CCSD(T)/ma-def2-TZVP gives 119.8 kJ/mol (28.6 kcal/mol), +4.8 kJ/mol (+1.1 kcal/mol) from experiment. The SMD toluene value at the same level, 118.1 kJ/mol (28.2 kcal/mol), is given for reference only, not as the comparison.
butadiene s-trans → butadiene s-cis, ΔE: 12.8 kJ/mol (3.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (gas phase).
Thermal corrections: working-level harmonic frequencies (SMD), quasi-RRHO entropy for ΔG‡; the s-trans reference is why ΔE‡ here exceeds the barrier from the s-cis conformer.

### Rate constant at 840 K

Transition-state theory, gas phase: k = (k_B T/h)(RT/p°) exp(−ΔG‡/RT), with ΔG‡ at 840 K for a 1 atm (101 325 Pa) standard state, no tunnelling correction, and the reaction-path degeneracy carried by the rotational symmetry numbers. ΔG‡(T) uses each level's gas-phase ΔE‡ and the working-level thermal corrections (quasi-RRHO entropy, checked against ORCA's own values at 298.15 K).

| Level | ΔG‡ (840 K, 1 atm) | k(TST) (L mol⁻¹ s⁻¹) | k(TST) (cm³ molecule⁻¹ s⁻¹) | k(TST) / k(experiment) |
|---|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP (working) | 221.1 kJ/mol (52.8 kcal/mol) | 21.7 | 3.60 × 10⁻²⁰ | **10.2** |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | 228.3 kJ/mol (54.6 kcal/mol) | 7.69 | 1.28 × 10⁻²⁰ | **3.63** |
| M06-2X/ma-def2-TZVP | 232.5 kJ/mol (55.6 kcal/mol) | 4.19 | 6.96 × 10⁻²¹ | **1.98** |
| DLPNO-CCSD(T)/ma-def2-TZVP | 246.1 kJ/mol (58.8 kcal/mol) | 0.604 | 1.00 × 10⁻²¹ | **0.29** |

Experiment at 840 K, from the measured Arrhenius parameters (A = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹ = 3.00 × 10⁷ L mol⁻¹ s⁻¹, Ea = 115 kJ/mol): k = **2.12 L mol⁻¹ s⁻¹** (3.52 × 10⁻²¹ cm³ molecule⁻¹ s⁻¹).

Experiment (gas phase): C₂H₄ + 1,3-butadiene → cyclohexene, second order, 760–921 K, Ea = 115 kJ/mol (27.5 kcal/mol), A = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹ (3.0 × 10⁷ L mol⁻¹ s⁻¹). Pressure 3.72 × 10⁴ to 1.01 × 10⁵ Pa; rate constant reported as pressure independent; static thermal experiment followed by pressure measurement; absolute values. Source: Rowley and Steiner 1951 (D. Rowley and H. Steiner, Discuss. Faraday Soc. 1951, 10, 198–213, doi:10.1039/DF9511000198), as compiled in the NIST Chemical Kinetics Database, squib 1951ROW/STE198-213. Activation energy and temperature range confirmed by J. Picione from the NIST entry (2026-11-10); the original paper has not yet been read. The measured activation energy refers to equilibrium butadiene (almost all s-trans) in the gas phase; this calculation is in SMD toluene and starts from s-cis butadiene, so the s-trans → s-cis energy is reported separately. Concerted, synchronous transition state: K. N. Houk, Y.-T. Lin and F. K. Brown, J. Am. Chem. Soc. 1986, 108, 554–556, doi:10.1021/ja00263a059; E. Goldstein, B. Beno and K. N. Houk, J. Am. Chem. Soc. 1996, 118, 6036–6043, doi:10.1021/ja9601494.

## Known artifacts and limitations

- B3LYP underestimates this barrier: ΔG‡ = 115.3 kJ/mol (27.6 kcal/mol) at the protocol level versus 134.3 kJ/mol (32.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (19 kJ/mol lower). The comparison with experiment is the gas-phase Ea table above.
- Reactant complex (a weak, dispersion-bound pair of the neutral reactants) lies 8.4 kJ/mol (2.0 kcal/mol) below separated reactants on the working-level curve but only 7.3 kJ/mol (1.7 kcal/mol) below at ma-def2-TZVP, and 5.0 kJ/mol (1.2 kcal/mol) below at DLPNO-CCSD(T)/ma-def2-TZVP; the barriers above are also given from separated reactants.
- Implicit solvent only: no explicit toluene molecules.
- The optimized reactant complex used for the energies (−8.4 kJ/mol) is a different local minimum on a flat surface from the first frame of the path (−2.8 kJ/mol, the IRC end shown in the viewer and video), 5.6 kJ/mol (1.3 kcal/mol) lower (working level, relative to separated reactants).
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- The canonical HOMO changes character abruptly along the path (near frame 58); canonical orbitals are delocalized and re-ordered by energy, which is why the localized IBOs are the default view.
- The path ends at the boat form of cyclohexene (mirror-symmetric, the new CH₂–CH₂ bond eclipsed): a synchronous path keeps the mirror plane of the TS. Optimized from the IRC end it stays in the boat form (lowest vibration 38 cm⁻¹, a very flat ring motion), 23.1 kJ/mol (5.5 kcal/mol) above the half-chair, the ground-state conformer. Reaction energies use half-chair cyclohexene; the boat is the end of the path shown in the viewer and video.
- M06-2X matches the rate better than coupled cluster only because its low activation energy cancels an entropy error. At 840 K every level shares the same thermal corrections (working-level frequencies), so the computed pre-exponential factor is the same at every level, 0.57 times the measured one (from the entropy of activation); the levels differ only in Ea. M06-2X/ma-def2-TZVP puts Ea 8.7 kJ/mol (2.1 kcal/mol) below experiment, which raises its rate by a factor of 3.5 and brings k/k(exp) to 2.0; DLPNO-CCSD(T)/ma-def2-TZVP puts Ea 4.8 kJ/mol (1.1 kcal/mol) above, which lowers its rate by a factor of 2.0, to 0.29 (gas phase, symmetry numbers 4 for ethylene, 2 for s-trans butadiene, 1 for the TS, as detected by ORCA).

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/da_butadiene_ethylene.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`da_butadiene_ethylene_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=da_butadiene_ethylene`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
