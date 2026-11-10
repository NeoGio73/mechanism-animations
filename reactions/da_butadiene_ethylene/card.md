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
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of the shown orbitals 0.83); canonical HOMO/LUMO optional in the viewer |
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

ΔG‡ adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic barrier.

## Comparison with the gas-phase experiment

Stationary points recomputed without the solvent model at every level (same geometries). Barriers from butadiene s-trans + ethylene (the measured reference) unless marked; Ea = ΔH‡(T) + 2RT for a gas-phase bimolecular reaction, at T = 800 K, provisional.

| Level | Phase | ΔE‡ from butadiene s-cis + ethylene | ΔE‡ | ΔH‡ (298.15 K) | ΔG‡ (298.15 K) | Ea (800 K, provisional) |
|---|---|---|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP (working) | SMD toluene | 58.0 kJ/mol (13.9 kcal/mol) | 71.8 kJ/mol (17.2 kcal/mol) | 76.3 kJ/mol (18.2 kcal/mol) | 119.1 kJ/mol (28.5 kcal/mol) (1 M) | **91.9 kJ/mol (22.0 kcal/mol)** |
| B3LYP-D3(BJ)/ma-def2-SVP (working) | gas phase | 59.8 kJ/mol (14.3 kcal/mol) | 73.7 kJ/mol (17.6 kcal/mol) | 78.2 kJ/mol (18.7 kcal/mol) | 128.9 kJ/mol (30.8 kcal/mol) (1 atm) | **93.8 kJ/mol (22.4 kcal/mol)** |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | SMD toluene | 65.2 kJ/mol (15.6 kcal/mol) | 79.5 kJ/mol (19.0 kcal/mol) | 84.0 kJ/mol (20.1 kcal/mol) | 126.8 kJ/mol (30.3 kcal/mol) (1 M) | **99.6 kJ/mol (23.8 kcal/mol)** |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | gas phase | 66.6 kJ/mol (15.9 kcal/mol) | 80.9 kJ/mol (19.3 kcal/mol) | 85.4 kJ/mol (20.4 kcal/mol) | 136.1 kJ/mol (32.5 kcal/mol) (1 atm) | **101.0 kJ/mol (24.1 kcal/mol)** |
| M06-2X/ma-def2-TZVP | SMD toluene | 71.4 kJ/mol (17.1 kcal/mol) | 84.0 kJ/mol (20.1 kcal/mol) | 88.5 kJ/mol (21.1 kcal/mol) | 131.2 kJ/mol (31.4 kcal/mol) (1 M) | **104.1 kJ/mol (24.9 kcal/mol)** |
| M06-2X/ma-def2-TZVP | gas phase | 72.5 kJ/mol (17.3 kcal/mol) | 85.2 kJ/mol (20.4 kcal/mol) | 89.7 kJ/mol (21.4 kcal/mol) | 140.4 kJ/mol (33.5 kcal/mol) (1 atm) | **105.3 kJ/mol (25.2 kcal/mol)** |
| DLPNO-CCSD(T)/ma-def2-TZVP | SMD toluene | 84.2 kJ/mol (20.1 kcal/mol) | 97.0 kJ/mol (23.2 kcal/mol) | 101.5 kJ/mol (24.3 kcal/mol) | 144.3 kJ/mol (34.5 kcal/mol) (1 M) | **117.1 kJ/mol (28.0 kcal/mol)** |
| DLPNO-CCSD(T)/ma-def2-TZVP | gas phase | 85.9 kJ/mol (20.5 kcal/mol) | 98.7 kJ/mol (23.6 kcal/mol) | 103.2 kJ/mol (24.7 kcal/mol) | 153.9 kJ/mol (36.8 kcal/mol) (1 atm) | **118.8 kJ/mol (28.4 kcal/mol)** |

Experiment (gas phase): Ea ≈ 115.0 kJ/mol (27.5 kcal/mol).
butadiene s-trans → butadiene s-cis, ΔE: 12.8 kJ/mol (3.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (gas phase).
Thermal corrections: working-level harmonic frequencies (SMD), quasi-RRHO entropy for ΔG‡; the s-trans reference is why ΔE‡ here exceeds the barrier from the s-cis conformer.

Experiment (gas phase): butadiene + ethylene → cyclohexene, Ea ≈ 115 kJ/mol (27.5 kcal/mol). Source: D. Rowley and H. Steiner, Discuss. Faraday Soc. 1951, 10, 198–213, doi:10.1039/DF9511000198. The measured activation energy refers to equilibrium butadiene (almost all s-trans) in the gas phase; this calculation is in SMD toluene and starts from s-cis butadiene, so the s-trans → s-cis energy is reported separately. Concerted, synchronous transition state: K. N. Houk, Y.-T. Lin and F. K. Brown, J. Am. Chem. Soc. 1986, 108, 554–556, doi:10.1021/ja00263a059; E. Goldstein, B. Beno and K. N. Houk, J. Am. Chem. Soc. 1996, 118, 6036–6043, doi:10.1021/ja9601494.

*Value to be confirmed against the primary paper by J. Picione.*

## Known artifacts and limitations

- B3LYP underestimates this barrier: ΔG‡ = 115.3 kJ/mol (27.6 kcal/mol) at the protocol level versus 134.3 kJ/mol (32.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (19 kJ/mol lower). The comparison with experiment is the gas-phase Ea table above.
- Reactant complex (a weak, dispersion-bound pair of the neutral reactants) lies 8.4 kJ/mol (2.0 kcal/mol) below separated reactants on the working-level curve but only 7.3 kJ/mol (1.7 kcal/mol) below at ma-def2-TZVP, and 5.0 kJ/mol (1.2 kcal/mol) below at DLPNO-CCSD(T)/ma-def2-TZVP; the barriers above are also given from separated reactants.
- Implicit solvent only: no explicit toluene molecules.
- The optimized reactant complex used for the energies (−8.4 kJ/mol) is a different local minimum on a flat surface from the first frame of the path (−2.8 kJ/mol, the IRC end shown in the viewer and video), 5.6 kJ/mol (1.3 kcal/mol) lower (working level, relative to separated reactants).
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- The canonical HOMO changes character abruptly along the path (near frame 58); canonical orbitals are delocalized and re-ordered by energy, which is why the localized IBOs are the default view.
- The path ends at the boat form of cyclohexene (mirror-symmetric, the new CH₂–CH₂ bond eclipsed): a synchronous path keeps the mirror plane of the TS. Optimized from the IRC end it stays in the boat form (lowest vibration 38 cm⁻¹, a very flat ring motion), 23.1 kJ/mol (5.5 kcal/mol) above the half-chair, the ground-state conformer. Reaction energies use half-chair cyclohexene; the boat is the end of the path shown in the viewer and video.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/da_butadiene_ethylene.mp4`: rendered video, energy inset at the working level (`da_butadiene_ethylene_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=da_butadiene_ethylene`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
