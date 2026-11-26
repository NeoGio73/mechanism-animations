# Reaction card: 1,2-Hydride shift, 3-methyl-2-butyl cation → 2-methyl-2-butyl cation

| | |
|---|---|
| Reaction | 1,2-Hydride shift, 3-methyl-2-butyl cation → 2-methyl-2-butyl cation |
| Model system | CH₃CH⁺CH(CH₃)₂ → CH₃CH₂C⁺(CH₃)₂ |
| Level of theory | M06-2X-D3(0)/ma-def2-SVP, RIJCOSX; barrier at M06-2X-D3(0)/ma-def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, dichloromethane |
| Barrier (ma-def2-TZVP, electronic) | **18.1 kJ/mol (4.3 kcal/mol)** for the path shown, from the eclipsed form (start of the path); 26.2 kJ/mol (6.3 kcal/mol) for the anti form, through its own rotation transition state |
| Barrier (ma-def2-SVP, electronic) | 16.0 kJ/mol (3.8 kcal/mol) for the path shown, from the eclipsed form (start of the path); 25.0 kJ/mol (6.0 kcal/mol) for the anti form, through its own rotation transition state |
| Reaction energy (ma-def2-TZVP, electronic) | -36.6 kJ/mol (-8.7 kcal/mol), from the anti form (lowest form of the secondary cation) |
| Imaginary frequency at TS | 92i cm⁻¹ (exactly one) |
| TS geometry | C2–H 2.01 Å (forming), C3–H 1.11 Å (breaking) |
| IRC frames | 407; reaction coordinate s from -9.57 to 10.81 bohr·amu^1/2 (TS at s = 0) |
| Path quality | bond-changing region IRC frames 134–341, median step there 0.050 bohr·amu½; steps longer than 2 × that median: none |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: σ(C3–H) → σ(C2–H) 1.000 (IRC frames 290→291); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.714, E2 1.000, Tertiary hydride shift 1.000, Cyanide + acetone 1.000, SN2 1.000); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | σ(C3–H) → σ(C2–H) |
| Also found by tracking (not arrows) | Tracking also finds two bonding pairs that are shared with C2 at the start of the path and end as ordinary bonds. They are not extra arrows. The C3–CH₃ pair starts as a three-centre bond: the methyl bridges C2 and C3 (C2–CH₃ 1.70 Å, C3–CH₃ 1.80 Å), so the starting ion is not the open secondary cation of the textbook drawing. A C1–H pair leans toward C2 by hyperconjugation. Both return to ordinary bonds during the C2–C3 rotation, before the transition state and before the hydrogen moves (Mulliken shares on C2: 0.36 and 0.17 at the start, 0.03 and 0.04 at the transition state). |
| Along the path | In this calculation (implicit dichloromethane) the transition state comes before the hydrogen moves. It is a twist of about 40° about the C2–C3 bond that turns the C3–H bond toward C2. At the transition state C2–H is still 2.01 Å (2.13 Å at the start); 88% of its shortening, to 1.12 Å, comes afterwards, downhill. The other form of the secondary cation (anti) does not use this path: it reaches the tertiary cation through its own rotation transition state. |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in by J. Picione) |

## Barrier at other levels (same geometries)

| Level | ΔE‡, path shown (from the eclipsed form (start of the path)) | ΔG‡ (298 K), path shown | ΔE‡, the anti form, through its own rotation transition state | ΔG‡ (298 K), the anti form, through its own rotation transition state | Reaction energy ΔE, from the anti form (lowest form of the secondary cation) | Reaction free energy ΔG (298 K) |
|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (working) | 16.0 kJ/mol (3.8 kcal/mol) | 13.0 kJ/mol (3.1 kcal/mol) | 25.0 kJ/mol (6.0 kcal/mol) | 20.9 kJ/mol (5.0 kcal/mol) | -39.4 kJ/mol (-9.4 kcal/mol) | -47.7 kJ/mol (-11.4 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP (protocol) | 18.1 kJ/mol (4.3 kcal/mol) | 15.1 kJ/mol (3.6 kcal/mol) | 26.2 kJ/mol (6.3 kcal/mol) | 22.2 kJ/mol (5.3 kcal/mol) | -36.6 kJ/mol (-8.7 kcal/mol) | -44.9 kJ/mol (-10.7 kcal/mol) |
| M06-2X/ma-def2-TZVP | 18.1 kJ/mol (4.3 kcal/mol) | 15.1 kJ/mol (3.6 kcal/mol) | 26.2 kJ/mol (6.3 kcal/mol) | 22.2 kJ/mol (5.3 kcal/mol) | -36.6 kJ/mol (-8.7 kcal/mol) | -44.9 kJ/mol (-10.7 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 14.1 kJ/mol (3.4 kcal/mol) | 11.1 kJ/mol (2.7 kcal/mol) | 22.6 kJ/mol (5.4 kcal/mol) | 18.6 kJ/mol (4.4 kcal/mol) | -37.1 kJ/mol (-8.9 kcal/mol) | -45.4 kJ/mol (-10.9 kcal/mol) |
| B3LYP/ma-def2-TZVP | 1.4 kJ/mol (0.3 kcal/mol) | -1.6 kJ/mol (-0.4 kcal/mol) | 9.3 kJ/mol (2.2 kcal/mol) | 5.2 kJ/mol (1.2 kcal/mol) | -55.0 kJ/mol (-13.1 kcal/mol) | -63.4 kJ/mol (-15.1 kcal/mol) |

ΔG adds the working-level thermal corrections (298.15 K) to each level's electronic energies; there is one molecule throughout, so no standard-state correction applies.

## Measured value

No measured value on this card. The 3-methyl-2-butyl cation has not been observed: it rearranges too fast. For context only, not as a comparison: for the 1,2-hydride shift in the 2-butyl cation (a different, degenerate shift between two secondary carbons), Saunders and Kates give an upper limit, ΔG‡ < 10 kJ/mol (2.4 kcal/mol), from the absence of ¹³C NMR line broadening at −140 °C (M. Saunders and M. R. Kates, J. Am. Chem. Soc. 1978, 100, 7082–7083; reviewed in Saunders and Jiménez-Vázquez 1991).

## Known artifacts and limitations

- M06-2X overestimates this barrier: ΔG‡ = 15.1 kJ/mol (3.6 kcal/mol) at M06-2X-D3(0)/ma-def2-TZVP versus 11.1 kJ/mol (2.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (4 kJ/mol higher).
- Implicit solvent only: no explicit dichloromethane molecules.
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = 10.81 bohr·amu½ (IRC step 217), after the bonds had finished changing; the product-end minimum was optimized from that last point.
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = -9.57 bohr·amu½ (IRC step 189), after the bonds had finished changing; the reactant-end minimum was optimized from that last point.
- Two forms of the secondary cation are minima at this level, anti and eclipsed, and both are partly methyl-bridged. A third arrangement, with C3–H in line with the empty orbital on C2, is not a minimum: it goes downhill to the tertiary cation. The path shown starts from the eclipsed form. The anti form is 2.8 kJ/mol (0.7 kcal/mol) lower in free energy at the working level and has its own, higher rotation barrier (18.6 against 11.1 kJ/mol at DLPNO-CCSD(T)/ma-def2-TZVP, 298 K).

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/hshift_3_methyl_2_butyl.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`hshift_3_methyl_2_butyl_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=hshift_3_methyl_2_butyl`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
