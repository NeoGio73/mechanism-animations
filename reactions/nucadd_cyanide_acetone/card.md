# Reaction card: Nucleophilic addition, cyanide + acetone

| | |
|---|---|
| Reaction | Nucleophilic addition, cyanide + acetone |
| Model system | N≡C⁻ + (CH₃)₂C=O → (CH₃)₂C(O⁻)C≡N |
| Level of theory | B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX; barrier at B3LYP-D3(BJ)/ma-def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, water |
| Barrier (ma-def2-TZVP, electronic) | **44.7 kJ/mol (10.7 kcal/mol)** relative to the reactant complex; 37.2 kJ/mol (8.9 kcal/mol) relative to separated reactants |
| Barrier (ma-def2-SVP, electronic, the path shown) | 39.2 kJ/mol (9.4 kcal/mol) relative to the reactant complex; 30.9 kJ/mol (7.4 kcal/mol) relative to separated reactants |
| Reaction energy (ma-def2-TZVP, electronic) | -3.6 kJ/mol (-0.9 kcal/mol) |
| Imaginary frequency at TS | 319i cm⁻¹ (exactly one) |
| TS geometry | C(N)–C 2.03 Å (forming), C–O 1.28 Å (becoming single), Bürgi–Dunitz angle C(N)···C=O 106.4° |
| IRC frames | 120; reaction coordinate s from -5.92 to 3.38 bohr·amu^1/2 (TS at s = 0) |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: C lone pair (cyanide) → σ(C–C) 0.992 (IRC frames 100→101); π(C=O) → O lone pair 0.994 (IRC frames 100→101); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.827, E2 0.990, SN2 0.993); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | C lone pair (cyanide) → σ(C–C); π(C=O) → O lone pair |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in by J. Picione) |

## Known artifacts and limitations

- Ion–dipole wells: the reactant complex lies 8.3 kJ/mol (2.0 kcal/mol) below separated reactants on the working-level curve but only 7.5 kJ/mol (1.8 kcal/mol) below at ma-def2-TZVP. In water the deep double well of gas-phase ion–molecule reactions is largely washed out; what remains here is a shallow, partly artificial well (a teaching point, not a failure).
- Implicit solvent only: no explicit water molecules, so no specific hydrogen bonds to the nucleophile or the product anion.
- The optimized reactant complex used for the energies (−8.3 kJ/mol) is a different local minimum on a flat surface from the first frame of the path (−2.3 kJ/mol, the IRC end shown in the viewer and video), 6.0 kJ/mol (1.4 kcal/mol) lower (working level, relative to separated reactants).
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = -5.92 bohr·amu½ (IRC step 99), after the bonds had finished changing; beyond it only the separating fragments drift apart. The reactant complex was optimized from that last point.
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- The canonical HOMO changes character abruptly along the path (near frame 81); canonical orbitals are delocalized and re-ordered by energy, which is why the localized IBOs are the default view.
- The path ends at the alkoxide, (CH₃)₂C(CN)O⁻, which lies above separated cyanide + acetone in free energy (ΔG = +34.4 kJ/mol (+8.2 kcal/mol) at the working level, 298 K, 1 M; electronically it is 14.5 kJ/mol (3.5 kcal/mol) below): the addition step on its own is uphill. The neutral cyanohydrin forms only after a later proton transfer (from HCN or water) that this path does not model.
- No measured rate is shown yet. Any measured rate depends on pH, because only free cyanide (CN⁻), not HCN, adds to the carbonyl; leads are listed in docs/references.md.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/nucadd_cyanide_acetone.mp4`: rendered video, energy inset at the working level (`nucadd_cyanide_acetone_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=nucadd_cyanide_acetone`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
