# Card details: 1,2-Hydride shift, 3-methyl-2-butyl cation → 2-methyl-2-butyl cation

## Energies at every level (same M06-2X-D3(0)/ma-def2-SVP geometries, SMD dichloromethane)

| Level | ΔE‡, path shown (from the eclipsed form (start of the path)) | ΔG‡ (298 K), path shown | ΔE‡, the anti form, through its own rotation transition state | ΔG‡ (298 K), the anti form, through its own rotation transition state | Reaction energy ΔrE, from the anti form (lowest form of the secondary cation) | Reaction free energy ΔrG (298 K) |
|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (geometries) | 16.0 kJ/mol (3.8 kcal/mol) | 13.0 kJ/mol (3.1 kcal/mol) | 25.0 kJ/mol (6.0 kcal/mol) | 20.9 kJ/mol (5.0 kcal/mol) | −39.4 kJ/mol (−9.4 kcal/mol) | −47.7 kJ/mol (−11.4 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP | 18.1 kJ/mol (4.3 kcal/mol) | 15.1 kJ/mol (3.6 kcal/mol) | 26.2 kJ/mol (6.3 kcal/mol) | 22.2 kJ/mol (5.3 kcal/mol) | −36.6 kJ/mol (−8.7 kcal/mol) | −44.9 kJ/mol (−10.7 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP (default curve) | 14.1 kJ/mol (3.4 kcal/mol) | 11.1 kJ/mol (2.7 kcal/mol) | 22.6 kJ/mol (5.4 kcal/mol) | 18.6 kJ/mol (4.4 kcal/mol) | −37.1 kJ/mol (−8.9 kcal/mol) | −45.4 kJ/mol (−10.9 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP | 1.4 kJ/mol (0.3 kcal/mol) | −1.6 kJ/mol (−0.4 kcal/mol) | 9.3 kJ/mol (2.2 kcal/mol) | 5.2 kJ/mol (1.2 kcal/mol) | −55.0 kJ/mol (−13.1 kcal/mol) | −63.4 kJ/mol (−15.1 kcal/mol) |

ΔG adds the M06-2X-D3(0)/ma-def2-SVP thermal corrections (298.15 K, quasi-RRHO) to each level's electronic energies; there is one molecule throughout, so no standard-state correction applies.

## Measured value

No measured value on this card. The 3-methyl-2-butyl cation has not been observed: it rearranges too fast. For context only, not as a comparison: for the 1,2-hydride shift in the 2-butyl cation (a different, degenerate shift between two secondary carbons), Saunders and Kates give an upper limit, ΔG‡ < 10 kJ/mol (2.4 kcal/mol), from the absence of ¹³C NMR line broadening at −140 °C (M. Saunders and M. R. Kates, J. Am. Chem. Soc. 1978, 100, 7082–7083; reviewed in Saunders and Jiménez-Vázquez 1991).

## Path and orbitals

- IRC frames: 407; reaction coordinate s from −9.57 to 10.81 bohr·amu½ (TS at s = 0).
- Path quality: bond-changing region IRC frames 134–341, median step there 0.050 bohr·amu½; steps longer than 2 × that median: none.
- Orbitals: IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^−3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: σ(C3–H) → σ(C2–H) 1.000 (IRC frames 290→291); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.714, E2 1.000, Tertiary hydride shift 1.000, Cyanide + acetone 1.000, SN2 1.000); canonical HOMO/LUMO optional in the viewer.

## Technical notes

- M06-2X overestimates this barrier: ΔG‡ = 15.1 kJ/mol (3.6 kcal/mol) at M06-2X-D3(0)/ma-def2-TZVP versus 11.1 kJ/mol (2.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (4 kJ/mol higher).
- Implicit solvent only: no explicit dichloromethane molecules.
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = 10.81 bohr·amu½ (IRC step 217), after the bonds had finished changing; the product-end minimum was optimized from that last point.
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = −9.57 bohr·amu½ (IRC step 189), after the bonds had finished changing; the reactant-end minimum was optimized from that last point.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu½), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/hshift_3_methyl_2_butyl.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`hshift_3_methyl_2_butyl_<level>.mp4`: the same video with the other energy curves)
- `web/viewer/index.html?r=hshift_3_methyl_2_butyl`: interactive viewer
- `run_log.md`: every ORCA job and rerun

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
