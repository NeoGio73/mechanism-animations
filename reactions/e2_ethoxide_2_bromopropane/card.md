# Reaction card: E2, ethoxide + 2-bromopropane

| | |
|---|---|
| Reaction | E2, ethoxide + 2-bromopropane |
| Model system | EtO⁻ + (CH₃)₂CHBr → CH₂=CHCH₃ + EtOH + Br⁻ |
| Level of theory | B3LYP-D3(BJ)/ma-def2-SVP, RIJCOSX; barrier at B3LYP-D3(BJ)/ma-def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, ethanol |
| Barrier (ma-def2-TZVP, electronic) | **51.1 kJ/mol (12.2 kcal/mol)** relative to the reactant complex; 44.2 kJ/mol (10.6 kcal/mol) relative to separated reactants |
| Barrier (ma-def2-SVP, electronic, the path shown) | 48.3 kJ/mol (11.6 kcal/mol) relative to the reactant complex; 39.9 kJ/mol (9.5 kcal/mol) relative to separated reactants |
| Reaction energy (ma-def2-TZVP, electronic) | -88.0 kJ/mol (-21.0 kcal/mol) |
| Imaginary frequency at TS | 888i cm⁻¹ (exactly one) |
| TS geometry | O–H 1.39 Å (forming), Cβ–H 1.28 Å (breaking), Cα–Br 2.46 Å (breaking), Cα–Cβ 1.43 Å (becoming double), H–Cβ–Cα–Br dihedral 176.3° |
| IRC frames | 146; reaction coordinate s from -6.38 to 20.58 bohr·amu^1/2 (TS at s = 0) |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of the shown orbitals 0.99); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | O lone pair → σ(O–H); σ(Cβ–H) → π(C=C); σ(Cα–Br) → Br lone pair |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in) |

## Barrier at other levels (same geometries)

| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |
|---|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP (working) | 39.9 kJ/mol (9.5 kcal/mol) | 48.3 kJ/mol (11.6 kcal/mol) | -84.9 kJ/mol (-20.3 kcal/mol) | 68.5 kJ/mol (16.4 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP (protocol) | 44.2 kJ/mol (10.6 kcal/mol) | 51.1 kJ/mol (12.2 kcal/mol) | -88.0 kJ/mol (-21.0 kcal/mol) | 72.9 kJ/mol (17.4 kcal/mol) |
| M06-2X/ma-def2-TZVP | 61.6 kJ/mol (14.7 kcal/mol) | 70.9 kJ/mol (16.9 kcal/mol) | -88.8 kJ/mol (-21.2 kcal/mol) | 90.2 kJ/mol (21.6 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 78.6 kJ/mol (18.8 kcal/mol) | 86.6 kJ/mol (20.7 kcal/mol) | -82.5 kJ/mol (-19.7 kcal/mol) | 107.2 kJ/mol (25.6 kcal/mol) |

ΔG‡ adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic barrier.

Experiment: 2-bromopropane + NaOEt in ethanol, 25 °C: overall second-order k ≈ 2.88 × 10⁻⁶ M⁻¹ s⁻¹ (mean of 19 points, 2.76–2.96 × 10⁻⁶ M⁻¹ s⁻¹), 80.3 % propene, the rest SN2 ether. Source: M. L. Dhar, E. D. Hughes and C. K. Ingold, J. Chem. Soc. 1948, 2058–2065 (Tables III and VII), as quoted by C. G. Jesudason, arXiv:1101.5060; the primary paper was not re-checked here. Eyring at 298.15 K (1 M standard state): E2 share k(E2) ≈ 2.31 × 10⁻⁶ M⁻¹ s⁻¹ gives ΔG‡(E2) ≈ 105.2 kJ/mol (25.1 kcal/mol); the SN2 share gives ΔG‡(SN2) ≈ 108.7 kJ/mol (26.0 kcal/mol). The two pathways differ by only about 3.5 kJ/mol (0.8 kcal/mol), which is why a secondary bromide with ethoxide gives mostly, not only, alkene.

## Known artifacts and limitations

- B3LYP underestimates this barrier: ΔG‡ = 72.9 kJ/mol (17.4 kcal/mol) at the protocol level versus 107.2 kJ/mol (25.6 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (34 kJ/mol lower). The experimental value is given above; the close agreement of the coupled-cluster number is partly fortuitous (implicit solvent, working-level geometries and thermal corrections).
- Ion–dipole wells: the reactant complex lies 8.5 kJ/mol (2.0 kcal/mol) below separated reactants on the working-level curve but only 6.9 kJ/mol (1.6 kcal/mol) below at ma-def2-TZVP, and 8.1 kJ/mol (1.9 kcal/mol) below at DLPNO-CCSD(T)/ma-def2-TZVP. In ethanol the deep double well of gas-phase ion–molecule reactions is largely washed out; what remains here is a shallow, partly artificial well (a teaching point, not a failure).
- Implicit solvent only: no explicit ethanol molecules, so no specific hydrogen bonds to the base/nucleophile or leaving group.
- The optimized reactant complex used for the energies (−8.5 kJ/mol) is a different local minimum on a flat surface from the first frame of the path (−2.1 kJ/mol, the IRC end shown in the viewer and video), 6.3 kJ/mol (1.5 kcal/mol) lower (working level, relative to separated reactants).
- The optimized product complex used for the energies (−103.3 kJ/mol) is a different local minimum on a flat surface from the last frame of the path (−91.3 kJ/mol, the IRC end shown in the viewer and video), 12.0 kJ/mol (2.9 kcal/mol) lower (working level, relative to separated reactants).
- The IRC point 1 step from the TS on the reactant side (s = -0.37 bohr·amu½) lies 1.3 kJ/mol (0.3 kcal/mol) above the TS: a straight-line first step off a curved hydrogen-transfer path. It is kept in the data files (02_irc/energies.csv) but left out of the viewer path and the video, so the displayed maximum is the TS.
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by decision at s = 20.58 bohr·amu½ (IRC step 115), after the bonds had finished changing; beyond it only the separating fragments drift apart. The product complex was optimized from that last point.
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- Only the anti-periplanar E2 path is computed; the competing SN2 path (about 20 % of the product in experiment) is not shown.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/e2_ethoxide_2_bromopropane.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`e2_ethoxide_2_bromopropane_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=e2_ethoxide_2_bromopropane`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
