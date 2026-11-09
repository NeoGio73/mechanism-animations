# Reaction card: SN2, hydroxide + bromomethane

| | |
|---|---|
| Reaction | SN2, hydroxide + bromomethane |
| Model system | HO⁻ + CH₃Br → CH₃OH + Br⁻ |
| Level of theory | B3LYP-D3(BJ)/def2-SVP, RIJCOSX; barrier at B3LYP-D3(BJ)/def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, water |
| Barrier (def2-TZVP, electronic) | **40.0 kJ/mol (9.6 kcal/mol)** relative to the reactant complex; 36.7 kJ/mol (8.8 kcal/mol) relative to separated reactants |
| Barrier (def2-SVP, electronic, the path shown) | 26.0 kJ/mol (6.2 kcal/mol) relative to the reactant complex; -10.9 kJ/mol (-2.6 kcal/mol) relative to separated reactants |
| Reaction energy (def2-TZVP, electronic) | -108.9 kJ/mol (-26.0 kcal/mol) |
| Imaginary frequency at TS | 441i cm⁻¹ (exactly one) |
| TS geometry | C–O 2.23 Å, C–Br 2.30 Å, angle 178.6° |
| IRC frames | 60; reaction coordinate s from -3.64 to 11.02 bohr·amu^1/2 (TS at s = 0) |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of the shown orbitals 0.99); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | σ(C–Br) → lone pair on Br; lone pair on O → σ(C–O) |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in by J. Picione) |

## Barrier at other levels (same geometries)

| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |
|---|---|---|---|---|
| B3LYP/def2-SVP (working) | -10.9 kJ/mol (-2.6 kcal/mol) | 26.0 kJ/mol (6.2 kcal/mol) | -175.3 kJ/mol (-41.9 kcal/mol) | 18.1 kJ/mol (4.3 kcal/mol) |
| B3LYP/def2-TZVP (protocol) | 36.7 kJ/mol (8.8 kcal/mol) | 40.0 kJ/mol (9.6 kcal/mol) | -108.9 kJ/mol (-26.0 kcal/mol) | 65.7 kJ/mol (15.7 kcal/mol) |
| B3LYP/ma-def2-TZVP | 46.7 kJ/mol (11.2 kcal/mol) | 42.0 kJ/mol (10.0 kcal/mol) | -101.8 kJ/mol (-24.3 kcal/mol) | 75.7 kJ/mol (18.1 kcal/mol) |
| M06-2X/ma-def2-TZVP | 60.3 kJ/mol (14.4 kcal/mol) | 59.5 kJ/mol (14.2 kcal/mol) | -114.0 kJ/mol (-27.2 kcal/mol) | 89.4 kJ/mol (21.4 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 66.6 kJ/mol (15.9 kcal/mol) | 61.4 kJ/mol (14.7 kcal/mol) | -105.7 kJ/mol (-25.3 kcal/mol) | 95.6 kJ/mol (22.9 kcal/mol) |

ΔG‡ adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic barrier.

Experiment: k(HO⁻ + CH₃Br, water, 25 °C) ≈ 1.4 × 10⁻⁴ M⁻¹ s⁻¹ (extrapolated from higher-temperature data; Ea = 96.3 kJ/mol), from the critical review W. Mabey and T. Mill, J. Phys. Chem. Ref. Data 1978, 7, 383–415, doi:10.1063/1.555572, Tables 4.2 and 5.1. Eyring at 298.15 K (1 M standard state) gives ΔG‡ ≈ 95.0 kJ/mol (22.7 kcal/mol). The primary measurements are Moelwyn-Hughes and co-workers (J. Chem. Soc. 1959), not re-checked here.

## Known artifacts and limitations

- On the working-level curve the TS lies 10.9 kJ/mol (2.6 kcal/mol) *below* separated reactants, which would suggest no barrier from separated reactants at all. This is mostly the small basis set describing HO⁻ poorly (at def2-TZVP the TS is 36.7 kJ/mol (8.8 kcal/mol) above). The viewer and the main video therefore show the same path with DLPNO-CCSD(T)/ma-def2-TZVP energies; the working-level curve stays selectable in the viewer.
- B3LYP underestimates this barrier: ΔG‡ = 65.7 kJ/mol (15.7 kcal/mol) at the protocol level versus 95.6 kJ/mol (22.9 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (30 kJ/mol lower). The experimental value is given above; the close agreement of the coupled-cluster number is partly fortuitous (implicit solvent, working-level geometries and thermal corrections).
- Ion–dipole wells: the reactant complex lies 36.9 kJ/mol (8.8 kcal/mol) below separated reactants on the working-level curve but only 3.3 kJ/mol (0.8 kcal/mol) below at def2-TZVP, and 5.2 kJ/mol (1.2 kcal/mol) above at DLPNO-CCSD(T)/ma-def2-TZVP. In water the double well of the gas-phase SN2 profile is largely washed out; what remains here is a shallow, partly artificial well (a teaching point, not a failure).
- The def2-SVP and def2-TZVP basis sets have no diffuse functions; small anions (HO⁻) are described less well than neutral molecules.
- Implicit solvent only: no explicit water molecules, so no specific hydrogen bonds to the nucleophile or leaving group.
- The canonical HOMO changes character abruptly along the path (near frame 54); canonical orbitals are delocalized and re-ordered by energy, which is why the localized IBOs are the default view.
- The pictures (geometries, IRC path, orbitals) are at the working level; only the energies change between the curves offered in the viewer. Single points along a path computed at another level can peak one or two frames away from the working-level TS.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/sn2_hydroxide_bromomethane.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`sn2_hydroxide_bromomethane_<level>.mp4`: same video with the other energy curves)
- `web/viewer/index.html?r=sn2_hydroxide_bromomethane`: interactive viewer
- `run_log.md`: every ORCA job, wall time, and reruns

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
