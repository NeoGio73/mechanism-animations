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
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: σ(C–Br) → lone pair on Br 0.995 (IRC frames 17→18); lone pair on O → σ(C–O) 0.993 (IRC frames 17→18); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.827, E2 0.990, Nucleophilic addition 0.992); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | σ(C–Br) → lone pair on Br; lone pair on O → σ(C–O) |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in by J. Picione) |

## Barrier at other levels (same geometries)

| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |
|---|---|---|---|---|
| B3LYP-D3(BJ)/def2-SVP (working) | -10.9 kJ/mol (-2.6 kcal/mol) | 26.0 kJ/mol (6.2 kcal/mol) | -175.3 kJ/mol (-41.9 kcal/mol) | 18.1 kJ/mol (4.3 kcal/mol) |
| B3LYP-D3(BJ)/def2-TZVP (protocol) | 36.7 kJ/mol (8.8 kcal/mol) | 40.0 kJ/mol (9.6 kcal/mol) | -108.9 kJ/mol (-26.0 kcal/mol) | 65.7 kJ/mol (15.7 kcal/mol) |
| B3LYP/ma-def2-TZVP | 46.7 kJ/mol (11.2 kcal/mol) | 42.0 kJ/mol (10.0 kcal/mol) | -101.8 kJ/mol (-24.3 kcal/mol) | 75.7 kJ/mol (18.1 kcal/mol) |
| M06-2X/ma-def2-TZVP | 60.3 kJ/mol (14.4 kcal/mol) | 59.5 kJ/mol (14.2 kcal/mol) | -114.0 kJ/mol (-27.2 kcal/mol) | 89.4 kJ/mol (21.4 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 66.6 kJ/mol (15.9 kcal/mol) | 61.4 kJ/mol (14.7 kcal/mol) | -105.7 kJ/mol (-25.3 kcal/mol) | 95.6 kJ/mol (22.9 kcal/mol) |

ΔG adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic energies.

### Barrier from the highest point of the DLPNO-CCSD(T)/ma-def2-TZVP curve

The barrier quoted first is the one at the working-level transition state (IRC frame 17), with its thermal correction. On the DLPNO-CCSD(T)/ma-def2-TZVP curve the highest point is IRC frame 18 (s = +0.74 bohr·amu½), 4.1 kJ/mol (1.0 kcal/mol) above the TS frame (70.7 kJ/mol (16.9 kcal/mol) against 66.6 kJ/mol (15.9 kcal/mol) relative to separated reactants, electronic). The barrier from that maximum adds this difference and keeps the thermal correction of the TS frame (there are no frequencies at the maximum). Taking the barrier from the highest point of high-level single points along a lower-level path is the IRCMax practice of Malick, Petersson and Montgomery (J. Chem. Phys. 1998, 108, 5704).

| ΔG‡ (298 K, 1 M), DLPNO-CCSD(T)/ma-def2-TZVP | Value | Computed − measured |
|---|---|---|
| From the working-level TS frame (IRC frame 17) | 95.6 kJ/mol (22.9 kcal/mol) | +0.6 kJ/mol (+0.2 kcal/mol) |
| From the curve maximum (IRC frame 18; thermal correction from the TS frame) | 99.7 kJ/mol (23.8 kcal/mol) | +4.7 kJ/mol (+1.1 kcal/mol) |
| **Experiment (water, 25 °C; Mabey and Mill 1978)** | 95.0 kJ/mol (22.7 kcal/mol) | — |

The offset of the coupled-cluster maximum from the TS frame is expected to shrink in the rerun of this reaction at the standard working level, B3LYP-D3(BJ)/ma-def2-SVP.

Experiment: k(HO⁻ + CH₃Br, water, 25 °C) ≈ 1.4 × 10⁻⁴ M⁻¹ s⁻¹ (Ea = 96.3 kJ/mol (23.0 kcal/mol)). Eyring at 298.15 K (1 M standard state) gives ΔG‡ ≈ 95.0 kJ/mol (22.7 kcal/mol). Sources: W. Mabey and T. Mill, J. Phys. Chem. Ref. Data 1978, 7, 383–415, doi:10.1063/1.555572 (critical review, Tables 4.2 and 5.1; free reprint at srd.nist.gov). Values checked against the review itself: Table 4.2 gives k = 0.352 M⁻¹ s⁻¹ at 373 K, Ea = 96.3 kJ/mol (23.0 kcal/mol) and log A = 13.017; Table 5.1 gives k = 1.41 × 10⁻⁴ M⁻¹ s⁻¹ at 298 K. The 298 K value is the review's extrapolation from the higher-temperature data, not a measurement at 298 K. The review's source for the methyl halide data (its ref. 7) is I. Fells and E. A. Moelwyn-Hughes, J. Chem. Soc. 1959, 398, doi:10.1039/JR9590000398. That paper's title concerns the chlorinated methanes, so the CH₃Br coefficients may come from a companion paper, R. H. Bathgate and E. A. Moelwyn-Hughes, J. Chem. Soc. 1959, 2642, doi:10.1039/JR9590002642.

*Value to be confirmed against the primary paper by J. Picione.*

## Known artifacts and limitations

- On the working-level curve the TS lies 10.9 kJ/mol (2.6 kcal/mol) *below* separated reactants, which would suggest no barrier from separated reactants at all. This is mostly the small basis set describing the anion poorly (at def2-TZVP the TS is 36.7 kJ/mol (8.8 kcal/mol) above). The viewer and the main video therefore show the same path with DLPNO-CCSD(T)/ma-def2-TZVP energies; the working-level curve stays selectable in the viewer.
- B3LYP underestimates this barrier: ΔG‡ = 65.7 kJ/mol (15.7 kcal/mol) at the protocol level versus 95.6 kJ/mol (22.9 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (30 kJ/mol lower). Against experiment the DLPNO-CCSD(T)/ma-def2-TZVP ΔG‡ is +0.6 kJ/mol (+0.2 kcal/mol) from the TS frame and +4.7 kJ/mol (+1.1 kcal/mol) from the curve maximum (IRC frame 18; table above). The offset of the coupled-cluster maximum from the TS frame is expected to shrink in the rerun of this reaction at the standard working level, B3LYP-D3(BJ)/ma-def2-SVP.
- Ion–dipole wells: the reactant complex lies 36.9 kJ/mol (8.8 kcal/mol) below separated reactants on the working-level curve but only 3.3 kJ/mol (0.8 kcal/mol) below at def2-TZVP, and 5.2 kJ/mol (1.2 kcal/mol) above at DLPNO-CCSD(T)/ma-def2-TZVP. In water the deep double well of gas-phase ion–molecule reactions is largely washed out; what remains here is a shallow, partly artificial well (a teaching point, not a failure).
- The def2-SVP and def2-TZVP basis sets have no diffuse functions; small anions (HO⁻, Br⁻) are described less well than neutral molecules.
- Implicit solvent only: no explicit water molecules, so no specific hydrogen bonds to the base/nucleophile or leaving group.
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
