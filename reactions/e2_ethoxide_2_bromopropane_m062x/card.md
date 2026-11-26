# Reaction card: E2, ethoxide + 2-bromopropane

| | |
|---|---|
| Reaction | E2, ethoxide + 2-bromopropane |
| Model system | EtO⁻ + (CH₃)₂CHBr → CH₂=CHCH₃ + EtOH + Br⁻ |
| Level of theory | M06-2X-D3(0)/ma-def2-SVP, RIJCOSX; barrier at M06-2X-D3(0)/ma-def2-TZVP, RIJCOSX (single points on the working-level geometries) |
| Solvent model | CPCM/SMD, ethanol |
| Barrier (ma-def2-TZVP, electronic) | **66.4 kJ/mol (15.9 kcal/mol)** relative to the reactant complex; 56.8 kJ/mol (13.6 kcal/mol) relative to separated reactants |
| Barrier (ma-def2-SVP, electronic, the path shown) | 59.7 kJ/mol (14.3 kcal/mol) relative to the reactant complex; 50.0 kJ/mol (12.0 kcal/mol) relative to separated reactants |
| Reaction energy (ma-def2-TZVP, electronic) | -88.4 kJ/mol (-21.1 kcal/mol) |
| Imaginary frequency at TS | 1389i cm⁻¹ (exactly one) |
| TS geometry | O–H 1.29 Å (forming), Cβ–H 1.35 Å (breaking), Cα–Br 2.26 Å (breaking), Cα–Cβ 1.44 Å (becoming double), H–Cβ–Cα–Br dihedral 176.9° |
| IRC frames | 532; reaction coordinate s from -6.59 to 19.95 bohr·amu^1/2 (TS at s = 0) |
| Path quality | bond-changing region IRC frames 18–497, median step there 0.050 bohr·amu½; steps longer than 2 × that median: none |
| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^-3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: O lone pair → σ(O–H) 1.000 (IRC frames 131→132); σ(Cβ–H) → π(C=C) 1.000 (IRC frames 136→137); σ(Cα–Br) → Br lone pair 1.000 (IRC frames 141→142); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.714, Tertiary hydride shift 1.000, Secondary hydride shift 1.000, Cyanide + acetone 1.000, SN2 1.000); canonical HOMO/LUMO optional in the viewer |
| Orbitals that change | O lone pair → σ(O–H); σ(Cβ–H) → π(C=C); σ(Cα–Br) → Br lone pair |
| Along the path | In this calculation (implicit ethanol) the reaction is one step with no intermediate: ethoxide takes the proton first, the C=C double bond forms as the C–H bond breaks, and bromide leaves last. The C–Br bond is already stretched at the transition state (2.26 Å, from 2.00 Å), but bromide moves away only after the proton has moved: just after the transition state (0.5 bohr·amu½ further along the path), O–H is 1.06 Å while C–Br is 2.28 Å. The energy falls at every frame from the transition state to the product end. |
| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, no explicit solvent molecules |
| Course mapping | McMurry chapter __ (to be filled in) |

## Barrier at other levels (same geometries)

| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ from separated reactants (298 K, 1 M) | ΔG‡ from reactant complex (298 K) | Reaction free energy ΔG (298 K, 1 M) |
|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (working) | 50.0 kJ/mol (12.0 kcal/mol) | 59.7 kJ/mol (14.3 kcal/mol) | -92.6 kJ/mol (-22.1 kcal/mol) | 77.6 kJ/mol (18.5 kcal/mol) | 48.4 kJ/mol (11.6 kcal/mol) | -119.3 kJ/mol (-28.5 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP (protocol) | 56.8 kJ/mol (13.6 kcal/mol) | 66.4 kJ/mol (15.9 kcal/mol) | -88.4 kJ/mol (-21.1 kcal/mol) | 84.4 kJ/mol (20.2 kcal/mol) | 55.1 kJ/mol (13.2 kcal/mol) | -115.1 kJ/mol (-27.5 kcal/mol) |
| M06-2X/ma-def2-TZVP | 56.8 kJ/mol (13.6 kcal/mol) | 66.4 kJ/mol (15.9 kcal/mol) | -88.4 kJ/mol (-21.1 kcal/mol) | 84.4 kJ/mol (20.2 kcal/mol) | 55.1 kJ/mol (13.2 kcal/mol) | -115.1 kJ/mol (-27.5 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 72.1 kJ/mol (17.2 kcal/mol) | 79.0 kJ/mol (18.9 kcal/mol) | -81.5 kJ/mol (-19.5 kcal/mol) | 99.6 kJ/mol (23.8 kcal/mol) | 67.7 kJ/mol (16.2 kcal/mol) | -108.3 kJ/mol (-25.9 kcal/mol) |
| B3LYP/ma-def2-TZVP | 48.6 kJ/mol (11.6 kcal/mol) | 54.3 kJ/mol (13.0 kcal/mol) | -88.4 kJ/mol (-21.1 kcal/mol) | 76.2 kJ/mol (18.2 kcal/mol) | 43.0 kJ/mol (10.3 kcal/mol) | -115.1 kJ/mol (-27.5 kcal/mol) |

ΔG adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M standard-state correction to each level's electronic energies; the reactant complex → TS step is unimolecular and needs no standard-state correction.

### Barrier from the highest point of the DLPNO-CCSD(T)/ma-def2-TZVP curve

The barrier quoted first is the one at the working-level transition state (IRC frame 132), with its thermal correction. On the DLPNO-CCSD(T)/ma-def2-TZVP curve the highest point is IRC frame 134 (s = +0.10 bohr·amu½), 1.1 kJ/mol (0.3 kcal/mol) above the TS frame (73.1 kJ/mol (17.5 kcal/mol) against 72.1 kJ/mol (17.2 kcal/mol) relative to separated reactants, electronic). The barrier from that maximum adds this difference and keeps the thermal correction of the TS frame (there are no frequencies at the maximum). Taking the barrier from the highest point of high-level single points along a lower-level path is the IRCMax practice of Malick, Petersson and Montgomery (J. Chem. Phys. 1998, 108, 5704).

| ΔG‡ (298 K, 1 M), DLPNO-CCSD(T)/ma-def2-TZVP | Value | Computed − measured |
|---|---|---|
| From the working-level TS frame (IRC frame 132) | 99.6 kJ/mol (23.8 kcal/mol) | — |
| From the curve maximum (IRC frame 134; thermal correction from the TS frame) | 100.7 kJ/mol (24.1 kcal/mol) | — |

Experiment: 2-bromopropane + 1 N NaOEt in dry ethanol, 25.0 °C (M. L. Dhar, E. D. Hughes and C. K. Ingold, J. Chem. Soc. 1948, 2058–2065, doi:10.1039/JR9480002058; Mechanism of elimination reactions, Part X). Table I: overall second-order 10⁵k₂ = 0.295 L mol⁻¹ s⁻¹ (k₂ = 2.95 × 10⁻⁶), corrected for solvent expansion, with 80.3 % olefin (two determinations, 82.7 and 78.0 %), the rest SN2 ether. Table VII's 10⁶k₂ = 2.88 is the uncorrected mean of a single run and is not used. Table II: 10⁵k(E2) = 0.237 L mol⁻¹ s⁻¹, i.e. k(E2) = 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹. Caveats stated in the paper: the first-order (solvolysis) processes contribute up to 2–3 % of the rate, and the olefin proportions are good to within 3 %, with a correction for loss of the volatile propene. Values checked against the original paper by J. Picione (2026-11-13; again 2026-11-20, with the table details above). See also E. D. Hughes, C. K. Ingold, S. Masterman and B. J. McNulty, J. Chem. Soc. 1940, 899–912, doi:10.1039/JR9400000899 (Part V), and M. L. Dhar, E. D. Hughes, C. K. Ingold, A. M. M. Mandour, G. A. Maw and L. I. Woolf, J. Chem. Soc. 1948, 2093, doi:10.1039/JR9480002093 (Part XVI). Eyring at 298.15 K (1 M standard state): k(E2) = 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹ gives ΔG‡(E2) = 105.1 kJ/mol (25.1 kcal/mol); the SN2 share (0.58 × 10⁻⁶ L mol⁻¹ s⁻¹) gives ΔG‡(SN2) = 108.6 kJ/mol (26.0 kcal/mol). The two pathways differ by only about 3.5 kJ/mol (0.8 kcal/mol), which is why a secondary bromide with ethoxide gives mostly, not only, alkene.

## Known artifacts and limitations

- M06-2X underestimates this barrier: ΔG‡ = 84.4 kJ/mol (20.2 kcal/mol) at the protocol level versus 99.6 kJ/mol (23.8 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (15 kJ/mol lower). The experimental value is given above; the close agreement of the coupled-cluster number is partly fortuitous (implicit solvent, working-level geometries and thermal corrections).
- Ion–dipole wells: the reactant complex lies 9.7 kJ/mol (2.3 kcal/mol) below separated reactants on the working-level curve but only 9.6 kJ/mol (2.3 kcal/mol) below at ma-def2-TZVP, and 6.9 kJ/mol (1.6 kcal/mol) below at DLPNO-CCSD(T)/ma-def2-TZVP. In ethanol the deep double well of gas-phase ion–molecule reactions is largely washed out; what remains here is a shallow, partly artificial well (a teaching point, not a failure).
- Implicit solvent only: no explicit ethanol molecules, so no specific hydrogen bonds to the base/nucleophile or leaving group.
- The optimized product complex used for the energies (−109.9 kJ/mol) is a different local minimum on a flat surface from the last frame of the path (−104.9 kJ/mol, the IRC end shown in the viewer and video), 5.0 kJ/mol (1.2 kcal/mol) lower (working level, relative to separated reactants).
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = 19.95 bohr·amu½ (IRC step 399), after the bonds had finished changing; beyond it only the separating fragments drift apart. The product complex was optimized from that last point.
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = -6.59 bohr·amu½ (IRC step 132), after the bonds had finished changing; beyond it only the separating fragments drift apart. The reactant complex was optimized from that last point.
- One conformer only: the path starts from a single hand-built TS guess; other rotamers of the reactants were not searched.
- The single-point curves peak just after the working-level transition state because the single points sit on a path optimized at another level: M06-2X/ma-def2-TZVP one frame after it (frame 133, s = +0.05 bohr·amu½, 0.05 kJ/mol (0.01 kcal/mol) above the TS frame) and DLPNO-CCSD(T)/ma-def2-TZVP two frames after it (frame 134, s = +0.10 bohr·amu½, O–H 1.243 Å, 1.1 kJ/mol (0.3 kcal/mol) above the TS frame). (Corrected 2026-11-25 from the final path; the earlier text, frame 28 and 0.6 kJ/mol with the coupled-cluster maximum on the TS frame, came from an earlier folder.)
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
