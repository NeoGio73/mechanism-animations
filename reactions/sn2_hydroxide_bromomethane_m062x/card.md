# SN2, hydroxide + bromomethane

| | |
|---|---|
| Reaction | HO⁻ + CH₃Br → CH₃OH + Br⁻, in water |
| Level of theory | Geometries, reaction path and orbitals: M06-2X-D3(0)/ma-def2-SVP, with the SMD model of water. Energies (the default curve in the viewer and video): DLPNO-CCSD(T)/ma-def2-TZVP single points on those geometries, SMD water; thermal corrections from M06-2X-D3(0)/ma-def2-SVP frequencies. |
| Free energy of activation, ΔG‡ (from separated reactants, 298 K, 1 M) | **100.4 kJ/mol (24.0 kcal/mol)** at DLPNO-CCSD(T)/ma-def2-TZVP; measured 94.9 kJ/mol (22.7 kcal/mol), from rates measured in water at 292.69–322.91 K (Eyring fit); computed − measured +5.5 kJ/mol (+1.3 kcal/mol) |
| Reaction free energy, ΔrG (from separated reactants, 298 K, 1 M) | −85.3 kJ/mol (−20.4 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; measured −100.8 kJ/mol (−24.1 kcal/mol), from tabulated standard formation free energies (water); computed − measured +15.5 kJ/mol (+3.7 kcal/mol) |
| Transition state | one imaginary frequency, 619i cm⁻¹; C–O 2.11 Å, C–Br 2.32 Å, angle 178.7° |
| Orbitals that change | lone pair on O → σ(C–O); σ(C–Br) → lone pair on Br (localized orbitals, IBOs, followed frame to frame) |
| What this is | The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, no solvent molecules. |
| Course mapping |  |

## Limitations

1. The solvent is a continuum (SMD water): there are no water molecules, so no hydrogen bonds to hydroxide or to the leaving bromide.
2. The coupled-cluster ΔG‡ is 5.5 kJ/mol (1.3 kcal/mol) above the measured value. That agreement should not be over-read: the solvent is a continuum, the geometries and thermal corrections come from M06-2X, and the choice of solvation scheme alone moves this barrier by about 3 kJ/mol. M06-2X/ma-def2-TZVP on the same geometries gives 94.1 kJ/mol (22.5 kcal/mol).
3. The solvation scheme is itself uncertain. Applying SMD directly at the coupled-cluster level (as here) and adding the M06-2X solvation energy to a gas-phase coupled-cluster energy are both reasonable choices; they differ by 3 to 7 kJ/mol (0.7 to 1.6 kcal/mol) for the reactions that start from hydroxide or ethoxide, mostly through the small anion. For this reaction the difference is 3.0 kJ/mol (0.7 kcal/mol) in ΔG‡ and 6.6 kJ/mol (1.6 kcal/mol) in ΔrG.
4. The computed reaction free energy is 15.5 kJ/mol (3.7 kcal/mol) less negative than the measured value. The measured value is not a single measurement: it is built from tabulated formation free energies, with bromomethane in water estimated from its gas-phase value and its solubility (± 0.4 kJ/mol).
5. In the continuum solvent, hydroxide and bromomethane form a loose complex 5.4 kJ/mol (1.3 kcal/mol) below the separated molecules (DLPNO-CCSD(T)/ma-def2-TZVP, electronic); this shallow well is partly an artifact of the model.
6. The pictures (geometries, path and orbitals) are at M06-2X/ma-def2-SVP; only the energies change between the curves in the viewer. A curve computed at another level can peak a frame or two away from the transition state shown.

Energies at every level, how the measured value was derived, and technical notes: `card_details.md`.

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
