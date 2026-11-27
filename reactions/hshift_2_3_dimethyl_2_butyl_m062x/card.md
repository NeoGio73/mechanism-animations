# 1,2-Hydride shift, 2,3-dimethyl-2-butyl cation (degenerate)

| | |
|---|---|
| Reaction | (CH₃)₂C⁺–CH(CH₃)₂ → (CH₃)₂CH–C⁺(CH₃)₂, in dichloromethane |
| Level of theory | Geometries, reaction path and orbitals: M06-2X-D3(0)/ma-def2-SVP, with the SMD model of dichloromethane. Energies (the default curve in the viewer and video): DLPNO-CCSD(T)/ma-def2-TZVP single points on those geometries, SMD dichloromethane; thermal corrections from M06-2X-D3(0)/ma-def2-SVP frequencies. |
| Free energy of activation at 135 K, ΔG‡ (from the lowest conformer of the cation) | **16.1 kJ/mol (3.8 kcal/mol)** at DLPNO-CCSD(T)/ma-def2-TZVP; measured 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol), from ¹³C NMR line broadening in a superacid solution (SO₂ClF/SO₂F₂), 135 K; computed − measured +3.1 kJ/mol (+0.7 kcal/mol) |
| Free energy of activation at 298 K, ΔG‡ (same zero) | 17.0 kJ/mol (4.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; computed only |
| Reaction free energy, ΔrG | 0.0 kJ/mol (0.0 kcal/mol): the product is the same ion (a degenerate shift) |
| Top of the barrier | The top of the barrier is flat: the transition state and the symmetric hydrogen-bridged structure lie within 0.06 kJ/mol (0.014 kcal/mol) of each other (M06-2X-D3(0)/ma-def2-SVP), so the hydrogen crosses a broad plateau, not a sharp peak. |
| Transition state | one imaginary frequency, 115i cm⁻¹; C2–H 1.37 Å (forming), C3–H 1.28 Å (breaking), C2–H–C3 angle 64.1° |
| Orbitals that change | σ(C3–H) → σ(C2–H) (three-centre at the TS). Localized orbitals (IBOs), followed frame to frame |
| What this is | The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, no solvent molecules. |
| Course mapping |  |

## Limitations

1. The measurement is in a superacid solution, which the solvent model cannot describe; dichloromethane (SMD, a continuum with no solvent molecules) is used as a polar, non-nucleophilic stand-in.
2. The coupled-cluster ΔG‡ at 135 K is 3.1 kJ/mol (0.7 kcal/mol) above the measured value, a factor of about 16 in rate, because RT is only 1.12 kJ/mol at 135 K. That difference should not be over-read either: the solvent is a continuum standing in for a superacid, the geometries and thermal corrections come from M06-2X, and the choice of solvation scheme alone moves this barrier by about 0.3 kJ/mol. M06-2X/ma-def2-TZVP on the same geometries gives 17.1 kJ/mol (4.1 kcal/mol).
3. Tunnelling of the moving hydrogen is not included. It is probably small here, because the top of the barrier is flat and wide: the Wigner estimate at 135 K is a factor of about 1.06 (a rough indicator only).
4. Only half of the path is computed: from the open cation to the symmetric hydrogen-bridged structure. The second half is its mirror image (the two carbons exchange roles); those frames are labelled as mirror images in the viewer and are not computed points.

Energies at every level, how the measured value was derived, and technical notes: `card_details.md`.

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
