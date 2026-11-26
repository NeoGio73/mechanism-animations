# E2, ethoxide + 2-bromopropane

| | |
|---|---|
| Reaction | EtO⁻ + (CH₃)₂CHBr → CH₂=CHCH₃ + EtOH + Br⁻, in ethanol |
| Level of theory | Geometries, reaction path and orbitals: M06-2X-D3(0)/ma-def2-SVP, with the SMD model of ethanol. Energies (the default curve in the viewer and video): DLPNO-CCSD(T)/ma-def2-TZVP single points on those geometries, SMD ethanol; thermal corrections from M06-2X-D3(0)/ma-def2-SVP frequencies. |
| Free energy of activation, ΔG‡ (from separated reactants, 298 K, 1 M) | **99.6 kJ/mol (23.8 kcal/mol)** at DLPNO-CCSD(T)/ma-def2-TZVP; measured 105.1 kJ/mol (25.1 kcal/mol), from the E2 part of the measured rate (1 M sodium ethoxide in ethanol, 25 °C); computed − measured −5.5 kJ/mol (−1.3 kcal/mol) |
| Reaction free energy, ΔrG (from separated reactants, 298 K, 1 M) | −108.3 kJ/mol (−25.9 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; no measured value |
| Transition state | one imaginary frequency, 1389i cm⁻¹; O–H 1.29 Å (forming), Cβ–H 1.35 Å (breaking), Cα–Br 2.26 Å (breaking), Cα–Cβ 1.44 Å (becoming double), H–Cβ–Cα–Br dihedral 176.9° |
| Orbitals that change | O lone pair → σ(O–H); σ(Cβ–H) → π(C=C); σ(Cα–Br) → Br lone pair (localized orbitals, IBOs, followed frame to frame) |
| Along the path | In this calculation (implicit ethanol) the reaction is one step with no intermediate: ethoxide takes the proton first, the C=C double bond forms as the C–H bond breaks, and bromide leaves last. The C–Br bond is already stretched at the transition state (2.26 Å, from 2.00 Å), but bromide moves away only after the proton has moved: just after the transition state (0.5 bohr·amu½ further along the path), O–H is 1.06 Å while C–Br is 2.28 Å. The energy falls at every frame from the transition state to the product end. |
| What this is | The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, no solvent molecules. |
| Course mapping |  |

## Limitations

1. The solvent is a continuum (SMD ethanol): there are no ethanol molecules, so no hydrogen bonds to ethoxide or to the leaving bromide.
2. The coupled-cluster ΔG‡ is 5.5 kJ/mol (1.3 kcal/mol) below the measured value. That agreement should not be over-read: the solvent is a continuum, the geometries and thermal corrections come from M06-2X, and the choice of solvation scheme alone moves this barrier by about 3 kJ/mol. M06-2X/ma-def2-TZVP on the same geometries gives 84.4 kJ/mol (20.2 kcal/mol).
3. The solvation scheme is itself uncertain. Applying SMD directly at the coupled-cluster level (as here) and adding the M06-2X solvation energy to a gas-phase coupled-cluster energy are both reasonable choices; they differ by 3 to 7 kJ/mol (0.7 to 1.6 kcal/mol) for the reactions that start from hydroxide or ethoxide, mostly through the small anion. For this reaction the difference is 3.2 kJ/mol (0.8 kcal/mol) in ΔG‡ and 4.0 kJ/mol (1.0 kcal/mol) in ΔrG.
4. Only the anti-periplanar E2 path is computed; the competing SN2 path (about 20 % of the product in experiment) is not shown.
5. In the continuum solvent, ethoxide and 2-bromopropane form a loose complex 6.9 kJ/mol (1.6 kcal/mol) below the separated molecules (DLPNO-CCSD(T)/ma-def2-TZVP, electronic); this shallow well is partly an artifact of the model.
6. The path stops where the energy flattens. On the product side the energies use a nearby minimum of the product complex, 5.0 kJ/mol (1.2 kcal/mol) below the last frame shown (M06-2X/ma-def2-SVP).

Energies at every level, how the measured value was derived, and technical notes: `card_details.md`.

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
