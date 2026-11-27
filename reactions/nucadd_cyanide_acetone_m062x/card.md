# Nucleophilic addition, cyanide + acetone

| | |
|---|---|
| Reaction | N≡C⁻ + (CH₃)₂C=O → (CH₃)₂C(O⁻)C≡N, in water |
| Level of theory | Geometries, reaction path and orbitals: M06-2X-D3(0)/ma-def2-SVP, with the SMD model of water. Energies (the default curve in the viewer and video): DLPNO-CCSD(T)/ma-def2-TZVP single points on those geometries, SMD water; thermal corrections from M06-2X-D3(0)/ma-def2-SVP frequencies. |
| Free energy of activation, ΔG‡ (from separated reactants, 298 K, 1 M) | **79.3 kJ/mol (19.0 kcal/mol)** at DLPNO-CCSD(T)/ma-def2-TZVP; measured 69.6 kJ/mol (16.6 kcal/mol), from the measured rate constant (water, 25 °C); computed − measured +9.7 kJ/mol (+2.3 kcal/mol) |
| Reaction free energy, ΔrG (from separated reactants, 298 K, 1 M) | +36.4 kJ/mol (+8.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; measured +10.5 kJ/mol (+2.5 kcal/mol), from the measured equilibrium constant (water, 25 °C); computed − measured +25.9 kJ/mol (+6.2 kcal/mol) |
| Transition state | one imaginary frequency, 367i cm⁻¹; C(N)–C 2.08 Å (forming), C–O 1.26 Å (becoming single), Bürgi–Dunitz angle C(N)···C=O 106.2° |
| Orbitals that change | C lone pair (cyanide) → σ(C–C); π(C=O) → O lone pair (localized orbitals, IBOs, followed frame to frame) |
| What this is | The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, no solvent molecules. |
| Course mapping |  |

## Limitations

1. The solvent is a continuum (SMD water): the hydrogen bonds that water gives to the alkoxide oxygen are missing, so the alkoxide is too high in energy. A test with one to three water molecules added around the ions supports this for the reaction free energy, which falls from +36.4 to +3.2 kJ/mol (measured +10.5 kJ/mol). The barrier does not respond as simply: it rises to about 85 kJ/mol with one or two waters and falls to 73.5 kJ/mol with three (measured 69.6 kJ/mol). The test shows the size of the effect, not a corrected value. The numbers on this card are without added waters.
2. The coupled-cluster ΔG‡ is 9.7 kJ/mol (2.3 kcal/mol) above the measured value. That difference should not be over-read either: besides the continuum solvent, the geometries and thermal corrections come from M06-2X, and the choice of solvation scheme alone moves this barrier by about 1 kJ/mol. M06-2X/ma-def2-TZVP on the same geometries gives 72.1 kJ/mol (17.2 kcal/mol).
3. The path ends at the alkoxide, (CH₃)₂C(CN)O⁻, which lies above cyanide + acetone in free energy: the addition step on its own is uphill. The neutral cyanohydrin forms only after a later proton transfer (from HCN or water), which is not part of this path.
4. The measured equilibrium constant rests on an estimated, not measured, acidity of the cyanohydrin O–H, so the measured ΔrG (+10.5 kJ/mol) is uncertain by a few kJ/mol. The measured rate constant does not depend on it.
5. The measured rate constant refers to free cyanide: only CN⁻, not HCN, adds to the carbonyl, so observed rates depend on pH.
6. In the continuum solvent, cyanide and acetone form a loose complex 8.1 kJ/mol (1.9 kcal/mol) below the separated molecules (DLPNO-CCSD(T)/ma-def2-TZVP, electronic); this shallow well is partly an artifact of the model.

Energies at every level, how the measured value was derived, and technical notes: `card_details.md`.

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
