# Diels–Alder, butadiene + ethylene

| | |
|---|---|
| Reaction | CH₂=CH–CH=CH₂ + CH₂=CH₂ → cyclohexene, in toluene |
| Level of theory | Geometries, reaction path and orbitals: M06-2X-D3(0)/ma-def2-SVP, with the SMD model of toluene. Energies (the default curve in the viewer and video): DLPNO-CCSD(T)/ma-def2-TZVP single points on those geometries, SMD toluene; thermal corrections from M06-2X-D3(0)/ma-def2-SVP frequencies. Comparison with the gas-phase measurement: separate gas-phase calculations at the same levels. |
| Activation energy, Ea (gas phase, 760–921 K; butadiene as its equilibrium conformer mixture, + ethylene) | **114.9 kJ/mol (27.5 kcal/mol)** at DLPNO-CCSD(T)/ma-def2-TZVP; measured 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol), from rates measured in the gas phase at 760–921 K; computed − measured −0.2 kJ/mol (−0.05 kcal/mol) |
| Free energy of activation in toluene, ΔG‡ (from separated s-trans butadiene + ethylene, 298 K, 1 M) | 144.4 kJ/mol (34.5 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; computed only (the measurement is in the gas phase) |
| Reaction free energy in toluene, ΔrG (from separated s-trans butadiene + ethylene, 298 K, 1 M) | −123.3 kJ/mol (−29.5 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP; no measured value |
| Transition state | one imaginary frequency, 512i cm⁻¹; C1–C6 2.26 Å (forming), C4–C5 2.26 Å (forming), C2–C3 1.41 Å (becoming double), C1–C2 1.38 Å (becoming single), C3–C4 1.38 Å (becoming single), C5–C6 1.38 Å (becoming single), C1=C2–C3=C4 dihedral 0.0° |
| Orbitals that change | π(C3=C4) → σ(C4–C5); π(C5=C6) → σ(C1–C6); π(C1=C2) → π(C2=C3). Localized orbitals (IBOs), followed frame to frame. The arrows could equally be drawn the other way round: the clockwise and counter-clockwise drawings are mirror images of the same mechanism, and the computed orbitals follow the clockwise one. At six frames just after the transition state the computed orbitals came out as the mirror image, so the mirror image is shown there to keep the drawing continuous. The viewer has a mirror switch that shows the other drawing at every frame. |
| What this is | The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, no solvent molecules. |
| Course mapping |  |

## Limitations

1. The measurement is in the gas phase at 760–921 K, so the comparison uses separate gas-phase calculations. The path, orbitals and toluene energies shown are in a continuum solvent (SMD toluene), with no toluene molecules; for these two non-polar molecules the solvent changes the barrier little.
2. The coupled-cluster Ea is 0.2 kJ/mol (0.05 kcal/mol) below the measured value. That agreement should not be over-read: the geometries and thermal corrections come from M06-2X, the measured value is uncertain by ± 3.4 kJ/mol, and counting the butadiene conformer mixture alone moves Ea by 3.7 kJ/mol (118.6 kJ/mol from s-trans butadiene only).
3. The rate is not reproduced as well as Ea: the computed pre-exponential factor is 0.17 times the measured one, an error almost entirely in the entropy of activation.
4. The path runs from s-cis butadiene to the boat form of cyclohexene, the forms the reaction passes through. The energies are measured from the more stable forms, s-trans butadiene and half-chair cyclohexene (23.0 kJ/mol (5.5 kcal/mol) lower than the boat), so they include the cost of turning s-trans into s-cis.
5. The ethylene π orbital shown before s = −2.0 bohr·amu½ is a display choice. The two orbitals of the C5=C6 double bond can be mixed freely, and followed frame to frame it is the σ orbital, not the π orbital, that ends as the new C1–C6 bond.

Energies at every level, how the measured value was derived, and technical notes: `card_details.md`.

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
