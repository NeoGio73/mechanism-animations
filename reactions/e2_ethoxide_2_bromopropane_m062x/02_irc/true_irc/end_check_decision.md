# End check, E2 production path: decision and new official end structures (J. Picione, 2026-11-16)

The automatic end check (`end_check.md`, 03:39) returned "decision": on each side, Opt Freq from the last path
frame reaches a minimum that differs slightly from the hand-built reference complex.

**Decision (J. Picione, 2026-11-16, 05:00):** the path is accepted with a note. The minima the path actually reaches
become the official end structures, because the path defines which minima the TS connects. The hand-built reference
complexes are kept on record with this comparison.

| | reactant side (backward) | product side (forward) |
|---|---|---|
| Official end structure (new) | `reactant_complex_irc` = `end_check/backward/opt` | `product_complex_irc` = `end_check/forward/opt` |
| Hand-built reference (kept on record) | `01_stationary/reactant_complex` (optimized from the B3LYP-level IRC end) | `01_stationary/product_complex` (same) |
| Key distance that differs | ethoxide O···Hβ: 2.623 Å (new) vs 2.688 Å (reference) | Cβ···H(ethanol O–H): 2.544 Å vs 2.443 Å |
| Other key distances | Cβ–H, Cα–Br, Cα–Cβ within 0.001 Å | O–H, Cα–Br, Cα–Cβ within 0.001 Å |
| Heavy-atom RMSD, new vs reference (aligned) | 0.226 Å (all atoms 0.277 Å) | 0.150 Å (all atoms 0.173 Å) |
| Largest atom-pair change | H8···O14 4.102 vs 4.707 Å (the ethoxide turns against the substrate) | H7···H18 2.922 vs 3.294 Å (the ethanol turns against propene) |
| E, new − reference | −0.267 kJ/mol (−0.064 kcal/mol) | −0.279 kJ/mol (−0.067 kcal/mol) |
| G (298.15 K), new − reference | +2.659 kJ/mol (+0.636 kcal/mol) | −0.306 kJ/mol (−0.073 kcal/mol) |
| Lowest frequency (new) | 20.4 cm⁻¹, no imaginary | 17.2 cm⁻¹, no imaginary |
| Last path frame: RMSD to new / to reference | 0.203 / 0.315 Å; +1.17 / +0.90 kJ/mol above it | 0.402 / 0.526 Å; +5.01 / +4.73 kJ/mol above it |

Effect on numbers measured from the complexes (none of these had been reported or used in a card or table):

| Quantity (M06-2X-D3(0)/ma-def2-SVP, SMD ethanol) | with the reference complexes | with the official (new) complexes | change |
|---|---|---|---|
| Barrier from the reactant complex, electronic | 59.40 kJ/mol (14.20 kcal/mol) | 59.67 kJ/mol (14.26 kcal/mol) | +0.27 kJ/mol |
| Barrier from the reactant complex, ΔG‡ (298.15 K) | 51.06 kJ/mol (12.20 kcal/mol) | 48.40 kJ/mol (11.57 kcal/mol) | −2.66 kJ/mol |
| Reactant complex → product complex, electronic | −100.23 kJ/mol (−23.96 kcal/mol) | −100.24 kJ/mol (−23.96 kcal/mol) | −0.01 kJ/mol |
| Reactant complex → product complex, ΔG (298.15 K) | −105.64 kJ/mol (−25.25 kcal/mol) | −108.61 kJ/mol (−25.96 kcal/mol) | −2.97 kJ/mol |

The free-energy changes are larger than the electronic ones because the two forms of each complex differ in
their softest vibrations (both below 25 cm⁻¹), which carry most of the thermal correction in a loose complex.
Barriers from separated reactants (the comparison with experiment, `docs/sensitivity_M06-2X.md`) do not use the
complexes and do not change.
