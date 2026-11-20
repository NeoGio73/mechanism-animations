kJ/mol (kcal/mol), 298.15 K, 1 M (gas) → 1 M (solution).

**Solvation free energies in water** (SMD, solution geometry, minus gas, gas geometry):

| Level | Hydroxide | Bromomethane | Methanol | Bromide | Bromide − hydroxide |
|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP | -409.5 (-97.9) | -7.9 (-1.9) | -22.5 (-5.4) | -224.0 (-53.5) | +185.5 (+44.3) |
| M06-2X/ma-def2-TZVP | -400.3 (-95.7) | -7.8 (-1.9) | -19.9 (-4.7) | -222.9 (-53.3) | +177.4 (+42.4) |
| B3LYP/ma-def2-TZVP | -388.1 (-92.8) | -7.7 (-1.8) | -19.0 (-4.6) | -221.9 (-53.0) | +166.1 (+39.7) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -393.8 (-94.1) | -7.0 (-1.7) | -19.5 (-4.7) | -222.9 (-53.3) | +170.8 (+40.8) |
| Measured | -438.1 (-104.7) | -3.4 (-0.8) | -21.3 (-5.1) | -285.8 (-68.3) | +152.3 (+36.4) |

**Reaction free energy and its error split** (measured: solution -100.8 (-24.1), gas -235.2 (-56.2) = solution minus the measured hydration terms):

| Level | ΔrG, solution | ΔrG, gas | Error, total | (a) gas-phase ΔrG | (b) hydroxide | (b) bromomethane | (b) methanol | (b) bromide |
|---|---|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP | -109.8 (-26.2) | -280.8 (-67.1) | -9.0 (-2.1) | -45.6 (-10.9) | -28.5 (-6.8) | +4.5 (+1.1) | -1.1 (-0.3) | +61.8 (+14.8) |
| M06-2X/ma-def2-TZVP | -93.7 (-22.4) | -259.0 (-61.9) | +7.1 (+1.7) | -23.8 (-5.7) | -37.7 (-9.0) | +4.3 (+1.0) | +1.5 (+0.4) | +62.8 (+15.0) |
| B3LYP/ma-def2-TZVP | -82.0 (-19.6) | -236.8 (-56.6) | +18.8 (+4.5) | -1.6 (-0.4) | -50.0 (-11.9) | +4.2 (+1.0) | +2.3 (+0.5) | +63.8 (+15.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -85.3 (-20.4) | -243.7 (-58.2) | +15.5 (+3.7) | -8.5 (-2.0) | -44.3 (-10.6) | +3.6 (+0.9) | +1.8 (+0.4) | +62.8 (+15.0) |

(b) is each species' contribution to the reaction error: ν × [computed − measured solvation], ν = −1 for reactants and +1 for products; (a) + (b) = total.

**Findings** (DLPNO-CCSD(T)/ma-def2-TZVP, the reference level, unless stated):

- The solvation of bromide carries most of the error: +62.8 (+15.0). SMD solvates bromide by -222.9 (-53.3) against -285.8 (-68.3) measured; the same at every level (+61.8 to +63.8 kJ/mol), because the solvation of a single ion barely depends on the level.
- Hydroxide is also under-solvated, by +44.3 (+10.6); as a reactant it pulls the error the other way (-44.3 (-10.6)). The two ions together give +18.5 (+4.4): SMD's bromide − hydroxide difference is +170.8 (+40.8) against +152.3 (+36.4) measured, so SMD favours hydroxide over bromide too much.
- That ion total does not depend on the absolute proton value behind the single-ion data (−265.9 kcal/mol, uncertain by at least 2 kcal/mol = 8.4 kJ/mol): one anion is consumed and one formed, so the proton value cancels. The split between the two ions does depend on it.
- The neutrals add +5.4 (+1.3) (bromomethane +3.6 (+0.9), methanol +1.8 (+0.4)).
- The gas-phase reaction free energy is off by -8.5 (-2.0). Without the solvation errors, DLPNO-CCSD(T) would give ΔrG in water within 8.5 kJ/mol of experiment.
- The total solvation error is M06-2X-D3(0)/ma-def2-SVP +36.6 (+8.7), M06-2X/ma-def2-TZVP +30.9 (+7.4), B3LYP/ma-def2-TZVP +20.4 (+4.9), DLPNO-CCSD(T)/ma-def2-TZVP +23.9 (+5.7). The bromide term hardly changes with the level; the hydroxide term does (-28.5, -37.7, -50.0, -44.3 kJ/mol in the same order). The gas-phase part: M06-2X-D3(0)/ma-def2-SVP -45.6 (-10.9), M06-2X/ma-def2-TZVP -23.8 (-5.7), B3LYP/ma-def2-TZVP -1.6 (-0.4), DLPNO-CCSD(T)/ma-def2-TZVP -8.5 (-2.0). So M06-2X/ma-def2-TZVP's small total error (+7.1 (+1.7)) is a cancellation of a gas-phase error against the solvation error, and B3LYP/ma-def2-TZVP has the gas-phase value closest to experiment but the largest total error.
- SMD cavity radius ORCA used for Br: 3.06 Å (ORCA 6.1.1, "element-dependent radii"). The electrostatic part of the computed bromide solvation, -221.8 kJ/mol, is close to the Born value for a point charge in a sphere of that radius (-224.1 kJ/mol). The measured value would need a Born sphere of about 2.4 Å. So the bromide error follows from the size of the cavity SMD gives Br; this radius was not checked against the SMD paper (Marenich, Cramer and Truhlar 2009).

**Sources and standard state.** Every solvation free energy is 1 M (gas) → 1 M (solution), 298.15 K. ORCA's G is for 1 atm in both phases, so the 1 atm → 1 M term (7.92 kJ/mol) cancels in G(solution) − G(gas). Neutrals: FreeSolv v0.52 (Mobley and Guthrie, J. Comput.-Aided Mol. Des. 2014, 28, 711), entries mobley_1636752 (methanol, −5.10 kcal/mol) and mobley_8983100 (bromomethane, −0.82 kcal/mol), experimental reference Rizzo et al., J. Chem. Theory Comput. 2006, 2, 128; FreeSolv's ±0.60 kcal/mol for these two is a default ("not presently available"), not a measured uncertainty, and FreeSolv does not state its standard state (its methanol value equals the 1 M → 1 M value of Kelly et al., Table 3). Ions: C. P. Kelly, C. J. Cramer and D. G. Truhlar, J. Phys. Chem. B 2006, 110, 16066, Table 8 ('this work'; p. 58 of the authors' final manuscript, comp.chem.umn.edu/Truhlar/docs/final744.pdf; standard state in Section 2, pp. 5–6): OH⁻ −104.7 and Br⁻ −68.3 kcal/mol with ΔG*solv(H⁺) = −265.9 kcal/mol. Cross-check, same table, column 'ref 23' (Pliego and Riveros 2002): OH⁻ −105.0, Br⁻ −68.6 kcal/mol, the same bromide − hydroxide difference.

**Caveats.** The measured solution ΔrG (−100.8 kJ/mol) built CH₃Br(aq) from a Henry constant (+4.70 kJ/mol, 1 bar → 1 mol/kg); FreeSolv's value used here for bromomethane corresponds to +4.53 kJ/mol for the same step. That moves the derived gas-phase value, and the split between (a) and the bromomethane term, by 0.17 kJ/mol; the total error is unchanged. Gas-phase G uses the working-level frequencies of the gas geometry; solution G those of the SMD geometry. No structure had an imaginary frequency.
