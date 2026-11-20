Sensitivity check only: every other number in this summary stays on standard SMD (Br radius 3.06 Å). SMD18 (Engelage, Mohamed and Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652) changed the SMD Coulomb radius of Br from 3.06 to 2.60 Å. The 5 bromine-containing structures (TS, Reactant complex, Product complex, Bromomethane, Bromide) were recomputed at all four levels with only that radius changed: single points on the existing M06-2X/ma-def2-SVP geometries, SMD water, working-level thermal corrections. Species without bromine do not change. kJ/mol (kcal/mol).

**Radius check.** In all 20 new outputs ORCA 6.1.1 printed "Radius for Br used is 4.9133 Bohr (= 2.6000 Ang.)", against 3.0600 Å in the standard outputs. The SMD-CDS term did not change (largest change 0.000 kJ/mol): only the electrostatic cavity uses this radius. Bromide's cavity area went from 420.2 to 303.4 Å². At the working level the standard single points reproduce the stored optimization energies within 0.062 kJ/mol.

**Energy change with the 2.60 Å radius, E(2.60 Å) − E(3.06 Å):**

| Structure | M06-2X/ma-def2-SVP (working level) | M06-2X/ma-def2-TZVP | B3LYP/ma-def2-TZVP | DLPNO-CCSD(T)/ma-def2-TZVP |
|---|---|---|---|---|
| TS | -10.7 (-2.6) | -10.0 (-2.4) | -9.9 (-2.4) | -9.6 (-2.3) |
| Reactant complex | -2.7 (-0.6) | -2.4 (-0.6) | -2.3 (-0.6) | -2.3 (-0.5) |
| Product complex | -31.2 (-7.5) | -30.4 (-7.3) | -29.4 (-7.0) | -30.5 (-7.3) |
| Bromomethane | -2.8 (-0.7) | -2.5 (-0.6) | -2.4 (-0.6) | -2.3 (-0.6) |
| Bromide | -33.2 (-7.9) | -32.1 (-7.7) | -31.0 (-7.4) | -32.2 (-7.7) |

**ΔG‡ from separated reactants (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 91.5 (21.9) | 83.6 (20.0) | -7.9 (-1.9) |
| M06-2X/ma-def2-TZVP | 94.1 (22.5) | 86.5 (20.7) | -7.6 (-1.8) |
| B3LYP/ma-def2-TZVP | 79.4 (19.0) | 71.9 (17.2) | -7.5 (-1.8) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 100.4 (24.0) | 93.1 (22.2) | -7.3 (-1.7) |
| Measured | 94.9 (22.7) | | |

**ΔG‡ from the reactant complex (298.15 K):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 74.9 (17.9) | 66.9 (16.0) | -8.0 (-1.9) |
| M06-2X/ma-def2-TZVP | 78.0 (18.6) | 70.3 (16.8) | -7.6 (-1.8) |
| B3LYP/ma-def2-TZVP | 61.3 (14.7) | 53.7 (12.8) | -7.6 (-1.8) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 82.8 (19.8) | 75.5 (18.0) | -7.3 (-1.8) |

**Reaction free energy ΔrG (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | -109.8 (-26.2) | -140.3 (-33.5) | -30.5 (-7.3) |
| M06-2X/ma-def2-TZVP | -93.7 (-22.4) | -123.4 (-29.5) | -29.7 (-7.1) |
| B3LYP/ma-def2-TZVP | -82.0 (-19.6) | -110.6 (-26.4) | -28.6 (-6.8) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -85.3 (-20.4) | -115.2 (-27.5) | -29.9 (-7.1) |
| Measured | -100.8 (-24.1) | | |

**Findings.**

- **Barrier.** With the SMD18 radius every barrier drops by 7.3–8.0 kJ/mol (1.7–1.9 kcal/mol). The TS gains about 30 % of what free bromide gains (−9.6 against −32.2 kJ/mol at DLPNO-CCSD(T)), while bromomethane and the reactant complex gain only 2.3 kJ/mol. DLPNO-CCSD(T)/ma-def2-TZVP ΔG‡ from separated reactants goes from 100.4 to 93.1 kJ/mol, against 94.9 measured: from 5.5 kJ/mol above to 1.8 kJ/mol below.
- **Reaction free energy.** ΔrG drops by 28.6–30.5 kJ/mol (6.8–7.3 kcal/mol), almost all of it from bromide. At DLPNO-CCSD(T) it goes from −85.3 to −115.2 kJ/mol, against −100.8 measured: the error changes from +15.5 to −14.4 kJ/mol. It is about as large as before, with the opposite sign.
- **Why the reaction energy does not improve (with section 9).** Bromide's computed solvation becomes −255.1 kJ/mol at DLPNO-CCSD(T) (was −222.9), still 30.7 kJ/mol short of the measured −285.8. Before, that shortfall was offset by hydroxide's over-solvation (−44.3 kJ/mol in the error split); with the smaller radius the offset is too large. Bromomethane's computed solvation also moves away from experiment (−9.4 against −3.4 measured, was −7.0).
- **With E2.** The same radius lowers the E2 barrier by 5.8 kJ/mol at DLPNO-CCSD(T) and moves it away from experiment (99.6 → 93.8 kJ/mol against 105.1 measured), so it does not improve both reactions at once.
- **Geometries** are the standard-SMD ones; no structure was reoptimized with the new radius.
