Sensitivity check only: every other number in this summary stays on standard SMD (Br radius 3.06 Å). SMD18 (Engelage, Mohamed and Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652) changed the SMD Coulomb radius of Br from 3.06 to 2.60 Å. The 5 bromine-containing structures (TS, Reactant complex, Product complex, 2-bromopropane, Bromide) were recomputed at all four levels with only that radius changed: single points on the existing M06-2X/ma-def2-SVP geometries, SMD ethanol, working-level thermal corrections. Species without bromine do not change. kJ/mol (kcal/mol).

**Radius check.** In all 20 new outputs ORCA 6.1.1 printed "Radius for Br used is 4.9133 Bohr (= 2.6000 Ang.)", against 3.0600 Å in the standard outputs. The SMD-CDS term did not change (largest change 0.000 kJ/mol): only the electrostatic cavity uses this radius. Bromide's cavity area went from 420.2 to 303.4 Å². At the working level the standard single points reproduce the stored optimization energies within 0.017 kJ/mol.

**Energy change with the 2.60 Å radius, E(2.60 Å) − E(3.06 Å):**

| Structure | M06-2X/ma-def2-SVP (working level) | M06-2X/ma-def2-TZVP | B3LYP/ma-def2-TZVP | DLPNO-CCSD(T)/ma-def2-TZVP |
|---|---|---|---|---|
| TS | -9.4 (-2.2) | -8.8 (-2.1) | -8.6 (-2.0) | -8.2 (-2.0) |
| Reactant complex | -2.9 (-0.7) | -2.6 (-0.6) | -2.5 (-0.6) | -2.4 (-0.6) |
| Product complex | -29.8 (-7.1) | -28.8 (-6.9) | -27.8 (-6.7) | -29.0 (-6.9) |
| 2-bromopropane | -2.8 (-0.7) | -2.6 (-0.6) | -2.5 (-0.6) | -2.4 (-0.6) |
| Bromide | -32.3 (-7.7) | -31.2 (-7.5) | -30.1 (-7.2) | -31.3 (-7.5) |

**ΔG‡ from separated reactants (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 77.6 (18.5) | 71.0 (17.0) | -6.6 (-1.6) |
| M06-2X/ma-def2-TZVP | 84.4 (20.2) | 78.1 (18.7) | -6.3 (-1.5) |
| B3LYP/ma-def2-TZVP | 76.2 (18.2) | 70.1 (16.8) | -6.1 (-1.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 99.6 (23.8) | 93.8 (22.4) | -5.8 (-1.4) |
| Measured | 105.1 (25.1) | | |

**ΔG‡ from the reactant complex (298.15 K):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | 48.4 (11.6) | 41.9 (10.0) | -6.5 (-1.6) |
| M06-2X/ma-def2-TZVP | 55.1 (13.2) | 48.9 (11.7) | -6.2 (-1.5) |
| B3LYP/ma-def2-TZVP | 43.0 (10.3) | 37.0 (8.8) | -6.0 (-1.4) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 67.7 (16.2) | 61.9 (14.8) | -5.8 (-1.4) |

**Reaction free energy ΔrG (298.15 K, 1 M):**

| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level) | -119.3 (-28.5) | -148.8 (-35.6) | -29.5 (-7.1) |
| M06-2X/ma-def2-TZVP | -115.1 (-27.5) | -143.7 (-34.4) | -28.7 (-6.8) |
| B3LYP/ma-def2-TZVP | -115.1 (-27.5) | -142.7 (-34.1) | -27.6 (-6.6) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -108.3 (-25.9) | -137.2 (-32.8) | -28.9 (-6.9) |

**Findings.**

- **Barrier.** With the SMD18 radius every barrier drops by 5.8–6.6 kJ/mol (1.4–1.6 kcal/mol), a little less than in the SN2 check (7.3–8.0 kJ/mol). The TS gains about 26 % of what free bromide gains (−8.2 against −31.3 kJ/mol at DLPNO-CCSD(T)), while 2-bromopropane and the reactant complex gain only 2.4 kJ/mol. DLPNO-CCSD(T)/ma-def2-TZVP ΔG‡ from separated reactants goes from 99.6 to 93.8 kJ/mol, against 105.1 measured: from 5.5 kJ/mol below to 11.3 kJ/mol below. For this reaction the new radius moves the barrier away from experiment.
- **Reaction free energy.** ΔrG drops by 27.6–29.5 kJ/mol (6.6–7.1 kcal/mol), almost all of it from bromide. There is no measured ΔrG for this reaction to compare with.
- **SN2 and E2 together.** The radius moves the two barriers by similar amounts in the same direction. It brings SN2 closer to experiment (+5.5 → −1.8 kJ/mol) and E2 further away (−5.5 → −11.3 kJ/mol), so it does not improve both at once.
- **Geometries** are the standard-SMD ones; no structure was reoptimized with the new radius.
