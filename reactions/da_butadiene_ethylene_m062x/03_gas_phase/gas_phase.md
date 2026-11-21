Gas-phase comparison with Rowley and Steiner (J. Picione, 2026-11-20). Gas phase, no solvent model, gas-phase geometries and frequencies at M06-2X-D3(0)/ma-def2-SVP; barrier from s-trans butadiene + ethylene; T = 840 K (midpoint of the measured 760–921 K).

| Level | ΔE‡ | ΔH‡ (840 K) | Ea = ΔH‡ + 2RT (840 K) | Ea − measured | A (L mol⁻¹ s⁻¹) | A / measured A | k(TST, 840 K) (L mol⁻¹ s⁻¹) | k / measured k |
|---|---|---|---|---|---|---|---|---|
| M06-2X_ma-def2-SVP (working) | 80.0 kJ/mol (19.1 kcal/mol) | 86.9 kJ/mol (20.8 kcal/mol) | 100.9 kJ/mol (24.1 kcal/mol) | −14.2 kJ/mol (−3.4 kcal/mol) | 1.36 × 10⁷ | 0.45 | 7.26 | 3.48 |
| M06-2X_ma-def2-TZVP | 85.6 kJ/mol (20.5 kcal/mol) | 92.5 kJ/mol (22.1 kcal/mol) | 106.5 kJ/mol (25.4 kcal/mol) | −8.6 kJ/mol (−2.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 3.27 | 1.57 |
| B3LYP_ma-def2-TZVP | 81.3 kJ/mol (19.4 kcal/mol) | 88.2 kJ/mol (21.1 kcal/mol) | 102.2 kJ/mol (24.4 kcal/mol) | −12.9 kJ/mol (−3.1 kcal/mol) | 1.36 × 10⁷ | 0.45 | 6.01 | 2.88 |
| DLPNO-CCSD(T)_ma-def2-TZVP | 98.4 kJ/mol (23.5 kcal/mol) | 105.3 kJ/mol (25.2 kcal/mol) | 119.3 kJ/mol (28.5 kcal/mol) | +4.2 kJ/mol (+1.0 kcal/mol) | 1.36 × 10⁷ | 0.45 | 0.524 | 0.25 |
| Measured (Rowley and Steiner 1951) | | | 115.1 ± 3.4 kJ/mol (27.5 ± 0.8 kcal/mol) | | 3.00 × 10⁷ | 1 | 2.09 | 1 |

**Units and standard state of A.** A and k(TST) are second-order rate constants in concentration units: L mol⁻¹ s⁻¹ (= 10³ cm³ mol⁻¹ s⁻¹, Rowley and Steiner's unit), also given per molecule (cm³ molecule⁻¹ s⁻¹). ΔG‡ and ΔS‡ use the 1 atm standard state; the factor RT/p° turns that into concentration units, so A and k do not depend on the choice of standard state and compare directly with the measured A. Measured A = 3.00 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 3.00 × 10⁷ L mol⁻¹ s⁻¹ = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹; computed A = 1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹ = 2.27 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹.

A is the same at every level: all levels share the working-level gas-phase frequencies, so ΔS‡ (840 K, 1 atm) = −168.8 J mol⁻¹ K⁻¹ at every level; the levels differ only in Ea.

Method: ΔH‡(840 K) = ΔE‡ + Δ[H(840 K) − E]; Ea = ΔH‡ + 2RT; ΔG‡(840 K, 1 atm) = ΔE‡ + Δ[G(840 K) − E] (quasi-RRHO entropy); k(TST) = (kB T/h)(RT/p°) exp(−ΔG‡/RT), κ = 1; A = k(TST) exp(Ea/RT). Thermal terms from the gas-phase working-level frequencies (rotational symmetry numbers from ORCA's point-group detection: ts σ = 1, butadiene_s-trans σ = 2, ethylene σ = 4). Same method as the earlier B3LYP-geometry folder, except that geometries and frequencies are now gas phase.

**Gas-phase TS:** one imaginary frequency, −514 cm⁻¹; overlap with the toluene TS mode 1.00; forming bonds C1–C6 2.265 Å and C4–C5 2.265 Å (spread 0.000 Å, stretch-rate ratio 1.00): synchronous, accepted by the folder's TS rules. Toluene TS: 2.263 Å.

**Butadiene conformers (gas phase, working level):** the gauche s-cis form lies 11.6 kJ/mol (2.8 kcal/mol) in E and 9.2 kJ/mol (2.2 kcal/mol) in G (298.15 K) above s-trans; the lowest conformer is s-trans.

**Findings (2026-11-20 20:52).**

- **Activation energy.** DLPNO-CCSD(T)/ma-def2-TZVP gives Ea = 119.3 kJ/mol (28.5 kcal/mol), 4.2 kJ/mol (1.0 kcal/mol) above the measured 115.1 ± 3.4 kJ/mol, just outside the fit error. M06-2X/ma-def2-TZVP is 8.6 kJ/mol (2.1 kcal/mol) below, B3LYP/ma-def2-TZVP 12.9 kJ/mol (3.1 kcal/mol) below, and the working level 14.2 kJ/mol (3.4 kcal/mol) below.
- **Pre-exponential factor.** A = 1.36 × 10⁷ L mol⁻¹ s⁻¹ (1.36 × 10¹⁰ cm³ mol⁻¹ s⁻¹), 0.45 times the measured 3.0 × 10⁷ L mol⁻¹ s⁻¹, at every level (ΔS‡ = −168.8 J mol⁻¹ K⁻¹ at 840 K, 1 atm). A and Ea from one Arrhenius fit over 760–921 K are correlated: the ±3.4 kJ/mol in Ea corresponds to about a factor of 1.6 in A at 840 K, so the computed A is low by somewhat more than the fit uncertainty.
- **Rate at 840 K.** DLPNO-CCSD(T) gives 0.25 times the measured k: 0.45 from A times 0.55 from the higher Ea. M06-2X/ma-def2-TZVP gives 1.57 times: 0.45 from A times 3.4 from its low Ea. As in the earlier folder, M06-2X matches the rate better only because its low Ea cancels the low A.
- **Against the withdrawn numbers** from the earlier B3LYP-geometry folder (solution geometries and SMD frequencies, gas-phase single points): Ea difference −8.8 → −8.6 kJ/mol (M06-2X/ma-def2-TZVP), +4.7 → +4.2 kJ/mol (DLPNO-CCSD(T)); A ratio 0.57 → 0.45.
- **Gas phase against toluene** (same levels, s-trans butadiene + ethylene): the DLPNO-CCSD(T) ΔE‡ is 98.4 kJ/mol in the gas phase and 96.7 kJ/mol in SMD toluene; the forming C–C bonds are 2.265 Å and 2.263 Å.
- **Checks.** Every job ended normally. The minima have no imaginary frequency. The TS has one (−514 cm⁻¹), synchronous, with overlap 1.00 with the toluene TS mode. The thermal terms reproduce ORCA's 298.15 K H and G within 0.001 kJ/mol. The largest DLPNO-CCSD(T) T1 is 0.0113.
- **Job time:** 2092 s in total (4 optimizations 820 s, 9 single points 1272 s, 4 cores each), 20:24–20:52, ahead of the path single points.
