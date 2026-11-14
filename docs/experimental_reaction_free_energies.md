# Experimental reaction free energies from tables (J. Picione, 2026-11-14; low priority)

Standard reaction free energies at 298.15 K, from tabulated standard formation free energies, compared
with the computed values. Sources and status are listed in `docs/references.md`. Further reactions are
added where tables exist.

## SN2: OH⁻(aq) + CH₃Br(aq) → CH₃OH(aq) + Br⁻(aq), water

**Tabulated values.** NBS tables (Wagman et al. 1982), 298.15 K, 1 bar; aqueous standard state is the ideal
solution at 1 mol/kg ("ao"). Each value was read in the scanned source.

| Species | State | ΔfG° | Source |
|---|---|---|---|
| OH⁻ | ao | −157.244 kJ/mol | Table 2:H, p. 2-38 |
| Br⁻ | ao | −103.96 kJ/mol | Table 11:Br, p. 2-50 |
| CH₃OH | ao | −175.31 kJ/mol | Table 23:C, p. 2-84 |
| CH₃Br | g | −25.9 kJ/mol | Table 23:C, p. 2-88 (no aqueous entry) |

**CH₃Br(aq), not tabulated.** It is the gas value plus the free energy of solution from 1 bar gas to the
1 mol/kg ideal solution:

ΔG_solv = −RT ln(Hscp · p° / ρ_water), with p° = 10⁵ Pa and ρ_water = 997.05 kg/m³.

| Henry constant Hscp (mol m⁻³ Pa⁻¹, 298.15 K) | Source (as listed by Sander 2023) | ΔG_solv |
|---|---|---|
| 1.5 × 10⁻³ (used) | Burkholder et al., JPL Publication 19-5 (2019), first-ranked entry | +4.70 kJ/mol |
| 1.6 × 10⁻³ | Glew and Moelwyn-Hughes 1953 (refit); Wilhelm, Battino and Wilcock 1977 | +4.54 kJ/mol |
| 1.3–1.8 × 10⁻³ | measured range: Hiatt 2013 to De Bruyn and Saltzman 1997 | +5.05 to +4.24 kJ/mol |

Cross-check: FreeSolv gives ΔG_hyd = −0.82 kcal/mol (1 M → 1 M), which is +4.53 kJ/mol for 1 bar → 1 M.

With the first-ranked constant, ΔfG°(CH₃Br, aq) = −25.9 + 4.70 = −21.20 kJ/mol.

**Reaction free energy:**

ΔrG° = (−175.31) + (−103.96) − (−157.244) − (−21.20) = **−100.8 kJ/mol (−24.1 kcal/mol)**

The spread of the measured Henry constants changes this by ±0.4 kJ/mol (−101.2 to −100.4 kJ/mol). Two
standard-state points:
- The tables use 1 mol/kg, the calculations 1 mol/L. Per species the two differ by RT ln(0.997) = −0.007 kJ/mol.
- The reaction is 2 → 2, so the standard-state term cancels in both.

**Against the computed values.** Single points on the M06-2X/ma-def2-SVP geometries; SMD water; thermal
corrections from the M06-2X/ma-def2-SVP frequencies; bromide from the standard formulas; 298.15 K, 1 M.

| Level | ΔrG° computed | Computed − table |
|---|---|---|
| M06-2X/ma-def2-TZVP | −93.7 kJ/mol (−22.4 kcal/mol) | +7.1 kJ/mol |
| DLPNO-CCSD(T)/ma-def2-TZVP | −85.3 kJ/mol (−20.4 kcal/mol) | +15.5 kJ/mol |
| B3LYP/ma-def2-TZVP | −82.0 kJ/mol (−19.6 kcal/mol) | +18.8 kJ/mol |
| M06-2X/ma-def2-SVP (working level) | −109.8 kJ/mol (−26.2 kcal/mol) | −9.0 kJ/mol |

All three single-point levels make the reaction less exergonic than the tables. Reading (not tested): the
reaction exchanges a small, strongly hydrated anion (OH⁻) for a large, weakly hydrated one (Br⁻), so the
error of the implicit solvent model for the two anions enters the reaction free energy directly, at every
level.