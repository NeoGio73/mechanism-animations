# Methyl-shift TS search from the anti form: reading (2026-11-23 14:40)

Jobs 08:34–10:22, all ended normally (`01_stationary/ts_methyl_shift/`; automated check in `rotation_ts.md`, "NOT confirmed").
Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane; zero = `secondary_anti`; 298.15 K; kJ/mol (kcal/mol).

## What was found

**One stationary point with one imaginary frequency, 165.5i cm⁻¹, slightly above anti. Its imaginary mode is not the methyl shift.**
- **Geometry:** the bridging methyl (C4) sits almost midway between the two carbons, C2–CH₃ 1.733 Å and C3–CH₃ 1.769 Å. In anti these are 1.804 and 1.674 Å; in the C2/C3-swapped anti they are 1.674 and 1.804 Å. C1–C2–C3–H is 8.4°.
- **The imaginary mode:** 96 % of its motion is on the three hydrogens of that methyl (H10, H11, H12): it is the internal rotation of the bridging CH₃ group. A 0.05 Å step along it changes C2–CH₃ and C3–CH₃ by only 0.006–0.007 Å.
- **Lowest real frequencies:** 218.5, 223.6, 238.0 cm⁻¹.
- **Energies from `secondary_anti`:**

| Level | ΔE | ΔG (working-level thermal correction) |
|---|---|---|
| M06-2X/ma-def2-SVP (working) | +1.7 (+0.4) | +3.1 (+0.7) |
| M06-2X/ma-def2-TZVP | +0.7 (+0.2) | +2.1 (+0.5) |
| B3LYP/ma-def2-TZVP | +0.8 (+0.2) | +2.2 (+0.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | +0.7 (+0.2) | +2.0 (+0.5) |

- **DLPNO-CCSD(T) T1:** 0.0095.

**Two-sided check (displaced 0.10 Å along the 165.5i mode and optimized):**
- **Minus side:** back to `secondary_anti`, 0.002 Å away, ΔE +0.01 kJ/mol.
- **Plus side:** stops at a second point that also has one imaginary frequency, 123.4i cm⁻¹.
  - It has the same character: 94 % on the same methyl hydrogens.
  - Its bridge is still near-symmetric, C2–CH₃ 1.718 Å and C3–CH₃ 1.781 Å.
  - It lies ΔE +1.55 kJ/mol (+0.37 kcal/mol) above anti.
  - It is not the swapped anti.

## Reading

- **The degenerate methyl shift region is flat.**
  - Every structure found with the methyl between C3 and C2 lies within 1.7 kJ/mol of anti in working-level electronic energy. That includes the scan points and both stationary points; the first scan's step past the midpoint is +1.6 kJ/mol.
  - No minimum was found other than anti, and none lies below anti in working-level free energy. The TS above lies +3.1 kJ/mol above anti in G.
  - No TS was found whose imaginary mode moves the methyl between the carbons. The only saddles found turn the methyl's hydrogens.
  - The methyl therefore exchanges C2 and C3 with at most a few kJ/mol, but this search did not locate a TS for that motion itself.
- **Per Section 2.4, anti stays the zero for the anti channel.** The anti channel's DLPNO-CCSD(T) ΔG‡ of 18.6 kJ/mol (4.4 kcal/mol) is unchanged, and so are its other levels.
- **Possible next step (not run):** a relaxed scan of the methyl-shift coordinate itself, with the methyl hydrogens free. The coordinate is C3–CH₃ minus C2–CH₃, from +0.13 Å (anti) to −0.13 Å (swapped anti). This would give the profile of the shift directly. About 20–30 min of job time.
