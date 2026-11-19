# Results summary: SN2, hydroxide + bromomethane

`HO⁻ + CH₃Br → CH₃OH + Br⁻`, in water (SMD), 298.15 K. Reaction folder `reactions/sn2_hydroxide_bromomethane_m062x`. Written 2026-11-19 08:43 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. "Separated reactants" means hydroxide + bromomethane, each optimized on its own (single conformers). Steps that change the number of molecules use the 1 M standard state (RT ln 24.47 = 7.93 kJ/mol per mole of solute: subtracted for each molecule lost, added for each gained); steps from the reactant complex involve no change in the number of molecules, so no standard-state term applies. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 342 (reactant side 153 frames, TS = frame 153, product side 188 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 342 of 342; M06-2X/ma-def2-TZVP 342 of 342; B3LYP/ma-def2-TZVP 342 of 342; DLPNO-CCSD(T)/ma-def2-TZVP 328 of 342. DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 17.05 bohr·amu½ in the mass-weighted coordinate s (from s = -7.65 at the reactant end to s = 9.40 at the product end; TS at s = 0); 3.15 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (C–O or C–Br between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 15 (s = -6.90, C–O 3.21 Å) to frame 327 (s = 8.70, 1.42 Å). At the TS the C–O distance is 2.112 Å, 63 % of the way from the reactant end (3.308 Å) to the product end (1.422 Å): most of the shortening has happened by the TS, and the remaining 37 % (0.69 Å) comes after it, over the product side (s = 0 to 9.40).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| C–O | forms | 3.308 | 2.112 | 1.422 | 63 % |
| C–Br | breaks | 1.949 | 2.323 | 3.490 | 24 % |

| Point | Frame | s (bohr·amu½) | C–O (Å) | C–Br (Å) | E, M06-2X/ma-def2-SVP, zero = separated reactants |
|---|---|---|---|---|---|
| reactant end | 0 | -7.65 | 3.308 | 1.949 | -6.2 (-1.5) |
| s ≈ -6 | 33 | -6.00 | 3.086 | 1.946 | -4.8 (-1.1) |
| s ≈ -4 | 73 | -4.00 | 2.812 | 1.945 | +2.5 (+0.6) |
| s ≈ -2 | 113 | -2.00 | 2.495 | 2.035 | +23.1 (+5.5) |
| s ≈ -1 | 133 | -1.00 | 2.304 | 2.177 | +46.4 (+11.1) |
| s ≈ -0.5 | 143 | -0.50 | 2.208 | 2.250 | +56.8 (+13.6) |
| TS | 153 | +0.00 | 2.112 | 2.323 | +61.2 (+14.6) |
| s ≈ +0.5 | 163 | +0.50 | 2.017 | 2.395 | +56.1 (+13.4) |
| s ≈ +1 | 173 | +1.00 | 1.921 | 2.467 | +40.7 (+9.7) |
| s ≈ +1.5 | 183 | +1.50 | 1.824 | 2.539 | +15.2 (+3.6) |
| s ≈ +2 | 193 | +2.00 | 1.728 | 2.611 | -17.6 (-4.2) |
| s ≈ +3 | 213 | +3.00 | 1.535 | 2.753 | -83.7 (-20.0) |
| s ≈ +6 | 273 | +6.00 | 1.421 | 3.114 | -126.7 (-30.3) |
| product end | 341 | +9.40 | 1.422 | 3.490 | -134.0 (-32.0) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD water)

- **TS:** 1 imaginary frequency, -619.2 cm⁻¹ (written 619.2i cm⁻¹). Stretch rates of the key bonds in the mode: C–O +1.27, C–Br -0.96; forming and breaking bonds move together. C–O 2.112 Å, C–Br 2.323 Å.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = separated reactants | G (298.15 K, 1 M), zero = separated reactants |
|---|---|---|---|---|---|
| hydroxide (`[OH-]`, job `hydroxide`) | separated reactant | O–H 0.964 Å | none | — | — |
| bromomethane (`CBr`, job `bromomethane`) | separated reactant | C–Br 1.946 Å | none | — | — |
| reactant complex (`reactant_complex_irc`) | reactant end of the path | C–O 3.309 Å, C–Br 1.949 Å | none | -6.2 (-1.5) | +16.7 (+4.0) |
| TS |  | C–O 2.112 Å, C–Br 2.323 Å | 1 | +61.1 (+14.6) | +91.5 (+21.9) |
| product complex | product end of the path | C–O 1.422 Å, C–Br 3.576 Å | none | -134.2 (-32.1) | -98.4 (-23.5) |
| methanol (`CO`, job `methanol`) | separated product | C–O 1.417 Å | none | — | — |
| bromide (`[Br-]`, job `bromide`) | separated product | single atom | none | — | — |
| methanol + bromide | separated products | | | -130.2 (-31.1) | -109.8 (-26.2) |

The products are 2 separate molecules (methanol + bromide); the product complex holds them together at the end of the path. Reaction energies use the separated products.
The reactant end of the path ended 0.001 Å (RMSD) from the reactant-end minimum, the product end 0.031 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = separated reactants (hydroxide + bromomethane), 298.15 K, ΔG at 1 M:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 61.1 (14.6) | 63.5 (15.2) | 91.5 (21.9) |
| M06-2X/ma-def2-TZVP | 63.7 (15.2) | 66.1 (15.8) | 94.1 (22.5) |
| B3LYP/ma-def2-TZVP | 49.0 (11.7) | 51.4 (12.3) | 79.4 (19.0) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 69.9 (16.7) | 72.3 (17.3) | 100.4 (24.0) |

**Barriers, zero = reactant complex, 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 67.4 (16.1) | 65.6 (15.7) | 74.9 (17.9) |
| M06-2X/ma-def2-TZVP | 70.4 (16.8) | 68.7 (16.4) | 78.0 (18.6) |
| B3LYP/ma-def2-TZVP | 53.8 (12.9) | 52.0 (12.4) | 61.3 (14.7) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 75.3 (18.0) | 73.5 (17.6) | 82.8 (19.8) |

**Reaction, zero = separated reactants; products = separated methanol + bromide, 298.15 K, ΔG at 1 M:**

| Level | ΔrE | ΔrH | ΔrG | Reactant complex, ΔE / ΔG | Reverse barrier from the product complex, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | -130.2 (-31.1) | -117.3 (-28.0) | -109.8 (-26.2) | -6.2 (-1.5) / +16.7 (+4.0) | 195.3 (46.7) / 189.9 (45.4) |
| M06-2X/ma-def2-TZVP | -114.1 (-27.3) | -101.2 (-24.2) | -93.7 (-22.4) | -6.8 (-1.6) / +16.1 (+3.9) | 182.0 (43.5) / 176.5 (42.2) |
| B3LYP/ma-def2-TZVP | -102.4 (-24.5) | -89.5 (-21.4) | -82.0 (-19.6) | -4.8 (-1.1) / +18.1 (+4.3) | 153.2 (36.6) / 147.7 (35.3) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -105.8 (-25.3) | -92.9 (-22.2) | -85.3 (-20.4) | -5.4 (-1.3) / +17.5 (+4.2) | 177.7 (42.5) / 172.3 (41.2) |

ΔrH was not asked for; it is given because it costs nothing. The reactant complex lies below separated reactants in E but above them in G at 1 M at every level: the binding energy does not pay for the entropy lost on association.

## 4. Highest point of the path at each level

Path energies at each level, zero = separated reactants (hydroxide + bromomethane) at that level. The TS is frame 153 (s = 0, C–O 2.112 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C–O (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 342 of 342 | 61.2 (14.6) | 153 | +0.00 | 2.112 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 342 of 342 | 63.7 (15.2) | 154 | +0.05 | 2.103 | +0.07 (+0.02) |
| B3LYP/ma-def2-TZVP | 342 of 342 | 49.0 (11.7) | 153 | +0.00 | 2.112 | +0.00 (+0.00) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 328 of 342 | 69.9 (16.7) | 155 | +0.10 | 2.093 | +0.39 (+0.09) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter C–O distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0147 | 0.051602 |
| reactant_complex | 0.0091 | 0.057079 |
| product_complex | 0.0069 | 0.034959 |
| hydroxide | 0.0132 | 0.037156 |
| bromomethane | 0.0070 | 0.055178 |
| methanol | 0.0091 | 0.029507 |
| bromide | 0.0041 | 0.033941 |

- **Path frames (328 computed):** T1 from 0.0069 to 0.0148, median 0.0092; largest at frame 159 (s = +0.30, C–O 2.055 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: passes** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 153 | 0 | 0.03° / 2.93° (20) | 0.97° / 3.22° (132) | yes |
| product side (forward) | 188 | 0 | 0.02° / 3.09° (20) | 0.34° / 1.45° (167) | yes |

- **End checks: ok** (`02_irc/true_irc/end_check_all.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame +0.099 kJ/mol from `reactant_complex` (RMSD 0.077 Å). Opt Freq from the last frame: E +0.036 kJ/mol from it, heavy-atom RMSD 0.055 Å (all-atom 0.077 Å), largest key-distance difference -0.116 (C–O) Å, 0 imaginary frequencies: a different form of reactant_complex with the same bonds (fails: key distances); now the official end structure reactant_complex_irc.
  - Product side (forward): last frame +0.232 kJ/mol from `product_complex` (RMSD 0.031 Å). Opt Freq from the last frame: E +0.011 kJ/mol from it, heavy-atom RMSD 0.002 Å (all-atom 0.004 Å), largest key-distance difference -0.001 (C–Br) Å, 0 imaginary frequencies: reaches the same minimum.
  - Official end structures now: reactant end `reactant_complex_irc`, product end `product_complex` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: 0.06 kJ/mol; product side: 0.23 kJ/mol (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 1 (stopping rule near a flat TS): the 20 flat steps on both sides began at least 10 kJ/mol below the TS, so the amended rule stops them at the same frame; not rerun.
  - Entry 4 (the minimum the path reaches is the official end): on the reactant side the reached form has C···O 3.309 Å instead of 3.425 Å (C–Br within 0.001 Å, heavy-atom RMSD 0.055 Å, same bonds), E +0.036 kJ/mol and G +0.589 kJ/mol from the earlier reference; it is now the official reactant end (`reactant_complex_irc`). This corrected the barriers from the reactant complex reported on 2026-11-14 (17:00 report); the values above are the corrected ones. Product side: same minimum.
  - Entry 5: the side-end row above.
  - Entry 6 (angle against gradient size): 2 of the 70 frames over 15° in all six paths are on this path, in the final approach to a minimum, where the gradient is nearly zero. Near this TS the angles stay small.
  - Entry 8 (experimental barrier): the measured ΔG‡ is 94.9 kJ/mol (Eyring fit of Moelwyn-Hughes 1949, Table 5), not the earlier 95.0 kJ/mol; every computed − measured difference rises by 0.1 kJ/mol. No calculation changes.
  - Entry 9 (lowest-free-energy conformers): hydroxide, bromomethane, methanol and bromide have single conformers; nothing changes.

## 7. Experiment

Experiment: k₂(HO⁻ + CH₃Br, water) measured at six temperatures, 292.69–322.91 K (E. A. Moelwyn-Hughes, Proc. R. Soc. Lond. A 1949, 196, 540–553, doi:10.1098/rspa.1949.0044, Table 5; NaOH in water, rates corrected for the competing reaction with water). An Eyring fit of the six points gives ΔH‡ = 93.7 kJ/mol (22.4 kcal/mol), ΔS‡ = −4.1 J mol⁻¹ K⁻¹ (−1.0 cal mol⁻¹ K⁻¹) and ΔG‡(298.15 K, 1 M) = 94.9 kJ/mol (22.7 kcal/mol). 298 K lies inside the measured range, so no extrapolation is involved. The paper's Table 9 gives k₂(298.1 K) = 1.435 × 10⁻⁴ L mol⁻¹ s⁻¹, A = 1.04 × 10¹³ L mol⁻¹ s⁻¹ and Ea = 96.2 ± 0.6 kJ/mol (23.0 ± 0.15 kcal/mol). The same values are summarized in I. Fells and E. A. Moelwyn-Hughes, J. Chem. Soc. 1959, 398, doi:10.1039/JR9590000398 (Table 11, p. 408) and R. H. Bathgate and E. A. Moelwyn-Hughes, J. Chem. Soc. 1959, 2642, doi:10.1039/JR9590002642 (Table 12, p. 2646). Checked against the original measurements by J. Picione, 2026-11-17. The earlier card value came from 1.41 × 10⁻⁴ L mol⁻¹ s⁻¹ (Mabey and Mill 1978, Table 5.1), calculated from the 1959 Arrhenius parameters at 298.0 K; our earlier records treated it as a 298.15 K value, giving 95.0 kJ/mol. Review: W. Mabey and T. Mill, J. Phys. Chem. Ref. Data 1978, 7, 383–415, doi:10.1063/1.555572. Reaction free energy in water (298.15 K; tables use 1 mol/kg, the calculation 1 mol/L; the 2 → 2 reaction makes the difference vanish): ΔrG° = −100.8 ± 0.4 kJ/mol (−24.1 kcal/mol), from standard formation free energies in the NBS tables (D. D. Wagman et al., J. Phys. Chem. Ref. Data 1982, 11, Suppl. 2): OH⁻(ao) −157.244 kJ/mol (p. 2-38), Br⁻(ao) −103.96 kJ/mol (p. 2-50), CH₃OH(ao) −175.31 kJ/mol (p. 2-84), CH₃Br(g) −25.9 kJ/mol (p. 2-88). CH₃Br(aq) is not tabulated: the gas value plus its hydration, +4.70 kJ/mol from the Henry constant 1.5 × 10⁻³ mol m⁻³ Pa⁻¹ (J. B. Burkholder et al., JPL Publication 19-5, 2019, as compiled by R. Sander, Atmos. Chem. Phys. 2023, 23, 10901–12440, doi:10.5194/acp-23-10901-2023); the measured Henry constants (1.3–1.8 × 10⁻³) give the ± 0.4 kJ/mol, and FreeSolv (D. L. Mobley and J. P. Guthrie, J. Comput.-Aided Mol. Des. 2014, 28, 711–720, doi:10.1007/s10822-014-9747-x) gives +4.53 kJ/mol for the same step. Derivation: docs/experimental_reaction_free_energies.md.

**Computed against measured** (computed − measured; zero = separated reactants (hydroxide + bromomethane), 298.15 K, 1 M):

| Level | ΔG‡ | ΔrG |
|---|---|---|
| Measured | 94.9 (22.7) | -100.8 (-24.1) |
| M06-2X/ma-def2-SVP (working level, geometries) | 91.5 (21.9), -3.4 (-0.8) | -109.8 (-26.2), -9.0 (-2.1) |
| M06-2X/ma-def2-TZVP | 94.1 (22.5), -0.8 (-0.2) | -93.7 (-22.4), +7.1 (+1.7) |
| B3LYP/ma-def2-TZVP | 79.4 (19.0), -15.5 (-3.7) | -82.0 (-19.6), +18.8 (+4.5) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 100.4 (24.0), +5.5 (+1.3) | -85.3 (-20.4), +15.5 (+3.7) |

Limitations recorded for this comparison (card notes, `00_inputs/card_notes.toml`):

- The pictures (geometries, IRC path, orbitals) are at the working level; only the energies change between the curves offered in the viewer. Single points along a path computed at another level can peak one or two frames away from the working-level TS.

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 342 | 3.2 | 34 | 35 | 38 (frame 38) | 2026-11-14 07:22 | 2026-11-14 11:20 |
| B3LYP/ma-def2-TZVP | 342 | 2.7 | 29 | 29 | 31 (frame 173) | 2026-11-14 11:20 | 2026-11-14 14:31 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 328 | 13.6 | 153 | 186 | 198 (frame 158) | 2026-11-14 14:32 | 2026-11-15 07:12 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 390 gradient jobs, 3.8 h job time (folders `backward`); product side (forward): 420 gradient jobs, 4.1 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.
