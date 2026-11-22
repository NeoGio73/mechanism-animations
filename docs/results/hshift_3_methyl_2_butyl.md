# Results summary: 1,2-Hydride shift, 3-methyl-2-butyl cation → 2-methyl-2-butyl cation

`CH₃CH⁺CH(CH₃)₂ → CH₃CH₂C⁺(CH₃)₂`, in dichloromethane (SMD), 298.15 K. Reaction folder `reactions/hshift_3_methyl_2_butyl_m062x`. Written 2026-11-22 05:29 by `pipeline/results_summary.py` from the finished jobs; no new calculation.

**Units and zeros.** Energies in kJ/mol with kcal/mol in brackets. Every number names its zero. The reaction has one molecule throughout, so there are no separated reactants, and every step keeps the number of molecules, so no standard-state term applies. The zero is the lowest-free-energy conformer of the secondary cation (`secondary_anti`, working-level G; rule of J. Picione, 2026-11-19, the same as for separated species), not the form the path reaches. The path starts from the secondary cation as the path reaches it (`secondary_eclip_irc`, frame 0), which lies above the zero by +2.77 (+0.66) kJ/mol in G (working level) and, in E, by M06-2X/ma-def2-SVP +3.06 (+0.73) kJ/mol; M06-2X/ma-def2-TZVP +2.41 (+0.58) kJ/mol; B3LYP/ma-def2-TZVP +2.32 (+0.55) kJ/mol; DLPNO-CCSD(T)/ma-def2-TZVP +3.16 (+0.76) kJ/mol. Barriers from the secondary cation as the path reaches it (`secondary_eclip_irc`, frame 0) are given as well. ΔE is electronic. ΔH and ΔG add the M06-2X/ma-def2-SVP thermal correction (quasi-RRHO, Grimme 2012, 100 cm⁻¹) to the electronic energy of the level named. ΔH needs no standard-state term.

## 1. Path

- **Frames:** 407 (reactant side 189 frames, TS = frame 189, product side 217 frames). Single points: M06-2X/ma-def2-SVP (working level, geometries) 407 of 407; M06-2X/ma-def2-TZVP 407 of 407; B3LYP/ma-def2-TZVP 407 of 407; DLPNO-CCSD(T)/ma-def2-TZVP 309 of 407. DLPNO-CCSD(T) uses your sampling (J. Picione, 2026-11-13): every frame in the bond-changing region and within |s| ≤ 1 bohr·amu½ of the TS, every second frame elsewhere, both ends.
- **Path length:** 20.37 bohr·amu½ in the mass-weighted coordinate s (from s = -9.57 at the reactant end to s = 10.81 at the product end; TS at s = 0); 7.64 Å as the summed all-atom Cartesian displacement between frames (not mass-weighted).
- **Where the bonds change:** the bond-changing region (C2–H or C3–H between 5 % and 95 % of its total change, the definition of `pipeline/path_quality.py`) runs from frame 134 (s = -2.75, C2–H 2.08 Å) to frame 341 (s = 7.57, 1.12 Å). At the TS the C2–H distance is 2.008 Å, 12 % of the way from the reactant end (2.127 Å) to the product end (1.117 Å): less than half of the shortening has happened by the TS, and the remaining 88 % (0.89 Å) comes after it, over the product side (s = 0 to 10.81).

| Bond | Role | Reactant end (Å) | TS (Å) | Product end (Å) | Share of the change done at the TS |
|---|---|---|---|---|---|
| C2–H | forms | 2.127 | 2.008 | 1.117 | 12 % |
| C3–H | breaks | 1.094 | 1.113 | 1.992 | 2 % |

| Point | Frame | s (bohr·amu½) | C2–H (Å) | C3–H (Å) | E, M06-2X/ma-def2-SVP, zero = secondary cation (`secondary_anti`) |
|---|---|---|---|---|---|
| reactant end | 0 | -9.57 | 2.127 | 1.094 | +3.1 (+0.7) |
| s ≈ -6 | 69 | -6.01 | 2.107 | 1.095 | +13.1 (+3.1) |
| s ≈ -4 | 109 | -4.01 | 2.092 | 1.098 | +15.5 (+3.7) |
| s ≈ -2 | 149 | -2.00 | 2.063 | 1.104 | +17.6 (+4.2) |
| s ≈ -1 | 169 | -1.00 | 2.040 | 1.108 | +18.6 (+4.5) |
| s ≈ -0.5 | 179 | -0.49 | 2.024 | 1.111 | +18.9 (+4.5) |
| TS | 189 | +0.00 | 2.008 | 1.113 | +19.0 (+4.5) |
| s ≈ +0.5 | 199 | +0.48 | 1.989 | 1.116 | +18.9 (+4.5) |
| s ≈ +1 | 209 | +0.98 | 1.962 | 1.120 | +18.7 (+4.5) |
| s ≈ +1.5 | 219 | +1.48 | 1.929 | 1.125 | +18.3 (+4.4) |
| s ≈ +2 | 229 | +1.98 | 1.890 | 1.131 | +17.7 (+4.2) |
| s ≈ +3 | 250 | +3.02 | 1.748 | 1.155 | +15.5 (+3.7) |
| s ≈ +6 | 310 | +6.02 | 1.139 | 1.789 | -24.5 (-5.9) |
| s ≈ +10 | 390 | +10.02 | 1.117 | 1.988 | -39.0 (-9.3) |
| product end | 406 | +10.81 | 1.117 | 1.992 | -39.4 (-9.4) |

## 2. Stationary points (M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane)

- **TS:** 1 imaginary frequency, -91.6 cm⁻¹ (written 91.6i cm⁻¹). Stretch rates of the key bonds in the mode: C2–H -0.12, C3–H +0.02; forming and breaking bonds move together. C2–H 2.008 Å, C3–H 1.113 Å.
- **References** (structures the barriers and reaction energies are measured from):

| Structure | Role | Key distances | Imaginary frequencies | E, zero = secondary cation (`secondary_anti`) | G (298.15 K, 1 M), zero = secondary cation (`secondary_anti`) |
|---|---|---|---|---|---|
| secondary cation, `secondary_anti` (job `01_stationary/secondary_anti/secondary_anti`) | lowest-free-energy conformer; the zero | C2–H 2.124 Å, C3–H 1.094 Å | none | 0 | 0 |
| secondary cation (`secondary_eclip_irc`) | reactant end of the path | C2–H 2.127 Å, C3–H 1.094 Å | none | +3.1 (+0.7) | +2.8 (+0.7) |
| TS |  | C2–H 2.008 Å, C3–H 1.113 Å | 1 | +19.0 (+4.5) | +15.8 (+3.8) |
| tertiary cation (`tertiary_irc`) | product end of the path | C2–H 1.117 Å, C3–H 1.992 Å | none | -39.4 (-9.4) | -47.7 (-11.4) |

The reaction has one molecule throughout. The zero is the secondary cation conformer lowest in free energy (`secondary_anti`); the path starts from another conformer, the secondary cation as the path reaches it (`secondary_eclip_irc`, frame 0), +3.06 (+0.73) kJ/mol above it in E and +2.77 (+0.66) kJ/mol in G. The product is the tertiary cation at the product end.
The reactant end of the path ended 0.000 Å (RMSD) from the reactant-end minimum, the product end 0.002 Å from the product-end minimum (section 6).

## 3. Barriers and reaction energies

All at M06-2X/ma-def2-SVP geometries. Electronic energies at the level named; ΔH and ΔG with the working-level thermal correction. kJ/mol (kcal/mol).

**Barriers, zero = the secondary cation, lowest-free-energy conformer (`secondary_anti`), 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 19.0 (4.5) | 13.1 (3.1) | 15.8 (3.8) |
| M06-2X/ma-def2-TZVP | 20.5 (4.9) | 14.6 (3.5) | 17.2 (4.1) |
| B3LYP/ma-def2-TZVP | 3.7 (0.9) | -2.2 (-0.5) | 0.4 (0.1) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 17.3 (4.1) | 11.4 (2.7) | 14.0 (3.3) |

**Barriers, zero = the secondary cation as the path reaches it (`secondary_eclip_irc`, frame 0), 298.15 K:**

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 16.0 (3.8) | 9.2 (2.2) | 13.0 (3.1) |
| M06-2X/ma-def2-TZVP | 18.1 (4.3) | 11.4 (2.7) | 15.1 (3.6) |
| B3LYP/ma-def2-TZVP | 1.4 (0.3) | -5.3 (-1.3) | -1.6 (-0.4) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 14.1 (3.4) | 7.4 (1.8) | 11.1 (2.7) |

**Reaction, zero = the secondary cation, lowest-free-energy conformer (`secondary_anti`); product = the tertiary cation (product-end minimum, `tertiary_irc`), 298.15 K:**

| Level | ΔrE | ΔrH | ΔrG | Secondary cation as the path reaches it, ΔE / ΔG | Reverse barrier from the tertiary cation, ΔE‡ / ΔG‡ |
|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | -39.4 (-9.4) | -44.1 (-10.6) | -47.7 (-11.4) | +3.06 (+0.73) / +2.77 (+0.66) | 58.4 (14.0) / 63.5 (15.2) |
| M06-2X/ma-def2-TZVP | -36.6 (-8.7) | -41.3 (-9.9) | -44.9 (-10.7) | +2.41 (+0.58) / +2.12 (+0.51) | 57.0 (13.6) / 62.1 (14.9) |
| B3LYP/ma-def2-TZVP | -55.0 (-13.1) | -59.8 (-14.3) | -63.4 (-15.1) | +2.32 (+0.55) / +2.03 (+0.49) | 58.7 (14.0) / 63.8 (15.2) |
| DLPNO-CCSD(T)/ma-def2-TZVP | -37.1 (-8.9) | -41.9 (-10.0) | -45.4 (-10.9) | +3.16 (+0.76) / +2.88 (+0.69) | 54.4 (13.0) / 59.5 (14.2) |

ΔrH was not asked for; it is given because it costs nothing.

## 4. Highest point of the path at each level

Path energies at each level, zero = the secondary cation, lowest-free-energy conformer (`secondary_anti`) at that level. The TS is frame 189 (s = 0, C2–H 2.008 Å).

| Level | Frames computed | TS frame energy | Highest frame | Its s (bohr·amu½) | Its C2–H (Å) | Highest − TS frame |
|---|---|---|---|---|---|---|
| M06-2X/ma-def2-SVP (working level, geometries) | 407 of 407 | 19.0 (4.5) | 189 | +0.00 | 2.008 | +0.00 (+0.00) |
| M06-2X/ma-def2-TZVP | 407 of 407 | 20.5 (4.9) | 210 | +1.03 | 1.959 | +0.29 (+0.07) |
| B3LYP/ma-def2-TZVP | 407 of 407 | 3.7 (0.9) | 200 | +0.53 | 1.986 | +0.10 (+0.02) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 309 of 407 | 17.3 (4.1) | 215 | +1.28 | 1.943 | +0.46 (+0.11) |

A highest frame after the TS frame means the level puts the top of the barrier later along the M06-2X/ma-def2-SVP path, at a shorter C2–H distance. The difference is the energy that a single point at the TS geometry misses.

## 5. DLPNO-CCSD(T) diagnostics

| Structure | T1 | Largest amplitude |
|---|---|---|
| ts | 0.0134 | 0.036676 |
| reactant_complex | 0.0098 | 0.043372 |
| product_complex | 0.0132 | 0.031442 |

- **Path frames (309 computed):** T1 from 0.0095 to 0.0135, median 0.0130; largest at frame 221 (s = +1.58, C2–H 1.922 Å).
- **Reading:** every T1 is below 0.02, the usual limit for closed-shell molecules (Lee and Taylor 1989), so a single-reference treatment is adequate for every structure and frame computed.
- **D1:** not available. ORCA prints T1 and the largest amplitudes for DLPNO-CCSD(T), not D1.

## 6. Validation

- **Angle check: fails (flagged)** (rule: in each zone, median < 5°, 90th percentile < 15°; mass-weighted angle between the path direction and the downhill gradient; `02_irc/true_irc/angles.md`):

| Side | Frames | Left out at the end minimum | |s| ≤ 1: median / 90th percentile (frames) | Beyond | Passes |
|---|---|---|---|---|---|
| reactant side (backward) | 189 | 0 | 7.63° / 25.37° (19) | 2.37° / 6.17° (169) | no |
| product side (forward) | 217 | 0 | 10.55° / 25.37° (20) | 0.98° / 7.50° (196) | no |

- **End checks: ok** (`02_irc/true_irc/end_check.md`; from 2026-11-16 each side is optimized from its last frame, and the minimum it reaches is the official end structure, SI entry 4):
  - Reactant side (backward): last frame -0.103 kJ/mol from `secondary_eclip` (RMSD 0.225 Å). Opt Freq from the last frame: E -0.103 kJ/mol from it, heavy-atom RMSD 0.055 Å (all-atom 0.225 Å), largest key-distance difference +0.013 (C2–H) Å, 0 imaginary frequencies: a different form of secondary_eclip with the same bonds (fails: energy); now the official end structure secondary_eclip_irc.
  - Product side (forward): last frame +0.185 kJ/mol from `tertiary` (RMSD 0.429 Å). Opt Freq from the last frame: E +0.181 kJ/mol from it, heavy-atom RMSD 0.125 Å (all-atom 0.428 Å), largest key-distance difference -0.091 (C3–H) Å, 0 imaginary frequencies: a different form of tertiary with the same bonds (fails: energy, key distances); now the official end structure tertiary_irc.
  - Official end structures now: reactant end `secondary_eclip_irc`, product end `tertiary_irc` (`02_irc/endpoints.json`).
- **Side-end table (SI entry 5):** reactant side: secondary cation: 0.0004 kJ/mol; product side: tertiary cation: 0.004 kJ/mol (last frame above the official end minimum, working level). The 0.5 kJ/mol target was dropped on 2026-11-16 (J. Picione): sides are not extended to it; the paper lists how far above its minimum each side ended.
- **SI corrections that touch this reaction** (`docs/si_corrections.md`):
  - Entry 1 (early stop at a flat TS): the first run stopped both sides 20 steps from this flat C2–C3 rotation TS (92i cm⁻¹); it is kept as not results in `02_irc/attempts/first_run_early_stop_flat_ts/`. Both sides were rerun from the TS under the amended rule.
  - Entry 5 (stopping rule near the end minimum): the rerun backward side stopped 12.15 kJ/mol (2.90 kcal/mol) above the secondary cation on a slow twist; it was continued from its last frame (`02_irc/true_irc/backward_part1` + `backward`, 189 frames) and ended 0.0004 kJ/mol above the form of the secondary cation it reaches. Both ends optimize to slightly different forms of the hand-built cations with the same bonds, now the official end structures (`secondary_eclip_irc`, `tertiary_irc`).
  - Entry 5 (angle check): within |s| ≤ 1 bohr·amu½ it fails on both sides (backward 7.63° / 25.37°, forward 10.55° / 25.37°; rule < 5° / < 15°); beyond it, it passes. Near the TS the gradient is nearly zero (entry 6). Kept as computed, both parts stated, no rerun (J. Picione, 2026-11-17).
  - Entry 9 (lowest-free-energy conformer as the zero, J. Picione, 2026-11-19): the zero is the conformer of the secondary cation lowest in working-level G, not automatically frame 0. Frame 0 is the eclip form (C1–C2–C3–H −161.9°, C2–C3–CH₃ 63.0°), not the lowest: the anti minimum (the partly methyl-bridged gas-phase minimum of `docs/protocol.md`, C2–C3–CH₃ 71.3° there), reoptimized at the working level (`01_stationary/secondary_anti`, C1–C2–C3–H −1.2°, C2–C3–CH₃ 71.1°, C3–H 1.094 Å, no imaginary frequency), is 2.77 kJ/mol (0.66 kcal/mol) lower in G. Frame 0 lies above it by 3.06 / 2.41 / 2.32 / 3.16 kJ/mol (0.73 / 0.58 / 0.55 / 0.76 kcal/mol) in E at M06-2X/ma-def2-SVP / M06-2X/ma-def2-TZVP / B3LYP/ma-def2-TZVP / DLPNO-CCSD(T)/ma-def2-TZVP. Nothing on the path was rerun; only the zero changed (`02_irc/zero.json`, `pipeline/unimolecular_zero.py`). Only the anti and eclip turned-away forms are minima; the third start (C3–H in line with the empty p orbital) goes downhill to the tertiary cation.
  - Entry 6 (angle against gradient size): 9 of the 70 frames over 15° in all six paths lie near this flat TS and 12 in the slow final approach to a minimum, all where the gradient is nearly zero.

## 7. Experiment

No measured value is recorded for this reaction (`00_inputs/card_notes.toml` has none).

## 8. Job times

ORCA job time (TOTAL RUN TIME in each output), 4 cores per job; several jobs ran side by side, so the wall-clock span of a level is longer than its job time divided by the number of lanes.

| Level | Jobs | Total job time (h) | Median per frame (s) | TS frame (s) | Slowest frame (s) | First job started | Last job finished |
|---|---|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | 407 | 4.1 | 36 | 35 | 73 (frame 6) | 2026-11-16 21:27 | 2026-11-17 11:25 |
| B3LYP/ma-def2-TZVP | 407 | 1.9 | 16 | 16 | 21 (frame 123) | 2026-11-17 11:27 | 2026-11-18 23:10 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 309 | 29.2 | 336 | 347 | 419 (frame 166) | 2026-11-18 23:16 | 2026-11-21 22:20 |

- **Path (IRC) gradient jobs, working level:** reactant side (backward): 622 gradient jobs, 4.6 h job time (folders `backward`, `backward_part1`); product side (forward): 621 gradient jobs, 4.8 h job time (folders `forward`). The Hessian jobs and the TS and end-minimum optimizations are not counted.
