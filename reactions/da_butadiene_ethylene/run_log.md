# Run log

Every ORCA job, its wall time, and what changed when a job failed. Times are wall-clock on the Omen (Ryzen 7 9800X3D).

| When | Job | Keywords | Wall time | Result | Note |
|---|---|---|---|---|---|
| 2026-11-09 21:04 | `01_ts/fragments/butadiene_s-cis/butadiene_s-cis` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) Opt Freq` | 243 s | E = -155.78653276 Eh; opt converged; 0 imaginary freq | reactant C=CC=C |
| 2026-11-09 21:05 | `01_ts/fragments/ethylene/ethylene` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) Opt Freq` | 61 s | E = -78.48227531 Eh; opt converged; 0 imaginary freq | reactant C=C |
| 2026-11-09 21:13 | `01_ts/fragments/cyclohexene/cyclohexene` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) Opt Freq` | 447 s | E = -234.34863635 Eh; opt converged; 0 imaginary freq | product C1CCC=CC1 |
| 2026-11-09 21:15 | `01_ts/fragments/butadiene_s-trans/butadiene_s-trans` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) Opt Freq` | 170 s | E = -155.79179810 Eh; opt converged; 0 imaginary freq | reference C=CC=C |
| 2026-11-09 21:21 | `01_ts/preopt/preopt` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) Opt` | 336 s | E = -234.24648384 Eh; opt converged | TS guess relaxed with forming/breaking bonds frozen |
| 2026-11-09 21:26 | `01_ts/ts/ts` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(toluene) OptTS Freq` | 319 s | E = -234.24672693 Eh; opt converged; 1 imaginary freq (-516.3 cm⁻¹) | TS search from guess |
| 2026-11-09 21:46 | — | — | — | — | STATUS 21:46: TS accepted as synchronous (C1–C6 = C4–C5 = 2.265 Å, 516i cm⁻¹, stretch-rate ratio 1.000). IRC running (forward direction first, about 1.5 min per step; expected to finish about 01:00–02:00 on 2026-11-10). Queued after a clean IRC (J. Picione: continue): endpoints, singlepoints, barriers, diagnostics, gas_phase (gas-phase single points, ΔH‡ and Ea = ΔH‡ + 2RT at 800 K, provisional), orbital_sps, orbital_track, orbital_cubes, meshes, video, card (console_rest.txt); then the coupled-cluster path pass (path_singlepoints, profile_levels; console_path_sp.txt). A failed stage stops the queue; nothing is restarted with changed settings. |
