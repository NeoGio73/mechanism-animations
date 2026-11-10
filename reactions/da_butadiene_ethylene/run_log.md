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
