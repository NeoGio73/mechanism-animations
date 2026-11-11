# Run log

Every ORCA job, its wall time, and what changed when a job failed. Times are wall-clock on the Omen (Ryzen 7 9800X3D).

| When | Job | Keywords | Wall time | Result | Note |
|---|---|---|---|---|---|
| 2026-11-11 11:42 | `01_ts/fragments/hydroxide/hydroxide` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 25 s | E = -75.86184068 Eh; opt converged; 0 imaginary freq | reactant [OH-] |
| 2026-11-11 11:44 | `01_ts/fragments/bromomethane/bromomethane` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 96 s | E = -2613.59128963 Eh; opt converged; 0 imaginary freq | reactant CBr |
| 2026-11-11 11:45 | `01_ts/fragments/methanol/methanol` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 73 s | E = -115.58964587 Eh; opt converged; 0 imaginary freq | product CO |
| 2026-11-11 11:45 | `01_ts/fragments/bromide/bromide` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water)` | 8 s | E = -2573.90543045 Eh | product [Br-] |
| 2026-11-11 11:48 | `01_ts/ts/ts` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water) OptTS Freq` | 201 s | E = -2689.43409938 Eh; opt converged; 1 imaginary freq (-486.9 cm⁻¹) | TS search from guess |
| 2026-11-11 12:32 | `02_irc/attempts/init_displ_test/test` | `B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water) IRC` | 124 s | E = -2689.43410109 Eh | test of Init_Displ length / Scale_Init_Displ (2 IRC iterations), not used |
| 2026-11-11 12:32 | — | — | — | — | IRC forward (started 12:2x with ORCA's default initial displacement, 2 mEh) stopped by Claude after J. Picione's decision on the initial displacement; partial run kept in 02_irc/attempts/run0_default_init_displ_stopped/. Test 02_irc/attempts/init_displ_test: 'Init_Displ length' with 'Scale_Init_Displ' accepted by ORCA 6.1.1, first point at s = 0.210 bohr·amu½ as requested. Rerun continues from the IRC with irc_init_step_s = 0.21. |
