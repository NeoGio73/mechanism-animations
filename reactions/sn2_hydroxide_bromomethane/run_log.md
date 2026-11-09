# Run log

Every ORCA job, its wall time, and what changed when a job failed. Times are wall-clock on the Omen (Ryzen 7 9800X3D).

| When | Job | Keywords | Wall time | Result | Note |
|---|---|---|---|---|---|
| 2026-11-08 23:10 | `01_ts/fragments/hydroxide/hydroxide` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 5 s | E = -75.79348904 Eh; opt converged; 0 imaginary freq | reactant [OH-] |
| 2026-11-08 23:11 | `01_ts/fragments/bromomethane/bromomethane` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 52 s | E = -2613.58484102 Eh; opt converged; 0 imaginary freq | reactant CBr |
| 2026-11-08 23:12 | `01_ts/fragments/methanol/methanol` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 32 s | E = -115.57501851 Eh; opt converged; 0 imaginary freq | product CO |
| 2026-11-08 23:12 | `01_ts/fragments/bromide/bromide` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water)` | 2 s | E = -2573.87006606 Eh | product [Br-] |
| 2026-11-08 23:13 | `01_ts/ts/ts` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) OptTS Freq` | 73 s | E = -2689.38248026 Eh; opt converged; 1 imaginary freq (-441.3 cm⁻¹) | TS search from guess |
| 2026-11-08 23:20 | `02_irc/irc` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) IRC` | 428 s | E = -2689.38246001 Eh | IRC both directions |
| 2026-11-08 23:21 | — | — | — | — | IRC with default step finished (428 s, 37 frames, path connects RC and PC) but only 7 frames span the steep central region. Archived to 02_irc/attempts/run1_default_step; rerunning with Scale_Displ_SD 0.05 for smoother animation. No other setting changed. |
| 2026-11-08 23:32 | `02_irc/irc` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) IRC` | 649 s | E = -2689.38246001 Eh | IRC both directions, Scale_Displ_SD 0.05 |
| 2026-11-08 23:33 | `02_irc/endpoints/reactant_complex/reactant_complex` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 60 s | E = -2689.39236718 Eh; opt converged; 0 imaginary freq | from IRC F end |
| 2026-11-08 23:34 | `02_irc/endpoints/product_complex/product_complex` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) Opt Freq` | 41 s | E = -2689.45133514 Eh; opt converged; 1 imaginary freq (-21.4 cm⁻¹) | from IRC B end |
| 2026-11-08 23:34 | `01_ts/sp_tzvp/ts/ts` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 7 s | E = -2689.87110965 Eh | def2-TZVP single point |
| 2026-11-08 23:34 | `01_ts/sp_tzvp/reactant_complex/reactant_complex` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 6 s | E = -2689.88633521 Eh | def2-TZVP single point |
| 2026-11-08 23:34 | `01_ts/sp_tzvp/product_complex/product_complex` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 6 s | E = -2689.92829127 Eh | def2-TZVP single point |
| 2026-11-08 23:34 | `01_ts/sp_tzvp/hydroxide/hydroxide` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 2 s | E = -75.93719672 Eh | def2-TZVP single point |
| 2026-11-08 23:34 | `01_ts/sp_tzvp/bromomethane/bromomethane` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 5 s | E = -2613.94789205 Eh | def2-TZVP single point |
| 2026-11-08 23:35 | `01_ts/sp_tzvp/methanol/methanol` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 3 s | E = -115.71927735 Eh | def2-TZVP single point |
| 2026-11-08 23:35 | `01_ts/sp_tzvp/bromide/bromide` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 3 s | E = -2574.20727991 Eh | def2-TZVP single point |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/ts/ts` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 9 s | E = -2689.87828806 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/reactant_complex/reactant_complex` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 8 s | E = -2689.89426958 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/product_complex/product_complex` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 8 s | E = -2689.93529517 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/hydroxide/hydroxide` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 2 s | E = -75.94758871 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/bromomethane/bromomethane` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 5 s | E = -2613.94847704 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/methanol/methanol` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 3 s | E = -115.72095661 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/bromide/bromide` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 3 s | E = -2574.21389614 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:35 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/ts/ts` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 9 s | E = -2690.01947927 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/reactant_complex/reactant_complex` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 7 s | E = -2690.04214202 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/product_complex/product_complex` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 8 s | E = -2690.08739601 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/hydroxide/hydroxide` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 2 s | E = -75.94657201 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/bromomethane/bromomethane` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 6 s | E = -2614.09589146 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/methanol/methanol` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 4 s | E = -115.72344746 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:36 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/bromide/bromide` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 3 s | E = -2574.36242952 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:37 | — | — | — | — | IBO single points at 60 IRC frames (`03_orbitals/frames/`, per-frame log in `03_orbitals/run_log_frames.md`): 248 s in total this run. |
| 2026-11-08 23:38 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/ts/ts` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 123 s | E = -2688.35291338 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:40 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/reactant_complex/reactant_complex` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 112 s | E = -2688.37630229 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:41 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/product_complex/product_complex` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 80 s | E = -2688.41902611 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:41 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/hydroxide/hydroxide` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 6 s | E = -75.84322129 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:42 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/bromomethane/bromomethane` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 53 s | E = -2612.53507153 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:42 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/methanol/methanol` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 18 s | E = -115.55016256 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:43 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/bromide/bromide` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 14 s | E = -2572.86837763 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:44 | `02_irc/endpoints/product_complex/product_complex_tight` | `B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water) TightOpt Freq` | 27 s | E = -2689.45135176 Eh; opt converged; 0 imaginary freq | rerun: product_complex had imaginary mode(s) [-21.4] cm-1; displaced 0.1 Å along it, TightOpt |
| 2026-11-08 23:44 | `01_ts/sp_tzvp/product_complex/product_complex` | `B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 4 s | E = -2689.92831483 Eh | def2-TZVP single point |
| 2026-11-08 23:44 | `01_ts/diagnostics/B3LYP_ma-def2-TZVP/product_complex/product_complex` | `B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 4 s | E = -2689.93531831 Eh | diagnostic B3LYP_ma-def2-TZVP |
| 2026-11-08 23:44 | `01_ts/diagnostics/M06-2X_ma-def2-TZVP/product_complex/product_complex` | `M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF CPCM(water)` | 5 s | E = -2690.08740767 Eh | diagnostic M06-2X_ma-def2-TZVP |
| 2026-11-08 23:46 | `01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/product_complex/product_complex` | `DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF CPCM(water)` | 80 s | E = -2688.41902804 Eh | diagnostic DLPNO-CCSD(T)_ma-def2-TZVP |
| 2026-11-08 23:51 | — | — | — | — | Single points at all 60 IRC frames at M06-2X_ma-def2-TZVP (`02_irc/path_sp/M06-2X_ma-def2-TZVP/`, log `02_irc/path_sp/run_log_path_sp.md`): 497 s of job time; highest frame 18, 64.0 kJ/mol above separated reactants. |
| 2026-11-09 00:48 | — | — | — | — | Single points at all 60 IRC frames at DLPNO-CCSD(T)_ma-def2-TZVP (`02_irc/path_sp/DLPNO-CCSD(T)_ma-def2-TZVP/`, log `02_irc/path_sp/run_log_path_sp.md`): 6820 s of job time; highest frame 18, 70.7 kJ/mol above separated reactants. |
| 2026-11-09 03:35 | — | — | — | — | Presentation only, no new calculations: default energy curve set to DLPNO-CCSD(T)/ma-def2-TZVP (J. Picione). Viewer data and videos rebuilt from the existing frames; <id>.mp4 now has the coupled-cluster inset, the working-level video is <id>_B3LYP.mp4. |
