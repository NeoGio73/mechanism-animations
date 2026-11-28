# IRC from ts_anti_rotation (J. Picione, 2026-11-27)

Started 2026-11-27 15:06:40; job time 6.11 h; 14 h limit not reached. Working level, SMD dichloromethane, 4 cores, %maxcore 2000; pysisyphus Gonzalez–Schlegel, step 0.05 bohr·amu½, exact TS Hessian from ts.hess; production stopping rule (no [irc_end_min]).

| Side | Ended | Points | Last frame − TS | Opt Freq: imaginary | Reached | Core RMSD / ΔE to each minimum |
|---|---|---|---|---|---|---|
| forward | own stopping rule | 172 | -64.45 kJ/mol | none | tertiary_irc | secondary_anti 0.785 Å / -39.49 kJ/mol; secondary_eclip_irc 0.939 Å / -42.55 kJ/mol; tertiary_irc 0.006 Å / -0.10 kJ/mol; secondary_anti_swapped 0.811 Å / -39.49 kJ/mol |
| backward | own stopping rule | 158 | -24.80 kJ/mol | none | secondary_anti | secondary_anti 0.001 Å / +0.17 kJ/mol; secondary_eclip_irc 1.083 Å / -2.89 kJ/mol; tertiary_irc 0.783 Å / +39.56 kJ/mol; secondary_anti_swapped 0.069 Å / +0.17 kJ/mol |

**Reading (rule fixed before the run):** secondary_anti + tertiary_irc: the present description stands.
