# End check, production path: sn2_hydroxide_bromomethane_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). From 2026-11-16 (J. Picione) every side is optimized from its last frame; a different form of the minimum with the same bonds becomes the official end structure (the path defines which minima the TS connects); different bonds stop the switch. Verdict: **ok**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | reactant_complex | +0.099 | 0.077 | yes | +0.036 | 0.055 | 0.077 | -0.116 (C–O) | 0 | a different form of reactant_complex with the same bonds (fails: key distances); now the official end structure reactant_complex_irc |
| forward | product_complex | +0.232 | 0.031 | yes | +0.011 | 0.002 | 0.004 | -0.001 (C–Br) | 0 | reaches the same minimum |
