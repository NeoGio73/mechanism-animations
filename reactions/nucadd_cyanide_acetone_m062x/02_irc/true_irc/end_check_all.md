# End check, production path: nucadd_cyanide_acetone_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). From 2026-11-16 (J. Picione) every side is optimized from its last frame; a different form of the minimum with the same bonds becomes the official end structure (the path defines which minima the TS connects); different bonds stop the switch. Verdict: **ok**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| forward | product_complex | -0.028 | 0.003 | yes | -0.031 | 0.003 | 0.003 | -0.000 (C(N)–C) | 0 | reaches the same minimum |
