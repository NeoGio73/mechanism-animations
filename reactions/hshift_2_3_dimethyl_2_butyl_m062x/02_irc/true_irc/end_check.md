# End check, production path: hshift_2_3_dimethyl_2_butyl_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). From 2026-11-16 (J. Picione) every side is optimized from its last frame; a different form of the minimum with the same bonds becomes the official end structure (the path defines which minima the TS connects); different bonds stop the switch. Verdict: **ok**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | open_cation | +1.468 | 0.260 | yes | -0.059 | 0.009 | 0.029 | +0.000 (C2–H) | 0 | reaches the same minimum |
| forward | bridged_c2 | -0.003 | 0.004 | yes | -0.002 | 0.001 | 0.004 | -0.005 (C3–H) | 0 | last kept frame matches the bridged reference; optimization lands nearest bridged_c2 (RMSD 0.004 Å, E − bridged reference -0.002 kJ/mol) |
