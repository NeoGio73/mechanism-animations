# End check, production path: hshift_3_methyl_2_butyl_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). From 2026-11-16 (J. Picione) every side is optimized from its last frame; a different form of the minimum with the same bonds becomes the official end structure (the path defines which minima the TS connects); different bonds stop the switch. Verdict: **ok**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | secondary_eclip | -0.103 | 0.225 | yes | -0.103 | 0.055 | 0.225 | +0.013 (C2–H) | 0 | a different form of secondary_eclip with the same bonds (fails: energy); now the official end structure secondary_eclip_irc |
| forward | tertiary | +0.185 | 0.429 | yes | +0.181 | 0.125 | 0.428 | -0.091 (C3–H) | 0 | a different form of tertiary with the same bonds (fails: energy, key distances); now the official end structure tertiary_irc |
