# End check, production path: hshift_3_methyl_2_butyl_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). Verdict: **decision**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | secondary_eclip | +15.433 | 0.597 | yes | +0.039 | 0.052 | 0.214 | +0.012 (C2–H) | 0 | reaches the same minimum |
| forward | secondary_eclip | +15.499 | 0.699 | yes | -42.613 | 0.628 | 1.152 | -0.998 (C2–H) | 0 | fails: energy, heavy-atom RMSD, key distances |
