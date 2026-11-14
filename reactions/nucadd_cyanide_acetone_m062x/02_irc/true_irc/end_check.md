# End check, production path: nucadd_cyanide_acetone_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). Verdict: **ok**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | reactant_complex | +0.941 | 0.245 | yes | -0.072 | 0.110 | 0.393 | +0.007 (C(N)–C) | 0 | reaches the same minimum |
| forward | product_complex | -0.028 | — | no | — | — | — | — | — | within the limit |
