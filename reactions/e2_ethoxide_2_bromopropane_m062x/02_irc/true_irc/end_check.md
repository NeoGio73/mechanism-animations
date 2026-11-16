# End check, production path: e2_ethoxide_2_bromopropane_m062x

A side ending more than 0.5 kJ/mol above its nearest optimized minimum is accepted with a note only if Opt Freq from its last frame reaches that minimum (|dE| <= 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD <= 0.2 Å, each key distance within 0.05 Å (all-atom RMSD reported for information)). Verdict: **decision**

| Side | minimum | end above it (kJ/mol) | end RMSD (Å) | check needed | optimized: E − minimum (kJ/mol) | heavy-atom RMSD (Å) | all-atom RMSD (Å, information) | largest key-distance difference (Å) | imaginary | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| backward | reactant_complex | +0.899 | 0.315 | yes | -0.267 | 0.226 | 0.277 | -0.065 (O–H) | 0 | fails: energy, heavy-atom RMSD, key distances |
| forward | product_complex | +4.730 | 0.526 | yes | -0.279 | 0.150 | 0.173 | +0.101 (Cβ–H) | 0 | fails: energy, key distances |
