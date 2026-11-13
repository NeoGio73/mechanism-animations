# Tail test: mass-weighted angle along the flat tails (downhill from the same frame of the 0.05 path)

| Run | n | mw median | mw 90 % | mw largest | frames > 10° | largest turn | ended | gradients / step | Hessians | s per step | end above min (kJ/mol) | ΔE vs 0.05 path (kJ/mol) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reactant: earlier pysisyphus s0100 (same stretch) | 56 | 6.23 | 34.16 | 67.56 | 18 | 4.8 | — | nan | — | nan | nan | nan |
| reactant: earlier pysisyphus s0050 (same stretch) | 115 | 1.55 | 8.12 | 55.92 | 9 | 4.2 | — | nan | — | nan | nan | nan |
| reactant: tail_eulerpc_reactant | not run |
| reactant: tail_eulerpc_r5_reactant | 53 | 1.21 | 73.02 | 82.97 | 12 | 55.1 | energy increased | 2.02 | 11 | 57.1 | 0.223 | 0.018 |
| reactant: tail_gs_reactant | 70 | 0.70 | 25.32 | 86.05 | 16 | 72.8 | own stopping rule | 5.34 | 1 | 116.8 | 0.103 | 0.013 |
| product: earlier pysisyphus s0100 (same stretch) | 45 | 7.98 | 45.84 | 60.66 | 22 | 8.8 | — | nan | — | nan | nan | nan |
| product: earlier pysisyphus s0050 (same stretch) | 89 | 0.12 | 26.76 | 52.61 | 27 | 6.6 | — | nan | — | nan | nan | nan |
| product: tail_eulerpc_product | not run |
| product: tail_eulerpc_r5_product | not run |
| product: tail_gs_product | 57 | 0.55 | 3.84 | 22.60 | 2 | 34.3 | own stopping rule | 2.95 | 1 | 64.2 | -0.009 | 0.003 |
