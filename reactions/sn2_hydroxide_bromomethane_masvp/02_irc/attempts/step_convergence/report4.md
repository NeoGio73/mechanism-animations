# Step-size test, four paths (SN2, B3LYP-D3(BJ)/ma-def2-SVP, SMD water)

| Path | First step set (bohr·amu½) | First steps measured (fwd, bwd) | Median step in region | Steps in region | Covers region |
|---|---|---|---|---|---|
| 0.10 (first step 0.10) | 0.100 | 0.100, 0.100 | 0.1000 | 148 | True |
| 0.05 (first step 0.05) | 0.050 | 0.050, 0.050 | 0.0500 | 297 | True |
| 0.05 (first step 0.10) | 0.100 | 0.100, 0.100 | 0.0503 | 293 | True |
| 0.025 (first step 0.10) | 0.100 | 0.100, 0.100 | 0.0249 | 576 | True |

## 0.10 vs 0.05, both first step 0.10

| Side (matched by) | Points | Largest energy difference | at distance | Largest other-bond difference | Within 0.2 kJ/mol and 0.005 Å |
|---|---|---|---|---|---|
| reactant (C–O) | 66 | -0.126 kJ/mol (-0.0301 kcal/mol) | C–O 2.461 Å | C–Br +0.0078 Å | no (5 points outside, s +2.32 to +2.65) |
| product (C–Br) | 80 | -0.288 kJ/mol (-0.0689 kcal/mol) | C–Br 2.896 Å | C–O +0.0067 Å | no (7 points outside, s -4.82 to -3.32) |

## 0.05 vs 0.025, both first step 0.10

| Side (matched by) | Points | Largest energy difference | at distance | Largest other-bond difference | Within 0.2 kJ/mol and 0.005 Å |
|---|---|---|---|---|---|
| reactant (C–O) | 135 | -0.128 kJ/mol (-0.0306 kcal/mol) | C–O 2.746 Å | C–Br +0.0024 Å | yes |
| product (C–Br) | 156 | +0.753 kJ/mol (+0.1800 kcal/mol) | C–Br 2.617 Å | C–O -0.0064 Å | no (65 points outside, s -4.91 to -1.06) |

## old 0.05 (first step 0.05) vs new 0.05 (first step 0.10)

| Side (matched by) | Points | Largest energy difference | at distance | Largest other-bond difference | Within 0.2 kJ/mol and 0.005 Å |
|---|---|---|---|---|---|
| reactant (C–O) | 135 | -0.145 kJ/mol (-0.0347 kcal/mol) | C–O 2.481 Å | C–Br +0.0030 Å | yes |
| product (C–Br) | 160 | +0.325 kJ/mol (+0.0776 kcal/mol) | C–Br 2.945 Å | C–O +0.0073 Å | no (18 points outside, s -4.99 to -4.09) |
