# Mass-weighted IRC test (pysisyphus 1.0.0, EulerPC) on SN2, against the ORCA step-test paths

Angles in degrees between the path direction (central chord) and the downhill gradient, median / largest.
plain: Δx vs −g. mass-weighted: M^½Δx vs −M^−½g. Cartesian form: Δx vs −M^−1g.

| Path | Zone | n | plain | mass-weighted | Cartesian form |
|---|---|---|---|---|---|
| s0100 | near_TS_|s|<=1 | 23 | 22.79 / 30.03 | 0.07 / 7.23 | 0.25 / 22.96 |
| s0100 | far_|s|>1 | 205 | 58.99 / 65.13 | 0.52 / 67.56 | 1.19 / 109.55 |
| s0050 | near_TS_|s|<=1 | 45 | 22.68 / 30.01 | 0.01 / 7.33 | 0.02 / 23.27 |
| s0050 | far_|s|>1 | 411 | 57.93 / 63.39 | 0.03 / 55.92 | 0.05 / 72.74 |
| orca_s010 | near_TS_|s|<=1 | 5 | 22.98 / 34.23 | 55.94 / 67.62 | 75.65 / 81.02 |
| orca_s010 | far_|s|>1 | 10 | 10.94 / 25.47 | 52.46 / 72.30 | 62.33 / 76.66 |
| orca_s005 | near_TS_|s|<=1 | 5 | 20.03 / 29.36 | 46.00 / 61.78 | 70.73 / 80.16 |
| orca_s005 | far_|s|>1 | 10 | 14.81 / 25.01 | 47.64 / 72.79 | 59.19 / 76.05 |
| orca_s005f10 | near_TS_|s|<=1 | 5 | 18.69 / 29.19 | 55.79 / 59.30 | 72.50 / 80.31 |
| orca_s005f10 | far_|s|>1 | 10 | 7.58 / 26.13 | 55.01 / 71.60 | 63.19 / 77.12 |
| orca_s0025 | near_TS_|s|<=1 | 5 | 1.05 / 18.87 | 37.58 / 51.86 | 43.70 / 71.97 |
| orca_s0025 | far_|s|>1 | 10 | 11.41 / 28.40 | 48.53 / 72.95 | 62.75 / 75.27 |

## Completeness (pysisyphus sides)

| Path | Side | Ended | Steps | Drop from TS (kJ/mol) | RMS grad last | Stopping rule met | Verdict |
|---|---|---|---|---|---|---|---|
| s0100 | reactant | rms(grad) converged! | 98 | 54.92 | 5.0e-04 | None | complete |
| s0100 | product | rms(grad) converged! | 132 | 161.74 | 4.8e-04 | None | complete |
| s0050 | reactant | rms(grad) converged! | 197 | 54.90 | 4.9e-04 | None | complete |
| s0050 | product | rms(grad) converged! | 261 | 161.68 | 5.0e-04 | 100 | complete |

## pysisyphus 0.10 vs 0.05

| Zone | Side | Points | Largest energy difference (kJ/mol) | at | Largest other-bond difference (Å) | Within limits |
|---|---|---|---|---|---|---|
| whole region | reactant (C–O) | 86 | +0.031 (+0.0074 kcal/mol) | C–O 2.768 Å | C–Br +0.0003 | yes (0 outside) |
| whole region | product (C–Br) | 120 | +0.016 (+0.0039 kcal/mol) | C–Br 3.276 Å | C–O -0.0001 | yes (0 outside) |
| |s| <= 1 | reactant (C–O) | 12 | -0.029 (-0.0070 kcal/mol) | C–O 2.328 Å | C–Br +0.0000 | yes (0 outside) |
| |s| <= 1 | product (C–Br) | 11 | +0.001 (+0.0003 kcal/mol) | C–Br 2.509 Å | C–O +0.0000 | yes (0 outside) |

## ORCA 0.10 vs pysisyphus 0.10

| Zone | Side | Points | Largest energy difference (kJ/mol) | at | Largest other-bond difference (Å) | Within limits |
|---|---|---|---|---|---|---|
| whole region | reactant (C–O) | 66 | -1.540 (-0.3681 kcal/mol) | C–O 2.445 Å | C–Br -0.0357 | no (61 outside) |
| whole region | product (C–Br) | 80 | +18.251 (+4.3620 kcal/mol) | C–Br 2.804 Å | C–O +0.0914 | no (52 outside) |
| |s| <= 1 | reactant (C–O) | 8 | -0.676 (-0.1615 kcal/mol) | C–O 2.305 Å | C–Br -0.0189 | no (6 outside) |
| |s| <= 1 | product (C–Br) | 8 | +2.334 (+0.5579 kcal/mol) | C–Br 2.504 Å | C–O +0.0253 | no (6 outside) |

## O–H on the product side (Å)

- s0100: O–H 0.9677–0.9702, largest distance from the relaxed value 0.0050
- s0050: O–H 0.9677–0.9702, largest distance from the relaxed value 0.0050

## Cost

- s0100: {"orca_calls": 460, "steps": 230, "calls_per_step": 2.0, "median_s_per_call": 21.6755, "s_per_step": 43.31536086956522}
- s0050: {"orca_calls": 916, "steps": 458, "calls_per_step": 2.0, "median_s_per_call": 21.47, "s_per_step": 42.88901091703057}
- orca_irc: {"s010": {"walls_s_2_cores": [4182], "steps": 162, "s_per_step": 25.814814814814813}, "s005f10": {"walls_s_2_cores": [7598, 8551], "steps": 318, "s_per_step": 50.783018867924525}}
