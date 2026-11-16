# Angle check, production path: hshift_2_3_dimethyl_2_butyl_m062x

Mass-weighted angle, M^1/2 Δx (central chord) vs −M^−1/2 g; median / 90th percentile (frames). Rule: each zone: median < 5°, 90th percentile < 15°. Left out: trailing frames at the end minimum: |dE| < 0.005 kJ/mol from the previous frame and RMS gradient < 5e-05 Eh/bohr (J. Picione, 2026-11-14). Path: **FAILS (flagged)**

| Side | ended | frames | left out at the end minimum | abs(s) ≤ 1 | beyond | passes | nearest optimized minimum at the end | RMSD (Å) | end above it (kJ/mol) |
|---|---|---|---|---|---|---|---|---|---|
| backward | stopped at the early comparison (J. Picione, 2026-11-16: compare once the rerun has 50 frames past the old side's last frame, then stop both runs and keep their frames) | 264 | 0 | 9.61° / 39.28° (21) | 5.74° / 10.18° (242) | **NO** | open_cation | 0.260 | +1.468 |
| forward | own stopping rule | 129 | 0 | 25.16° / 44.82° (20) | 4.27° / 8.38° (108) | **NO** | bridged_c2 | 0.427 | -11.415 |
