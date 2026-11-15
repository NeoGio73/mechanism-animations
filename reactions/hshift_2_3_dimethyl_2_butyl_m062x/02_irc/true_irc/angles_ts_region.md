# Angle check near the TS, tertiary hydride shift — both ways (J. Picione, 2026-11-15)

Decision (J. Picione, 2026-11-15): the end-minimum exclusion is NOT extended to the TS region. The angle check
stays as computed (`angles.md`); for the bridged side both numbers are kept for the paper.

Zone |s| ≤ 1 bohr·amu½, median / 90th percentile (number of frames); rule: median < 5°, 90th percentile < 15°.
End-minimum test: |ΔE| < 0.005 kJ/mol from the previous frame and RMS gradient < 5×10⁻⁵ Eh/bohr.

| Side | as computed (the check) | frames with angle ≥ 5° | of those meeting the end-minimum test | without those frames (information only) |
|---|---|---|---|---|
| backward (open side) | 9.61° / 39.28° (21) | 13 (frames 0–12, 54.2° → 5.5°) | 0 (energy falls 0.006–0.16 kJ/mol per step) | same as computed |
| forward (bridged side) | 25.16° / 44.82° (20) | 14 | 12 (frames 0–7, 9–11, 13); not frames 8 (40.0°, ΔE −0.0052 kJ/mol) and 12 (17.4°) | 4.25° / 24.19° (8) |

The 7 bridged-side frames used in the path (frames 0–6) all meet the end-minimum test (as computed: 36.74° / 62.94°).
The open side's TS region is rechecked at half the step after the secondary rerun (J. Picione, 2026-11-15).
