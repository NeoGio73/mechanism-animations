# Relaxed scan of the degenerate methyl shift (J. Picione, 2026-11-23)

q = d(C3–CH₃) − d(C2–CH₃); anti q = −0.130 Å, swapped anti q = +0.130 Å. Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane; ΔE electronic, zero = secondary_anti. Each point optimized with q fixed, started from the previous point.

| Point | q (Å) | C3–CH₃ (Å) | C2–CH₃ (Å) | H10–C4–C2–C3 (°) | ΔE kJ/mol (kcal/mol) | Converged | Flag |
|---|---|---|---|---|---|---|---|
| p00 | -0.130 | 1.673 | 1.803 | -168.7 | -0.00 (-0.00) | True |  |
| p01 | -0.110 | 1.682 | 1.792 | -168.8 | +0.04 (+0.01) | True |  |
| p02 | -0.090 | 1.695 | 1.785 | -167.7 | +0.16 (+0.04) | True |  |
| p03 | -0.070 | 1.708 | 1.778 | -165.6 | +0.35 (+0.08) | True |  |
| p04 | -0.050 | 1.720 | 1.770 | -162.6 | +0.57 (+0.14) | True |  |
| p05 | -0.030 | 1.730 | 1.760 | -162.5 | +0.82 (+0.20) | True |  |
| p06 | -0.010 | 1.744 | 1.754 | -157.6 | +1.07 (+0.26) | True |  |
| p07 | +0.010 | 1.756 | 1.746 | -154.5 | +1.34 (+0.32) | True |  |
| p08 | +0.030 | 1.768 | 1.738 | -150.0 | +1.62 (+0.39) | True |  |
| p09_from_next | +0.050 | 1.768 | 1.718 | -108.4 | +0.86 (+0.21) | True |  |
| p10 | +0.070 | 1.779 | 1.709 | -108.2 | +0.64 (+0.15) | True |  |
| p11 | +0.090 | 1.785 | 1.695 | -109.8 | +0.45 (+0.11) | True |  |
| p12 | +0.110 | 1.794 | 1.684 | -109.7 | +0.33 (+0.08) | True |  |
| p13 | +0.130 | 1.803 | 1.673 | -110.3 | +0.28 (+0.07) | True |  |

Flagged (rotamer jumps, not part of the profile): p09 dihedral -108.5°, ΔE +0.86 kJ/mol; p09_retry dihedral -108.4°, ΔE +0.86 kJ/mol

Highest accepted point: p08.
- Freq: {'imaginary_cm': [-106.09194708361268], 'lowest_real_cm': [216.9, 227.3, 234.8, 281.3], 'share_methyl_H': 0.946, 'share_methyl_C': 0.011, 'dq_for_0p05_step': 0.0158, 'moves_methyl_between_C3_and_C2': False}
- DLPNO-CCSD(T)/ma-def2-TZVP: {'E_eh': -196.537471954604, 'dE_kJ_vs_anti': 0.552, 'T1': 0.009519365}
