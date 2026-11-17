# Core-count check: SN2 TS single points

Identical input except `%pal nprocs` (4 → 2); clean folders; SCF noise floor 0.034 kJ/mol (VeryTightSCF).

| Level | E, original (Eh) | E, check (Eh) | ΔE (Eh) | ΔE (kJ/mol) | wall time, check (s) |
|---|---|---|---|---|---|
| M06-2X/ma-def2-TZVP | -2690.018420385918 | -2690.018420385922 | -4.09e-12 | -1.07e-08 | 21 |
| DLPNO-CCSD(T)/ma-def2-TZVP | -2688.351838728033 | -2688.351838728029 | 3.64e-12 | 9.55e-09 | 187 |
