# End check, open side (backward), tertiary hydride shift — 2026-11-15

The automatic end check compared the open-side end only with `bridged_c2`, because plain RMSD does not match it
to `open_cation` (atom numbering differs by the C2 ⇄ C3 swap and methyl exchanges). Checked here by hand with a
permutation-aware heavy-atom RMSD (every relabelling of the carbons that keeps the C–C bond graph), same rule as
the end check (docs/protocol.md): a side ending more than 0.5 kJ/mol above its minimum is accepted with a note only
if Opt Freq from its last frame reaches that minimum (|dE| ≤ 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD
≤ 0.2 Å, key distances within 0.05 Å).

| | E − open_cation (kJ/mol) | heavy-atom RMSD, permutation-aware (Å) | C2–H (Å) | C3–H (Å) | imaginary |
|---|---|---|---|---|---|
| last frame of the side | +4.106 | 0.485 | 1.9273 | 1.1253 | — |
| Opt Freq from it (`end_check/backward/opt`) | −0.015 | 0.002 | 2.0858 | 1.1017 | 0 |
| open_cation (reference) | 0 | 0 | 2.0862 | 1.1016 | 0 |

Verdict: **accepted with a note** — the side stopped 4.1 kJ/mol (1.0 kcal/mol) above the open cation, on its
slow final descent; optimizing from its last frame reaches the open cation. Script: `end_check/backward/open_side_check.py`.
