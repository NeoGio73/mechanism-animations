# Orbital tracking log

543 IRC frames, 18 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.
Bond type from the s-character of each two-atom orbital: π below 0.10; when every orbital on a multiple bond is above it, the bond was localized as bent bonds (σ and π mixed; the two need not be equivalent). Atom numbers in brackets are 0-based.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | σ(C–C) (C0–C1) | σ(C–C) (C1–C0) | no | 1.000 (272) |
| 1 | σ(C–H) (C1–H2) | σ(C–H) (C0–H2) | **yes** | 1.000 (272) |
| 2 | σ(C–C) (C0–C3) | σ(C–C) (C0–C3) | no | 1.000 (272) |
| 3 | σ(C–C) (C0–C4) | σ(C–C) (C0–C4) | no | 1.000 (272) |
| 4 | σ(C–C) (C1–C5) | σ(C–C) (C1–C5) | no | 1.000 (272) |
| 5 | σ(C–C) (C1–C6) | σ(C–C) (C1–C6) | no | 1.000 (272) |
| 6 | σ(C–H) (C3–H7) | σ(C–H) (C3–H7) | no | 1.000 (272) |
| 7 | σ(C–H) (C3–H8) | σ(C–H) (C3–H8) | no | 1.000 (217) |
| 8 | σ(C–H) (C3–H9) | σ(C–H) (C3–H9) | no | 1.000 (272) |
| 9 | σ(C–H) (C4–H10) | σ(C–H) (C4–H10) | no | 1.000 (272) |
| 10 | σ(C–H) (C4–H11) | σ(C–H) (C4–H11) | no | 1.000 (175) |
| 11 | σ(C–H) (C4–H12) | σ(C–H) (C4–H12) | no | 1.000 (272) |
| 12 | σ(C–H) (C5–H13) | σ(C–H) (C5–H13) | no | 1.000 (351) |
| 13 | σ(C–H) (C5–H14) | σ(C–H) (C5–H14) | no | 1.000 (368) |
| 14 | σ(C–H) (C5–H15) | σ(C–H) (C5–H15) | no | 1.000 (347) |
| 15 | σ(C–H) (C6–H16) | σ(C–H) (C6–H16) | no | 1.000 (272) |
| 16 | σ(C–H) (C6–H17) | σ(C–H) (C6–H17) | no | 1.000 (326) |
| 17 | σ(C–H) (C6–H18) | σ(C–H) (C6–H18) | no | 1.000 (272) |

Frames with a match below 0.7: 0

## Expected orbital changes (reaction.toml)

| Expected | Tracked orbital | Lowest overlap | Found? |
|---|---|---|---|
| σ(C3–H) → σ(C2–H) (three-centre at the TS) | 1 | 1.0 | yes |

Orbitals that change but were not expected: none
Frames where an unexpected orbital changes its atoms: none

Canonical orbitals (sign tracking only): lowest overlap HOMO 1.000, LUMO 0.997
