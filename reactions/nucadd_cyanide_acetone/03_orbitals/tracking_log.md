# Orbital tracking log

120 IRC frames, 17 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.
Bond type from the s-character of each two-atom orbital: π below 0.10; when every orbital on a multiple bond is above it, the bond was localized as bent bonds (σ and π mixed; the two need not be equivalent). Atom numbers in brackets are 0-based.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | lone pair on O (O1) | lone pair on O (O1) | no | 1.000 (101) |
| 1 | lone pair on O (O1) | lone pair on O (O1) | no | 1.000 (100) |
| 2 | lone pair on C (C4) | σ(C–C) (C4–C0) | **yes** | 0.992 (101) |
| 3 | lone pair on N (N5) | lone pair on N (N5) | no | 1.000 (101) |
| 4 | π(C=O) (O1–C0) | lone pair on O (O1) | **yes** | 0.994 (101) |
| 5 | σ(C–O) (O1–C0) | σ(C–O) (O1–C0) | no | 0.999 (101) |
| 6 | σ(C–C) (C0–C2) | σ(C–C) (C0–C2) | no | 0.999 (101) |
| 7 | σ(C–C) (C0–C3) | σ(C–C) (C0–C3) | no | 0.999 (101) |
| 8 | σ(C–N) (N5–C4) | σ(C–N) (C4–N5) | no | 0.993 (101) |
| 9 | π(C=N) (N5–C4) | π(C=N) (N5–C4) | no | 0.993 (101) |
| 10 | π(C=N) (N5–C4) | π(C=N) (N5–C4) | no | 1.000 (103) |
| 11 | σ(C–H) (C2–H6) | σ(C–H) (C2–H6) | no | 0.999 (101) |
| 12 | σ(C–H) (C2–H7) | σ(C–H) (C2–H7) | no | 1.000 (100) |
| 13 | σ(C–H) (C2–H8) | σ(C–H) (C2–H8) | no | 1.000 (97) |
| 14 | σ(C–H) (C3–H9) | σ(C–H) (C3–H9) | no | 0.999 (101) |
| 15 | σ(C–H) (C3–H10) | σ(C–H) (C3–H10) | no | 1.000 (97) |
| 16 | σ(C–H) (C3–H11) | σ(C–H) (C3–H11) | no | 1.000 (100) |

Frames with a match below 0.7: 0

## Expected orbital changes (reaction.toml)

| Expected | Tracked orbital | Lowest overlap | Found? |
|---|---|---|---|
| C lone pair (cyanide) → σ(C–C) | 2 | 0.992 | yes |
| π(C=O) → O lone pair | 4 | 0.994 | yes |

Orbitals that change but were not expected: none
Frames where an unexpected orbital changes its atoms: none

Canonical orbitals (sign tracking only): lowest overlap HOMO 0.002, LUMO 0.968
