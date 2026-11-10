# Orbital tracking log

110 IRC frames, 17 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.
Bond type from the s-character of each two-atom orbital: π below 0.10; when every orbital on a multiple bond is above it, the bond was localized as equivalent bent bonds. Atom numbers in brackets are 0-based.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | σ(C–C) (C0–C1) | σ(C–C) (C1–C0) | no | 0.995 (58) |
| 1 | π(C=C) (C0–C1) | σ(C–C) (C5–C0) | **yes** | 0.835 (60) |
| 2 | σ(C–C) (C2–C1) | σ(C–C) (C2–C1) | no | 1.000 (58) |
| 3 | σ(C–C) (C3–C2) | σ(C–C) (C2–C3) | no | 0.999 (60) |
| 4 | π(C=C) (C3–C2) | π(C=C) (C1–C2) | **yes** | 0.827 (60) |
| 5 | bent bond (C=C) (C5–C4) | σ(C–C) (C5–C4) | no | 0.999 (3) |
| 6 | bent bond (C=C) (C4–C5) | σ(C–C) (C4–C3) | **yes** | 0.838 (60) |
| 7 | σ(C–H) (H6–C0) | σ(C–H) (H6–C0) | no | 0.997 (59) |
| 8 | σ(C–H) (H7–C0) | σ(C–H) (H7–C0) | no | 1.000 (59) |
| 9 | σ(C–H) (H8–C1) | σ(C–H) (H8–C1) | no | 1.000 (59) |
| 10 | σ(C–H) (H9–C2) | σ(C–H) (H9–C2) | no | 1.000 (59) |
| 11 | σ(C–H) (H10–C3) | σ(C–H) (H10–C3) | no | 1.000 (59) |
| 12 | σ(C–H) (H11–C3) | σ(C–H) (H11–C3) | no | 0.998 (59) |
| 13 | σ(C–H) (C4–H12) | σ(C–H) (H12–C4) | no | 0.998 (59) |
| 14 | σ(C–H) (C4–H13) | σ(C–H) (C4–H13) | no | 0.998 (59) |
| 15 | σ(C–H) (C5–H14) | σ(C–H) (C5–H14) | no | 0.999 (59) |
| 16 | σ(C–H) (C5–H15) | σ(C–H) (H15–C5) | no | 0.999 (59) |

Frames with a match below 0.7: 0

## Expected orbital changes (reaction.toml)

| Expected | Tracked orbital | Lowest overlap | Found? |
|---|---|---|---|
| π(C1=C2) → σ(C1–C6) | 1 | 0.835 | yes |
| π(C5=C6) → σ(C4–C5) | 6 | 0.838 | yes |
| π(C3=C4) → π(C2=C3) | 4 | 0.827 | yes |

Orbitals that change but were not expected: none
Frames where an unexpected orbital changes its atoms: [58, 59]

Canonical orbitals (sign tracking only): lowest overlap HOMO 0.000, LUMO 0.006
