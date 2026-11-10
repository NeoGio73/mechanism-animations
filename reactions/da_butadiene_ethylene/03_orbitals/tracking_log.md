# Orbital tracking log

110 IRC frames, 17 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | σ(C–C) (C0–C1) | σ(C–C) (C1–C0) | no | 0.995 (58) |
| 1 | σ(C–C) (C0–C1) | σ(C–C) (C5–C0) | **yes** | 0.835 (60) |
| 2 | σ(C–C) (C2–C1) | σ(C–C) (C2–C1) | no | 1.000 (58) |
| 3 | σ(C–C) (C3–C2) | σ(C–C) (C2–C3) | no | 0.999 (60) |
| 4 | σ(C–C) (C3–C2) | σ(C–C) (C1–C2) | **yes** | 0.827 (60) |
| 5 | σ(C–C) (C5–C4) | σ(C–C) (C5–C4) | no | 0.999 (3) |
| 6 | σ(C–C) (C4–C5) | σ(C–C) (C4–C3) | **yes** | 0.838 (60) |
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
| π(C1=C2) → π(C2=C3) | — | — | **NO** |
| π(C3=C4) → σ(C4–C5) | — | — | **NO** |
| π(C5=C6) → σ(C1–C6) | — | — | **NO** |

Orbitals that change but were not expected: [1, 4, 6]
Frames where an unexpected orbital changes its atoms: [58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109]

Canonical orbitals (sign tracking only): lowest overlap HOMO 0.000, LUMO 0.006
