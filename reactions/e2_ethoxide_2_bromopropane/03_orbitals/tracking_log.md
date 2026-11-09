# Orbital tracking log

146 IRC frames, 28 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.991 (88) |
| 1 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.604 (9) |
| 2 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.624 (3) |
| 3 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.721 (7) |
| 4 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.491 (2) |
| 5 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.503 (10) |
| 6 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.988 (34) |
| 7 | lone pair on Br (Br3) | lone pair on Br (Br3) | no | 0.988 (34) |
| 8 | lone pair on O (O13) | lone pair on O (O13) | no | 0.999 (30) |
| 9 | lone pair on O (O13) | lone pair on O (O13) | no | 0.998 (30) |
| 10 | lone pair on O (O13) | σ(O–H) (O13–H5) | **yes** | 0.990 (30) |
| 11 | σ(C–C) (C1–C0) | σ(C–C) (C1–C0) | no | 1.000 (29) |
| 12 | σ(C–C) (C1–C2) | σ(C–C) (C1–C2) | no | 1.000 (30) |
| 13 | σ(C–Br) (Br3–C1) | lone pair on Br (Br3) | **yes** | 0.997 (30) |
| 14 | σ(C–H) (C0–H4) | σ(C–H) (C0–H4) | no | 0.999 (30) |
| 15 | σ(C–H) (C0–H5) | σ(C–C) (C0–C1) | **yes** | 0.991 (30) |
| 16 | σ(C–H) (C0–H6) | σ(C–H) (C0–H6) | no | 0.999 (30) |
| 17 | σ(C–H) (C1–H7) | σ(C–H) (C1–H7) | no | 1.000 (30) |
| 18 | σ(C–H) (C2–H8) | σ(C–H) (C2–H8) | no | 1.000 (30) |
| 19 | σ(C–H) (C2–H9) | σ(C–H) (C2–H9) | no | 1.000 (30) |
| 20 | σ(C–H) (C2–H10) | σ(C–H) (C2–H10) | no | 1.000 (30) |
| 21 | σ(C–C) (C11–C12) | σ(C–C) (C12–C11) | no | 1.000 (29) |
| 22 | σ(C–O) (O13–C12) | σ(C–O) (O13–C12) | no | 1.000 (29) |
| 23 | σ(C–H) (C11–H14) | σ(C–H) (C11–H14) | no | 1.000 (30) |
| 24 | σ(C–H) (C11–H15) | σ(C–H) (C11–H15) | no | 1.000 (2) |
| 25 | σ(C–H) (C11–H16) | σ(C–H) (C11–H16) | no | 1.000 (29) |
| 26 | σ(C–H) (C12–H17) | σ(C–H) (C12–H17) | no | 1.000 (30) |
| 27 | σ(C–H) (H18–C12) | σ(C–H) (C12–H18) | no | 1.000 (30) |

Frames with a match below 0.7: 32
- frame 1, orbital 4: overlap 0.534
- frame 2, orbital 1: overlap 0.619
- frame 2, orbital 2: overlap 0.682
- frame 2, orbital 4: overlap 0.491
- frame 2, orbital 5: overlap 0.686
- frame 3, orbital 2: overlap 0.624
- frame 3, orbital 4: overlap 0.568
- frame 3, orbital 5: overlap 0.674
- frame 4, orbital 1: overlap 0.641
- frame 5, orbital 4: overlap 0.651
- frame 7, orbital 5: overlap 0.647
- frame 9, orbital 1: overlap 0.604
- frame 9, orbital 4: overlap 0.681
- frame 10, orbital 1: overlap 0.630
- frame 10, orbital 5: overlap 0.503
- frame 14, orbital 2: overlap 0.656
- frame 14, orbital 4: overlap 0.497
- frame 15, orbital 4: overlap 0.584
- frame 15, orbital 5: overlap 0.609
- frame 17, orbital 4: overlap 0.492
- frame 24, orbital 2: overlap 0.649
- frame 27, orbital 1: overlap 0.650
- frame 27, orbital 4: overlap 0.694
- frame 33, orbital 1: overlap 0.682
- frame 34, orbital 1: overlap 0.689
- frame 35, orbital 4: overlap 0.683
- frame 36, orbital 4: overlap 0.692
- frame 39, orbital 4: overlap 0.698
- frame 43, orbital 1: overlap 0.611
- frame 46, orbital 1: overlap 0.690
- frame 46, orbital 4: overlap 0.670
- frame 51, orbital 1: overlap 0.696

## Expected orbital changes (reaction.toml)

| Expected | Tracked orbital | Lowest overlap | Found? |
|---|---|---|---|
| O lone pair → σ(O–H) | 10 | 0.99 | yes |
| σ(Cβ–H) → π(C=C) | 15 | 0.991 | yes |
| σ(Cα–Br) → Br lone pair | 13 | 0.997 | yes |

Orbitals that change but were not expected: none
Frames where an unexpected orbital changes its atoms: none

Canonical orbitals (sign tracking only): lowest overlap HOMO 0.953, LUMO 0.992
