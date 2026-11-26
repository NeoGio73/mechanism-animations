# Orbital tracking log

407 IRC frames, 15 localized valence IBOs per frame, matched by maximum |overlap| between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are positive. Frames where a best match falls below 0.7 are flagged.
Bond type from the s-character of each two-atom orbital: π below 0.10; when every orbital on a multiple bond is above it, the bond was localized as bent bonds (σ and π mixed; the two need not be equivalent). Atom numbers in brackets are 0-based.

| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |
|---|---|---|---|---|
| 0 | σ(C–C) (C1–C0) | σ(C–C) (C1–C0) | no | 1.000 (292) |
| 1 | σ(C–H) (C1–H2) | σ(C–H) (C1–H2) | no | 1.000 (302) |
| 2 | σ(C–C) (C3–C1) | σ(C–C) (C3–C1) | no | 1.000 (291) |
| 3 | σ(C–C) (C3–C4) | σ(C–C) (C3–C4) | no | 1.000 (273) |
| 4 | σ(C–H) (C0–H6) | σ(C–H) (C0–H6) | no | 1.000 (343) |
| 5 | σ(C–H) (C0–H7) | σ(C–H) (C0–H7) | no | 1.000 (364) |
| 6 | delocalized over C–C–H (C0–H8–C1) | σ(C–H) (C0–H8) | **yes** | 1.000 (339) |
| 7 | σ(C–H) (C3–H9) | σ(C–H) (C1–H9) | **yes** | 1.000 (291) |
| 8 | σ(C–H) (C4–H10) | σ(C–H) (C4–H10) | no | 1.000 (395) |
| 9 | σ(C–H) (C4–H11) | σ(C–H) (C4–H11) | no | 1.000 (397) |
| 10 | σ(C–H) (C4–H12) | σ(C–H) (C4–H12) | no | 1.000 (393) |
| 11 | σ(C–H) (C5–H13) | σ(C–H) (C5–H13) | no | 1.000 (15) |
| 12 | σ(C–H) (C5–H14) | σ(C–H) (C5–H14) | no | 1.000 (15) |
| 13 | σ(C–H) (C5–H15) | σ(C–H) (C5–H15) | no | 1.000 (287) |
| 14 | delocalized over C–C–C (C5–C1–C3) | σ(C–C) (C3–C5) | **yes** | 1.000 (25) |

Frames with a match below 0.7: 0

## Expected orbital changes (reaction.toml)

| Expected | Tracked orbital | Lowest overlap | Found? |
|---|---|---|---|
| σ(C3–H) → σ(C2–H) | 7 | 1.0 | yes |

## Further changes found by tracking, documented (not in the textbook drawing)

| Documented change | Tracked orbital | Lowest overlap | Found? | Note |
|---|---|---|---|---|
| σ(C1–H), hyperconjugating into C2 → ordinary σ(C1–H) | 6 | 1.0 | yes | hyperconjugation into the cationic centre C2 withdrawn; not an extra arrow |
| σ(C3–CH₃), three-centre bond with C2 at the start → ordinary σ(C3–CH₃) | 14 | 1.0 | yes | at the path start a three-centre bond: the methyl bridges C2 and C3 (J. Picione, 2026-11-26: not hyperconjugation); it opens to an ordinary bond before the TS; not an extra arrow |

Orbitals that change but were neither expected nor documented: none
Frames where any other orbital's atoms (share above 0.15) differ from its reactant-end atoms: 40–144

Canonical orbitals (sign tracking only): lowest overlap HOMO 0.969, LUMO 0.999
