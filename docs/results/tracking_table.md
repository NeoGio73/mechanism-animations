# IBO tracking along the six final paths

Frame-to-frame overlap = |<orbital at frame k−1 | matched orbital at frame k>| (03_orbitals/tracking.json). A swap is a frame where an orbital's largest overlap with the next frame is with an orbital other than the one it was matched to. s in bohr·amu½.

| Reaction | Tracked orbitals | Active orbital | Lowest frame-to-frame overlap | Frame | s | Swaps |
|---|---|---|---|---|---|---|
| sn2_hydroxide_bromomethane_m062x | 16 | lone pair on O → σ(C–O) | 1.000 | 84 | -3.45 | 1 |
|  |  | σ(C–Br) → lone pair on Br | 1.000 | 172 | +0.95 |  |
| e2_ethoxide_2_bromopropane_m062x | 28 | O lone pair → σ(O–H) | 1.000 | 132 | +0.00 | 52 |
|  |  | σ(Cβ–H) → π(C=C) | 1.000 | 137 | +0.25 |  |
|  |  | σ(Cα–Br) → Br lone pair | 1.000 | 142 | +0.50 |  |
| nucadd_cyanide_acetone_m062x | 17 | C lone pair (cyanide) → σ(C–C) | 1.000 | 180 | +1.10 | 0 |
|  |  | π(C=O) → O lone pair | 1.000 | 180 | +1.10 |  |
| da_butadiene_ethylene_m062x | 17 | π(C3=C4) → σ(C4–C5) | 0.714 | 189 | +0.50 | 0 |
|  |  | π(C5=C6) → σ(C1–C6) | 0.725 | 189 | +0.50 |  |
|  |  | π(C1=C2) → π(C2=C3) | 0.716 | 189 | +0.50 |  |
| hshift_3_methyl_2_butyl_m062x | 15 | σ(C3–H) → σ(C2–H) | 1.000 | 291 | +5.08 | 0 |
|  |  | σ(C1–H), hyperconjugating into C2 → ordinary σ(C1–H) | 1.000 | 339 | +7.47 |  |
|  |  | σ(C3–CH₃), three-centre bond with C2 at the start → ordinary σ(C3–CH₃) | 1.000 | 25 | -8.21 |  |
| hshift_2_3_dimethyl_2_butyl_m062x | 18 | σ(C3–H) → σ(C2–H) (three-centre at the TS) | 1.000 | 272 | +0.44 | 0 |

Lowest overlap over all tracked orbitals: sn2_hydroxide_bromomethane_m062x: 0.615 (frame 87); e2_ethoxide_2_bromopropane_m062x: 0.444 (frame 62); nucadd_cyanide_acetone_m062x: 1.000 (frame 234); da_butadiene_ethylene_m062x: 0.714 (frame 189); hshift_3_methyl_2_butyl_m062x: 1.000 (frame 291); hshift_2_3_dimethyl_2_butyl_m062x: 1.000 (frame 272).

Near-ties (another orbital overlapping above 0.5): sn2_hydroxide_bromomethane_m062x: 51; e2_ethoxide_2_bromopropane_m062x: 424; nucadd_cyanide_acetone_m062x: 0; da_butadiene_ethylene_m062x: 7; hshift_3_methyl_2_butyl_m062x: 0; hshift_2_3_dimethyl_2_butyl_m062x: 0.
