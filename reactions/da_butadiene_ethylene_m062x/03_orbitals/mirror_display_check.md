# Diels–Alder: is the mirror image of the computed orbitals exact? (2026-11-25)

**Result: no. The reflected set fails at 27 of the 302 frames, so it is not displayed.** J. Picione
decides what is shown instead (computed clockwise set with the drawing flipped, or re-localization started
from the reflected orbitals).

**Background.** The production path (M06-2X-D3(0)/ma-def2-SVP, SMD) localized the three active IBOs in the
clockwise pattern:
- π(C3=C4) → σ(C4–C5)
- π(C5=C6) → σ(C1–C6)
- π(C1=C2) → π(C2=C3)

The pilot path (B3LYP-D3(BJ)/ma-def2-SVP, `reactions/da_butadiene_ethylene`) localized them counterclockwise, the
pattern written in reaction.toml. The two are mirror images through the TS mirror plane (C1↔C4, C2↔C3, C5↔C6)
and describe the same mechanism.

**Tracking itself is continuous.** The lowest frame-to-frame overlap of the three active orbitals is **0.714**,
at frame 189 (s = +0.50 bohr·amu½). The same value is the lowest over all 17 tracked orbitals. Overlaps are below 0.9 at only 6 frames, all between frames 183 and 198 (s = +0.20 to +0.95). There the three
active orbitals rotate quickly into one another. The pilot's lowest active-orbital overlap was 0.827.

## The two tests, at every frame (limits set by J. Picione)

| Test | Limit | Worst value | Where | Frames that fail |
|---|---|---|---|---|
| Geometry: largest atom displacement between the frame reflected through the viewer's mirror plane and the frame itself (atoms swapped) | 0.001 Å | 0.0040 Å | frame 299, s = +5.96 | 0–2, 4, 6 (s = −8.94 to −8.65; 0.0010–0.0035 Å) and 299 (0.0040 Å) |
| Orbitals: smallest singular value of the 3×3 overlap matrix between the reflected and the computed active orbitals | 0.999 | 0.9924 (determinant 0.9848) | frame 181, s = +0.10 | 179–199 (s = 0.00 to +1.00, from the TS onwards) |

**At the TS (frame 179)** the geometry is mirror-symmetric to 0.0000015 Å, but the orbital value is 0.9953.

**Interpretation.** The reflection is an exact symmetry of the TS. From the TS to s = +1.0, however, the clockwise
and counterclockwise localizations are two different sets of orbitals. They do not differ only by a relabelling
inside the active three, so the mirror image of the computed set is not the computed set there.

**At the 275 frames that pass:**
- geometry within 0.0008 Å;
- orbital value at least 0.9990 (median 0.9999).

**Product end** (last 30 frames, s = +4.65 to +5.97):
- geometry 0.0040 Å at frame 299 alone; every other frame is within 0.0008 Å;
- orbital value at least 0.9990.

So the path stays mirror-symmetric (within 0.0008 Å, apart from the single frame 299) to its end. Its end is
therefore not cyclohexene's own minimum, which is not mirror-symmetric.

**Method.** Each frame is reflected through its own mirror plane, fitted to the atom pairs as the viewer does.
The orbitals are compared on the cube files: the reflected orbital is interpolated trilinearly onto the grid of
the computed one, and each orbital is normalized on its grid. Code: `pipeline/mirror_display_check.py`. Data: `mirror_display_check.json`, which
also holds the 3×3 matrix at every frame.
