# Hydride shifts: orbital check against the mechanism, and hyperconjugation into the cationic carbons (2026-11-26)

J. Picione, 2026-11-26. All numbers come from the IBOs already on disk (`03_orbitals/frames/fNNN/f_loc.npz`); no ORCA
job was run. Code: `pipeline/hshift_hyperconj.py` (writes `03_orbitals/hyperconjugation.json` in both folders) and
the comparison in `pipeline/orbital_steps.py`. s is in bohr·amu½; energies are working-level electronic energies
(M06-2X/ma-def2-SVP, SMD dichloromethane), in kJ/mol (kcal/mol).

**Numbering.** Secondary shift: the 3-methyl-2-butyl cation, CH₃(C1)–C2H⁺–C3H(CH₃)–C4H₃, kept along the whole path.
At the product end the cationic carbon is C3 (C2 in the 2-methyl-2-butyl cation's own name). "CH₃ on C3" is the methyl
that leans over C2. Tertiary shift: the 2,3-dimethyl-2-butyl cation, (CH₃)₂C2⁺–C3H(CH₃)₂.

**Caveat (as for E2).** A "share" is the Mulliken population of one atom in one normalized IBO (the shares of an orbital
sum to 1). It is computed in the full ma-def2-SVP basis, with no cut-off. These are not IAO charges of single
orbitals, which ORCA's orca_loc cannot print (`docs/results/e2_orbital_measures.md`). They can be slightly negative.

**Threshold.** Tracking treats an atom as carrying an orbital when its share is above **0.15**
(`POP_CUT`, `pipeline/orbital_steps.py`). An orbital has **changed** when its set of carrying atoms differs between the
two ends of the path. Since today the comparison also flags every frame where an orbital that is neither expected nor
documented has a carrier set different from its reactant-end set.

## 1. Secondary shift: comparison with the mechanism

The expected list is now one change, the arrow the textbook draws: **σ(C3–H) → σ(C2–H)**. It was written from the
mechanism, with atoms from `[key_atoms]`, set on 2026-11-13 from the TS before any orbital calculation. The other two
changes are recorded under a separate key, `[[documented_orbital_changes]]`.

| Category | Orbital | Found | Lowest frame-to-frame overlap |
|---|---|---|---|
| Expected (mechanism) | σ(C3–H) → σ(C2–H) | yes | 1.000 |
| Documented (tracking; not in the drawing) | σ(C1–H), hyperconjugating into C2 → ordinary σ(C1–H) | yes | 1.000 |
| Documented (tracking; not in the drawing) | σ(C3–CH₃), hyperconjugating into C2 → ordinary σ(C3–CH₃) | yes | 1.000 |

- **Result:** 1 expected change, found; 2 further changes, documented.
- **Unexpected changes:** none.
- **Flagged frames: 40–144** (s −7.46 to −2.25). This is new, because before today nothing was compared. Two C–H
  orbitals with the same atoms at both ends cross the 0.15 threshold on C3 for part of the C2–C3 rotation:
  - σ(C4–H): up to 0.210 on C3 at s −5.76 (frames 40–144); 0.082 at the reactant end, 0.049 at the product end;
  - σ(C–H) of the CH₃ on C3: up to 0.158 on C3 at s −3.81 (frames 92–136); −0.030 at the reactant end, 0.110 at the
    product end.

  These are transient crossings, not changes. They show that the 0.15 threshold is close to the share a C–H orbital
  can have on a neighbouring carbon.

## 2. The two documented orbitals: share on the cationic carbon C2 along the path

| Point | s | σ(C1–H): share on C2 | σ(C3–CH₃): share on C2 | σ(C3–CH₃): share on C3 | C2–CH₃ | C3–CH₃ | C2–C3–CH₃ angle |
|---|---|---|---|---|---|---|---|
| Reactant end (frame 0) | −9.57 | 0.171 | 0.359 | 0.220 | 1.70 Å | 1.80 Å | 63° |
| TS (frame 189; the C2–C3 rotation TS) | 0.00 | 0.040 | 0.033 | 0.543 | 2.42 Å | 1.56 Å | 107° |
| Product end (frame 406) | +10.81 | 0.003 | −0.001 | 0.621 | 2.52 Å | 1.46 Å | 119° |

Where the share on C2 changes fastest:

| Orbital | Fastest change, s | Half of the change done by s | Rate there (per bohr·amu½) |
|---|---|---|---|
| σ(C1–H) | −7.46 | −7.26 | 0.072 |
| σ(C3–CH₃) | −8.06 | −6.91 | 0.128 |

For comparison, the hydrogen moves later. C2–H changes fastest at s +3.87 and C3–H at s +5.47. The migrating σ(C3–H)
has 0.067 of its share on C2 at the reactant end and 0.050 at the TS; it then becomes σ(C2–H), with 0.534 on C2 at the
product end.

- **The start of this path is partly methyl-bridged.** At the reactant end the C3–CH₃ bonding pair has more of its
  share on C2 (0.36) than on C3 (0.22). The methyl sits 1.70 Å from C2 and 1.80 Å from C3, and C2–C3 is 1.40 Å. This
  is the flat methyl-shift region found on 2026-11-23. It is the strong end of σ(C–C) hyperconjugation, not ordinary
  hyperconjugation.
- **Answer to the open question (does hyperconjugation weaken near the TS?): on this path, yes.** Both orbitals lose
  almost all of their share on C2 before the TS: from 0.36 and 0.17 to 0.03 and 0.04. Half of each loss is done by
  s ≈ −7, well before the TS at s = 0 and long before the hydrogen moves. Over the same stretch the energy rises by
  16.0 kJ/mol (3.8 kcal/mol), from the reactant end to the TS. No other orbital visibly takes over until the hydrogen
  itself moves: the migrating C–H has only 0.05–0.08 on C2 up to the TS.
- **Limits:**
  - This path's TS is the C2–C3 rotation TS; the hydride moves downhill after it. So "near the TS" here means "during
    the rotation that brings the C3–H into line with the empty orbital".
  - These are Mulliken shares, with the caveat above.

## 3. The new cationic carbon C3 (product end)

At the product end the positive charge is on C3. The orbitals next to it are the bonds of the atoms attached to C3:
C2, C4 and the CH₃ on C3. Their shares on C3, ordered by gain:

| Orbital | Share on C3, reactant end | At the TS | At the product end | Largest anywhere |
|---|---|---|---|---|
| σ(C–H) of the CH₃ on C3 (one of three) | −0.030 | 0.115 | **0.110** | 0.158 (s −3.81) |
| σ(C1–C2) | −0.063 | −0.030 | −0.004 | −0.002 |
| σ(C4–H) (one of three) | 0.079 | 0.061 | 0.085 | 0.086 |
| σ(C–H) of the CH₃ on C3, second | 0.051 | 0.002 | 0.039 | 0.059 |
| σ(C2–H) (the H that was on C2 from the start) | 0.103 | 0.083 | 0.080 | 0.118 |
| σ(C–H) of the CH₃ on C3, third | 0.052 | 0.005 | 0.029 | 0.052 |
| σ(C4–H), second | 0.082 | 0.112 | 0.049 | 0.210 (s −5.76) |
| σ(C4–H), third | 0.073 | −0.003 | 0.034 | 0.082 |
| σ(C2–H) of the migrated H (the expected orbital) | 0.534 (as σ(C3–H)) | 0.583 | 0.076 | — |

- **No orbital next to C3 crosses 0.15 at the product end.** The largest share there is 0.110. Only one orbital gains
  noticeably (+0.14). The next, σ(C1–C2), goes from −0.063 to −0.004 (+0.06), and the rest gain less than 0.01 or
  fall.
- **The same kind of share appears without any charge.** At the reactant end, C–H orbitals of C4 have 0.07–0.08 on
  C3 while C3 is not cationic, and mid-path one reaches 0.21. So in this basis, Mulliken shares do not separate
  hyperconjugation into C3 from ordinary sharing with a neighbouring carbon. In an ordinary open tertiary cation, the
  hyperconjugation from each bond is too small for this measure. The large shares at the reactant end come from its
  partly bridged structure.

## 4. Tertiary shift (2,3-dimethyl-2-butyl cation): the same check

- **Written from the mechanism: yes.** The expected list, σ(C3–H) → σ(C2–H) (three-centre at the TS), was added on
  2026-11-10 in commit 673e45a ("inputs prepared (no calculations yet)"), before any calculation.
- **Delocalized orbitals not listed: none.** Tracking finds no delocalized orbital at either end. No orbital other than
  the expected one changes its carrying atoms at any of the 543 frames, computed or mirror.
- **Closest approach to the threshold:** 0.147, which is a σ(C–H) of a CH₃ on C3 with its share on C3 at the reactant
  end, where C3 is not cationic.
- **Methyl C–H shares on the cationic carbon C2:**
  - the three largest at the reactant end are 0.087, 0.110 and 0.121;
  - at the TS the C2-methyl C–H shares on C2 are 0.037–0.066;
  - the product end mirrors the reactant end on C3.

  So they about halve at the TS, the direction Saunders and Kates infer. But the same kind of share on the
  non-cationic C3 reaches 0.125 and 0.147 at the reactant end. By this measure, a weaker hyperconjugation at the
  tertiary TS is therefore not resolved.
