# Questions from three papers (Knizia 2013, Shiner 1952, Farcasiu and Norton 1997), 2026-11-27

J. Picione's five questions of 2026-11-27. Answered from files on disk and the ORCA 6.1 manual; no ORCA job and no new
geometry. Energies in kJ/mol (kcal/mol).

## 1. IBO localization exponent

**What the inputs, output and manual say.** Every IBO job uses `%loc LocMet IAOIBO end` and no other %loc option
(e.g. `reactions/sn2_hydroxide_bromomethane_m062x/03_orbitals/frames/f000/f.inp`).
- **ORCA output:** it prints the method but not the exponent:
  - `Localization criterion                   ... IAO-IBO`
  - `Your calculation utilizes the IAOIBO localization method` (citing Knizia, JCTC 2013, 9, 4834)
- **ORCA 6.1 manual** (section 9.2.5, orca_loc): the IBO step is described as "(ii) The PM functional is used to
  localize these IAOs into IBOs". The %loc block has no exponent keyword.

**What the printed numbers show.** The output does print the converged localization sum, and that sum tells the
exponents apart.
- **Printed sum:** `ITERATION   2 : L=   11.8237569491`, SN2 frame 0, 16 orbitals.
- **The test:** a two-centre bond orbital with charges near ½ / ½ adds about 0.5 to the sum with exponent 2, but only
  about 0.13 with exponent 4.
- **Estimate from the printed compositions** (11 orbitals on one atom; five bonds with populations 0.43/0.59 three
  times, 0.35/0.65 and 0.60/0.40):
  - exponent 4 gives about 11.9;
  - exponent 2 gives about 13.7.

The printed 11.82 matches exponent 4, Knizia's functional. Mulliken populations are not IAO charges, but the gap
between the two estimates (about 1.9) is far larger than that difference.

**Reading:** ORCA 6.1.1 ran the exponent-4 IBO functional, although its manual calls the step "the PM functional".

## 2. σ/π display rotation

| Reaction | Display rotation |
|---|---|
| SN2 | none |
| E2 | none |
| Cyanide + acetone | none |
| Tertiary hydride shift | none |
| Secondary hydride shift | none |
| Diels–Alder | yes: tracked IBO 6, π(C5=C6) → σ(C1–C6), IRC frames 0–139 |

Only `reactions/da_butadiene_ethylene_m062x/03_orbitals/display_rotation.json` exists among the production folders,
set by reaction.toml `[display_rotation]` (atoms 4, 5; max_angle_deg 0.5).

**Diels–Alder details:**
- **What is rotated:** the displayed orbital is the least-s-character combination of IBO 6 and its partner IBO 5, the
  other C5=C6 orbital.
- **Frames:** IRC frames 0–139 (s from −8.94 to −2.00 bohr·amu½; TS at frame 179). From frame 140 on, as computed;
  overlap across the join 0.99998.
- **Rotation angle along the path:**
  - 86.8° at frame 0, 74.2° at frame 70, 39.8° at frame 90, 11.2° at frame 110 and 0.51° at frame 139;
  - s-character of IBO 6 is 0.43 at frame 0, 0.17 at frame 90 and −0.02 at frame 139.
- **What the computed pair is:**
  - At the reactant end, the computed pair is already σ + π (s-character 0.43 and 0.00), with the tracked orbital
    the σ.
  - The mixed (bent-bond) character appears only on the way, at about frames 70–110.
  - So on the production path the rotation mostly shows the π partner instead of the tracked σ, rather than splitting
    two bent bonds.
- **Protocol text describes the pilot:** docs/protocol.md (display-rotation paragraph) gives the pilot path's numbers
  (two bent bonds with s-character 0.18 and 0.24, frames 0–47), not the production path's.

## 3. Secondary hydride shift: ts_anti_rotation

**How the two sides were followed:** not an IRC.
- **Procedure:** `00_inputs/rotation_ts_m062x.py`, `STEP_A = 0.10`, function `displaced`.
  - The TS (OptTS Freq, one imaginary frequency, 162.4i cm⁻¹) was displaced both ways along its imaginary normal mode
    from `ts.hess`, scaled so that the largest atom displacement is 0.10 Å.
  - Each side was then optimized with Opt Freq.
- **Result** (`01_stationary/rotation_ts.md`):
  - The plus side reached secondary_anti (RMSD 0.001 Å, ΔE +0.00 kJ/mol).
  - The minus side reached tertiary_irc (RMSD 0.007 Å, ΔE −0.11 kJ/mol).
- **The match test:**
  - core atoms within 0.10 Å RMSD and energies within 1.0 kJ/mol of the reference minimum;
  - the mirror image and the two C3 methyls swapped are allowed.

**Geometry** (`01_stationary/ts_anti_rotation/ts/ts.xyz`; H on C2 = atom 2, shifting H = atom 9):

| Structure | H–C2–C3–H | C1–C2–C3–H | C3–H | C2–H |
|---|---|---|---|---|
| ts_anti_rotation | 127.2° | −55.4° | 1.113 Å | 2.020 Å |
| TS of the animated path (`01_stationary/ts/ts.xyz`), for comparison | 53.6° | −122.7° | 1.113 Å | 2.008 Å |

## 4. Secondary hydride shift: reverse step

**What is computed:** from the tertiary cation to the TS of the animated path, 298.15 K.

**Zero:** `tertiary_irc`, the minimum the path's product end reaches.
- Source: `02_irc/endpoints.json` product_complex; the Opt Freq is in `02_irc/true_irc/end_check/forward/opt.out`.
- A separate tight optimization of the tertiary cation (`01_stationary/tertiary/tertiary_tight.out`) lies 0.18
  kJ/mol lower in E and 0.46 kJ/mol higher in G at the working level.

**TS:** `01_stationary/ts/ts.out`, 91.6i cm⁻¹.

**Thermal terms:** from the M06-2X-D3(0)/ma-def2-SVP frequencies, SMD dichloromethane.
- **Enthalpy:** ORCA's H = E + ZPE + thermal vibrational, rotational and translational energy + kT.
  - Vibrational energy is from the harmonic oscillator.
  - ORCA leaves the TS's imaginary mode out of the thermal sums.
- **Free energy:** G = H − TS, with the vibrational entropy from Grimme's quasi-RRHO ("Quasi RRHO True" in the
  output).
- **Standard state:** one molecule on both sides, so no standard-state term.
- **Corrections at the other levels:** E at each level + the working-level (H − E) and (G − E).
- **Size of the thermal differences, TS − tertiary_irc:**
  - (H − E): −1.15 kJ/mol;
  - (G − E): +5.09 kJ/mol.

| Level | ΔE‡ | ΔH‡ | ΔG‡ |
|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (working) | 58.4 (14.0) | 57.2 (13.7) | 63.5 (15.2) |
| M06-2X/ma-def2-TZVP | 57.1 (13.6) | 55.9 (13.4) | 62.1 (14.9) |
| DLPNO-CCSD(T)/ma-def2-TZVP | 54.4 (13.0) | 53.2 (12.7) | 59.5 (14.2) |
| B3LYP/ma-def2-TZVP | 58.7 (14.0) | 57.6 (13.8) | 63.8 (15.2) |

**Check:** the same terms give the forward ΔG‡ from the eclipsed form at DLPNO-CCSD(T) as 11.13 kJ/mol, the card's
11.1.

## 5. Tertiary hydride shift: the 108 forward frames of the angle check

**Which frames they are:** computed forward IRC steps 21–128 from the TS (`02_irc/true_irc/forward/path.xyz`, 129
steps).

**What the computed forward side does:**
- **Steps 1–7:** it reaches the C₂-symmetric bridged structure at step 7 (RMSD 0.0036 Å to bridged_c2, 0.015 kJ/mol
  below the TS).
- **From step 8 on:** it continues past the bridged structure, back down toward an open cation with H on C3.
  - step 21: C2–H 1.43 Å, C3–H 1.24 Å, 0.15 kJ/mol below the TS;
  - step 51: 1.81 / 1.14 Å, −5.9 kJ/mol;
  - step 129: 1.93 / 1.125 Å, −11.4 kJ/mol, RMSD 0.43 Å from bridged_c2.

**What the displayed path uses:** `pipeline/mirror_path.py` cuts the side at the frame closest to the bridged
reference.
- The 543-frame path uses forward steps 1–7 only (frames 265–271), then the mirror half.
- Steps 8–129 are kept in `path.xyz` but are not IRC points of the displayed path.

**What the angle-check numbers count:** both forward zones were computed from all 128 interior frames of the computed
side.
- The 108 frames beyond |s| = 1 (median 4.27°, 90th percentile 8.38°, passes) all lie past the bridged structure, on
  the stretch the displayed path does not contain.
- **Near zone:**
  - The near zone's 20 frames are steps 1–20. Steps 1–6 are on the displayed path; their angles use only steps 0–7.
  - Their angles: 75.4°, 54.6°, 43.7°, 36.7°, 31.5° and 26.2°; median 40.2°, 90th percentile 65.0°.
  - Step 7's angle (10.9°) uses step 8, which is not on the displayed path.

**End check:** the forward end check (bridged form, −0.003 kJ/mol, 0.004 Å) is for step 7, the last frame kept.
