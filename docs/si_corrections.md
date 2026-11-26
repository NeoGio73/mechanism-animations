# Corrections list for the Supporting Information

Problems found in the computational workflow, what they affected, and how they were fixed. Each entry names
the files that hold the withdrawn output (kept, labelled not results) and the corrected output. Started
2026-11-15 at J. Picione's request.

## 1. Early stop of an IRC side at a flat TS (secondary hydride shift, 3-methyl-2-butyl cation)

- **Found:** 2026-11-15, at the join of the production path (M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane).
- **What happened:** the TS is a flat C2–C3 rotation TS (92i cm⁻¹). On the first 20 steps from it, each step
  changed the energy by less than 0.05 kJ/mol, and no key bond changed. The project stopping rule (20
  consecutive steps below 0.05 kJ/mol after the last bond-forming or bond-breaking event) was therefore met at
  step 20 on both sides, while the path was still leaving the TS: 0.32 kJ/mol (forward) and 0.39 kJ/mol
  (backward) below the TS (the last stored points, one step further, 0.36 and 0.42 kJ/mol), about 15.4 kJ/mol
  above the secondary cation.
- **Effect:** the joined path covered only the TS region. It was not presented as the hydride shift path.
  Optimizations from the two ends reached the tertiary cation (forward) and the eclipsed secondary cation
  (backward), so the TS connects the expected minima.
- **Fix (J. Picione, 2026-11-15):** the 20 flat steps count only once the energy is at least 1 kJ/mol below the
  TS; a side whose minimum lies less than 1 kJ/mol below its TS stops when the end-minimum test (|ΔE| < 0.005
  kJ/mol per step, RMS gradient < 5×10⁻⁵ Eh/bohr) holds on 40 consecutive frames (a single hit occurs right next
  to a flat TS). Last backstop: 3000 steps per side. Both sides are rerun from the TS in clean folders, after
  E2. The first run is kept as not results in `02_irc/attempts/first_run_early_stop_flat_ts/`.
- **Other paths:** for every other finished side (SN2, cyanide + acetone, Diels–Alder, both tertiary hydride
  shift sides) the 20 flat steps began at least 10 kJ/mol below the TS, so the new rule stops them at the same
  frame; none is rerun. The Diels–Alder product side ended by the energy-rise stop and was accepted by
  J. Picione (2026-11-14); the rule change does not apply to it.
- **Files:** `reactions/hshift_3_methyl_2_butyl_m062x/02_irc/attempts/first_run_early_stop_flat_ts/` (first run,
  WHY.txt), `run_log.md` there; rule replay test `tests/test_stopping_rule.py`.

## 2. Join that picked the wrong bridged side (tertiary hydride shift, 2,3-dimethyl-2-butyl cation)

- **Found:** 2026-11-15, at the first production join (12:13).
- **What happened:** the path is built as the computed half (open cation → TS → bridged structure) plus its
  mirror image under the C2 ⇄ C3 swap. The code chose the side leading to the bridged structure by the RMSD of
  each side's last frame to the optimized bridged structure (`bridged_c2`). On the noise-limited plateau (the
  bridged structure lies 0.012 kJ/mol below the TS), the forward side reached the bridged structure (frame 6,
  RMSD 0.004 Å) and then slid back to the open cation it started from. Both sides ended at the same open
  cation (last-frame RMSD 0.424 Å and 0.427 Å), and the backward side was chosen by 0.003 Å.
- **Effect:** the mirror path was built from the wrong half, and the bridged-side check ran on the frame one
  step from the TS, which failed the key C–H limit (+0.064 / −0.051 Å against 0.05 Å). No IRC calculation
  was affected; only the join (path assembly, angle check, end checks).
- **Fix:** the bridged side is the side that comes closest to `bridged_c2` at any frame
  (`pipeline/mirror_path.py`). The join was redone with the same settings at 12:18: bridged end 0.004 Å from
  `bridged_c2`, end check accepted; path 241 frames (121 computed, 120 mirror).
- **Files:** withdrawn output in `reactions/hshift_2_3_dimethyl_2_butyl_m062x/02_irc/attempts/join_wrong_bridged_side_1213/`
  (WHY.txt); corrected output in `02_irc/` there.

## 3. E2 end structures: the minima the path reaches replace the hand-built complexes (E2, ethoxide + 2-bromopropane)

- **Found:** 2026-11-16, at the end check of the production path (join 03:39; M06-2X-D3(0)/ma-def2-SVP, SMD ethanol).
- **What happened:** the reactant and product complexes had been optimized beforehand from the B3LYP-level IRC end
  points (hand-built references). On each side, Opt Freq from the last frame of the M06-2X path reached a slightly
  different form of the same complex: covalent key distances within 0.001 Å; one non-covalent contact different
  (reactant side ethoxide O···Hβ 2.623 vs 2.688 Å; product side Cβ···H(ethanol) 2.544 vs 2.443 Å); heavy-atom RMSD
  0.226 Å and 0.150 Å; E lower by 0.267 and 0.279 kJ/mol; no imaginary frequency.
- **Decision (J. Picione, 2026-11-16):** the path is accepted with a note; the minima the path actually reaches are
  the official end structures, because the path defines which minima the TS connects. The references are kept on
  record with the comparison.
- **Effect:** the energy zero and every quantity measured from the complexes use the new structures. Electronic
  barrier from the reactant complex 59.67 instead of 59.40 kJ/mol (+0.27); ΔG‡ from the reactant complex
  48.40 instead of 51.06 kJ/mol (−2.66; the two forms differ in their softest vibrations, below 25 cm⁻¹); ΔG from
  complex to complex −108.61 instead of −105.64 kJ/mol. None of the old values had been reported, carded or
  tabulated. Barriers from separated reactants (comparison with experiment) are unaffected.
- **Files:** `reactions/e2_ethoxide_2_bromopropane_m062x/02_irc/true_irc/end_check_decision.md`;
  new structures `02_irc/true_irc/end_check/{backward,forward}/opt`; references `01_stationary/{reactant,product}_complex`;
  `docs/stationary_m062x.json` (`superseded_by`, `official_end`); `02_irc/endpoints.json`.


## 4. Every path: the minimum the path reaches is the official end structure (all six reactions)

- **Decision (J. Picione, 2026-11-16):** extends entry 3 to every path, the finished ones included, so all six are
  treated alike. Every side is optimized from its last frame (Opt Freq as its reference minimum); where that reaches
  a different form of the reference with the same bonds, the reached form replaces it, and every quantity measured
  from it is redone. Different bonds: nothing is switched and the case goes to J. Picione.
- **Switches so far** (filled in as the optimizations finish; references kept, marked `superseded_by`):
  - SN2, reactant complex (2026-11-16): the reached form has C···O 3.309 Å instead of 3.425 Å (C–Br within
    0.001 Å, heavy-atom RMSD 0.055 Å, same bonds), E +0.036 kJ/mol, G +0.589 kJ/mol relative to the reference.
    Single points and diagnostics on it redone in clean folders (old ones in `01_ts/superseded_20261116/`).
    **Corrects values reported on 2026-11-14 (17:00 report)**, ΔG‡ from the reactant complex: M06-2X/ma-def2-TZVP
    78.35 → 77.96 kJ/mol; DLPNO-CCSD(T)/ma-def2-TZVP 83.99 → 82.83 kJ/mol; B3LYP/ma-def2-TZVP 62.50 → 61.31 kJ/mol;
    working level 75.47 → 74.88 kJ/mol. Electronic barriers from the complex: working 67.40 → 67.36, M06-2X/TZVP
    70.28 → 70.44, DLPNO 75.92 → 75.31, B3LYP 54.43 → 53.79 kJ/mol. Barriers from separated reactants and reaction
    energies unchanged.
  - Diels–Alder, reactant complex (2026-11-16): the reached form (Opt, then TightOpt after an imaginary mode) has the
    C4···C5 contact 0.185 Å longer, heavy-atom RMSD 0.250 Å, same bonds; E −0.297 kJ/mol, G +2.245 kJ/mol relative to
    the reference (lowest frequency 37.8 cm⁻¹). Electronic barrier from the complex 75.34 → 75.63 kJ/mol, ΔG‡ from the
    complex 98.40 → 96.16 kJ/mol, ΔG complex to complex −158.04 → −160.28 kJ/mol (working level). No Diels–Alder
    after-path numbers existed yet; its single points start on the new structure. Product complex: same minimum.
  - Cyanide product side (2026-11-16): same minimum (E −0.031 kJ/mol, heavy-atom RMSD 0.003 Å); no change.
  - E2, both complexes: entry 3.
- **Free energies** (checked 2026-11-16): every free energy uses ORCA's quasi-RRHO vibrational entropy (Grimme 2012,
  ω₀ = 100 cm⁻¹), recomputed independently for all 21 stationary points (`pipeline/thermo.py`, within 0.03 kJ/mol);
  plain harmonic entropy would change G − E by up to 13.2 kJ/mol. No change needed.
- **Files:** `02_irc/true_irc/end_check_all.{json,md}` per reaction; `docs/stationary_m062x.json`;
  `02_irc/endpoints.json`; `pipeline/path_minimum.py`.

## 5. Stopping rule near the end minimum; secondary and tertiary shift sides rerun

- **Found:** 2026-11-16 05:46. The secondary shift's rerun backward side (rule of 2026-11-15) stopped 3.71 kJ/mol
  below its flat TS and 12.15 kJ/mol (2.90 kcal/mol) above the secondary cation, while the energy still fell about
  0.035 kJ/mol per step. The tertiary shift's open side had stopped 4.1 kJ/mol above the open cation (accepted with
  a note on 2026-11-15).
- **Decision (J. Picione, 2026-11-16):** on these sides the 20 flat steps count only within 0.5 kJ/mol of the end
  minimum's energy; end-minimum test (40 frames) and 3000-step backstop unchanged. Both sides rerun from the TS in
  clean folders (secondary first, then tertiary), both joins redone. The earlier sides are kept as not results in
  `02_irc/attempts/backward_stopped_short_0546` (secondary) and `02_irc/attempts/open_side_stopped_short_1218`
  (tertiary). (Changed the same day after the continuation test, below: both sides are continued from their last frame
  instead, so the earlier frames are the first part of each side: secondary `02_irc/true_irc/backward_part1`, the
  earlier join's records in `02_irc/attempts/backward_stopped_short_0546`.) The tertiary half-step check covers only the region near the TS and stays valid (checked 2026-11-16:
  heavy-atom RMSD ≤ 0.0007 Å against the production side up to s = 2.0 bohr·amu½).
- **Pitfall (for the SI): the stopping rule fired early on a slow twist.** The secondary shift's backward side leaves
  a flat rotation TS (92i cm⁻¹) along the C2–C3 twist. Along that twist the energy falls by only about 0.035 kJ/mol
  per 0.05 bohr·amu½ step, below the 0.05 kJ/mol flat-step threshold, so a rule that counts "flat" steps can end the
  side in the middle of the twist, with most of the descent (12.15 kJ/mol, 2.90 kcal/mol) still missing. Flat steps
  must therefore not be read as "near a minimum" on a side that crosses a soft torsion; for this side they count only
  within 0.5 kJ/mol of the secondary cation's energy (rule (c)), and the side is extended to it.
- **Revised the same afternoon (J. Picione, 2026-11-16):** a test showed that the final approach to a minimum along
  soft modes is far slower than expected (the tertiary open side lost about 0.01 kJ/mol per step 3 kJ/mol above
  the open cation) and changes nothing visible. The target of ending every side within 0.5 kJ/mol of its minimum is
  dropped; standard practice applies: each end is confirmed by optimization from its last frame, and the table of side
  ends (below) gives how far above its minimum each side ended. Only the secondary backward side is extended (above).
  The tertiary open side ends where the continuation test was stopped (below).
- **Continuation test result (2026-11-16 16:13):** the rerun from the TS was stopped once it had 54 frames past the old
  side's last frame, and the continuation at the same moment. At equal arc length the 54 compared frames agree to
  heavy-atom RMSD ≤ 0.000013 Å (median 0.000004 Å) and |ΔE| ≤ 0.00007 kJ/mol, against criteria of 0.05 Å and
  0.1 kJ/mol fixed before either run finished; the 112 earlier frames are identical. Continuing a side downhill from
  its last frame (exact Hessian there) is therefore equivalent to rerunning it from the TS. The tertiary open side is
  the earlier side's 113 frames plus 151 continued frames (264 frames, s = 13.29 bohr·amu½), ending 1.47 kJ/mol
  (0.35 kcal/mol) above the open cation; the rerun's frames are kept as the test record
  (`02_irc/attempts/rerun_from_ts_test_1614/`, comparison in `02_irc/attempts/continuation_test/compare_backward.md`).
  The secondary backward side is continued the same way from its last frame (started 2026-11-16 16:17). It ended by
  rule (c) at 19:01 (104 new frames, 189 in all), 0.0004 kJ/mol above the form of the secondary cation it reaches;
  joined 19:07, end check ok.
- **Angle check on the continued frames:** in the slow final approach the mass-weighted angle between step and
  gradient grows (frames past the earlier end: median 6.83°, 90th percentile 10.92°, 150 frames; earlier frames
  beyond |s| = 1: 4.41° / 7.81°). The rerun from the TS has the same frames, so this belongs to the region, not to the
  continuation. The open side's zone beyond |s| = 1 now has median 5.74° (rule < 5°), 90th percentile 10.18° (rule
  < 15°); kept as computed, both parts stated (J. Picione, 2026-11-16 evening). All 151 continued frames are kept
  (every frame computed by the time of the comparison, not only the 54 that overlap the rerun; J. Picione, same
  evening); cutting at the overlap would have ended the side about 3.1 kJ/mol (0.75 kcal/mol) above the open cation.
- **Angle check, secondary hydride shift:** within |s| ≤ 1 bohr·amu½ it fails on both sides (backward median 7.63°,
  90th percentile 25.37°; forward 10.55° / 25.37°; rule < 5° / < 15°); beyond it, it passes (backward 2.37° / 6.17°,
  forward 0.98° / 7.50°). Same pattern as the tertiary path: near the TS the gradient is nearly zero (entry 6). Kept
  as computed, both parts stated, no rerun (J. Picione, 2026-11-17).
- **Side ends** (last frame above the official end minimum, working level; filled in as sides finish):

  | Path | Reactant side | Product side |
  |---|---|---|
  | SN2 | 0.06 kJ/mol | 0.23 kJ/mol |
  | Diels–Alder | 0.33 kJ/mol | 0.01 kJ/mol |
  | Cyanide + acetone | 0.94 kJ/mol | −0.03 kJ/mol |
  | E2 | 1.17 kJ/mol | 5.01 kJ/mol |
  | Tertiary hydride shift | open side: 1.47 kJ/mol (end check ok) | bridged side: at the bridged structure |
  | Secondary hydride shift | secondary cation: 0.0004 kJ/mol | tertiary cation: 0.004 kJ/mol |

  Both secondary-shift ends optimize to a slightly different form of the hand-built cation with the same bonds
  (backward: E 0.10 kJ/mol lower, heavy-atom RMSD 0.055 Å; forward: E 0.18 kJ/mol higher, C3···H 0.091 Å shorter);
  by the rule of entry 4 these are the official end structures (`secondary_eclip_irc`, `tertiary_irc`). Working-level
  ΔG‡ from the secondary cation 12.59 → 13.00 kJ/mol (3.01 → 3.11 kcal/mol); ΔG_r −50.45 → −50.49 kJ/mol
  (−12.06 → −12.07 kcal/mol). No after-path value existed for this reaction.

## 6. Angle check against gradient size (all six paths; for the SI figure)

- **Requested:** J. Picione, 2026-11-16 evening, after the tertiary open side's angle check (entry 5). Existing data only,
  no new calculation: `pipeline/angle_vs_gradient.py`; figure `docs/si_figures/angle_vs_gradient.png`, per-frame values
  `docs/si_figures/angle_vs_gradient.json`.
- **What is plotted:** for every computed frame of every production side (2212 frames; mirror frames of the tertiary
  shift have no jobs and are left out), the mass-weighted angle between the path direction and the downhill gradient
  (the angle of the angle check) against the RMS Cartesian gradient at that frame (Eh/bohr, log scale). Frames are marked
  by zone: near the TS (|s| ≤ 1 bohr·amu½, as in the angle check); final approach (from the frame after which the RMS
  gradient stays below 5 × 10⁻⁴ Eh/bohr to the side's end, ten times the end-minimum test's 5 × 10⁻⁵); middle (the rest).
  Redrawn 2026-11-16 20:10 with the secondary shift's finished backward side (189 frames).
- **Result: the large angles occur only where the gradient is nearly zero.** All 70 frames with an angle above 15° have
  an RMS gradient of at most 2.0 × 10⁻⁴ Eh/bohr: 30 lie near the TS of the two hydride shifts, whose TSs are very flat
  (tertiary 21, secondary 9), and 40 in the slow final approach to a minimum (cyanide product side 18, secondary 12,
  Diels–Alder 5, tertiary open side 3, SN2 2). None is in the middle of a side. Where the RMS gradient exceeds
  5 × 10⁻⁴ Eh/bohr (1241 frames) the largest angle is 5.65° (median 0.10°), and above 1 × 10⁻³ Eh/bohr (927 frames)
  4.07° (median 0.04°). Near the TS of the four reactions with sharp TSs (SN2, cyanide, Diels–Alder, E2) the angles
  stay small (largest 5.65°). The largest angles (up to 117°) are on the cyanide product side, which ends at its
  minimum (−0.03 kJ/mol), where the gradient is about 1.5 × 10⁻⁵ Eh/bohr and its direction carries little information.
- **Reading for the paper:** a large angle in these paths marks a region where the energy surface is nearly flat (a very
  flat TS, or the last approach to a minimum along soft modes), not a path that leaves the valley: wherever the gradient
  is appreciable, path and gradient agree within a few degrees. The angle-check results are kept as computed (J. Picione,
  2026-11-15 and 2026-11-16).

## 7. Tertiary hydride shift: mirror frames take their partner's energy at B3LYP and DLPNO-CCSD(T) (J. Picione, 2026-11-17)

- **What:** frames 272–542 of the tertiary path (2,3-dimethyl-2-butyl cation) are the C2 ⇄ C3 mirror images of the
  computed frames 270 … 0 (frame m is the image of frame 542 − m). Their energies equal their partners' at every level
  by symmetry.
- **Check:** M06-2X/ma-def2-TZVP, the first path level, computes all 543 frames, mirror frames included, and each
  mirror frame is compared with its partner (`02_irc/path_sp/mirror_check.{json,md}`). Tolerance: the SCF noise floor,
  0.034 kJ/mol. Result (2026-11-17 20:52): **pass**. All 271 mirror frames agree with their partners; the largest
  |ΔE| is 1.40 × 10⁻³ kJ/mol (3.3 × 10⁻⁴ kcal/mol), at frame 460 against frame 82, and the median is
  3.2 × 10⁻⁵ kJ/mol (7.7 × 10⁻⁶ kcal/mol); no frame is over the tolerance.
- **If every mirror frame agrees:** B3LYP/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP copy the partner's energy for
  each mirror frame (about 271 and 244 frames) instead of computing it. Copied values are marked in the data files
  (`computed` = 0, `copied_from_frame` in `02_irc/energies_<level>.csv`; `copied_frames` in `energies_<level>.json`).
  This saves about 270 B3LYP and 245 DLPNO-CCSD(T) jobs.
- **If any differs:** every frame is computed at every level, and J. Picione is told.
- **SI wording (draft):** "For the degenerate tertiary shift the second half of the path is the mirror image of the
  first. Single points at M06-2X/ma-def2-TZVP were computed on every frame; mirror frames agreed with their partners
  to within [max] kJ/mol, so B3LYP and DLPNO-CCSD(T) energies of mirror frames were taken from their partners."

## 8. SN2 experimental barrier: 95.0 → 94.9 kJ/mol from the original measurements (J. Picione, 2026-11-17)

- **What:** the measured ΔG‡ for HO⁻ + CH₃Br in water was 95.0 kJ/mol (22.7 kcal/mol), from
  1.41 × 10⁻⁴ L mol⁻¹ s⁻¹ (Mabey and Mill 1978, Table 5.1), calculated from the 1959 Arrhenius parameters at 298.0 K;
  our earlier records treated it as a 298.15 K value, giving 95.0 kJ/mol.
  Check: A = 1.04 × 10¹³ L mol⁻¹ s⁻¹ and Ea = 23.0 kcal/mol give 1.4107 × 10⁻⁴ L mol⁻¹ s⁻¹ at 298.0 K; Eyring with
  1.41 × 10⁻⁴ gives 95.00 kJ/mol taken at 298.15 K (94.95 kJ/mol at 298.0 K). Neither 1959 table prints a 298 K rate
  constant (Fells and Moelwyn-Hughes Table 11: k₂ only at 373.16 K, 3.52 × 10⁻¹ L mol⁻¹ s⁻¹, with A and Ea;
  Bathgate and Moelwyn-Hughes Table 12: log A = 13.0 and Ea only). J. Picione checked it against the original measurements.
- **Source now:** E. A. Moelwyn-Hughes, Proc. R. Soc. Lond. A 1949, 196, 540–553, doi:10.1098/rspa.1949.0044 (NaOH in
  water; rates corrected for the competing reaction with water). Table 5, k₂ (10⁻⁴ L mol⁻¹ s⁻¹): 292.69 K, 0.708;
  297.62 K, 1.37; 303.41 K, 2.87; 308.46 K, 5.33; 313.29 K, 9.31; 322.91 K, 29.1. Table 9: k₂(298.1 K) =
  1.435 × 10⁻⁴ L mol⁻¹ s⁻¹, A = 1.04 × 10¹³ L mol⁻¹ s⁻¹, Ea = 96.2 ± 0.6 kJ/mol (23,000 ± 150 cal/mol). The same values
  are summarized in Fells and Moelwyn-Hughes, J. Chem. Soc. 1959, 398 (Table 11, p. 408) and Bathgate and
  Moelwyn-Hughes, J. Chem. Soc. 1959, 2642 (Table 12, p. 2646).
- **Fit:** Eyring, least squares of ln(k₂/T) against 1/T over the six Table 5 points (`pipeline/eyring_fit.py`):
  ΔH‡ = 93.7 ± 0.5 kJ/mol (22.4 kcal/mol), ΔS‡ = −4.1 ± 1.5 J mol⁻¹ K⁻¹, ΔG‡(298.15 K, 1 M) = 94.9 kJ/mol
  (22.7 kcal/mol). 298 K lies inside the measured range, so no extrapolation is involved.
- **Effect:** every computed − measured SN2 difference rises by 0.1 kJ/mol. Def2-SVP pilot card, DLPNO-CCSD(T)/ma-def2-TZVP:
  +0.6 → +0.7 kJ/mol from the TS frame, +4.7 → +4.8 kJ/mol from the curve maximum. M06-2X sensitivity test:
  +5.6 → +5.7 kJ/mol (B3LYP geometry), +5.3 → +5.4 kJ/mol (M06-2X geometry). No calculation changes.
- **Files:** `reactions/sn2_hydroxide_bromomethane{,_masvp,_m062x}/00_inputs/card_notes.toml`,
  `reactions/sn2_hydroxide_bromomethane/card.md`, `pipeline/level_test.py`, `docs/sensitivity_M06-2X.{json,md}`,
  `docs/protocol.md`, `docs/references.md`. Reports written before 2026-11-17 keep 95.0 as written.

## 9. Reference conformers: lowest free energy for every separated species (J. Picione, 2026-11-17)

- **What:** every separated reactant and product now uses its lowest-energy conformer (working-level free energy) as
  the reference; the path and the animation keep the form the path reaches. Before, the Diels–Alder zero was gauche
  butadiene (C=C–C=C 32°, the form the path reaches) + ethylene.
- **Diels–Alder:** s-trans butadiene is 9.48 kJ/mol (2.27 kcal/mol) lower in G, so it is now the zero. Barriers from
  separated reactants rise by that gap; M06-2X/ma-def2-TZVP ΔG‡(298 K, 1 M) 121.88 → 131.89 kJ/mol (29.13 →
  31.52 kcal/mol), DLPNO-CCSD(T)/ma-def2-TZVP 134.18 → 144.36 kJ/mol (32.07 → 34.50 kcal/mol). Path profiles move up
  by the electronic gap at each level (M06-2X 12.23, B3LYP 14.27, DLPNO-CCSD(T) 12.40 kJ/mol). Barriers from the
  reactant complex and the gas-phase comparison with experiment (already s-trans) are unchanged.
- **E2:** trans ethanol computed; 0.90 kJ/mol (0.22 kcal/mol) above gauche in G, so gauche stays; nothing changes.
- **Secondary hydride shift (J. Picione, 2026-11-19):** the same rule for the one-molecule reaction: the zero is the
  secondary-cation conformer lowest in G, not automatically frame 0. Frame 0 is the eclip form; the anti minimum,
  reoptimized at the working level (`01_stationary/secondary_anti`), is 2.77 kJ/mol (0.66 kcal/mol) lower in G and is the
  zero. Frame 0 lies above it in E by 3.06 / 2.41 / 2.32 / 3.16 kJ/mol (0.73 / 0.58 / 0.55 / 0.76 kcal/mol) at
  M06-2X/ma-def2-SVP / M06-2X/ma-def2-TZVP / B3LYP/ma-def2-TZVP / DLPNO-CCSD(T)/ma-def2-TZVP; path profiles move up by
  that much at each level. Barriers from the form the path reaches are unchanged and are kept. No path rerun.
- **Tertiary hydride shift, zero (J. Picione, 2026-11-19):** the same rule as for the secondary shift. The other
  C2–C3 rotamers of the tertiary cation were optimized at the working level (Opt Freq, `01_ts/fragments/rotamer_h_perp`,
  `rotamer_h_plane`) from two starts turned about C2–C3: C3–H along the empty p orbital (Me–C2–C3–H 90°) and C3–H in
  the cation plane (180°). Both relaxed, with no imaginary frequency, to the rotamer already computed
  (`dimethylbutyl_cation`: one C3 methyl 14° from the p axis, C3–H 62° from it; all-atom RMSD 0.018 and 0.012 Å, heavy
  atoms 0.005 Å). At this level the cation has one C2–C3 rotamer. The four optimizations of it differ in G by up to
  1.19 kJ/mol (0.28 kcal/mol) at 298.15 K and 0.94 kJ/mol (0.22 kcal/mol) at 135 K, from the lowest frequencies
  (48–72 cm⁻¹) of the same flat minimum, not from different conformers. The lowest in G is `dimethylbutyl_cation` at
  both 298.15 K and 135 K (the temperature of the Saunders–Kates measurement), so the zero is unchanged. Frame 0 lies
  above it in E by 1.27 kJ/mol (0.30 kcal/mol) at M06-2X/ma-def2-SVP, 1.86 (0.44) at M06-2X/ma-def2-TZVP, 1.67 (0.40)
  at B3LYP/ma-def2-TZVP and 2.22 (0.53) at DLPNO-CCSD(T)/ma-def2-TZVP. ΔG‡ from it, same order: 14.1 / 18.0 / 27.3 /
  17.0 kJ/mol at 298.15 K and 13.2 / 17.1 / 26.4 / 16.1 kJ/mol (3.2 / 4.1 / 6.3 / 3.8 kcal/mol) at 135 K; measured
  13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol) at 135 K. Script `00_inputs/rotamer_zero.py`, output `01_ts/rotamer_zero.md`.
  No path rerun.
- **Working-level zero chosen by G (fix, 2026-11-19):** since 2026-11-17 the other levels chose the zero of the path
  profile by working-level G, but the working-level profile (`pipeline/steps.profile`) still compared E. For the
  tertiary shift alone this made the working-level zero the path-end form (0.20 kJ/mol lower in E than
  `dimethylbutyl_cation`, 1.19 kJ/mol higher in G). Now G at every level. The working-level tertiary profile moved
  by −0.195 kJ/mol (TS frame 15.5 → 15.3 kJ/mol, 3.7 kcal/mol); Diels–Alder and cyanide zeros unchanged. Earlier
  reports that quote the working-level tertiary TS frame (15.5 kJ/mol) used the path-end form as zero.
- **Others:** single conformers, already compliant (details in `docs/protocol.md`).
- **Files:** `pipeline/steps.py` (`reference_names`), `reactions/{da_butadiene_ethylene,e2_ethoxide_2_bromopropane}_m062x/
  00_inputs/reaction.toml`, `00_inputs/ethanol_trans_start.xyz`; rewritten from existing jobs: `01_ts/diagnostics.json`,
  `02_irc/barriers.json`, `02_irc/energies*.csv/json` of both reactions. Reports before 2026-11-17 21:00 give the
  Diels–Alder values from gauche butadiene.

## 10. E2 experimental rate constant: 2.88 × 10⁻⁶ (uncorrected) → 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹ (J. Picione, 2026-11-13; rechecked 2026-11-20)

- **What:** from 2026-11-09 to 2026-11-13 the E2 records used k₂ = 2.88 × 10⁻⁶ L mol⁻¹ s⁻¹ as the overall
  second-order rate constant of 2-bromopropane + NaOEt in ethanol at 25 °C ("mean of 19 points, 2.76–2.96 × 10⁻⁶",
  read at first from a later reproduction of the tables). That is Table VII's 10⁶k₂ = 2.88, the uncorrected mean of a
  single run. With 80.3 % olefin it gave k(E2) ≈ 2.31 × 10⁻⁶ L mol⁻¹ s⁻¹ and ΔG‡(E2) ≈ 105.2 kJ/mol (25.1 kcal/mol),
  ΔG‡(SN2) ≈ 108.7 kJ/mol.
- **Source now:** M. L. Dhar, E. D. Hughes and C. K. Ingold, J. Chem. Soc. 1948, 2058–2065, doi:10.1039/JR9480002058;
  2-bromopropane + 1 N NaOEt in dry ethanol at 25.0 °C. Table I: 10⁵k₂ = 0.295 L mol⁻¹ s⁻¹ (corrected for solvent
  expansion), 80.3 % olefin (82.7, 78.0). Table II: 10⁵k(E2) = 0.237 L mol⁻¹ s⁻¹. Eyring (298.15 K, 1 M):
  k(E2) = 2.37 × 10⁻⁶ L mol⁻¹ s⁻¹ gives ΔG‡(E2) = 105.1 kJ/mol (25.1 kcal/mol); the SN2 share
  (0.58 × 10⁻⁶ L mol⁻¹ s⁻¹) gives ΔG‡(SN2) = 108.6 kJ/mol (26.0 kcal/mol). Checked against the original paper by
  J. Picione on 2026-11-13 and again, with the table details, on 2026-11-20.
- **Caveats stated in the paper:** first-order (solvolysis) processes contribute up to 2–3 % of the rate; olefin
  proportions are good to within 3 %, with a correction for loss of the volatile propene. A 3 % change in k(E2) moves
  ΔG‡ by 0.08 kJ/mol (0.02 kcal/mol).
- **Effect:** the measured E2 ΔG‡ falls by 0.1 kJ/mol (105.2 → 105.1); every computed − measured E2 difference rises
  by 0.1 kJ/mol. M06-2X sensitivity test (B3LYP geometries, M06-2X geometries): +2.0 → +2.1 and −5.4 → −5.3 kJ/mol.
  No calculation changes. Corrected on 2026-11-13 (commit `be06378`), before any production E2 number was compared.
- **Files that carried 2.88 × 10⁻⁶ or 105.2 kJ/mol:** `reactions/e2_ethoxide_2_bromopropane/00_inputs/card_notes.toml`,
  `reactions/e2_ethoxide_2_bromopropane/card.md`, `docs/references.md`, `pipeline/level_test.py`,
  `docs/sensitivity_M06-2X.{json,md}`. All now use 2.37 × 10⁻⁶ and 105.1 kJ/mol, as do
  `reactions/e2_ethoxide_2_bromopropane_m062x/00_inputs/card_notes.toml` and `summary_notes.toml` (from which the E2
  results summary takes its measured row). Reports written before 2026-11-13 keep the old value as written.


## 11. Diels–Alder gas-phase comparison: equilibrium butadiene conformer mixture as the reactant (J. Picione, 2026-11-21)

- **What:** the gas-phase comparison with Rowley and Steiner (760–921 K) of 2026-11-20 took s-trans butadiene alone
  as the reactant and gave Ea = ΔH‡ + 2RT at 840 K. At these temperatures about 30 % of the butadiene is in the two
  mirror-image gauche (s-cis) forms, and the measured value refers to that equilibrium mixture. The summary and card
  also said the measured value refers to butadiene that is "almost all s-trans", true at room temperature only.
- **Now:** gauche butadiene single points at the three single-point levels (05:06–05:09); conformer-weighted
  butadiene free energy G = −RT ln[exp(−G(s-trans)/RT) + 2 exp(−G(gauche)/RT)] (s-trans C2h σ = 2; gauche C2, σ = 2
  per enantiomer, × 2 for the two enantiomers); k(TST) at 760, 800, 840, 880 and 921 K and an ln k against 1/T fit,
  as in the experiment (`pipeline/gas_conf.py`, `reactions/da_butadiene_ethylene_m062x/03_gas_phase/conformers.json`).
- **Effect (DLPNO-CCSD(T)/ma-def2-TZVP):** Ea − measured +4.2 kJ/mol (s-trans, 840 K point value; +3.5 from the same
  fit) → −0.2 kJ/mol (−0.05 kcal/mol) with the mixture; A / measured A 0.45 → 0.17; k(840 K) / measured k
  0.25 → 0.17. M06-2X/ma-def2-TZVP: Ea −8.6 → −13.0 kJ/mol, k ratio 1.57 → 1.08. Gauche fraction 1.8 % at 298 K,
  30.4 % at 840 K.
- **Conclusion changed:** the coupled-cluster activation energy now agrees with experiment within the fit error; the
  remaining disagreement is the pre-exponential factor (0.17 times the measured one, beyond the factor of about 1.6
  that the ±3.4 kJ/mol Ea error allows), i.e. the entropy of activation. M06-2X still matches the rate only through
  cancellation of a low Ea against a low A.
- **Files:** `reactions/da_butadiene_ethylene_m062x/00_inputs/card_notes.toml` (experimental text and gas-phase
  limitation), `reactions/da_butadiene_ethylene_m062x/03_gas_phase/gas_phase.md` (conformer section appended),
  `docs/results/da_butadiene_ethylene.md` (section 9). The 2026-11-20 s-trans-only numbers stay in the gas-phase
  section as the s-trans rows. The earlier B3LYP-geometry folder (`reactions/da_butadiene_ethylene`) keeps its old
  wording; it is superseded.

## 12. Diels–Alder orbitals: the production path localized clockwise, the pilot path counterclockwise (2026-11-25)

- **What:** the pilot path (B3LYP-D3(BJ)/ma-def2-SVP, `reactions/da_butadiene_ethylene`) localized the three active
  IBOs counterclockwise round the ring, the pattern written in reaction.toml (π(C1=C2) → σ(C1–C6),
  π(C5=C6) → σ(C4–C5), π(C3=C4) → π(C2=C3)). The production path (M06-2X-D3(0)/ma-def2-SVP, SMD) localized them
  clockwise: π(C3=C4) → σ(C4–C5), π(C5=C6) → σ(C1–C6), π(C1=C2) → π(C2=C3). The two are mirror images through the
  TS mirror plane (C1↔C4, C2↔C3, C5↔C6). Both are equally valid descriptions of the same mechanism.
- **Tracking:** continuous along all 302 frames. The lowest frame-to-frame overlap of the active orbitals is 0.714,
  at frame 189 (s = +0.50 bohr·amu½); it is the lowest of all 17 tracked orbitals. The pilot's lowest was 0.827.
- **Can the viewer show the mirror image of the computed set? Checked at every frame**
  (`reactions/da_butadiene_ethylene_m062x/03_orbitals/mirror_display_check.md`), with J. Picione's limits of
  0.001 Å for the geometry and 0.999 for the orbitals:
  - **Geometry** (largest displacement of the reflected frame from itself, atoms swapped): within 0.0008 Å at 296
    frames. It fails at frames 0–2, 4 and 6 (s = −8.94 to −8.65; up to 0.0035 Å) and at frame 299
    (s = +5.96; 0.0040 Å).
  - **Orbitals** (smallest singular value of the 3×3 overlap matrix, reflected against computed active orbitals):
    at least 0.9990 elsewhere. It fails from the TS to s = +1.00 (frames 179–199), worst 0.9924 (determinant
    0.9848) at frame 181 (s = +0.10). At the TS it is 0.9953, although the TS geometry is mirror-symmetric to
    0.0000015 Å.
  - **Meaning:** between the TS and s = +1.0 the clockwise and counterclockwise localizations are different sets
    of orbitals, not one set relabelled.
- **Display (J. Picione, 2026-11-25, option 1):** the viewer and video show the computed clockwise set, and the
  card and this SI describe the clockwise arrows. No re-localization. The viewer's mirror toggle is off until a
  3×17 check (each reflected active orbital against all 17 computed valence orbitals, frames 179–199) shows that
  the reflected set lies in the computed occupied space.
- **Why the 3×3 test was stricter than needed (J. Picione):** where the geometry is mirror-symmetric, the mirror
  image of a valid localized set is itself a valid localized set of the same wavefunction. The three active
  orbitals need not span the same space as their mirror images, because the other 14 orbitals also differ
  slightly between the two localizations.
- **Flicker between the TS and s = +1.0 (J. Picione, 2026-11-25 evening):** at 6 frames the localizer landed on
  the counterclockwise set instead of the clockwise one: frames 181, 182, 185, 187, 188 and 198 (s = +0.10, +0.15,
  +0.30, +0.40, +0.45, +0.95 bohr·amu½). This is the cause of the 0.714 overlap at frame 189. Two tests, frames
  179–199 (`reactions/da_butadiene_ethylene_m062x/03_orbitals/mirror_continuity.json`):
  - **(a) 3×17:** each reflected active orbital against all 17 computed valence IBOs of the same frame. The
    smallest singular value is at least 0.99982 at every frame, so the reflected orbitals lie in the computed
    occupied space.
  - **(b) Continuity:** this frame's three active orbitals, as computed and mirror-reflected, against the
    previous displayed frame. The continuing choice overlaps at least 0.99834 at every frame (frames 180–199). The
    other choice is 0.67–0.96 (computed at the six frames above, mirror at the rest).
  - **Display:** at those 6 frames the viewer and video show the mirror image of the computed active orbitals,
    which continues the previous frame. This is valid because the geometry there is mirror-symmetric within
    0.0008 Å, so the mirror image of a localized set of this wavefunction is itself a localized set of it. The
    displayed frame-to-frame overlap of the three orbitals is now at least 0.9957 between frames 178 and 201. The
    computed cubes, .loc files and tracking.json are unchanged; the displayed copies are
    `03_orbitals/cubes/fNNN/ibo<id>_mirror.npz`.
- **Mirror toggle back on (J. Picione, 2026-11-25 evening):** the viewer's "mirror direction" switch reflects the
  displayed set through the TS mirror plane at every frame. The path is mirror-symmetric within 0.0008 Å except at
  frames 0–2, 4, 6 and 299 (up to 0.0040 Å). The 3×17 check holds at frames 179–199 (0.99982 or higher).
- **Tracking swaps (all six reactions, `docs/results/tracking_table.md`):** the only swaps are among the five
  equivalent bromine lone pairs (52 frames in E2, 1 in SN2), which do not change their atoms. No active orbital is
  involved in a swap.
- **Code:** `pipeline/orbital_steps._check_expected` now also accepts the mirror image of the written pattern
  (all items, never mixed). `pipeline/mirror_display_check.py` performs the check.

## 13. Solvation scheme for the DLPNO-CCSD(T) energies: scheme A kept; B − A given as an uncertainty (J. Picione, 2026-11-26)

Not a correction: no number in any table changes. A check of how SMD enters the DLPNO-CCSD(T)/ma-def2-TZVP energies.

- **Scheme A (in every table, card and video):** DLPNO-CCSD(T)/ma-def2-TZVP with SMD.
- **Scheme B (check only):** gas-phase DLPNO-CCSD(T)/ma-def2-TZVP + [E(SMD) − E(gas)] at M06-2X/ma-def2-TZVP, at the
  same geometries. Both schemes use the same working-level thermal corrections, reference species and 1 M standard
  state. 80 gas-phase single points, 2026-11-25 (`pipeline/smd_scheme.py`; jobs in
  `reactions/<reaction>/01_ts/smd_scheme/gas/`).
- **B − A, kJ/mol (kcal/mol), 298.15 K:**

| Reaction | ΔG‡ | ΔrG |
|---|---|---|
| SN2, HO⁻ + CH₃Br (water) | +3.0 (+0.7) | +6.6 (+1.6) |
| E2, EtO⁻ + 2-bromopropane (ethanol) | +3.2 (+0.8) | +4.0 (+1.0) |
| CN⁻ + acetone (water) | −1.2 (−0.3) | −0.1 (0.0) |
| Diels–Alder, butadiene + ethylene | +0.5 (+0.1) | +0.6 (+0.2) |
| 1,2-H shift, 3-methyl-2-butyl cation (DCM), path shown, from the eclipsed form | +0.4 (+0.1) | +0.3 (+0.1) |
| 1,2-H shift, 2,3-dimethyl-2-butyl cation (DCM) | +0.3 (+0.1) | 0.0 (0.0) |

(Secondary hydride shift row corrected 2026-11-26: ΔG‡ is now the barrier of the path shown, from the eclipsed form, 11.1 kJ/mol in scheme A; it was the TS shown measured from the anti form, 14.0 kJ/mol, which is not the barrier of any path (J. Picione, 2026-11-26). The anti form's own rotation TS has no gas-phase single point, so it has no B − A.)

- **Where the difference sits (per-species table in `docs/results/smd_scheme_check.md`):** the hydroxide and
  ethoxide ions lie 6.6 and 4.9 kJ/mol (1.6 and 1.2 kcal/mol) lower in scheme B; bromide is the same in both
  schemes (0.0 kJ/mol), and the neutral reactants and products shift by 1.1 kJ/mol or less. The SN2 and E2 transition
  states shift by 4.3 and 2.7 kJ/mol, so the barriers change by less than the reaction free energies. In the cyanide
  reaction the species shift by 1.2 (cyanide) to 6.6 kJ/mol (TS), including 4.2 kJ/mol for acetone, but the shifts
  cancel to within 1.2 kJ/mol.
- **Decision:** scheme A is kept, for consistency with the paths, cards and videos. The two reactions that start from
  hydroxide or ethoxide carry an extra 3 to 7 kJ/mol (0.7 to 1.6 kcal/mol) of uncertainty from the solvation scheme;
  the other four agree within 1.2 kJ/mol (0.3 kcal/mol).
- **Scheme B is not uniformly closer to experiment:** it moves the SN2 barrier and reaction free energy away from the
  measured values (by 3.0 and 6.6 kJ/mol), and the E2 barrier (by 3.2 kJ/mol) and the cyanide barrier (by 1.2 kJ/mol)
  toward them; the cyanide reaction free energy moves by 0.1 kJ/mol. There is no measured E2 reaction free energy.
- **Files:** `docs/results/smd_scheme_check.md` and `.json` (`pipeline/smd_scheme_report.py`).

## 14. Secondary hydride shift: orbital check written from the mechanism; a three-centre bond and a hyperconjugation change documented (J. Picione, 2026-11-26)

A correction. Until 2026-11-26 the orbital check of the secondary hydride shift (3-methyl-2-butyl cation) compared
nothing, because its input had no expected list. The IBO output check passed (energies within 5 × 10⁻¹¹ Eh of the
path), but tracking found three changing orbitals and none of them was checked against the mechanism.

- **Expected list, from the mechanism:** one change, the arrow the textbook draws, σ(C3–H) → σ(C2–H). Its atoms come
  from `[key_atoms]`, set on 2026-11-13 from the TS before any orbital calculation. It was not written from the
  tracking output.
- **Documented, a separate category:** two changes found by tracking that the textbook drawing does not show. They
  are not extra arrows:
  - σ(C1–H), hyperconjugating into C2 at the start → ordinary σ(C1–H): hyperconjugation into the cationic centre
    withdrawn;
  - σ(C3–CH₃), a three-centre bond at the start → ordinary σ(C3–CH₃): the methyl bridges C2 and C3 (C2–CH₃ 1.70 Å,
    C3–CH₃ 1.80 Å), so the starting ion is not the open secondary cation of the textbook drawing.

  (First recorded the same day as two hyperconjugation changes; the C3–CH₃ pair was named a three-centre bond by
  J. Picione, 2026-11-26.)
- **Result:** 1 expected change, found; 2 further changes, documented; no unexpected change. Lowest frame-to-frame
  overlap 1.000 for all three. The comparison now has a "documented" category beside "expected" and "unexpected"
  (`pipeline/orbital_steps.py`).
- **Also new:** frames 40–144 are flagged. Two C–H orbitals with the same atoms at both ends briefly cross the 0.15
  threshold on C3 during the C2–C3 rotation: σ(C4–H) up to 0.21, and σ(C–H) of the CH₃ on C3 up to 0.16.
- **Numbering:** atoms are named in the chemical numbering of the 3-methyl-2-butyl cation,
  CH₃(C1)–C2H⁺–C3H(CH₃)–C4H₃, everywhere a reader sees them (viewer labels, tracking table, card). The 0-based atom
  indices stay only in the tracking log.
- **Shares on C2 (Mulliken populations of single IBOs in ma-def2-SVP; same caveat as for E2, not IAO charges):**
  - σ(C1–H), hyperconjugation: 0.17 at the reactant end, 0.04 at the TS, 0.00 at the product end;
  - σ(C3–CH₃), three-centre bond at the start: 0.36, 0.03 and 0.00.

  Half of each loss is done by s ≈ −7 bohr·amu½, before the TS (s = 0, the C2–C3 rotation TS) and before the hydrogen
  moves (s +3.9 to +5.5).
- **New cationic carbon C3:** no orbital next to C3 crosses the 0.15 threshold at the product end. The largest share
  there is 0.11. Shares of this size also occur on carbons that carry no charge, so this measure does not resolve
  hyperconjugation into C3.
- **Tertiary shift, same check:** its expected list was written from the mechanism on 2026-11-10, before any
  calculation. Tracking finds no delocalized orbital and no other change at any frame.
- **Files:** `docs/results/hshift_hyperconjugation.md`; `reactions/hshift_*_m062x/03_orbitals/hyperconjugation.json`,
  `tracking.json` and `tracking_log.md`.
