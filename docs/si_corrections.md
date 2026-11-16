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
  The secondary backward side is continued the same way from its last frame (started 2026-11-16 16:17).
- **Angle check on the continued frames:** in the slow final approach the mass-weighted angle between step and
  gradient grows (frames past the earlier end: median 6.83°, 90th percentile 10.92°, 150 frames; earlier frames
  beyond |s| = 1: 4.41° / 7.81°). The rerun from the TS has the same frames, so this belongs to the region, not to the
  continuation. The open side's zone beyond |s| = 1 now has median 5.74° (rule < 5°), 90th percentile 10.18° (rule
  < 15°); kept as computed, both parts stated (J. Picione, 2026-11-16 evening). All 151 continued frames are kept
  (every frame computed by the time of the comparison, not only the 54 that overlap the rerun; J. Picione, same
  evening); cutting at the overlap would have ended the side about 3.1 kJ/mol (0.75 kcal/mol) above the open cation.
- **Side ends** (last frame above the official end minimum, working level; filled in as sides finish):

  | Path | Reactant side | Product side |
  |---|---|---|
  | SN2 | 0.06 kJ/mol | 0.23 kJ/mol |
  | Diels–Alder | 0.33 kJ/mol | 0.01 kJ/mol |
  | Cyanide + acetone | 0.94 kJ/mol | −0.03 kJ/mol |
  | E2 | 1.17 kJ/mol | 5.01 kJ/mol |
  | Tertiary hydride shift | open side: 1.47 kJ/mol (end check pending) | bridged side: at the bridged structure |
  | Secondary hydride shift | pending (continued from 2026-11-16 16:17) | pending (forward side, join waits) |

## 6. Angle check against gradient size (all six paths; for the SI figure)

- **Requested:** J. Picione, 2026-11-16 evening, after the tertiary open side's angle check (entry 5). Existing data only,
  no new calculation: `pipeline/angle_vs_gradient.py`; figure `docs/si_figures/angle_vs_gradient.png`, per-frame values
  `docs/si_figures/angle_vs_gradient.json`.
- **What is plotted:** for every computed frame of every production side (2108 frames; mirror frames of the tertiary
  shift have no jobs and are left out), the mass-weighted angle between the path direction and the downhill gradient
  (the angle of the angle check) against the RMS Cartesian gradient at that frame (Eh/bohr, log scale). Frames are marked
  by zone: near the TS (|s| ≤ 1 bohr·amu½, as in the angle check); final approach (from the frame after which the RMS
  gradient stays below 5 × 10⁻⁴ Eh/bohr to the side's end, ten times the end-minimum test's 5 × 10⁻⁵); middle (the rest).
  The secondary shift's backward side is still being extended; its first part (85 frames) is used, and the figure is
  redrawn when the side ends.
- **Result: the large angles occur only where the gradient is nearly zero.** All 59 frames with an angle above 15° have
  an RMS gradient of at most 2.0 × 10⁻⁴ Eh/bohr: 30 lie near the TS of the two hydride shifts, whose TSs are very flat
  (tertiary 21, secondary 9), and 29 in the slow final approach to a minimum (cyanide product side 18, Diels–Alder 5,
  tertiary open side 3, SN2 2, secondary 1). None is in the middle of a side. Where the RMS gradient exceeds
  5 × 10⁻⁴ Eh/bohr (1196 frames) the largest angle is 5.65° (median 0.09°), and above 1 × 10⁻³ Eh/bohr (927 frames)
  4.07° (median 0.04°). Near the TS of the four reactions with sharp TSs (SN2, cyanide, Diels–Alder, E2) the angles
  stay small (largest 5.65°). The largest angles (up to 117°) are on the cyanide product side, which ends at its
  minimum (−0.03 kJ/mol), where the gradient is about 1.5 × 10⁻⁵ Eh/bohr and its direction carries little information.
- **Reading for the paper:** a large angle in these paths marks a region where the energy surface is nearly flat (a very
  flat TS, or the last approach to a minimum along soft modes), not a path that leaves the valley: wherever the gradient
  is appreciable, path and gradient agree within a few degrees. The angle-check results are kept as computed (J. Picione,
  2026-11-15 and 2026-11-16).
