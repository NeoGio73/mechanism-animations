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

