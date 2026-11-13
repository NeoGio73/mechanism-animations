# Computational protocol (draft methods text)

## Working level for all six reactions (decided by J. Picione, 2026-11-11)

The working level for geometries, transition states, IRCs and orbitals of all six reactions is
**M06-2X-D3(0)/ma-def2-SVP** on the DefGrid3 integration grid, RIJCOSX (def2/J), VeryTightSCF (from
2026-11-12; TightSCF before, see "Step, first step and SCF for production"), in each
reaction's SMD solvent model. This replaces the B3LYP-D3(BJ)/ma-def2-SVP standard below and the
3-methyl-2-butyl exception, both kept below as the record of earlier work.

Why:
- B3LYP fails for the two cations. It understates the methyl bridging of the 3-methyl-2-butyl
  cation, and gives a rotation barrier 13–14 kJ/mol below M06-2X. For the 2,3-dimethyl-2-butyl
  cation it makes the hydrogen-bridged ion a transition state, where M06-2X finds a shallow minimum.
- In a functional sensitivity test (`docs/sensitivity_M06-2X.md`, kept as it is for the supporting
  information), the DLPNO-CCSD(T)/ma-def2-TZVP barriers of SN2, cyanide + acetone and Diels–Alder
  changed by 0.5 kJ/mol or less between B3LYP and M06-2X geometries.
- E2 differed by 7.4 kJ/mol and waits for an arbitration (ωB97X-D3(BJ) and RI-MP2 geometries,
  coupled cluster along both paths near the TS).

SN2, Diels–Alder and cyanide + acetone are recomputed in full at this level in new folders once the
IRC step-size test has passed (single points along the path: see "Single points along the path",
2026-11-13); the earlier results stay. The
B3LYP-D3(BJ)/ma-def2-SVP SN2 transition state and its single points
(`reactions/sn2_hydroxide_bromomethane_masvp`) are kept for the sensitivity table.

### A crashed or incomplete job is never a result (J. Picione, 2026-11-11)

Every automatic decision (pipeline steps, queue tasks, test scripts) treats a job that did not end
with ORCA's normal termination, or a path that did not reach the range it was meant to cover, as
incomplete: no value is taken from it and no follow-up decision is made on it. It is logged in
`run_log.md` and rerun with the same settings, or reported.

### Core and memory budget (J. Picione, 2026-11-11)

All ORCA jobs running at the same time, across all lanes, use at most the 8 physical cores of the
machine and about 28 GB of memory (31 GB installed). Memory is counted as declared: nprocs ×
%maxcore per job. The budget is enforced by the queue itself (`pipeline/cores.py`): every ORCA job
asks for its cores and memory before it starts and waits, first come first served, until they are
free. Each lane start, job start, wait and job end is logged with the allocation in
`docs/level_test_queue/allocation.log`. Lanes are run one behind another when they would not fit
side by side; the E2 arbitration runs behind the step-size test. Paused jobs count against the
budget with the memory they actually hold. A lane is not paused when that would take the total above
28 GB or make Windows page; it is stopped and later rerun from scratch, and any job that does not end
normally after a pause is rerun from scratch with no wall time reported.

### One ORCA job per folder, always (J. Picione, 2026-11-11)

A job starts in a folder with no files left from a different job (one outside its own production
chain), or from a crashed run (J. Picione, 2026-11-11). Jobs that run one after another in one production chain may share a
folder, as before (J. Picione, 2026-11-12): for example the backward IRC after the forward IRC in
02_irc, or the single points of one path chain. A rerun after a crash starts in a new folder; the
crashed run's folder is kept whole under crashed_<time>_<direction>/. Before a test job starts,
anything in its folder besides its own input and the Hessian is moved aside
(leftover_<time>_<direction>/).

No two ORCA jobs ever run in the same folder at the same time. ORCA writes fixed-name files in the
job folder (redirect.tmp, the text output of sub-steps such as the dispersion printout, from
orca_tools/qcmsg.cpp), so two jobs running in one folder can overwrite each other's files. The queue
enforces this: a job whose folder is busy waits. Jobs that run strictly one after the other in one
folder (for example the forward and then the backward IRC in 02_irc) cannot overlap; the folder
layout of the pipeline stays as it is (J. Picione, 2026-11-11: no layout change mid-project). The
only exception is the deliberate reproduction test of this effect, in a scratch folder; it is low
priority, runs after the step-size verdict and the E2 lane, and takes the whole machine for each
round so it never runs beside a production job. If it does not crash, the cause of the step-test
crashes is recorded as unknown.

### An IRC is a result only if both sides reached their end point (J. Picione, 2026-11-12)

A side has reached its end point when ORCA ends it with "THE IRC HAS CONVERGED" and it is not a false
stop next to the TS (at least 10 steps, end at least 1 kJ/mol below the TS), or when the IRC stopping
rule (section below) is met within it. A side ended by MaxIter before the rule, by a crash, or by the gradient
criterion while still next to the TS (the SN2 0.025 runs stopped after one step that way: near the TS
the gradient was already below ORCA's thresholds) is incomplete, never finished. pipeline/path_quality.py
applies this to every IRC (docs/irc_completeness.json); the card marks an incomplete path, and the
step-size report is not written from one.

## Standard for the reaction set (decided by J. Picione, 2026-11-09; replaced 2026-11-11, see above)

From E2 onward, every reaction uses one working level for geometries, the TS, the IRC and the
orbitals: **B3LYP-D3(BJ)/ma-def2-SVP**, RIJCOSX (def2/J), TightSCF, SMD solvent. The minimally
augmented basis adds diffuse functions to the heavy atoms, which the anions in most of these
reactions need. Barrier single points use B3LYP-D3(BJ)/ma-def2-TZVP. Energies at every IRC frame are
recomputed at DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO), and that curve is the default in the viewer and
the main video; M06-2X-D3(0)/ma-def2-TZVP is also computed along the path for comparison.

Why: in the SN2 pilot the def2-SVP working level (no diffuse functions) put the TS 10.9 kJ/mol
(2.6 kcal/mol) below separated reactants. The barrier needed rescuing with coupled-cluster single
points, which gave ΔG‡ = 95.6 kJ/mol (22.9 kcal/mol) against 95.0 kJ/mol (22.7 kcal/mol) measured.
The SN2 pilot below was computed at def2-SVP. It will be rerun at ma-def2-SVP so the whole set shares
one working level (`docs/tasks.md`); the def2-SVP results stay until then.

### Points next to the TS that lie above it (decided by J. Picione, 2026-11-09)

ORCA's first IRC step leaves the TS along a straight line. On a strongly curved path (a hydrogen
transfer) that point can lie slightly above the TS. Any IRC point within two steps of s = 0 that
lies above the TS energy is flagged by `profile`, kept in the data files, left out of the viewer
path and the video so the displayed maximum is the TS, and stated on the card with its energy.
E2: the first reactant-side point is 1.3 kJ/mol (0.3 kcal/mol) above the TS. On E2 the M06-2X/ma-def2-TZVP
curve peaks one frame from the working-level TS (that same point, 0.6 kJ/mol (0.1 kcal/mol) above the TS
frame) because the single points sit on a path optimized at another level, while the
DLPNO-CCSD(T)/ma-def2-TZVP maximum falls on the TS frame.

### Exception: 3-methyl-2-butyl cation at M06-2X (decided by J. Picione, 2026-11-11; since then M06-2X is the working level for all reactions)

The working level for this one reaction is M06-2X-D3(0)/ma-def2-SVP in SMD dichloromethane, for
the minimum, the TS, the IRC and the orbitals, instead of the standard B3LYP-D3(BJ)/ma-def2-SVP.
Reason: B3LYP understates the methyl bridging of the secondary cation (gas-phase "anti" minimum:
C2–C3–CH₃ 81.6° at B3LYP against 71.3° at M06-2X and 71.6° at RI-MP2/ma-def2-SVP; M06-2X and RI-MP2
agree with each other, not with B3LYP), and it gives a rotation barrier 13–14 kJ/mol lower than
M06-2X (gas phase: "eclip" TS 6.9 against 20.0 kJ/mol, "anti" TS 16.0 against 29.8 kJ/mol, ΔE).
The path animated is the eclip-side one (lower at every level and in both phases).

### M06-2X integration grid

M06-2X geometry optimizations and frequencies from 2026-11-11 on use DefGrid3 (ORCA 6 default:
DefGrid2), because meta-GGA functionals such as M06-2X are sensitive to the integration grid.
M06-2X/ma-def2-TZVP single points stay on the default grid, as for the single points already on
the cards.

### First IRC step from the TS (decided by J. Picione, 2026-11-11)

ORCA's default first step is an initial displacement along the TS mode chosen for an expected
energy change of 2 mEh (ORCA 6.1 manual, %irc block: `DE_Init_Displ 2.0`, default; `Init_Displ`
DE or length; `Scale_Init_Displ`, default 0.1 a.u., used with length). On the finished paths it
gave one long step on each side of the TS (cyanide + acetone 1.020 against a median step of 0.042
bohr·amu½; Diels–Alder 0.630 against 0.093; SN2 (def2-SVP) 0.737 against 0.210; E2 0.366 against
0.163). The first point ORCA stores is the displaced geometry itself, so its distance from the TS is
s₁ = Scale_Init_Displ · √(Σ m_i |v_i|²) for the normalized Cartesian TS mode v (reproduces all four
first steps). New IRCs set `irc_init_step_s` in reaction.toml (about the median step of the path);
the pipeline writes `Init_Displ length` and the matching `Scale_Init_Displ` (tested on the SN2 rerun
TS with ORCA 6.1.1: first point at s = 0.210 bohr·amu½ as requested).

### IRC step length (decided by J. Picione, 2026-11-11)

Every path computed from 2026-11-11 on uses ORCA's non-adaptive step (`Adapt_Scale_Displ false`;
`irc_adapt_scale_displ = false` in reaction.toml), a first step of 0.10 bohr·amu^1/2 set by length,
and `Scale_Displ_SD` calibrated so that the median step in the bond-changing region is about
0.10 bohr·amu^1/2. ORCA's correction step after each descent step stays on (`Do_SD_Corr` unchanged).
The steps are not equal and are not required to be.

Why: a test on the SN2 rerun TS (`02_irc/attempts/constant_step/`) showed that switching adaptation
off does not give equal steps (first step 0.100, later steps 0.065–0.158 bohr·amu^1/2; switching
off the parabolic fit as well changed nothing). The adaptive path, with steps of 0.21–0.34 there,
lay up to 1.3 kJ/mol (0.3 kcal/mol) above the 0.10 path in the steep product-side stretch at
s = −1.7 to −3.7 bohr·amu^1/2, 1.55 kJ/mol at matching C–O distance. The larger steps are probably
the less accurate ones, so agreement with the adaptive path was not used as a criterion.

Step-size convergence test, on SN2 before any path is recomputed: the same non-adaptive IRC with
first step and median at 0.05 bohr·amu^1/2, compared with the 0.10 path at matching C–O distance
on the reactant side and at matching C–Br distance on the product side (energy and the other key
bond). Criterion: within 0.2 kJ/mol and 0.005 Å throughout the bond-changing region. If it passes,
0.10 is the step for the whole set.

Result (2026-11-12): four paths (0.10; 0.05 with first step 0.05; 0.05 and 0.025 with first step
0.10; `02_irc/attempts/step_convergence/report4.md`). Near the TS (|s| ≤ 1 bohr·amu^1/2) every pair is
within 0.154 kJ/mol and 0.0028 Å, and the DLPNO-CCSD(T)/ma-def2-TZVP maximum lies at the same place on
all four paths (s −0.175 to −0.187 bohr·amu^1/2, C–Br 2.3967–2.3978 Å, height within 0.06 kJ/mol;
`cc_near_ts.md`). Further down the product side the paths differ by up to 0.75 kJ/mol and 0.008 Å,
and the differences do not shrink with the step: they come from the O–H bond, which drifts by up
to 0.03 Å along every path (`oh_check.json`, `angle_check.md`), not from the step size. SCF noise
contributes at most 0.094 kJ/mol to a difference of two paths.

### Step, first step and SCF for production (decided by J. Picione, 2026-11-12)

For all six reactions: IRC step **0.10 bohr·amu^1/2** (ORCA sets it in plain bohr through
`Scale_Displ_SD`, calibrated per reaction; `Adapt_Scale_Displ false`), first step **0.10
bohr·amu^1/2** (`Init_Displ length`, `Scale_Init_Displ` from the TS mode), and **VeryTightSCF** for
every IRC gradient and every per-frame single point at all three levels (M06-2X/ma-def2-TZVP,
B3LYP/ma-def2-TZVP, DLPNO-CCSD(T)/ma-def2-TZVP).

Why VeryTightSCF: each IRC step starts its SCF from the previous step's orbitals, so the same
geometry gives slightly different energies depending on the starting orbitals. Measured at the
M06-2X working level (DefGrid3, SMD water, six SN2 geometries; `check1/scf_noise.json`): TightSCF
up to 3.9e-5 Eh (0.102 kJ/mol) and 1.4e-5 Eh/bohr; VeryTightSCF 1.3e-5 Eh (0.034 kJ/mol) and
2.6e-6 Eh/bohr, for about 2 s (9 %) more per evaluation. The solvent model adds no scatter (gas
phase is not quieter); DefGrid3 halves it at B3LYP.

The path method is set by the rule below (2026-11-13); ORCA's IRC is no longer the production path.

### Production path method (by J. Picione's rule of 2026-11-13; validated 2026-11-13, 15:58)

Rule (J. Picione, 2026-11-13): take the pysisyphus tail option that completes both sides with a
mass-weighted angle in the tails of median below 5° and 90th percentile below 15° (the cheaper one if
both pass); if neither passes at 0.10, test at 0.05; then the full SN2 path must end both sides
complete, each end within 0.5 kJ/mol of its optimized minimum.

Applied: at 0.10 neither option passed. At 0.05 Gonzalez–Schlegel passed (reactant tail 1.25° / 4.29°,
product tail 0.00° / 1.29°) and EulerPC with a Hessian every 5 steps failed (reactant 90th percentile
59.0°) (`true_irc2/tails.md`). Full SN2 validation path `true_irc2/final_gs_s005` (`final_final_gs_s005.md`):
**PASS**. Reactant side 179 steps, s = 8.94 bohr·amu^1/2, end −0.005 kJ/mol from the minimum (RMSD
0.027 Å); product side 203 steps, s = 10.15, end +0.322 kJ/mol (0.077 kcal/mol) above the minimum
(RMSD 0.080 Å); both ended by the stopping rule; mass-weighted angle median / 90th percentile
0.28° / 5.44° (reactant) and 0.14° / 5.22° (product) within |s| ≤ 1. 939 gradients, 3.0 h wall.

Production settings (all reactions, working level M06-2X-D3(0)/ma-def2-SVP, DefGrid3, SMD as in each
reaction.toml, VeryTightSCF; `pipeline/true_irc_prod.py`):
- Integrator: pysisyphus 1.0.0 Gonzalez–Schlegel in mass-weighted Cartesian coordinates, with the
  project's 2026-11-13 fixes (Bofill update; translations and rotations projected out in the current
  frame; bracketed λ search; a failed search stops the side as a failure, never a result).
- Step 0.05 bohr·amu^1/2, measured in mass-weighted coordinates; first step 0.05 along the
  mass-weighted TS mode.
- Hessian: ORCA's exact TS Hessian at the working level at the start (`01_stationary/ts/ts.hess`),
  Bofill updates after that, no recomputed Hessians.
- End of a side: the project stopping rule only (|ΔE| < 0.05 kJ/mol per step for 20 consecutive steps
  after the last bond event); pysisyphus's own gradient and energy stops are switched off; 3000 steps
  at most. A side that ends any other way makes the path incomplete; nothing is joined.
- One ORCA job per folder (`02_irc/true_irc/<direction>/g####`), through the core ledger.
- Joined path: backward (reversed) + TS + forward, all energies from ORCA EnGrad jobs at the working
  level (the TS energy from the path's own gradient job at the TS).

References:
- C. Gonzalez, H. B. Schlegel, "An improved algorithm for reaction path following", J. Chem. Phys.
  1989, 90, 2154–2161. doi:10.1063/1.456010
- C. Gonzalez, H. B. Schlegel, "Reaction path following in mass-weighted internal coordinates",
  J. Phys. Chem. 1990, 94, 5523–5527. doi:10.1021/j100377a021
(Metadata checked against Crossref, 2026-11-13.)

Modifications to pysisyphus: every change is kept as a patch file against pysisyphus 1.0.0 in
`patches/pysisyphus/`, with one line each on what it changes and the test that shows it works
(J. Picione, 2026-11-13: the paper states exactly what was modified).

Order of production (J. Picione, 2026-11-13): SN2, cyanide + acetone, Diels–Alder, tertiary hydride
shift, secondary hydride shift, then E2 (decided 2026-11-13 after the arbitration, same method).

### Angle check on every production path (J. Picione, 2026-11-13)

Every production path gets the mass-weighted angle check: at each frame the angle between the step
M^1/2 Δx and the mass-weighted downhill direction −M^−1/2 g (central chord), reported as median and
90th percentile within |s| ≤ 1 bohr·amu^1/2 and beyond. A path fails the rule if, in either zone, the
median is 5° or more or the 90th percentile is 15° or more (the tail rule of 2026-11-13); a failing
path is flagged in the report and on its card. The SN2 validation path: 0.28° / 5.44° and 0.14° /
5.22° within |s| ≤ 1, 1.04° / 8.01° and 0.66° / 1.53° beyond (reactant / product side).

### Single points along the path (J. Picione, 2026-11-13; replaces the sampling of 2026-11-09)

All six reactions, on the working-level production path:
- M06-2X/ma-def2-TZVP and B3LYP/ma-def2-TZVP at every frame.
- DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO) at every frame within the bond-changing region (frames at
  which a forming or breaking key bond has covered more than 5 % and less than 95 % of its change,
  as in the path-quality rule) and at every frame within |s| ≤ 1 bohr·amu^1/2; at every second frame
  elsewhere (counted from the TS); both end frames always.
- No interpolation: curves show only computed points, and no value between computed frames is
  reported or used.

### Secondary hydride shift: key bonds (accepted by J. Picione, 2026-11-13)

3-methyl-2-butyl cation at M06-2X: the TS is the 92i cm⁻¹ C2–C3 rotation TS of the eclipsed
secondary cation (H9 1.113 Å from C3, 2.008 Å from C2). Key bonds: C2–H9 forms, C3–H9 breaks. They
label the path and set the bond-changing region; they do not change the path. If the production
path ends in another conformer instead of the tertiary cation, that is reported plainly and the path
is not presented as the hydride shift.

### Mass-weighted IRC test with pysisyphus (2026-11-12; a test, not yet the method)

Software: pysisyphus 1.0.0 (GPL-3), in its own environment `.venv-pysis` (CPython 3.12.15 from uv's
user folder; numpy 2.5.4, scipy 1.18.1, h5py 3.16.0; full list `docs/pysisyphus_env_freeze.txt`),
with ORCA 6.1.1 for every energy and gradient through `pipeline/true_irc.py` (one ORCA job per call,
alone in its own folder, through the core ledger). The pipeline's own environment (Python 3.15) cannot
host pysisyphus: h5py has no package for Python 3.15.

Integrator and settings: EulerPC, the Euler predictor–corrector of Hratchian, Frisch and Schlegel, in
mass-weighted Cartesian coordinates; Bofill Hessian updates; start Hessian = ORCA's exact TS Hessian at
the same level; corrector on the distance-weighted interpolant (pysisyphus `corr_func = "mbs"`), up to
5000 Euler predictor substeps per step (pysisyphus default 500; with 500 the predictor fell short
of the step in the flat tail and the path stalled), 3 loose cycles. pysisyphus measures
`step_length` in plain (un-mass-weighted) Cartesian bohr, like ORCA's `Scale_Displ_SD`; only the first
displacement is mass-weighted. Steps: 0.02359 bohr for the 0.10 bohr·amu^1/2 path and 0.01179 bohr
for the 0.05 path (the ORCA step-test calibration, so both methods take the same Cartesian steps;
mass-weighted steps are measured and reported); first step 0.10 bohr·amu^1/2 along the mass-weighted
TS mode (`displ = "length"`). A first attempt with `step_length` read as bohr·amu^1/2 took steps of
0.35 and 0.18 bohr·amu^1/2 and is kept as not a result (`true_irc/attempt1_step_in_plain_bohr`).
A direction ends when the RMS of the (non-mass-weighted) gradient falls below 5e-4 Eh/bohr after the
inflection point, when the energy rises, or after 1500 steps. Test level: B3LYP-D3(BJ)/ma-def2-SVP,
RIJCOSX, SMD water, VeryTightSCF.

Result on SN2 (2026-11-12, 21:11; `true_irc/report.md`, `detail.json`; a test result, the method is
not decided):
- Angle check validated. On the pysisyphus paths the mass-weighted angle, angle(M^1/2 Δx, −M^−1/2 g),
  has median 0.07° (0.10 path) and 0.01° (0.05 path) within |s| ≤ 1, against 22.8° to the plain
  gradient; the Cartesian form angle(Δx, −M^−1 g) gives 0.25° and 0.02°. The same check on the ORCA
  paths gives 38–56° (mass-weighted) and 1–23° (plain), so ORCA's path is a Cartesian
  steepest-descent path, not a mass-weighted IRC.
- Off-gradient stretches on the pysisyphus path itself, in the flat tails: reactant side from C–O
  about 2.9 Å, product side from C–Br about 3.26 Å, RMS gradient 0.6–2e-3 Eh/bohr. There the steps run
  straight (turn between steps at most 10°) but up to 30–67° off the mass-weighted gradient; the
  off-path gradient sits on the hydroxide H (reactant side), and on the methanol O–H H and the methyl
  H atoms (product side). Energies and key bonds are unaffected at the 0.03 kJ/mol, 0.0003 Å level
  (0.10 vs 0.05 below).
- Measured mass-weighted steps: 0.10 path 0.082–0.084 bohr·amu^1/2 near the TS, median 0.064–0.069
  beyond; 0.05 path 0.041–0.042 near the TS, 0.032–0.034 beyond (first step 0.10 on both).
- Step convergence (pysisyphus 0.10 vs 0.05, report4 matching and limits): whole region 0.031 kJ/mol
  and 0.0003 Å, |s| ≤ 1 0.029 kJ/mol and 0.0000 Å; converged at 0.10, below the 0.094 kJ/mol noise floor.
- O–H: 0.962–0.970 Å along every pysisyphus side, within 0.005 Å of its relaxed value (ORCA paths
  0.959–0.992 Å, up to 0.016 Å from it).
- ORCA 0.10 vs pysisyphus 0.10, matched on the key bond: whole region up to 18.3 kJ/mol and 0.091 Å
  (product side, C–Br 2.80 Å, ORCA higher), |s| ≤ 1 up to 2.3 kJ/mol and 0.025 Å.
- Completeness: all four sides end on the RMS-gradient criterion (5e-4 Eh/bohr), none on an energy
  rise; ends agree between 0.10 and 0.05 (RMSD 0.028 Å reactant, 0.002 Å product) and with the ORCA
  ends (0.032 Å, 0.018 Å). The 0.10 reactant side was still falling 0.10–0.27 kJ/mol per step at its
  end; the distance of each end from its minimum is measured by `pipeline/true_irc_ends.py`.
- Cost: 2 ORCA gradients per step, 21.5 s each at 2 cores, 43 s per step (ORCA IRC at 0.10: 26 s per step).

References to cite if adopted:
- J. Steinmetzer, S. Kupfer, S. Gräfe, "pysisyphus: Exploring potential energy surfaces in ground and
  excited states", Int. J. Quantum Chem. 2021, 121, e26390. doi:10.1002/qua.26390
- H. P. Hratchian, M. J. Frisch, H. B. Schlegel, "Steepest descent reaction path integration using a
  first-order predictor–corrector method", J. Chem. Phys. 2010, 133, 224101. doi:10.1063/1.3514202
- H. P. Hratchian, H. B. Schlegel, "Accurate reaction paths using a Hessian based predictor–corrector
  integrator", J. Chem. Phys. 2004, 120, 9918–9924. doi:10.1063/1.1724823
- J. M. Bofill, "Updated Hessian matrix and the restricted step method for locating transition
  structures", J. Comput. Chem. 1994, 15, 1–11. doi:10.1002/jcc.540150102 (Hessian update)
(Metadata checked against Crossref, 2026-11-12.)

Each card states the first step, the median, the shortest and the longest step in the
bond-changing region, the largest change of a key bond in one step, and one sentence on the
convergence test.

### Path-quality rule (decided by J. Picione, 2026-11-11)

Within the bond-changing region (frames at which a forming or breaking bond has covered more than
5 % and less than 95 % of its change), every IRC step longer than 2.0 times the median step of that
region (not of the whole path, whose long flat approach lowers the median) is flagged and stated on
the card (pipeline/path_quality.py; full list in `02_irc/path_quality.json`). The factor was first
set at 1.5 and raised to 2.0 with the non-adaptive step, whose steps vary by about ±50 %
(J. Picione, 2026-11-11).

### Barrier from the TS frame and from the curve maximum (decided by J. Picione, 2026-11-11)

The barrier quoted first is the one at the working-level transition state, with its thermal
correction. Single points at a higher level along the working-level path can peak on another
frame. Whenever the default-level curve has its highest point on another frame and more than
1 kJ/mol above the TS frame, the card also gives the barrier from that maximum (the electronic
difference added to the TS-frame barrier; the thermal correction stays that of the TS frame,
stated as such, since there are no frequencies at the maximum), the frame number, and both
differences from experiment. This is the IRCMax practice of D. K. Malick, G. A. Petersson and
J. A. Montgomery Jr. (J. Chem. Phys. 1998, 108, 5704–5713). The card builder applies the rule
automatically (pipeline/card.py). At present it applies to SN2 only: the DLPNO-CCSD(T)/ma-def2-TZVP
maximum is IRC frame 18, 4.1 kJ/mol (1.0 kcal/mol) above the TS frame 17, so ΔG‡ is 95.6 kJ/mol
(22.9 kcal/mol) from the TS frame and 99.7 kJ/mol (23.8 kcal/mol) from the maximum, +0.6 and
+4.7 kJ/mol from experiment (95.0 kJ/mol); the offset is expected to shrink in the ma-def2-SVP
rerun. E2, Diels–Alder and cyanide + acetone have their maximum on the TS frame.

A gas-phase measurement is compared with the gas-phase computed value; the solution value is
given for reference only (Diels–Alder: Ea 119.8 kJ/mol (28.6 kcal/mol) in the gas phase, +4.8 kJ/mol
from experiment; 118.1 kJ/mol (28.2 kcal/mol) in SMD toluene, for reference).

### Implicit solvent and anionic products (noted by J. Picione, 2026-11-10)

SMD represents water as a continuum, with no explicit water molecules. It therefore misses the specific hydrogen bonds that water donates to a localized anion, and it under-stabilizes an alkoxide product. A computed step that forms an alkoxide (cyanide + acetone) is therefore likely too endergonic. The error comes from the solvent model, not the electronic-structure level: for cyanide + acetone the alkoxide lies +34.4 kJ/mol (+8.2 kcal/mol) above separated reactants in free energy at the working level and +36.3 kJ/mol (+8.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (298 K, 1 M). No explicit water molecules are added; the card states the limitation.

### Result: the 3-methyl-2-butyl cation has no classical secondary minimum to shift from (accepted by J. Picione, 2026-11-10)

The first choice for the carbocation 1,2-hydride shift was 3-methyl-2-butyl → 2-methyl-2-butyl. The rule set beforehand was that the secondary cation must be a true minimum (all real frequencies, not hydrogen-bridged) with a barrier to the shift; otherwise switch to the degenerate shift in the 2,3-dimethyl-2-butyl cation.

The secondary cation was optimized at the working level, in SMD dichloromethane and in the gas phase, from three rotamers about C2–C3 (`reactions/hshift_3_methyl_2_butyl`, run log there). Two findings came out:

- **No barrier for the shift itself.** With the C3–H in line with the empty p orbital on C2, the optimization runs downhill without a barrier to the tertiary 2-methyl-2-butyl cation, in both phases.
- **The minimum that does exist is partly methyl-bridged, not classical.** With the C3–H turned 64–78° away from the p axis, the optimization finds a true minimum: all frequencies real, the hydrogen still on C3. But one C3–CH₃ bond is stretched to 1.60–1.61 Å in DCM and 1.66–1.67 Å in the gas phase. That methyl leans toward C2, with a C2–C3–CH₃ angle of 92–95° in DCM and 82–85° in the gas phase. This is strong C–C hyperconjugation. Only rotation about C2–C3 keeps the hydrogen from shifting.

The tertiary cation lies 56.9 kJ/mol (13.6 kcal/mol) below this minimum in DCM and 54.2 kJ/mol (13.0 kcal/mol) below in the gas phase (ΔE).

This is a chemical result, not a failure of the calculation. A secondary cation next to a tertiary centre does not survive as a classical ion at this level; it is held only by its conformation. The degenerate tertiary → tertiary shift in the 2,3-dimethyl-2-butyl cation is used instead.

The literature on the simplest secondary case, the 2-butyl cation, points the same way. Experimentally, its ¹³C NMR spectrum shows no exchange broadening at −140 °C, which puts the barrier to hydrogen scrambling below 2.4 kcal/mol (10 kJ/mol) (Saunders and Kates, *J. Am. Chem. Soc.* 1978, 100, 7082–7083, p. 7083); the scrambling is not frozen out even at −190 °C in the solid (Myhre and Yannoni, *J. Am. Chem. Soc.* 1981, 103, 230–232, not read; cited through the review and Carneiro et al.); and isotope experiments indicate two forms in rapid equilibrium, a hydrogen-bridged form about 0.4 kcal/mol (1.7 kJ/mol) more stable (Saunders and Jiménez-Vázquez, *Chem. Rev.* 1991, 91, 375–397, which also reviews the 2.4 kcal/mol limit). Ab initio, Carneiro, Schleyer, Koch and Raghavachari (*J. Am. Chem. Soc.* 1990, 112, 4064–4066) found for the 2-butyl cation:

- the open classical form is a transition structure (one imaginary frequency), 2.8 kcal/mol (11.7 kJ/mol) above the methyl-bridged form at MP4(sdtq)/6-311G**//MP2(fu)/6-31G** with zero-point energy;
- the methyl-bridged and hydrogen-bridged forms are both minima, within 0.8 kcal/mol (3.3 kJ/mol); with zero-point energy the hydrogen-bridged form is lower by 0.3 kcal/mol (1.3 kJ/mol).

The partially methyl-bridged minimum found here for the 3-methyl-2-butyl cation is the same kind of structure reported for the 2-butyl cation. The bridging angle depends on the method in the same way: Hartree–Fock and B3LYP give a wider angle, MP2 and M06-2X a narrower one.

| Cation | Level | C⁺–C–CH₃ angle of the methyl-bridged form |
|---|---|---|
| 2-butyl (Carneiro et al. 1990) | HF/6-31G* | 102.4° |
| 2-butyl (Carneiro et al. 1990) | MP2(fu)/6-31G* | 77.4° |
| 2-butyl (Carneiro et al. 1990) | MP2/6-31G** | 76.9° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | B3LYP-D3(BJ)/ma-def2-SVP | 81.6° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | M06-2X/ma-def2-SVP | 71.3° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | RI-MP2/ma-def2-SVP | 71.6° |

One difference: in the 2-butyl cation the hydrogen-bridged form is a minimum; in the 3-methyl-2-butyl cation it is not, because the hydride shift from it leads downhill to the tertiary 2-methyl-2-butyl cation. Full references in docs/references.md.

Two follow-ups were asked for:
- whether the partial methyl bridging survives at M06-2X and at a wavefunction level;
- whether a C2–C3 rotation transition state connects this minimum to the tertiary cation.

### IRC stopping rule (decided by J. Picione, 2026-11-09)

Each IRC direction is ended when the energy has changed by less than **0.05 kJ/mol
(0.012 kcal/mol) per step for 20 consecutive steps** after the last bond-forming or
bond-breaking event. A forming or breaking key bond counts as bonded at ≤ 1.25 × the sum of the
covalent radii and as not bonded at ≥ 1.5 ×. The last event is the first point from which every
such bond keeps its final state to the end. Past that point the fragments only drift apart
through the solvent, and ORCA's gradient criteria can take hundreds of steps to be met. The pipeline
runs the two directions as separate jobs, checks the rule every 30 s, cuts the path at the step
where the rule was met, and logs the stop in `run_log.md`. The end complex is then optimized from
that last point (Opt Freq, TightOpt if an imaginary mode remains), and the card states the s value
where the path was cut.

With the median IRC step of about 0.10 bohr·amu^1/2 used for paths computed from 2026-11-11 on,
the rule is about the same for every reaction: an energy change below 0.5 kJ/mol (0.12 kcal/mol)
per bohr·amu^1/2, sustained over about 2.0 bohr·amu^1/2 of path (J. Picione, 2026-11-11; the steps
are not exactly equal, so these are approximate).

First use: E2 (2026-11-09). The product-side IRC was stopped by hand at step 115, before the rule
existed. Products had formed by step 56 (O–H 0.98 Å, C=C 1.34 Å). Each step changed the energy by
less than 0.05 kJ/mol only from step 79 on, and the rule would have stopped that direction at
step 98.

## SN2 pilot

Settings as actually used for the SN2 pilot (`reactions/sn2_hydroxide_bromomethane`). Every
job, its wall time and every rerun are in that folder's `run_log.md`.

## Electronic structure

All calculations used ORCA 6.1.1. The working level, used for all geometries, the transition
state (TS), the intrinsic reaction coordinate (IRC) and the orbitals, was
B3LYP-D3(BJ)/def2-SVP with the RIJCOSX approximation (def2/J auxiliary basis) and TightSCF
convergence, in implicit water: C-PCM with SMD parameters (`! CPCM(water)` with
`%cpcm smd true SMDsolvent "WATER" end`). ORCA's final energies include the SMD CDS term
and the D3(BJ) correction.

Separated reactants and products were built from SMILES (RDKit ETKDG + MMFF) and optimized
with frequency calculations at the working level (bromide: single point).

The TS was optimized (`OptTS`) from a collinear guess (C–O 2.0 Å, C–Br 2.3 Å, planar CH₃,
C–O–H 105°) with the exact Hessian computed at the start and every 5 steps
(`Calc_Hess true`, `Recalc_Hess 5`), followed by a frequency calculation. The TS was
accepted only with exactly one imaginary frequency whose normal mode is dominated by the
central carbon moving along the O···Br axis (share of the mode on C and alignment with the
axis are recorded in `01_ts/ts.json`).

The IRC was followed in both directions from the TS using the TS Hessian (`InitHess read`),
with `Scale_Displ_SD 0.05` to obtain 60 frames (the default step gave 37, with only 7 across
the steep central region). Frames are ordered reactant → product, and the mass-weighted path
length s (bohr·amu^1/2) is accumulated from consecutive frames, with s = 0 at the TS. Both IRC
end geometries were reoptimized (Opt Freq) into the reactant and product complexes; a
minimum showing a small imaginary mode was displaced 0.1 Å along that mode and reoptimized
with TightOpt.

Barrier heights were recomputed by single points on the working-level geometries at
B3LYP-D3(BJ)/def2-TZVP (protocol), B3LYP-D3(BJ)/ma-def2-TZVP, M06-2X-D3(0)/ma-def2-TZVP and
DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO, AutoAux), all in SMD water. Free-energy barriers add
the working-level thermal corrections (298.15 K, 1 atm, ideal gas, rigid rotor/harmonic
oscillator with Grimme's quasi-RRHO vibrational entropy, reference frequency 100 cm⁻¹, ORCA's
default) and the 1 atm → 1 M standard-state correction (RT ln 24.47 =
7.92 kJ/mol per species). The M06-2X and DLPNO-CCSD(T) energies were also computed at every IRC
frame, giving energy curves at those levels along the working-level path. The viewer and the
main video show the DLPNO-CCSD(T) curve by default; the working-level and M06-2X curves are
selectable in the viewer and rendered as separate videos.

## Orbitals

At every IRC frame a working-level single point was followed by intrinsic bond orbital
localization of the occupied valence orbitals (`%loc LocMet IAOIBO end`; ORCA localizes
the orbitals above the frozen core). Coefficients and the AO overlap matrix were exported
with `orca_2json`.

Tracking: orbitals of neighbouring frames k and k+1 are matched by maximizing |O| with
O = C_k^T S C_{k+1}, S = (S_k + S_{k+1})/2 (Hungarian assignment), propagating from the
reactant end. Each match's sign is fixed so that the overlap is positive, so phases never
flip between frames. Matches with |O| < 0.7 are flagged (none occurred for SN2; the two
orbitals shown never fell below 0.99). Each tracked orbital is labelled from its Mulliken
atomic populations (atoms above 0.15 e) at both ends of the path; a two-atom orbital with
less than 0.10 of its population in s-type basis functions is labelled π, otherwise σ, and a
multiple bond whose orbitals are all above 0.10 is labelled as bent bonds (IBO localization
sometimes mixes σ and π into two bent bonds, not necessarily equivalent). The orbitals whose set of
atoms changes between the ends are the ones shown by default: for SN2, an O lone pair →
σ(C–O), and σ(C–Br) → a Br lone pair. The overall sign of each orbital is chosen once, at
the first frame, so that its larger lobe is the + phase. Canonical HOMO and LUMO are tracked
for sign only (their identity is fixed by energy order).

Abrupt switch of the localization (Diels–Alder, tested 2026-11-10): the three shown orbitals
have their lowest frame-to-frame overlap at IRC step 59→60, one to two steps after the TS
(0.827–0.838). IBO single points at ¼ and ½ of that step (Cartesian interpolation between the two
frames, same input as the frames; `03_orbitals/step_test_59_60/`) show that the localized
orbitals switch abruptly within the first quarter of the step: overlap 0.79 from frame 59 to the ¼
point, then 0.996–0.997 on to frame 60. More IRC frames would therefore not spread the change, and
none were added. SN2 and E2 show no such switch: every frame-to-frame overlap of their shown
orbitals is at least 0.99 (lowest 0.993 and 0.990); nor does the cyanide + acetone addition
(lowest 0.992).

Display rotation of bent bonds (decided by J. Picione, 2026-11-10; first use: Diels–Alder): where
IBO localization gives an isolated double bond as two bent bonds (for ethylene at the Diels–Alder
reactant end, two non-equivalent ones, s-character 0.18 and 0.24), the pair is rotated within itself,
for display only, into the σ + π pair that spans the same two orbitals, and the π orbital is shown.
This applies from the reactant end at every frame where the computed orbital is more than 0.5° from
that split (π = the pair's combination with the least s-character), so the rotated frames join the
unrotated ones smoothly (Diels–Alder: IRC frames 0–47, join overlap 0.9998); a rotation inside an
occupied pair changes no energy or density, and the unrotated orbitals stay in the data files
(`03_orbitals/display_rotation.json` records the frames, angles and overlaps).

Arrow direction in pericyclic reactions (decided by J. Picione, 2026-11-10): in a pericyclic
reaction the clockwise and counter-clockwise curved-arrow drawings are mirror images of the
same mechanism, and a mirror-symmetric transition state favours neither. The direction shown
is the one the orbital localization settled on, not a physical preference; `reaction.toml`
lists that assignment, and where `[mirror]` is set the viewer can show the mirror-image set,
made by reflecting the computed orbital surfaces through the TS mirror plane at every frame
(no new calculation).

Cube files (70 points per axis, orca_plot) were contoured at ±0.08 e^1/2 bohr^-3/2 by marching
cubes (scikit-image) and decimated to at most 900 triangles per lobe (fast-simplification).
Vertex positions are stored as 16-bit integers in a common box, with 16-bit triangle indices:
about 15 kB per orbital per frame, 3.7 MB for the whole SN2 page.

## Rendering

Web: three.js r160, static files; ball-and-stick model linearly interpolated between IRC
frames, orbitals shown per nearest frame. The forming and breaking bonds are drawn solid
below 1.15 × their equilibrium length (taken from the path ends) and dashed up to 1.6 ×.

Video: Blender 5.2 (EEVEE, headless Python), 1920 × 1080, 30 fps, 20 s. EEVEE renders on the GPU
(NVIDIA GeForce RTX 5080, OpenGL; E2: 0.33 s per frame), and every render log records the engine,
graphics backend and GPU (`RENDER_DEVICE` line in `05_render/blender.log`). Playback moves at
constant speed in s for the SN2 pilot. From E2 on, playback time is allotted in proportion to the
summed change in the key bond lengths (for E2: Cβ–H, O–H, Cα–Br, Cα–Cβ), with a floor of about
0.3 s per bohr·amu^1/2 of path and the tail after the last bond change compressed to 1 s (decided
by J. Picione, 2026-11-09), so the bond-changing region gets the time and a long drifting-fragment
tail passes quickly. Timing by summed atomic displacement is kept as a selectable alternative
(`video_timing` in `reaction.toml`); the interactive viewer keeps every frame. Start and stop are eased; orbital surfaces cross-fade between
neighbouring IRC frames. The energy inset and labels are composited in Python (Pillow,
matplotlib) and encoded with H.264 (ffmpeg, CRF 18). Rendered frames carry a SHA-256 fingerprint
of everything that decides their look (orbital meshes, geometries, timeline, render settings and
the render script; not the energy curves, which are drawn by the compositor). The renderer
deletes all earlier frames when the fingerprint changes, and the compositor refuses frames whose
fingerprint does not match (added 2026-11-10 after a Diels–Alder render reused frames made before
the orbitals were localized; SN2 and E2 were checked and were not affected).
