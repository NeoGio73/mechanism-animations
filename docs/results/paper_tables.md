# Tables for the manuscript: wall times, integral settings, path checks

These tables were collected on 2026-11-27 from the run logs and the outputs on disk; no new calculation was run. Each
cell names its source in the notes below its table.

## 1. Wall time per reaction

**Hardware:** one workstation.
- **CPU:** AMD Ryzen 7 9800X3D, 8 cores / 16 threads.
- **RAM:** 31.1 GB.
- **GPU (rendering only):** NVIDIA GeForce RTX 5080.
- **Software:** Windows 11 Pro, ORCA 6.1.1, pysisyphus 1.0.0 (IRC driver), Blender 5.2 (EEVEE).

**How the times are counted:**
- Times are summed job wall times. Each ORCA job ran on 4 cores (`%pal nprocs 4`) unless noted, and two jobs usually
  ran at the same time, so the sums are not elapsed calendar time.
- h = hours. Frame counts are given where a level runs over path frames.

| Reaction | TS (OptTS + Freq) | Path (IRC gradients + end-point Opt Freq) | Single points, M06-2X/ma-def2-TZVP | Single points, B3LYP-D3(BJ)/ma-def2-TZVP | Single points, DLPNO-CCSD(T)/ma-def2-TZVP | Orbitals (IBO single points, M06-2X/ma-def2-SVP) | Rendering (viewer meshes + 4 videos, GPU) |
|---|---|---|---|---|---|---|---|
| SN2, HO⁻ + CH₃Br | 0.06 h | 8.21 h (812 jobs) | 3.33 h, 342 frames | 2.83 h, 342 frames | 13.74 h, 328 frames | 0.70 h, 342 frames | 9.6 min |
| E2, EtO⁻ + 2-bromopropane | 0.35 h | 18.86 h (1239 jobs) | 9.78 h, 532 frames | 6.08 h, 532 frames | 149.61 h, 507 frames | 4.48 h, 532 frames | 15.9 min |
| CN⁻ + acetone | 0.15 h | 10.11 h (659 jobs) | 2.41 h, 253 frames | 1.18 h, 253 frames | 26.31 h, 227 frames | 1.07 h, 253 frames | 7.2 min |
| Tertiary 1,2-H shift | 0.59 h | 26.08 h (1891 jobs); earlier attempts 8.76 h (785 jobs) | 8.85 h, 543 frames | 2.30 h, 272 computed + 271 mirror | 46.56 h, 245 computed + 244 mirror | 3.64 h, 543 frames | 10.3 min |
| Secondary 1,2-H shift | 0.23 h | 11.34 h (1389 jobs); separate check 6.11 h (946 jobs), not counted | 4.11 h, 407 frames | 1.88 h, 407 frames | 29.21 h, 309 frames | 1.86 h, 407 frames | 12.3 min + 4.3 min re-render |
| Diels–Alder | 0.24 h | 13.65 h (759 jobs) | 3.80 h, 302 frames | 1.79 h, 302 frames | 61.21 h, 275 frames | 1.53 h, 302 frames | 10.0 min |

**Notes, column by column:**
- **TS:** the production OptTS Freq job `01_stationary/ts/ts` in each reaction's `run_log.md`. It started from the
  level-test TS. The earlier TS searches (level tests, pilot folders) are not counted.
- **Path:**
  - All `02_irc/true_irc` job rows in `run_log.md`: the EnGrad steps of the IRC plus the Opt Freq end checks
    (0.04–1.28 h of each total).
  - Tertiary shift: the rerun attempts in `02_irc/attempts` are listed separately and were not used for the path.
  - Secondary shift, separate check (J. Picione, 2026-11-28): the IRC from the anti form's rotation TS
    (`02_irc/ts_anti_rotation_irc`, 2026-11-27; 944 gradient jobs and one Opt Freq at each end, 4 cores), 6.11 h of
    job time. It is not part of the path shown and is not in the path column or the total.
- **Single points:**
  - Source: `02_irc/path_sp/run_log_path_sp.md`; frame counts from `02_irc/energies_<level>.json` (`n_computed`).
  - The "TOTAL RUN TIME" lines in the ORCA outputs sum 0.3–4 % lower.
  - DLPNO-CCSD(T) covers only the frames J. Picione sampled.
  - Tertiary shift: the mirror half of the degenerate path takes the energies of its computed partner frames.
  - One failed DLPNO-CCSD(T) attempt each for SN2 (frame 302) and the tertiary shift (frame 71) is not counted. With
    them the totals are 13.77 h and 46.61 h.
  - **Cores per job, where not 4:**
    - E2: 21 M06-2X jobs at 2 cores and 1 at 3; 8 B3LYP jobs at 3 cores and 1 at 2.
    - Cyanide: 1 DLPNO-CCSD(T) job at 3 cores.
    - Tertiary shift: 2 DLPNO-CCSD(T) jobs at 8 cores.
    - Secondary shift: 8 M06-2X jobs at 2 cores; 10 DLPNO-CCSD(T) jobs at 3.
    - Diels–Alder: 1 B3LYP job at 3 cores; 5 DLPNO-CCSD(T) jobs at 3 and 1 at 2.
- **Orbitals:** `03_orbitals/run_log_frames.md`. All IBO jobs ran on 4 cores, on 2026-11-25. Tracking, cube files
  and meshes took seconds to minutes of CPU and are not counted.
- **Rendering:**
  - Elapsed time of the final "meshes + video" run in `docs/console_renders.txt`: the viewer meshes, Blender frames
    (600 frames at 1920 × 1080) and four mp4 files, one per energy curve.
  - Secondary shift: the 2026-11-25 run (12.3 min), plus the video-only re-render of 2026-11-27 for the new orbital
    labels (4.3 min).
  - Diels–Alder: the second run of 2026-11-25 (mirror display frames); the first run took 9.9 min.
- **Totals:** summed over the columns above (rendering, the tertiary attempts and the secondary separate check not counted), the six reactions took
  28.9, 189.2, 41.2, 88.0, 48.6 and 82.2 h of job time, in table order. DLPNO-CCSD(T) is 48–79 % of each total;
  E2's 149.6 h is the largest single item.

## 2. Integral approximations and auxiliary basis sets, as in the inputs

Every input also has `%pal nprocs 4 end` (except as noted in Table 1), `%maxcore 2000` and
`%cpcm smd true SMDsolvent "<solvent>" end`, with `CPCM(<solvent>)` on the keyword line. There is no `%method`,
`%cosx` or `%basis` block in any production input, so COSX uses ORCA's default grids for the chosen DefGrid.

| Level (use) | Keyword line, as in the inputs | Auxiliary basis (from the outputs) | Integration grid |
|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (geometries, TS, IRC gradients, frequencies, IBO single points) | `! M062X D3ZERO ma-def2-SVP def2/J RIJCOSX VeryTightSCF DefGrid3 CPCM(<solvent>)`, with `EnGrad NoPop` for IRC steps, `OptTS Freq` / `Opt Freq` for stationary points, and `%loc LocMet IAOIBO end` for the IBO single points | AuxJ def2/J | DefGrid3 (Lebedev-590) |
| M06-2X-D3(0)/ma-def2-TZVP (path single points) | `! M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX VeryTightSCF DefGrid3 CPCM(<solvent>)` | AuxJ def2/J | DefGrid3 (Lebedev-590) |
| B3LYP-D3(BJ)/ma-def2-TZVP (path single points) | `! B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX VeryTightSCF CPCM(<solvent>)` | AuxJ def2/J | no grid keyword: ORCA default (Lebedev-302 in the output) |
| DLPNO-CCSD(T)/ma-def2-TZVP (path single points) | `! DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX VeryTightSCF CPCM(<solvent>)` | AutoAux for AuxJ, AuxJK and AuxC | no grid keyword (the COSX grid of the HF reference) |

**Notes:**
- The solvents are water (SN2, cyanide), ethanol (E2), toluene (Diels–Alder) and dichloromethane (both hydride
  shifts).
- The gas-phase twins of the SMD-scheme check use the same lines without `CPCM(...)` (`01_ts/smd_scheme`).
- Source: the `.inp` files of every production job, with the explicit-water test and the discarded attempts left
  out. Each keyword line is the same in all six reactions apart from the solvent.
- The auxiliary-basis and grid lines were read from one output per level (`02_irc/path_sp/<level>/f*/f.out`,
  cyanide + acetone).
- B3LYP was run without `DefGrid3`, unlike M06-2X; this is recorded here as found. M06-2X is the functional that is
  sensitive to the grid; ORCA prints that warning for it.

## 3. Path checks of the six production paths

**Angle check** (`pipeline/path_angles.py`):
- **What it measures:** at each interior frame, the angle between the path step (mass-weighted central chord) and the
  downhill direction −M^−½g.
- **Pass rule:** frames are split into two zones, |s| ≤ 1 bohr·amu½ and beyond. In each zone the median must be below
  5° and the 90th percentile below 15°.
- **Frames left out:** trailing frames at the end minimum (|ΔE| < 0.005 kJ/mol and RMS gradient < 5 × 10⁻⁵ Eh/bohr).

**End-point agreement** (`pipeline/end_check.py`):
- **What it compares:** each path end against the optimized minimum, both directly and through an Opt Freq started
  from the last frame.
- **Same-minimum test:** |ΔE| ≤ 0.1 kJ/mol, no imaginary frequency, heavy-atom RMSD ≤ 0.2 Å, and each key distance
  within 0.05 Å.
- **When the test is not met:** the minimum the path actually reaches, with the same bonds, became the official end
  structure (decisions of J. Picione, 2026-11-16).
- **Energies:** electronic, at M06-2X-D3(0)/ma-def2-SVP with SMD.
- **Stopping rule:** |ΔE| < 0.05 kJ/mol per step for 20 consecutive steps after the last bond event.

| Reaction | Angle check | Largest angle, backward / forward side (frame from the TS) | Reactant end: last frame vs minimum; Opt Freq from it | Product end: last frame vs minimum; Opt Freq from it |
|---|---|---|---|---|
| SN2 | passes | 20.4° (frame 152, beyond zone) / 5.7° (frame 1) | reactant complex: +0.10 kJ/mol, 0.08 Å; opt +0.04 kJ/mol, 0.06 Å heavy atoms. The C–O distance differs by 0.12 Å, so the optimized structure became the official end | product complex: +0.23 kJ/mol, 0.03 Å; opt +0.01 kJ/mol, 0.002 Å: same minimum |
| E2 | passes | 8.5° (frame 129) / 2.6° (frame 264) | reactant complex: +0.90 kJ/mol, 0.32 Å; opt −0.27 kJ/mol, 0.23 Å. The official end is the reached minimum; the last frame is +1.17 kJ/mol and 0.20 Å from it | product complex: +4.73 kJ/mol, 0.53 Å; opt −0.28 kJ/mol, 0.15 Å. The official end is the reached minimum; the last frame is +5.01 kJ/mol and 0.40 Å from it |
| CN⁻ + acetone | passes | 11.4° (frame 157) / 15.1° (frame 76, beyond zone) | reactant complex: +0.94 kJ/mol, 0.25 Å; opt −0.07 kJ/mol, 0.11 Å: same minimum | product complex: −0.03 kJ/mol, 0.003 Å; opt −0.03 kJ/mol, 0.003 Å: same minimum |
| Tertiary 1,2-H shift | flagged: the near-TS zone fails (median / 90th percentile 9.6° / 39.3° on the open side, 21 frames; 40.2° / 65.0° on the bridged side, steps 1–6, 75.4° down to 26.2°); beyond \|s\| = 1 the open side fails the median narrowly, 5.74° / 10.18° (242 frames). The bridged side has no frame beyond \|s\| = 1 on the path, and the product half is the mirror image of the open side, so neither has an entry of its own | 54.2° (frame 1) / 75.4° (frame 1) | open cation: +1.47 kJ/mol, 0.26 Å; opt −0.06 kJ/mol, 0.009 Å: same minimum | bridged form: −0.003 kJ/mol, 0.004 Å; opt −0.002 kJ/mol, 0.001 Å: same minimum |
| Secondary 1,2-H shift | flagged: the near-TS zone fails (7.6° / 25.4° backward, 10.6° / 25.4° forward) | 61.9° (frame 1) / 66.9° (frame 1) | eclipsed form: −0.10 kJ/mol, 0.23 Å; opt −0.10 kJ/mol, 0.06 Å. Outside the energy limit, so the optimized structure became the official end | tertiary cation: +0.19 kJ/mol, 0.43 Å; opt +0.18 kJ/mol, 0.13 Å. The optimized structure became the official end |
| Diels–Alder | passes | 17.2° (frame 177) / 39.8° (frame 119) | reactant complex: +0.33 kJ/mol, 0.18 Å; opt −0.30 kJ/mol, 0.25 Å. The optimized structure became the official end | product (cyclohexene): +0.01 kJ/mol, 0.002 Å; opt +0.004 kJ/mol, 0.001 Å: same minimum |

**Notes:**
- **Sources:**
  - Angles: `02_irc/true_irc/angles.md/.json` (tertiary near-TS detail in `angles_ts_region.md`).
  - End checks: `02_irc/true_irc/end_check.md`, `end_check_all.md` and `end_check_decision.md` (E2), with
    `02_irc/endpoints.json` and `02_irc/irc_meta.json`.
  - RMSDs: in the "last frame vs minimum" column, all atoms; in the "opt" column, heavy atoms.
- **The two flagged paths are the hydride shifts.** Their TS regions are very flat: 115i and 92i cm⁻¹. In the SI
  analysis, all 70 frames of the six paths with an angle above 15° have an RMS gradient of 2.0 × 10⁻⁴ Eh/bohr or less
  (`docs/si_figures/angle_vs_gradient.md`): 30 lie near the flat hydride-shift TSs and 40 in the final approach to a
  minimum. The largest angles sit at the first frame from these TSs, where the gradient is near zero.
- **Frames the angle check counts:** every interior frame of each side (the last frame of a side has no angle,
  since the angle uses the steps on both sides of the frame), split at |s| = 1 bohr·amu½, except trailing frames at
  the end minimum (|ΔE| < 0.005 kJ/mol from the previous frame and RMS gradient < 5 × 10⁻⁵ Eh/bohr; rule of
  2026-11-14). Only two sides have such frames:
  - Cyanide product side: frames 77–93 from the TS (17 frames, angles up to 117.0° at frame 77) are left out; the
    counted zone beyond |s| = 1 is frames 21–76 (56 frames), largest 15.1° (frame 76). The 117° in SI entry 6 comes
    from the angle-against-gradient plot, which shows every frame, these 17 included.
  - Diels–Alder product side: frames 120–121 from the TS (2 frames, up to 89.1°) are left out.
  - The tertiary backward side has none left out: its 242 frames beyond |s| = 1 run to the side's last interior frame
    and include all 151 continued frames.
- **Tertiary, frames on the shown path only (J. Picione, 2026-11-27):** the bridged side of the path is computed
  forward steps 0–7 (the cut at the bridged structure, `pipeline/mirror_path.py`). The row counts its steps 1–6; step
  7's angle needs step 8 and is not counted. Angles from `02_irc/true_irc/angles.json` (forward frames 1–6).
  - Forward steps 8–129 are kept on disk (`02_irc/true_irc/forward/path.xyz`) and are not part of the path. Past the
    bridged structure they run back down toward an open cation with the hydrogen on C3, the same carbon as at the
    reactant end (step 129: C2–H 1.93 Å, C3–H 1.125 Å, 4.08 kJ/mol above the open-cation minimum). The side ended by
    the stopping rule at step 129, not at a minimum, and no optimization was run from it.
- **Exceptions recorded in `irc_meta.json`:**
  - Tertiary, backward side: not ended by the stopping rule. It was accepted where the comparison runs were stopped
    (2026-11-16).
  - Diels–Alder, forward side: ended on "energy increased" and accepted as complete (2026-11-14).
- **Tertiary product end (closed, J. Picione 2026-11-27):**
  - The two files do not conflict. The computed half of the degenerate path ends at the bridged form (7 frames on
    that side, frames 0–6), and `end_check.md` tests that end; this is the comparison in the table.
  - The open cation in `endpoints.json` (0.71 Å from the last frame) is the end of the full degenerate path, reached
    through the mirror half: the far end is the mirror image of the reactant end by construction. `endpoints.json`
    now says so in a note field; no value in it changed.
  - The manuscript quotes the end check as the table shows it.
- **Diels–Alder reactant end:** the first Opt Freq had an imaginary frequency (92.7i cm⁻¹) and was redone with
  TightOpt (`run_log.md`). `endpoints.json` gives 0.60 Å from the last frame to the official end structure.
