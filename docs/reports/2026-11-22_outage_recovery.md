# Power-outage recovery, 2026-11-22 (checked 05:22–05:40)

**Work stopped at 01:29.**
- Windows logged an unexpected shutdown at 01:29:07. The last ORCA output was written at 01:29:23.
- The machine started again at 05:15:33.
- The Claude app restarted at 05:19:53. On that restart it installed the pending update (6.28899.6, at 05:19:44).
- Calculation time lost: about 4 h (01:29 to 05:31).

## 1. Jobs that were running (not results)

| Job | Started | What I did |
|---|---|---|
| Tertiary hydride shift, DLPNO-CCSD(T) frame 222 | 01:24:30 | All files moved to `f222/outage_20261122_0129_f/` (not deleted); rerun started 05:31 in the clean folder |
| Tertiary hydride shift, DLPNO-CCSD(T) frame 223 | 01:25:25 | All files moved to `f223/outage_20261122_0129_f/`; queued, runs next in the clean folder |

- E2 had no job running; it was waiting for cores.
- The secondary hydride shift had already finished at 22:20 yesterday.
- **Every other path single-point output was checked (all reactions, all levels).** Each one ends with "ORCA TERMINATED NORMALLY", and its energy matches the run log and, where the level is finished, the CSV.
  - E2 DLPNO-CCSD(T): 74 done.
  - Tertiary DLPNO-CCSD(T): 195 done.
  - Secondary DLPNO-CCSD(T): 309 of 309 done.
- No other ORCA output written since yesterday noon is incomplete.

## 2. Files cut short by the outage

- **One file was damaged: the core scheduler's state (`docs/level_test_queue/allocation.json`).** It was 2249 bytes of zeros.
  - **Rebuilt:** every job and waiter it listed was a process that ended in the outage, and no lane was paused. So the correct state is empty, and I wrote it that way.
  - The damaged copy is kept in `docs/outage_2026-11-22/`. The rebuild is logged in `allocation.log`.
- **All other files are intact:**
  - Of the text files written since 2026-11-20 12:00 (1602 checked), none has zero bytes.
  - Every JSON file in the project parses.
  - The energy CSVs and run logs end on complete lines.
- **Git:** `git fsck --full` is clean (only the usual dangling blobs). HEAD equals GitHub (`3c343fe`), and the working-tree changes are the same as before the outage.

## 3. Secondary hydride shift (finished before the outage)

- **Zero check run** (`pipeline.unimolecular_zero`; nothing recomputed, E_eh unchanged).
  - The DLPNO-CCSD(T) CSV had been written by a lane started on 2026-11-16, before your zero rule of 2026-11-19. So its kJ/mol columns used 0 Eh as the zero.
  - They are now relative to `secondary_anti`, the lowest-free-energy conformer. Frame 0 lies +3.16 kJ/mol (+0.76 kcal/mol) above it.
- **Results summary written:** `docs/results/hshift_3_methyl_2_butyl.md`. At DLPNO-CCSD(T)/ma-def2-TZVP, zero = `secondary_anti`:

| Quantity | Value |
|---|---|
| ΔG‡ | 14.0 kJ/mol (3.3 kcal/mol) |
| ΔE‡ | 17.3 kJ/mol (4.1 kcal/mol) |
| ΔrG | −45.4 kJ/mol (−10.9 kcal/mol) |
| Highest path frame | frame 215, 0.46 kJ/mol (0.11 kcal/mol) above the TS frame, C2–H 1.943 Å |
| Largest T1 | 0.0135 |

- **B3LYP is out of line:** ΔG‡ is 0.4 kJ/mol (0.1 kcal/mol), against 17.2 kJ/mol at M06-2X/ma-def2-TZVP.
- The card notes have no numbers that needed fixing.

## 4. Restarted (05:31–05:33, all hidden through WMI)

- **foreign_orca check:** no ORCA process from outside this project.
- **Lanes:** the single points for the tertiary hydride shift and for E2. The explicit-water run is waiting behind them again.
- **Stall watchdog:** restarted.
- **Keep-awake:** now in `pipeline/keep_awake.ps1`, the same request as before, ending Wed 2026-12-02 12:00.
- **Restart guard:** started at 05:33 (pid 27588, until Wed 2026-12-02 12:00).
  - My first start of the guard failed: its encoded command is now over Windows' length limit. I started it from the script file instead, and noted this in the guard's header.
- **Report schedule:** both jobs are still scheduled. The next report runs at 06:39 today, and the renewal is Fri 2026-11-27 12:13.
- I changed no system, power or McAfee setting.

## 5. Revised finish times

Same ratios as yesterday: tertiary 1.04–1.06 ×, E2 1.08 ×.

| Item | Expected |
|---|---|
| Tertiary hydride shift, last single point (50 left) | About Sun 2026-11-22 21:10 (range Sun 18:00 to Mon 2026-11-23 01:50) |
| **Last path single point (E2, 433 left)** | **About Thu 2026-11-26 04:00 (range Wed 2026-11-25 09:05 to Fri 2026-11-27 08:15)** |
| Cyanide explicit-water test | Right after that, about Thu 2026-11-26 05:00 |

- The outage moved the E2 finish about 2 h later than yesterday's 17:00 estimate (Thu 01:50).
- If E2 stays at its measured 0.94 ×, it finishes about Wed 2026-11-25 12:15.
