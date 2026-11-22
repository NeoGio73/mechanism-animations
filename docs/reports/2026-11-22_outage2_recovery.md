# Second power blip, 2026-11-22 (recovery checked 13:33–13:40)

**Work stopped at 13:01.**
- Windows logged an unexpected shutdown at 13:00:58. The last ORCA output was written at 13:01:45.
- The machine started again at 13:08:37, and the app at 13:17:07 (still version 6.28899.6).
- Everything was restarted at 13:37, so about 36 min of calculation time was lost.

## 1. Jobs that were running (not results)

| Job | Started | What I did |
|---|---|---|
| Tertiary hydride shift, DLPNO-CCSD(T) frame 258 | 12:53:18 | All files moved to `f258/outage_20261122_1301_f/` (not deleted); rerun started 13:37 in the clean folder |
| E2, DLPNO-CCSD(T) frame 108 | 13:00:59 | All files moved to `f108/outage_20261122_1301_f/`; queued to rerun in the clean folder |

- **Every other path single-point output is complete** (all reactions, all levels). Each ends with "ORCA TERMINATED NORMALLY", and its energy matches the run log and, where the level is finished, the CSV.
  - E2 DLPNO-CCSD(T): 99 done.
  - Tertiary DLPNO-CCSD(T): 231 done.
- No other ORCA output written today is incomplete.

## 2. Damaged files: none

- **Text files:** of those written today (285 checked), none is zero bytes or NUL-filled. No file of any kind written today is zero bytes.
- **JSON:** every JSON file in the project parses.
- **Core scheduler state (`allocation.json`):** intact this time. All four of its entries belonged to processes that ended in the reboot. I cleared them before restarting, so that a reused process ID could not hold cores. The old copy is in `docs/outage_2026-11-22b/`.
- **Git:** `git fsck --full` is clean, the local copy matches GitHub, and the working-tree changes are the same as before.

## 3. Restarted (13:37–13:38, all hidden through WMI)

- **foreign_orca check:** no ORCA process from outside this project.
- **Lanes:** the tertiary hydride shift and E2 single points.
- **Waiting jobs:**
  - The explicit-water run is waiting behind the path single points.
  - The rotation-TS job for the secondary shift had been killed by the reboot, before any job started. It is restarted with the same start condition: after the tertiary lane's last single point.
- **Helpers:** stall watchdog; keep-awake (until Wed 2026-12-02 12:00).
- **Restart guard:** started from its script file at 13:38 (pid 28280, until Wed 2026-12-02 12:00). It watches the watchdog and both lanes.
- **Report schedule:** both jobs are still scheduled (next report 16:46 today; renewal Fri 2026-11-27 12:13).
- **Late-report rule:** a report that comes due while I am busy is written as soon as the task ends, marked late. I check for an overdue report at the end of every task. This rule is saved in my standing notes.
- **Settings:** I changed no system, power or McAfee setting.

## 4. Revised finish times

The E2 rows include about 1.5 h for the rotation-TS jobs.

| Item | Expected |
|---|---|
| Tertiary hydride shift, last single point (14 left) | About Sun 2026-11-22 18:00 (range 17:10 to 19:20) |
| Rotation-TS jobs (scans first) | Start then; done about Sun 21:00 (late case about Sun 22:30), for Monday's 07:00 report |
| **Last path single point (E2, 408 left)** | **About Thu 2026-11-26 04:20 (range Wed 2026-11-25 11:20 to Fri 2026-11-27 05:55)** |
| Cyanide explicit-water test | Right after that, about Thu 2026-11-26 05:30 |

- E2's finish is about the same as this morning's estimate (Thu 04:45).
- If E2 keeps its measured 0.97 ×, it finishes about Wed 2026-11-25 18:15.
