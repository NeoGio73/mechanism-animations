#!/usr/bin/env bash
# After-path single points at normal priority on all 8 cores (J. Picione 2026-11-16 evening: "now that no path is
# running, the single points no longer need low priority. Let them use all 8 cores at normal priority, within the
# memory limit"). CORES_FIT=1: first come first served in the core queue; a job runs on the cores free for it
# (pipeline/cores._grant; only ORCA's parallel core count changes, never the method). Finished jobs are reused.
# Replaces pipeline/run_after_path_low.sh for the lanes switched by pipeline/sp_switch.py. Started with no window.
# Usage: pipeline/run_singlepoints.sh <reaction>
R="$(cd "$(dirname "$0")/.." && pwd)"
cd "$R" || exit 1
rx="$1"; L="$R/docs/end_check_all.log"
if [ -f docs/END_BONDS_STOP ] && grep -q "$rx" docs/END_BONDS_STOP; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points not started (END_BONDS_STOP)" >> "$L"; exit 2
fi
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points start (normal priority, all 8 cores)" >> "$L"
CORES_FIT=1 .venv/Scripts/python.exe -X utf8 -m pipeline "reactions/$rx" from_stationary fragments singlepoints diagnostics barriers profile path_singlepoints \
    >> "reactions/$rx/console_after_path.txt" 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points exit $?" >> "$L"
