#!/usr/bin/env bash
# Single points on a final path at low priority (J. Picione 2026-11-16: "keep free cores busy with single points on
# paths that are final"; CORES_LOW=1, pipeline/cores._acquire_low): the same after-path steps as SN2 (2026-11-14).
# Waits first for the reaction's end check of 2026-11-16 (docs/end_check_all.log) and stops if docs/END_BONDS_STOP
# names the reaction. Started with no window through pipeline/launch_hidden.py.
# Usage: pipeline/run_after_path_low.sh <reaction>
R="$(cd "$(dirname "$0")/.." && pwd)"
cd "$R" || exit 1
rx="$1"; L="$R/docs/end_check_all.log"
# "nowait": the reaction's end check ran at its own join (E2, 2026-11-16 03:39; end_check_decision.md)
[ "$2" = "nowait" ] || until grep -q "$rx [a-z,]*: exit" "$L"; do sleep 60; done
if [ -f docs/END_BONDS_STOP ] && grep -q "$rx" docs/END_BONDS_STOP; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points not started (END_BONDS_STOP)" >> "$L"; exit 2
fi
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points start (low priority)" >> "$L"
CORES_LOW=1 .venv/Scripts/python.exe -X utf8 -m pipeline "reactions/$rx" from_stationary fragments singlepoints diagnostics barriers profile path_singlepoints \
    > "reactions/$rx/console_after_path.txt" 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: after-path single points exit $?" >> "$L"
