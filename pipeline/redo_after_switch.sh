#!/usr/bin/env bash
# After an end structure is switched to the minimum the path reaches (J. Picione, 2026-11-16): the after-path results
# computed on the old structure are moved (not deleted) to 01_ts/superseded_<label>_<date>/ with a copy of the old
# summary files, and the affected steps rerun in clean folders: single points, diagnostics, barriers, profile.
# Usage: pipeline/redo_after_switch.sh <reaction> <label> [<label> ...]   (label: reactant_complex | product_complex)
R="$(cd "$(dirname "$0")/.." && pwd)"
cd "$R" || exit 1
rx="$1"; shift
X="$R/reactions/$rx"; A="$X/01_ts/superseded_$(date +%Y%m%d)"
mkdir -p "$A/summaries"
for f in 01_ts/sp_tzvp.json 01_ts/diagnostics.json 01_ts/diagnostics.md 02_irc/barriers.json 02_irc/energy_profile.png; do
  [ -f "$X/$f" ] && cp -p "$X/$f" "$A/summaries/"
done
for label in "$@"; do
  [ -d "$X/01_ts/sp_tzvp/$label" ] && mkdir -p "$A/sp_tzvp" && mv "$X/01_ts/sp_tzvp/$label" "$A/sp_tzvp/"
  for d in "$X"/01_ts/diagnostics/*/; do
    lv="$(basename "$d")"
    [ -d "$d$label" ] && mkdir -p "$A/diagnostics/$lv" && mv "$d$label" "$A/diagnostics/$lv/"
  done
done
echo "Superseded end structure(s) $*: results on the old structure, kept on record (J. Picione 2026-11-16: the minimum the path reaches is the official end structure). Not used." > "$A/WHY.txt"
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: start redo after switch ($*); old results in 01_ts/$(basename "$A")" >> "$R/docs/end_check_all.log"
.venv/Scripts/python.exe -X utf8 -m pipeline "$X" singlepoints diagnostics barriers profile > "$X/console_redo_after_switch.txt" 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx: redo after switch exit $?" >> "$R/docs/end_check_all.log"
