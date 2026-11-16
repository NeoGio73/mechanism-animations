#!/usr/bin/env bash
# J. Picione, 2026-11-16: the minimum a path reaches becomes the official end structure for every path, finished ones
# included. Sides of finished paths not yet optimized from their last frame are optimized now (02_irc/true_irc/
# end_check_all/<side>, Opt Freq as their reference minimum, through the core queue); a different form with the same
# bonds is switched in (pipeline/path_minimum.py); different bonds: the lane stops and nothing is switched (docs/END_BONDS_STOP).
# Already optimized: E2 both sides, cyanide backward, tertiary bridged side. Secondary and tertiary open sides: at their rejoin.
# Started with no window through pipeline/launch_hidden.py.
R="$(cd "$(dirname "$0")/.." && pwd)"
cd "$R" || exit 1
L="$R/docs/end_check_all.log"
for job in "sn2_hydroxide_bromomethane_m062x backward,forward" "da_butadiene_ethylene_m062x backward,forward" \
           "nucadd_cyanide_acetone_m062x forward"; do
  set -- $job
  if [ -f "$R/docs/END_BONDS_STOP" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') $1: not started, docs/END_BONDS_STOP exists (waiting for J. Picione)" >> "$L"; exit 2
  fi
  echo "$(date '+%Y-%m-%d %H:%M:%S') $1 $2: start" >> "$L"
  .venv/Scripts/python.exe -X utf8 -m pipeline.end_check "$1" sides=$2 sub=end_check_all \
      > "$R/reactions/$1/02_irc/console_end_check_all.txt" 2>&1
  echo "$(date '+%Y-%m-%d %H:%M:%S') $1 $2: exit $?" >> "$L"
done
