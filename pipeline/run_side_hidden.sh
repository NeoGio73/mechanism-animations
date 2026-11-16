#!/usr/bin/env bash
# One production IRC side (or check/continuation run) as its own lane, started with no window through
# pipeline/launch_hidden.py (J. Picione 2026-11-15/16).
# Usage: pipeline/run_side_hidden.sh <reaction> <direction> <log tag> [key=value ...]   (passed to pipeline.true_irc_prod)
R="$(cd "$(dirname "$0")/.." && pwd)"
rx="$1"; d="$2"; tag="$3"; shift 3
D="$R/reactions/$rx/02_irc"
mkdir -p "$D/logs/$tag"
cd "$D/logs/$tag" || exit 1
PYTHONPATH="$R" "$R/.venv-pysis/Scripts/python.exe" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" "$d" "$@" \
    > "$D/console_$tag.txt" 2>&1
rc=$?
echo "end $(date +%H:%M) exit $rc" >> "$D/console_$tag.txt"
echo "$(date '+%Y-%m-%d %H:%M:%S') $rx $d ($tag): exit $rc" >> "$R/docs/side_lanes.log"
