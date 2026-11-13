#!/usr/bin/env bash
# Launch pipeline/true_irc2.py runs, each from its own log folder (pysisyphus writes irc.log to the
# working directory at import). Usage: pipeline/run_true_irc2.sh <run> <direction> [step] [fifo]
R="$(cd "$(dirname "$0")/.." && pwd)"
T="$R/reactions/sn2_hydroxide_bromomethane_masvp/02_irc/attempts/step_convergence/true_irc2"
L="$T/logs/$1_$2"; mkdir -p "$L"
cd "$L" && PYTHONPATH="$R" "$R/.venv-pysis/Scripts/python.exe" -X utf8 -W ignore -u -m pipeline.true_irc2 "$@" > "$T/console_$1_$2.txt" 2>&1
echo "end $(date +%H:%M) exit $?" >> "$T/console_$1_$2.txt"
