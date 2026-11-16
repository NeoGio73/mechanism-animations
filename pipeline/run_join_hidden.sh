#!/usr/bin/env bash
# Join one reaction's production path as its own lane (no window, through pipeline/launch_hidden.py), with the same
# production-log lines as the chain (pipeline/resume_true_irc_prod.sh), so the restart guard reads it as joined.
# First use: the tertiary hydride shift after the continuation test (J. Picione 2026-11-16 afternoon): its open side
# ends where the test runs were stopped; the join (with its end check) need not wait for the secondary shift.
# Usage: pipeline/run_join_hidden.sh <reaction>
R="$(cd "$(dirname "$0")/.." && pwd)"
rx="$1"; L="$R/docs/production_true_irc.log"; D="$R/reactions/$rx/02_irc"
mkdir -p "$D/logs"
echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: join start (J. Picione 2026-11-16: side set after the continuation test)" >> "$L"
( cd "$D/logs" && PYTHONPATH="$R" "$R/.venv-pysis/Scripts/python.exe" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" join \
    > "$D/console_true_irc_join.txt" 2>&1 ); rc=$?
echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: join exit $rc ($(tail -1 "$D/console_true_irc_join.txt" | cut -c1-200))" >> "$L"
