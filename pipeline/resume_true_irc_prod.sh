#!/usr/bin/env bash
# Resume the production chain after a stop (used by the restart guard, pipeline/restart_guard.ps1; J. Picione
# approved 2026-11-15). Same steps and settings as run_true_irc_prod.sh, with two differences:
#  - a reaction whose latest "join exit" or "join withdrawn" line in docs/production_true_irc.log is a join exit is skipped;
#  - a side that already finished (true_irc/<side>/path.json) is kept, not rerun.
# A killed side must be archived first (the guard moves it to 02_irc/attempts/killed_<time>/), because
# pipeline.true_irc_prod starts a side only in a new, clean folder.
# Usage: pipeline/resume_true_irc_prod.sh reaction [reaction ...]
R="$(cd "$(dirname "$0")/.." && pwd)"
PY="$R/.venv-pysis/Scripts/python.exe"
L="$R/docs/production_true_irc.log"
for rx in "$@"; do
  # joined = its latest "join exit" / "join withdrawn" line is a join exit (a withdrawn join is rerun, 2026-11-15)
  last="$(grep -E "production true IRC $rx: join (exit|withdrawn)" "$L" | tail -1)"
  case "$last" in *"join exit"*) continue ;; esac
  D="$R/reactions/$rx/02_irc"; mkdir -p "$D/logs"
  echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: start (resumed by the restart guard)" >> "$L"
  for d in forward backward; do
    if [ -f "$D/true_irc/$d/path.json" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: $d side finished earlier, kept" >> "$L"; continue
    fi
    mkdir -p "$D/logs/$d"
    ( cd "$D/logs/$d" && PYTHONPATH="$R" "$PY" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" $d \
        > "$D/console_true_irc_$d.txt" 2>&1; rc=$?; echo "end $(date +%H:%M) exit $rc" >> "$D/console_true_irc_$d.txt" ) &
    sleep 5
  done
  wait
  ( cd "$D/logs" && PYTHONPATH="$R" "$PY" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" join \
      > "$D/console_true_irc_join.txt" 2>&1 ); rc=$?
  echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: join exit $rc ($(tail -1 "$D/console_true_irc_join.txt" | cut -c1-200))" \
      >> "$L"
done
