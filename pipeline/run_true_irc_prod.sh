#!/usr/bin/env bash
# Production true IRCs at M06-2X in J. Picione's order (2026-11-13): SN2, cyanide + acetone,
# Diels–Alder, tertiary hydride shift, secondary hydride shift. Both sides of a reaction run at the
# same time (through the core queue), then the path is joined. A reaction whose side does not end by
# the stopping rule is logged as incomplete (not a result) and the chain goes on to the next one.
# Started only after the SN2 validation path passed. Usage: pipeline/run_true_irc_prod.sh [reaction ...]
R="$(cd "$(dirname "$0")/.." && pwd)"
PY="$R/.venv-pysis/Scripts/python.exe"
LIST=("$@")
[ ${#LIST[@]} -eq 0 ] && LIST=(sn2_hydroxide_bromomethane_m062x nucadd_cyanide_acetone_m062x da_butadiene_ethylene_m062x
                               hshift_2_3_dimethyl_2_butyl_m062x hshift_3_methyl_2_butyl_m062x)
for rx in "${LIST[@]}"; do
  D="$R/reactions/$rx/02_irc"; mkdir -p "$D/logs"
  echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: start" >> "$R/docs/production_true_irc.log"
  for d in forward backward; do
    mkdir -p "$D/logs/$d"
    ( cd "$D/logs/$d" && PYTHONPATH="$R" "$PY" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" $d \
        > "$D/console_true_irc_$d.txt" 2>&1; rc=$?; echo "end $(date +%H:%M) exit $rc" >> "$D/console_true_irc_$d.txt" ) &
    sleep 5
  done
  wait
  ( cd "$D/logs" && PYTHONPATH="$R" "$PY" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$rx" join \
      > "$D/console_true_irc_join.txt" 2>&1 ); rc=$?
  echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $rx: join exit $rc ($(tail -1 "$D/console_true_irc_join.txt" | cut -c1-200))" \
      >> "$R/docs/production_true_irc.log"
done
