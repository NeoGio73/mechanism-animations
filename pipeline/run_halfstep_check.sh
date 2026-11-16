#!/usr/bin/env bash
# Tertiary hydride shift, open side at half the step (J. Picione, 2026-11-15): a check run, never a production path.
# Same settings as production except step 0.025 bohr·amu½ and 80 steps; output in 02_irc/attempts/halfstep_check/backward.
# Started with no window through pipeline/launch_hidden.py. Afterwards: python -m pipeline.halfstep_compare.
R="$(cd "$(dirname "$0")/.." && pwd)"
PY="$R/.venv-pysis/Scripts/python.exe"
RX="$R/reactions/hshift_2_3_dimethyl_2_butyl_m062x"
D="$RX/02_irc"
mkdir -p "$D/logs/halfstep_backward"
cd "$D/logs/halfstep_backward" || exit 1
PYTHONPATH="$R" "$PY" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$RX" backward step=0.025 check_dir=attempts/halfstep_check max_cycles=80 \
    > "$D/console_halfstep_check_backward.txt" 2>&1
rc=$?
echo "end $(date +%H:%M) exit $rc" >> "$D/console_halfstep_check_backward.txt"
echo "$(date '+%Y-%m-%d %H:%M:%S') half-step check hshift_2_3_dimethyl_2_butyl_m062x backward: exit $rc" >> "$R/docs/production_true_irc.log"
