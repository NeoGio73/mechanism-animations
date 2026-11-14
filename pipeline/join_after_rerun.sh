#!/usr/bin/env bash
# Join a production path whose side was rerun outside the chain (2026-11-14: tertiary hydride shift forward,
# after the energy-increase stop fix). Waits for both sides; skips if the chain already joined it successfully.
R="$(cd "$(dirname "$0")/.." && pwd)"; RX="$1"; D="$R/reactions/$RX/02_irc"; L="$R/docs/production_true_irc.log"
until [ -f "$D/true_irc/forward/path.json" ] && [ -f "$D/true_irc/backward/path.json" ]; do sleep 120; done
if grep -q "$RX: join exit 0" "$L"; then echo "already joined by the chain"; exit 0; fi
( cd "$D/logs" && PYTHONPATH="$R" "$R/.venv-pysis/Scripts/python.exe" -X utf8 -W ignore -u -m pipeline.true_irc_prod "$R/reactions/$RX" join \
    > "$D/console_true_irc_join_rerun.txt" 2>&1 ); rc=$?
echo "$(date '+%Y-%m-%d %H:%M:%S') production true IRC $RX: join after rerun, exit $rc ($(tail -1 "$D/console_true_irc_join_rerun.txt" | cut -c1-200))" >> "$L"
