#!/usr/bin/env bash
# E2 production true IRC at M06-2X (J. Picione, 2026-11-13: "yes, after the five reactions").
# Waits until the five-reaction chain has written its last join line, then runs E2 the same way.
R="$(cd "$(dirname "$0")/.." && pwd)"
L="$R/docs/production_true_irc.log"
until grep -q "hshift_3_methyl_2_butyl_m062x: join" "$L" 2>/dev/null; do sleep 120; done
echo "$(date '+%Y-%m-%d %H:%M:%S') E2 queued after the five: starting" >> "$L"
bash "$R/pipeline/run_true_irc_prod.sh" e2_ethoxide_2_bromopropane_m062x
