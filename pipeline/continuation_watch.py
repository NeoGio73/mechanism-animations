"""Re-estimate each continuation after its first 50 new points; pause it if the estimate passes its limit
(J. Picione, 2026-11-16: "After the first 50 points of each continuation, re-estimate from that side's own rate. If the
new estimate passes its limit, pause it and tell me.").

Limits (his answers of 2026-11-16): secondary shift backward 3 days; E2 product side 2 days; cyanide reactant side and
E2 reactant side were started without a stated limit, so the limit from his earlier checkpoint applies (about a day,
24 h). The limit is compared with the side's whole continuation: hours already spent + hours still estimated.

Own-rate estimate: the side writes <side folder>/progress.json after every step (energies of the new points, step
times; pipeline/true_irc2.OwnStop). From the first 50 new points: drop per point d_i (10-point windows) against the gap
g_i above the end minimum, fitted as d = a * g^p (log-log least squares, p kept within 0..2); points still needed to
reach 0.5 kJ/mol = integral dg / (a g^p) from 0.5 to the current gap, + 20 flat points. If the drops do not fall with
the gap (no fit), the current drop is used as constant. Hours per point: the side's own mean over its last 30 points.
A side over its limit gets <side folder>/PAUSE (it then ends cleanly, ended "paused", nothing joined) and a line in
docs/CONTINUATION_PAUSED. Every check is logged in docs/continuation_watch.log; later checks every 50 points are for
information only.
Config: docs/continuation_watch.json [{"reaction", "side", "folder" (under 02_irc), "minimum", "limit_h"}]
Usage: python -m pipeline.continuation_watch        (runs until every listed side has ended)
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HARTREE_KJ = 2625.4996394799
NEAR, FLAT, FIRST, WIN = 0.5, 20, 50, 10
CFG = ROOT / "docs/continuation_watch.json"
LOG = ROOT / "docs/continuation_watch.log"
FLAG = ROOT / "docs/CONTINUATION_PAUSED"


def log(t: str) -> None:
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} {t}\n")


def own_estimate(E_kj: list[float], times: list[float], e_min_kj: float) -> dict:
    spent = float((times[-1] - times[0]) / 3600) if len(times) > 1 else 0.0
    E = np.asarray(E_kj)[-(FIRST + 1):]          # the side's own recent rate: its last 50 points
    times = list(times)[-(FIRST + 1):]
    g = E - e_min_kj
    d = (E[:-WIN] - E[WIN:]) / WIN                     # drop per point over 10-point windows
    gm = g[WIN // 2: WIN // 2 + len(d)]
    ok = (d > 1e-6) & (gm > 0)
    gap, d_now = float(g[-1]), float(d[-1]) if len(d) else 0.0
    fit = None
    if ok.sum() >= 5:
        p, la = np.polyfit(np.log(gm[ok]), np.log(d[ok]), 1)
        if p > 0:
            p = float(min(p, 2.0))      # clamped: a refitted with p fixed
            fit = (float(np.exp(np.mean(np.log(d[ok]) - p * np.log(gm[ok])))), p)
    if gap <= NEAR:
        n = 0.0
    elif fit:
        a, p = fit
        n = (math.log(gap / NEAR) / a) if abs(p - 1) < 1e-9 else (gap ** (1 - p) - NEAR ** (1 - p)) / (a * (1 - p))
    else:
        n = (gap - NEAR) / d_now if d_now > 1e-6 else math.inf
    n += FLAT
    t = np.asarray(times)
    k = min(30, len(t) - 1)
    s_pt = float((t[-1] - t[-1 - k]) / k) if k > 0 else math.nan
    return {"gap_kJ": gap, "drop_now_kJ": d_now, "fit_a": fit[0] if fit else None, "fit_p": fit[1] if fit else None,
            "points_left": n, "s_per_point": s_pt, "hours_left": n * s_pt / 3600,
            "hours_spent": spent}


def main() -> None:
    st = json.loads((ROOT / "docs/stationary_m062x.json").read_text(encoding="utf-8"))
    sides = json.loads(CFG.read_text(encoding="utf-8"))
    state = {i: {"checked": 0, "done": False} for i in range(len(sides))}
    names = ["%s %s (limit %g h)" % (s["reaction"], s["side"], s["limit_h"]) for s in sides]
    log("watch started: " + ", ".join(names))
    while not all(v["done"] for v in state.values()):
        for i, s in enumerate(sides):
            v = state[i]
            if v["done"]:
                continue
            wd = ROOT / "reactions" / s["reaction"] / "02_irc" / s["folder"]
            if (wd / "path.json").exists():
                j = json.loads((wd / "path.json").read_text(encoding="utf-8"))
                log(f"{s['reaction']} {s['side']}: side ended ({j['ended']}), {j['n_points']} points")
                v["done"] = True
                continue
            pj = wd / "progress.json"
            if not pj.exists():
                continue
            try:
                pr = json.loads(pj.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            n = pr["n_new"]
            if n < FIRST + 1 or n // FIRST <= v["checked"]:
                continue
            v["checked"] = n // FIRST
            e_min = st[s["minimum"]]["E_eh"] * HARTREE_KJ
            est = own_estimate([e * HARTREE_KJ for e in pr["energies"]], pr["times"], e_min)
            total = est["hours_spent"] + est["hours_left"]
            first = v["checked"] == 1
            over = total > s["limit_h"]
            log(f"{s['reaction']} {s['side']}: {n - 1} new points, {est['gap_kJ']:.2f} kJ/mol above {s['minimum']}, drop "
                f"{est['drop_now_kJ']:.4f} kJ/mol/point, fit p = {est['fit_p']}, {est['points_left']:.0f} points left at "
                f"{est['s_per_point']:.0f} s/point = {est['hours_left']:.1f} h; spent {est['hours_spent']:.1f} h; total "
                f"{total:.1f} h vs limit {s['limit_h']:g} h" + (" -> OVER THE LIMIT" if over else "")
                + ("" if first else " (information only)"))
            (wd / f"estimate_{n - 1:04d}.json").write_text(json.dumps({**s, **est, "total_h": total}, indent=1), encoding="utf-8")
            if first and over:
                why = (f"paused by pipeline.continuation_watch {datetime.now():%Y-%m-%d %H:%M}: re-estimate after {n - 1} new "
                       f"points {total:.1f} h > limit {s['limit_h']:g} h (J. Picione 2026-11-16)")
                (wd / "PAUSE").write_text(why, encoding="utf-8")
                with FLAG.open("a", encoding="utf-8") as fh:
                    fh.write(f"{s['reaction']} {s['side']}: {why}\n")
                log(f"{s['reaction']} {s['side']}: PAUSE written")
        time.sleep(120)
    log("watch ended: every listed side has ended")


if __name__ == "__main__":
    main()
