"""Recover path.xyz / path.json of a pipeline/true_irc2.py run that stopped with an exception before
writing them (runs started before the integrator-failure handling of 2026-11-12 23:20). Points come
from pysisyphus's own dump (<dir>/pysis/<direction>_irc_data.h5, written every 5 steps, so up to 4
steps before the failure can be missing), ORCA call timings from the job folders. The side is marked
FAILED with the reason given; it is never treated as complete.
Runs in .venv-pysis (h5py). Usage: python -m pipeline.true_irc2_recover <run> <direction> "<reason>"
"""
from __future__ import annotations

import json
import re
import sys

import h5py
import numpy as np

from .true_irc2 import BOHR, OUT, RUNS

def main(run, direction, reason):
    wd = OUT / run / direction
    h = h5py.File(wd / "pysis" / f"{direction}_irc_data.h5", "r")
    C, E, G = np.array(h["coords"]), list(np.array(h["energies"])), np.array(h["gradients"])
    sym = [a.decode().capitalize() if isinstance(a, bytes) else str(a).capitalize() for a in h["atoms"]]
    calls = []
    for d in sorted(p for p in wd.iterdir() if p.is_dir() and re.fullmatch(r"[gh]\d{4}", p.name)):
        out = (d / "j.out").read_text(errors="replace") if (d / "j.out").exists() else ""
        m = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", out)
        if "ORCA TERMINATED NORMALLY" in out and m:
            calls.append({"folder": d.name, "kind": "EnGrad" if d.name[0] == "g" else "Freq",
                          "wall_s": int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]) + int(m[5]) / 1000})
    with open(wd / "path.xyz", "w") as fh:
        for c, e in zip(C, E):
            fh.write(f"{len(sym)}\ntrue IRC {run} {direction} E {e:.12f}\n")
            fh.write("".join(f"{a:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n" for a, (x, y, z) in zip(sym, c.reshape(-1, 3) * BOHR)))
    res = {"run": run, "direction": direction, "settings": RUNS[run], "step_bohr_amu05": 0.10,
           "ended": f"FAILED: {reason}", "own_stop": None, "n_points": len(C),
           "n_engrad": sum(c["kind"] == "EnGrad" for c in calls), "n_hessian": sum(c["kind"] == "Freq" for c in calls),
           "calls": calls, "energies": [float(e) for e in E], "gradients": G.tolist(), "recovered_from": "pysisyphus h5 dump"}
    (wd / "path.json").write_text(json.dumps(res), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k not in ("gradients", "energies", "calls")}))


if __name__ == "__main__":
    main(*sys.argv[1:4])
