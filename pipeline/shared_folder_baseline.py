"""Check 1, Part A baseline (2026-11-12): how closely does a recompute alone reproduce the energy and
gradient of an evaluation that did NOT overlap another job? Same level, same recompute as Part A
(shared_folder_check.engrad), on gradient evaluations 7-12 of the s005 forward rerun snapshot, whose
only overlapping evaluations were the first 6. Without this baseline, Part A's A2 differences cannot
be told apart from ordinary recompute scatter.
Result: <step_convergence>/check1/baseline.json
Usage: python -m pipeline.shared_folder_baseline
"""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from . import shared_folder_check as sfc


def main() -> dict:
    prop = sfc.CK / "snapshots/s005_forward_rerun_IRC_F.property.txt"
    ev = sfc.read_property(prop)
    gev = [e for e in ev if e["grad"] is not None and "final" in e]
    base = gev[6:12]                 # the six gradient evaluations right after the overlapping ones
    with ThreadPoolExecutor(2) as ex:
        eg = list(ex.map(lambda a: sfc.engrad("baseline_s005_forward_rerun", a[0] + 6, a[1]), enumerate(base)))
    dE = [abs(x["E"] - e["final"]) for x, e in zip(eg, base)]
    dG = [float(np.max(np.abs(x["grad"] - e["grad"]))) for x, e in zip(eg, base)]
    out = {"property_file": str(prop.relative_to(sfc.OUT)), "evaluations": [e["idx"] for e in base],
           "overlapped_another_job": False, "dE_Eh": dE, "dgrad_max_Eh_per_bohr": dG,
           "max_abs_dE_Eh": max(dE), "max_abs_dgrad_Eh_per_bohr": max(dG)}
    (sfc.CK / "baseline.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1), flush=True)
    return out


if __name__ == "__main__":
    main()
