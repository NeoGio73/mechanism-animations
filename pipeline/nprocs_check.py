"""Does ORCA's parallel core count change the energy? (J. Picione, 2026-11-16 evening, after the switch to
CORES_FIT, which may run a job on fewer or more cores than it asked for.)

Reruns two finished SN2 single points on the TS, M06-2X/ma-def2-TZVP (01_ts/sp_tzvp/ts) and DLPNO-CCSD(T)/
ma-def2-TZVP (01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/ts), with the identical input except %pal nprocs
(4 -> 2), each in a clean folder under 01_ts/nprocs_check/. Through the core queue with fit mode off, so the job
really runs on 2 cores. Compares FINAL SINGLE POINT ENERGY against the SCF noise floor (VeryTightSCF, 1.3e-5 Eh =
0.034 kJ/mol; docs/protocol.md). Writes 01_ts/nprocs_check/compare.{md,json}.
Usage: .venv/Scripts/python.exe -X utf8 -m pipeline.nprocs_check   (started hidden)
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

from . import orca

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions/sn2_hydroxide_bromomethane_m062x"
OUT = RX / "01_ts/nprocs_check"
JOBS = [("M06-2X/ma-def2-TZVP", RX / "01_ts/sp_tzvp/ts", "ts", "m062x_tzvp"),
        ("DLPNO-CCSD(T)/ma-def2-TZVP", RX / "01_ts/diagnostics/DLPNO-CCSD(T)_ma-def2-TZVP/ts", "ts", "dlpno")]
NEW = 2
EH_KJ = 2625.4996
FLOOR_KJ = 0.034


def energy(out: Path) -> float:
    t = out.read_text(errors="replace")
    assert "ORCA TERMINATED NORMALLY" in t, out
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])


def main() -> dict:
    os.environ.pop("CORES_FIT", None)
    OUT.mkdir(parents=True, exist_ok=True)
    log = orca.RunLog(RX / "run_log.md")
    rows = []
    for level, d, name, sub in JOBS:
        inp = (d / f"{name}.inp").read_text()
        old_n = int(re.search(r"%pal\s+nprocs\s+(\d+)", inp).group(1))
        new_inp = re.sub(r"%pal\s+nprocs\s+\d+\s+end", f"%pal nprocs {NEW} end", inp, count=1)
        assert orca._same_input(inp, new_inp) and new_inp != inp
        wd = OUT / sub
        assert not (wd / f"{name}.out").exists(), f"{wd} not clean"
        r = orca.run(wd, name, new_inp, log, note=f"core-count check: same input as {d.relative_to(RX).as_posix()} "
                     f"except %pal nprocs {old_n} -> {NEW} (J. Picione 2026-11-16 evening)", reuse=False)
        e_old, e_new = energy(d / f"{name}.out"), energy(wd / f"{name}.out")
        rows.append({"level": level, "original": d.relative_to(RX).as_posix(), "check": wd.relative_to(RX).as_posix(),
                     "nprocs_original": old_n, "nprocs_check": NEW, "E_original_eh": e_old, "E_check_eh": e_new,
                     "dE_eh": e_new - e_old, "dE_kJ": (e_new - e_old) * EH_KJ, "wall_s_check": r.wall_s})
    res = {"floor_kJ": FLOOR_KJ, "rows": rows}
    (OUT / "compare.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    md = ["# Core-count check: SN2 TS single points", "",
          f"Identical input except `%pal nprocs` ({rows[0]['nprocs_original']} → {NEW}); clean folders; "
          f"SCF noise floor {FLOOR_KJ} kJ/mol (VeryTightSCF).", "",
          "| Level | E, original (Eh) | E, check (Eh) | ΔE (Eh) | ΔE (kJ/mol) | wall time, check (s) |",
          "|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['level']} | {r['E_original_eh']:.12f} | {r['E_check_eh']:.12f} | {r['dE_eh']:.2e} | "
                  f"{r['dE_kJ']:.2e} | {r['wall_s_check']:.0f} |")
    (OUT / "compare.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return res


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
