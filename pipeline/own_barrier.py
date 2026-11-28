"""Barrier of a reactant conformer through its own transition state (one-molecule reactions).

For the secondary hydride shift the anti form of the 3-methyl-2-butyl cation does not reach the transition state of
the path shown; it reaches the tertiary cation through its own rotation TS (01_stationary/rotation_ts.json). An energy
"from the anti form" to the TS shown is therefore not the barrier of any path (J. Picione, 2026-11-26). This module
gives each form's own barrier instead.

reaction.toml, under the [[zero_conformers]] entry:  own_ts = {result = "01_stationary/rotation_ts.json",
key = "ts_anti_rotation"}. The TS entry must be confirmed (it connects the conformer with the product).
Optional irc = "<result.json of an IRC from that TS>" (J. Picione, 2026-11-27: the IRC is the main evidence for the
connection, the displaced optimizations of rotation_ts.json agree with it); `irc_sentence` words it from that file.

`own_barriers(rx)` -> {name, ts_key, ts_folder, rows: {level: {dE, dH, dG}}}, level None = working level; ΔH and ΔG
add the working-level thermal corrections (ORCA quasi-RRHO block of each .out, 298.15 K) to the electronic energy of
the level, as everywhere else. One molecule throughout: no standard-state term. Reads finished jobs only.
"""
from __future__ import annotations

import json

from .steps import Reaction
from .thermo import orca_block

EH = 2625.499639


def own_barriers(rx: Reaction) -> dict | None:
    c = next((c for c in rx.cfg.get("zero_conformers", []) if c.get("own_ts")), None)
    if c is None:
        return None
    R = rx.root
    form = json.loads((R / c["result"]).read_text(encoding="utf-8"))
    ts = json.loads((R / c["own_ts"]["result"]).read_text(encoding="utf-8"))["ts"][c["own_ts"]["key"]]
    if not ts.get("confirmed") or c["name"] not in ts.get("connects", []):
        raise RuntimeError(f"{c['own_ts']['key']} is not a confirmed TS of {c['name']}")
    ts_out = R / ts["folder"] / "ts" / "ts.out"
    form_out = R / f"{form['job']}.out"
    bt, bf = orca_block(ts_out), orca_block(form_out)
    gcorr = bt["G_minus_E_kJ"] - bf["G_minus_E_kJ"]
    hcorr = bt["H_minus_E_kJ"] - bf["H_minus_E_kJ"]
    assert abs((ts["G_eh"] - ts["E_eh"]) * EH - bt["G_minus_E_kJ"]) < 0.01      # the .out is the job the JSON reports
    assert abs((form["G_eh"] - form["E_eh"]) * EH - bf["G_minus_E_kJ"]) < 0.01
    E = {None: (ts["E_eh"], form["E_eh"]), **{lv: (ts["sp"][lv], form["sp"][lv]) for lv in ts["sp"]}}
    rows = {}
    for lv, (et, ef) in E.items():
        de = (et - ef) * EH
        rows[lv] = {"dE": de, "dH": de + hcorr, "dG": de + gcorr}
    irc = c["own_ts"].get("irc")
    return {"name": c["name"], "ts_key": c["own_ts"]["key"], "ts_folder": ts["folder"],
            "imaginary_cm": ts["ts_job"]["imaginary_cm"], "rows": rows,
            "irc": irc_sentence(R / irc, irc) if irc else None}


def irc_sentence(path, shown: str) -> str:
    """One sentence on where an IRC from a rotation TS ends (result.json of a run like
    00_inputs/irc_ts_anti_rotation.py). Raises if either end is missing, not optimized or has an imaginary frequency."""
    from pathlib import Path
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    parts = []
    for side in ("backward", "forward"):
        e = d["ends"][side]
        o = e["opt"]
        if not isinstance(o, dict) or o["imaginary_cm"] or not o["reached"]:
            raise RuntimeError(f"{shown}: the {side} end is not an optimized minimum")
        m = o["matches"][o["reached"]]
        parts.append(f"`{o['reached']}` on one side ({e['n_points']} points; optimized end within "
                     f"{m['core_rmsd_A']:.3f} Å and {abs(m['dE_kJ']):.2f} kJ/mol)")
    return (f"an IRC from this TS ends at {parts[0]} and at {parts[1].replace('on one side', 'on the other')}, "
            f"with no imaginary frequency at either end (`{shown.replace('result.json', 'result.md')}`)")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    r = own_barriers(Reaction.load(Path(sys.argv[1])))
    for lv, q in r["rows"].items():
        print(lv or "working", {k: round(v, 2) for k, v in q.items()})
