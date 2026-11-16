"""The minimum a path reaches is the official end structure (J. Picione: E2 on 2026-11-16, every path the same day).

Every side of every production path is optimized from its last frame (pipeline/end_check.py). If the optimized
structure is not the same minimum as the current reference (same-minimum criteria in end_check.py) but has the same
bonds, it becomes the official end structure: a new entry <prefix>/<name>_irc in docs/stationary_m062x.json
(official_end true), the reference kept on record with superseded_by (left out by path_angles._minima), and
02_irc/endpoints.json rewritten (steps.from_stationary) when the path's summary files exist.
Safeguard (J. Picione, 2026-11-16): if the optimized structure has different bonds from the expected end, nothing is
switched; the side is reported with "STOP: different bonds" and a line goes to docs/END_BONDS_STOP for him.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import orca
from .steps import COV_RADII_A

ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / "docs" / "stationary_m062x.json"
BOND_FACTOR = 1.25          # bonded: distance <= 1.25 x sum of covalent radii (as steps._bond_states)


def bonds(sym, x) -> set[tuple[int, int]]:
    x = np.asarray(x, dtype=float)
    out = set()
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            if np.linalg.norm(x[i] - x[j]) <= BOND_FACTOR * (COV_RADII_A[sym[i]] + COV_RADII_A[sym[j]]):
                out.add((i, j))
    return out


def bond_difference(sym, x_opt, x_ref) -> dict:
    a, b = bonds(sym, x_opt), bonds(sym, x_ref)
    lab = lambda p: f"{sym[p[0]]}{p[0] + 1}–{sym[p[1]]}{p[1] + 1}"
    return {"same": a == b, "only_in_optimized": [lab(p) for p in sorted(a - b)],
            "only_in_reference": [lab(p) for p in sorted(b - a)]}


def stop_flag(rx_name: str, side: str, minimum: str, diff: dict) -> None:
    with open(ROOT / "docs" / "END_BONDS_STOP", "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M')} {rx_name} {side}: optimization from the last frame has different bonds "
                 f"from {minimum} (only in the optimized structure: {diff['only_in_optimized'] or 'none'}; only in the "
                 f"reference: {diff['only_in_reference'] or 'none'}). Nothing switched; waiting for J. Picione.\n")


def _key_of(js: dict, rx_name: str, name: str) -> str:
    for k, v in js.items():
        if v.get("name") == name and v.get("reaction_folder", v["folder"].split("/")[0]) == rx_name and not v.get("superseded_by"):
            return k
    raise KeyError(f"{name} of {rx_name} not in {ST.name}")


def adopt(rx_root: Path, side: str, minimum: str, jd: Path, job: str, start: str, res: dict) -> str:
    """Make the structure optimized in jd/<job> the official end structure in place of `minimum`. Returns its name."""
    js = json.loads(ST.read_text(encoding="utf-8"))
    old_key = _key_of(js, rx_root.name, minimum)
    prefix = old_key.split("/")[0]
    name = f"{minimum}_irc"
    while f"{prefix}/{name}" in js:
        name += "2"
    r = orca.parse(jd, job, float("nan"))
    h = orca.read_hess(jd / f"{job}.hess")
    fr = np.sort(h["freqs"])
    lowest = float(fr[np.abs(fr) > 1e-3][0])
    t = time.strftime("%Y-%m-%d")
    js[f"{prefix}/{name}"] = {
        "reaction": prefix, "folder": jd.relative_to(ROOT / "reactions").as_posix(), "reaction_folder": rx_root.name,
        "name": name, "kind": "min", "job": job, "ok": bool(r.ok), "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
        "imaginary_cm": list(r.imaginary), "lowest_freq_cm": round(lowest, 2), "accepted": True, "official_end": True,
        "start": start,
        "comparison_with_reference": {k: res[k] for k in ("opt_E_minus_min_kJ", "opt_rmsd_heavy_A", "opt_rmsd_all_atoms_A",
                                                          "key_distances", "criteria") if k in res},
        "decision": (f"J. Picione 2026-11-16: the minimum the path reaches is the official end structure for every path "
                     f"(the path defines which minima the TS connects); replaces {minimum}, kept on record (switched {t})")}
    js[old_key].update(superseded_by=name, official_end=False,
                       note=(js[old_key].get("note", "") + " " if js[old_key].get("note") else "")
                       + f"superseded as the official end structure by {name}, the minimum the path reaches ({t}; "
                         "J. Picione 2026-11-16); kept on record")
    ST.write_text(json.dumps(js, indent=1), encoding="utf-8")
    return name


def rewrite_endpoints(rx_root: Path) -> bool:
    """02_irc/endpoints.json from the official minima, if the path's summary files exist."""
    from .steps import Reaction, from_stationary
    if not ((rx_root / "02_irc/irc_IRC_F.xyz").exists() and (rx_root / "01_stationary/ts/ts.out").exists()):
        return False
    from_stationary(Reaction.load(rx_root))
    return True
