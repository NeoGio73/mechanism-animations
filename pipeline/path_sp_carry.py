"""Carry path single points over to a re-joined path (J. Picione, 2026-11-16: E2's single points start on the existing
frames while its sides are extended).

Path single points live in 02_irc/path_sp/<label>/f<frame>, numbered along the joined path. Extending the backward
side adds frames before the TS, so after the re-join every frame number changes. This step, run after the re-join and
before path_singlepoints:
  1. moves 02_irc/path_sp to 02_irc/path_sp_part1 (kept as computed; its frame numbers refer to the earlier join,
     whose trajectory is copied there as irc_full_trj_part1.xyz);
  2. for each frame of the new path whose geometry equals a frame of the earlier join (largest coordinate difference
     <= 1e-6 Å: the extension keeps the earlier frames verbatim), copies that frame's finished job folder to the new
     frame number. path_singlepoints then finds an identical ORCA input there and reuses it (orca.run, reuse=True);
     every other frame is computed in a clean folder.
Usage: python -m pipeline.path_sp_carry <reaction> <earlier joined trajectory>
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

from .geometry import read_xyz_frames

ROOT = Path(__file__).resolve().parents[1]


def carry(rx: str, earlier_trj: str) -> dict:
    R = ROOT / "reactions" / rx / "02_irc"
    old = [np.asarray(f.coords) for f in read_xyz_frames(Path(earlier_trj))]
    new = [np.asarray(f.coords) for f in read_xyz_frames(R / "irc_IRC_Full_trj.xyz")]
    sp, part1 = R / "path_sp", R / "path_sp_part1"
    if part1.exists():
        raise RuntimeError(f"{part1} exists already: carry runs once per re-join")
    if not sp.exists():
        return {"reaction": rx, "carried": 0, "note": "no earlier path single points"}
    shutil.move(str(sp), str(part1))
    shutil.copy2(earlier_trj, part1 / "irc_full_trj_part1.xyz")
    match = {}
    for k, x in enumerate(new):
        for j, y in enumerate(old):
            if x.shape == y.shape and float(np.max(np.abs(x - y))) <= 1e-6:
                match[k] = j
                break
    carried = []
    for lab in [p for p in part1.iterdir() if p.is_dir()]:
        for k, j in match.items():
            src = lab / f"f{j:03d}"
            if (src / "f.out").exists():
                dst = sp / lab.name / f"f{k:03d}"
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copytree(src, dst)
                carried.append((lab.name, j, k))
    out = {"reaction": rx, "earlier_frames": len(old), "new_frames": len(new), "matched_frames": len(match),
           "carried_jobs": len(carried), "map_new_to_earlier": match}
    (part1 / "carry.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    o = carry(sys.argv[1], sys.argv[2])
    print(json.dumps({k: v for k, v in o.items() if k != "map_new_to_earlier"}, indent=1))
