"""Is the mirror image of the computed orbital set exact at every frame? (J. Picione, 2026-11-25, Diels–Alder.)

The production Diels–Alder path localized the three active IBOs in the clockwise pattern; the viewer shows their
mirror image through the TS mirror plane (reaction.toml [mirror]) so that it opens on the counterclockwise drawing.
That is only allowed where the reflection is exact. For every frame:

1. geometry: the frame reflected through the plane the viewer uses (pipeline/meshes.mirror_data: fitted to the atom
   pairs of that frame) against the frame itself with the atom swap; largest atom displacement in Å;
2. orbitals: the 3×3 overlap matrix M_ij = <σφ_i|φ_j> between the reflected active orbitals σφ_i and the computed
   active orbitals φ_j (cube files, the reflected cube interpolated trilinearly onto the computed grid, each
   orbital normalized on its grid). Smallest singular value and |det M|: 1 means the reflected set spans exactly
   the space of the computed set; each reflected orbital's best match is listed too.

Limits (J. Picione): geometry within 0.001 Å and orbital agreement at least 0.999 at every frame; otherwise the
reflected set is not to be displayed and he decides (show the computed set and flip the drawing, or re-localize
from the reflected orbitals).

Also reported: the lowest frame-to-frame overlap of every tracked orbital (tracking.json), the active ones apart.

Usage: python -X utf8 -m pipeline.mirror_display_check <reaction folder name> ibo1,ibo4,ibo6
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

from .geometry import read_xyz_frames
from .ibo_check import BOHR, _interp
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
GEOM_LIMIT_A = 0.001
ORB_LIMIT = 0.999


def _cube(p: Path):
    lines = gzip.open(p, "rt").read().splitlines()
    nat = abs(int(lines[2].split()[0]))
    origin = np.array([float(v) for v in lines[2].split()[1:4]])
    dims, axes = [], []
    for i in range(3):
        q = lines[3 + i].split()
        dims.append(int(q[0]))
        axes.append([float(v) for v in q[1:4]])
    atoms = np.array([[float(v) for v in lines[6 + a].split()[2:5]] for a in range(nat)])
    start = 6 + nat + (1 if int(lines[2].split()[0]) < 0 else 0)
    vals = np.array(" ".join(lines[start:]).split(), dtype=np.float64).reshape(dims)
    return vals, origin, np.array(axes), atoms


def mirror_perm(rx: Reaction, X: np.ndarray, symbols: list[str]) -> list[int]:
    """Same permutation as pipeline/meshes.mirror_data (heavy pairs from [mirror]; hydrogens matched at the TS)."""
    perm = list(range(len(X)))
    for a, b in rx.cfg["mirror"]["pairs"]:
        perm[a], perm[b] = b, a
    n, d = plane(X, perm)
    R = reflect(X, n, d)
    heavy = {i for pr in rx.cfg["mirror"]["pairs"] for i in pr}
    for i, sym in enumerate(symbols):
        if i in heavy:
            continue
        cand = [j for j, s in enumerate(symbols) if s == sym and j not in heavy]
        perm[i] = min(cand, key=lambda j: np.linalg.norm(R[i] - X[j]))
    assert sorted(perm) == list(range(len(X))) and all(perm[perm[i]] == i for i in range(len(X))), perm
    return perm


def plane(X, perm):
    """Identical to pipeline/meshes.mirror_data.plane (the plane the viewer reflects through)."""
    diffs = [X[perm[i]] - X[i] for i in range(len(X)) if perm[i] > i]
    ref = max(diffs, key=np.linalg.norm)
    n = sum(d if d @ ref >= 0 else -d for d in diffs)
    n = n / np.linalg.norm(n)
    d = float(np.mean([(X[i] + X[perm[i]]) @ n / 2 for i in range(len(X))]))
    return n, d


def reflect(X, n, d):
    return X - 2 * ((X @ n) - d)[:, None] * n


def run(name: str, keys: list[str]) -> dict:
    rx = Reaction.load(ROOT / "reactions" / name)
    frames = read_xyz_frames(rx.root / "02_irc/irc_full_trj.xyz")
    syms = frames[0].symbols
    meta = json.loads((rx.root / "02_irc/irc_meta.json").read_text(encoding="utf-8"))
    ts = meta["ts_frame"]
    perm = mirror_perm(rx, np.asarray(frames[ts].coords), syms)
    rows = []
    for k, fr in enumerate(frames):
        X = np.asarray(fr.coords)
        n, d = plane(X, perm)
        dev = np.linalg.norm(reflect(X, n, d) - X[perm], axis=1)
        cubes = [_cube(rx.root / f"03_orbitals/cubes/f{k:03d}/{key}.cube.gz") for key in keys]
        cube_xyz_A = cubes[0][3] * BOHR
        geom_cube_A = float(np.abs(cube_xyz_A - X).max())        # cube written at this frame's geometry
        M = np.zeros((len(keys), len(keys)))
        for j, (vj, oj, aj, _) in enumerate(cubes):
            g = np.stack(np.meshgrid(*[np.arange(s) for s in vj.shape], indexing="ij"), -1).reshape(-1, 3)
            pts_A = (oj + g * np.diag(aj)) * BOHR
            back_A = reflect(pts_A, n, d)                          # σ is its own inverse
            b = vj.reshape(-1)
            for i, (vi, oi, ai, _) in enumerate(cubes):
                a = _interp(vi, oi, ai, back_A / BOHR)            # (σφ_i)(r) = φ_i(σ r)
                M[i, j] = float(a @ b / np.sqrt((vi.reshape(-1) @ vi.reshape(-1)) * (b @ b)))
        sv = np.linalg.svd(M, compute_uv=False)
        best = [(keys[int(np.argmax(np.abs(M[i])))], round(float(np.abs(M[i]).max()), 5)) for i in range(len(keys))]
        rows.append({"frame": k, "geom_max_dev_A": float(dev.max()), "geom_rms_A": float(np.sqrt((dev ** 2).mean())),
                     "cube_vs_frame_geometry_A": geom_cube_A,
                     "smallest_singular_value": float(sv.min()), "abs_det": float(abs(np.linalg.det(M))),
                     "best_match_of_reflected": dict(zip(keys, best)), "M": np.round(M, 5).tolist()})
        if k % 25 == 0:
            print(f"frame {k}: geometry {dev.max():.2e} Å, smallest singular value {sv.min():.5f}, "
                  f"|det| {abs(np.linalg.det(M)):.5f}", flush=True)
    tr = json.loads((rx.root / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    active = [int(k[3:]) for k in keys]
    ov = [(min(pf["overlap_with_previous"][i] for i in range(len(pf["overlap_with_previous"]))), pf["frame"])
          for pf in tr["per_frame"][1:]]
    ova = [(min(pf["overlap_with_previous"][i] for i in active), pf["frame"]) for pf in tr["per_frame"][1:]]
    worst_all = min(ov)
    worst_act = min(ova)
    gw = max(rows, key=lambda r: r["geom_max_dev_A"])
    ow = min(rows, key=lambda r: r["smallest_singular_value"])
    dw = min(rows, key=lambda r: r["abs_det"])
    passed = gw["geom_max_dev_A"] <= GEOM_LIMIT_A and ow["smallest_singular_value"] >= ORB_LIMIT
    fail = [r["frame"] for r in rows if r["geom_max_dev_A"] > GEOM_LIMIT_A or r["smallest_singular_value"] < ORB_LIMIT]
    n_end = 30
    tail = rows[-n_end:]
    res = {"reaction": name, "orbitals": keys, "perm": perm, "frames": len(rows), "ts_frame": ts,
           "limits": {"geometry_A": GEOM_LIMIT_A, "orbital_smallest_singular_value": ORB_LIMIT},
           "passed": passed, "failing_frames": fail,
           "worst_geometry": {"frame": gw["frame"], "max_dev_A": gw["geom_max_dev_A"], "rms_A": gw["geom_rms_A"]},
           "worst_smallest_singular_value": {"frame": ow["frame"], "value": ow["smallest_singular_value"]},
           "worst_abs_det": {"frame": dw["frame"], "value": dw["abs_det"]},
           "product_end_last_30_frames": {
               "max_geom_dev_A": max(r["geom_max_dev_A"] for r in tail),
               "min_smallest_singular_value": min(r["smallest_singular_value"] for r in tail),
               "min_abs_det": min(r["abs_det"] for r in tail)},
           "lowest_frame_to_frame_overlap_all_tracked": {"value": worst_all[0], "frame": worst_all[1]},
           "lowest_frame_to_frame_overlap_active": {"value": worst_act[0], "frame": worst_act[1]},
           "max_cube_vs_frame_geometry_A": max(r["cube_vs_frame_geometry_A"] for r in rows),
           "method": "plane fitted to the [mirror] atom pairs at each frame (as the viewer); cube overlaps on the "
                     "computed orbital's grid, reflected cube interpolated trilinearly",
           "rows": rows}
    (rx.root / "03_orbitals/mirror_display_check.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))
    return res


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2].split(","))
