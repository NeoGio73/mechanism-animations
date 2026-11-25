"""Diels–Alder flicker test (J. Picione, 2026-11-25 evening). From files on disk; orca_plot only (no ORCA calculation).

His reading: between the TS and s ≈ +1 the localizer lands on the clockwise or the counterclockwise set frame by
frame. Two tests at frames 179–199:

(a) 3×17: each reflected active orbital (through that frame's mirror plane, as the viewer) against all 17 computed
    valence IBOs of the same frame. Smallest singular value of the 3×17 overlap matrix: 1 means the reflected
    active orbitals lie inside the computed occupied valence space. Cubes of all 17 orbitals are written to a
    scratch folder (not the repository) from the frame's f.loc.
(b) continuity, frames 180–199: the three displayed orbitals of the previous frame against this frame's three, once
    as computed and once mirror-reflected (each option matched to the previous orbitals by the best assignment);
    keep the option with the larger lowest overlap. Previous orbitals are interpolated onto this frame's grid.

Limits: (a) ≥ 0.999 at every frame and (b) ≥ 0.99 at every frame -> the choice may be used for the display.
Writes <rx>/03_orbitals/mirror_continuity.json.

Usage: python -X utf8 -m pipeline.mirror_continuity <scratch folder> [first last] [--write]
  --write: for the frames where the mirror image continues the previous frame, write cubes/fNNN/ibo<id>_mirror.npz
  (display only; the computed cubes, .loc files and tracking.json are not changed)
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from .geometry import read_xyz_frames
from .ibo_check import BOHR, _interp
from .mirror_display_check import _cube, mirror_perm, plane, reflect
from .meshes import read_cube_gz
from .orbitals import plot_cubes, read_cube
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
WRITE = "--write" in sys.argv
NAME = "da_butadiene_ethylene_m062x"


def _grid_pts_A(v, o, a):
    g = np.stack(np.meshgrid(*[np.arange(s) for s in v.shape], indexing="ij"), -1).reshape(-1, 3)
    return (o + g * np.diag(a)) * BOHR


def main(scratch: Path, first: int = 179, last: int = 199) -> dict:
    rx = Reaction.load(ROOT / "reactions" / NAME)
    frames = read_xyz_frames(rx.root / "02_irc/irc_full_trj.xyz")
    syms = frames[0].symbols
    tr = json.loads((rx.root / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    meta = json.loads((rx.root / "02_irc/irc_meta.json").read_text(encoding="utf-8"))
    perm = mirror_perm(rx, np.asarray(frames[meta["ts_frame"]].coords), syms)
    act = [e["tracked_id"] for e in tr["expected_changes"]["items"]]
    m = len(tr["tracked"])
    s = [float(x.split(",")[1]) for x in (rx.root / "02_irc/energies.csv").read_text(encoding="utf-8").splitlines()[1:]]
    rows, prev = [], None
    for k in range(first, last + 1):
        wd = scratch / f"f{k:03d}"
        wd.mkdir(parents=True, exist_ok=True)
        cubes = {}
        for j in range(m):
            p = wd / f"ibo{j}.cube"
            if p.exists():
                cubes[j] = p
        if len(cubes) < m:
            for ext in ("loc", "densities", "densitiesinfo"):       # orca_plot needs the density files beside it
                shutil.copy(rx.root / f"03_orbitals/frames/f{k:03d}/f.{ext}", wd / f"f.{ext}")
            idx = {tr["per_frame"][k]["orbital_index"][j]: j for j in range(m)}
            made = plot_cubes(wd / "f.loc", sorted(idx))
            for p in made:
                o = int(re.search(r"\.mo(\d+)a\.cube$", p.name).group(1))
                p.rename(wd / f"ibo{idx[o]}.cube")
            cubes = {j: wd / f"ibo{j}.cube" for j in range(m)}
        sign = tr["per_frame"][k]["sign"]
        V = {}
        for j in range(m):
            v, o, a = read_cube(cubes[j])[:3]
            V[j] = v.astype(np.float64) * sign[j]
        v0, o0, a0 = read_cube(cubes[act[0]])[:3]
        pts = _grid_pts_A(v0, o0, a0)
        X = np.asarray(frames[k].coords)
        n, d = plane(X, perm)
        back = reflect(pts, n, d) / BOHR
        flat = {j: V[j].reshape(-1) for j in range(m)}
        norm = {j: np.sqrt(flat[j] @ flat[j]) for j in range(m)}
        refl = {i: _interp(V[i], o0, a0, back) for i in act}          # σφ_i on this frame's grid
        rnorm = {i: np.sqrt(refl[i] @ refl[i]) for i in act}
        M = np.array([[refl[i] @ flat[j] / (rnorm[i] * norm[j]) for j in range(m)] for i in act])
        sv17 = np.linalg.svd(M, compute_uv=False)
        row = {"frame": k, "s": round(s[k], 3), "smallest_sv_3x17": round(float(sv17.min()), 5),
               "sv_3x17": np.round(sv17, 5).tolist(),
               "captured_per_reflected_orbital": np.round((M ** 2).sum(1), 5).tolist()}
        comp = {i: flat[i] / norm[i] for i in act}
        mirr = {i: refl[i] / rnorm[i] for i in act}
        if prev is not None:
            pv = {i: _interp(prev["vals"][i], prev["o"], prev["a"], pts / BOHR) for i in act}
            pv = {i: x / np.sqrt(x @ x) for i, x in pv.items()}
            best = None
            for opt, cur in (("computed", comp), ("mirror", mirr)):
                O = np.array([[pv[i] @ cur[j] for j in act] for i in act])
                r, c = linear_sum_assignment(-np.abs(O))
                low = float(np.abs(O[r, c]).min())
                row[f"lowest_overlap_{opt}"] = round(low, 5)
                row[f"assignment_{opt}"] = [act[x] for x in c]
                if best is None or low > best[1]:
                    best = (opt, low, {act[a_]: cur[act[b_]] * np.sign(O[a_, b_]) for a_, b_ in zip(r, c)},
                            [(act[a_], act[b_], float(np.sign(O[a_, b_]))) for a_, b_ in zip(r, c)])
            row["chosen"], row["continuity_overlap"] = best[0], round(best[1], 5)
            disp = best[2]
            if best[0] == "mirror":
                row["shown_as"] = [{"label_of": a_, "reflection_of": b_, "sign": sg} for a_, b_, sg in best[3]]
                if WRITE:     # display file in the raw phase of the label's tracked orbital (meshes applies it)
                    cd = rx.root / f"03_orbitals/cubes/f{k:03d}"
                    _, oc, ac = read_cube_gz(cd / f"ibo{act[0]}.cube.gz")
                    assert np.allclose(oc, o0) and np.allclose(ac, a0), f"frame {k}: grid differs"
                    for a_, b_, sg in best[3]:
                        vals = sign[a_] * sg * refl[b_].reshape(v0.shape)
                        np.savez_compressed(cd / f"ibo{a_}_mirror.npz", values=vals.astype(np.float32),
                                            origin=o0, axes=a0)
        else:
            row["chosen"], row["continuity_overlap"] = "computed", None
            disp = comp
        prev = {"vals": {i: disp[i].reshape(v0.shape) for i in act}, "o": o0, "a": a0}
        rows.append(row)
        print({x: y for x, y in row.items() if x not in ("sv_3x17",)}, flush=True)
    a_ok = min(r["smallest_sv_3x17"] for r in rows)
    b_vals = [r["continuity_overlap"] for r in rows if r["continuity_overlap"] is not None]
    res = {"reaction": NAME, "frames": [first, last], "active": act,
           "worst_3x17": a_ok, "worst_continuity": min(b_vals),
           "passes": a_ok >= 0.999 and min(b_vals) >= 0.99,
           "mirror_frames": [r["frame"] for r in rows if r["chosen"] == "mirror"], "rows": rows,
           "note": "cubes of all 17 orbitals in a scratch folder, not committed; previous displayed orbitals "
                   "interpolated trilinearly onto this frame's grid"}
    (rx.root / "03_orbitals/mirror_continuity.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({x: y for x, y in res.items() if x != "rows"}, indent=1))
    return res


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if x != "--write"]
    main(Path(a[0]), *(int(x) for x in a[1:3]))
