"""Isosurface meshes for the viewer and the renderer.

Output (reactions/<id>/04_meshes/):
  manifest.json   frames (atoms, s, energy), orbital list, byte offsets into meshes.bin
  meshes.bin      for each orbital/frame/lobe: int16 xyz positions (quantized in a global
                  box) followed by uint16 triangle indices; little-endian.
"""
from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import numpy as np

from .geometry import read_xyz_frames
from .steps import Reaction, _read_json, _write_json

BOHR_TO_A = 0.529177210903
QUANT = 32767


def read_cube_gz(path: Path):
    from .orbitals import read_cube
    import tempfile, os
    data = gzip.open(path, "rb").read()
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cube") as fh:
        fh.write(data)
        tmp = fh.name
    try:
        return read_cube(Path(tmp))
    finally:
        os.unlink(tmp)


def lobes(values, origin, axes, iso: float, max_tris: int):
    """Marching cubes at +iso and -iso, decimated; returns [(verts Å, faces), (verts, faces)]."""
    from skimage.measure import marching_cubes
    import fast_simplification

    out = []
    for sgn in (1, -1):
        v = sgn * values
        if v.max() <= iso:
            out.append((np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64)))
            continue
        # pad so surfaces close at the box edge
        vp = np.pad(v, 1, constant_values=v.min())
        verts, faces, _, _ = marching_cubes(vp, level=iso)
        verts -= 1.0
        xyz = (origin + verts @ axes) * BOHR_TO_A
        if len(faces) > max_tris:
            xyz, faces = fast_simplification.simplify(xyz.astype(np.float32), faces.astype(np.int32),
                                                      target_reduction=1 - max_tris / len(faces))
        out.append((np.asarray(xyz, dtype=np.float64), np.asarray(faces, dtype=np.int64)))
    return out


def _pretty(label: str) -> str:
    """'M06-2X_ma-def2-TZVP' -> 'M06-2X/ma-def2-TZVP'"""
    return label.replace("_", "/")


def _curves(rx: Reaction, rows, meta, keep) -> list[dict]:
    """Energy curves along the displayed path (IRC frames in `keep`): working level first, then
    any path_singlepoints levels."""
    curves = [{"key": "working", "label": f"{rx.level()} (working level, IRC)",
               "E_kJ": [round(float(rows[k]["E_kJ_per_mol"]), 3) for k in keep],
               "E_separated_products_kJ": meta["E_separated_products_kJ"]}]
    for label in rx.cfg.get("path_levels", []):
        csv = rx.root / "02_irc" / f"energies_{label}.csv"
        if not csv.exists():
            continue
        lines = csv.read_text().splitlines()[1:]
        info = _read_json(rx.root / "02_irc" / f"energies_{label}.json")
        vals = [lines[k].split(",")[2] for k in keep]
        c = {"key": label, "label": f"{_pretty(label)} single points on the same path",
             # null where the frame was not computed: no interpolation (J. Picione, 2026-11-13)
             "E_kJ": [round(float(v), 3) if v else None for v in vals],
             "E_separated_products_kJ": info["E_separated_products_kJ"]}
        if any(not v for v in vals):
            c["label"] += f" ({sum(1 for v in vals if v)} computed frames)"
            c["sparse"] = True
        curves.append(c)
    return curves


def default_curve(rx: Reaction, curves: list[dict]) -> str:
    """Curve shown first (reaction.toml `default_curve`); falls back to the working level."""
    want = rx.cfg.get("default_curve", "working")
    if want not in {c["key"] for c in curves}:
        print(f"WARNING default_curve '{want}' not available (run path_singlepoints); using working level")
        return "working"
    return want


def _mirror_label(label: str, perm: list[int]) -> str:
    """Relabel atom numbers (1-based, as in 'π(C1=C2) → σ(C1–C6)') through the mirror permutation,
    keeping each bond's atoms in ascending order: 'π(C3=C4) → σ(C4–C5)'."""
    import re

    def bond(m):
        parts = re.split(r"([–=])", m.group(1))
        atoms = [re.sub(r"\d+", lambda d: str(perm[int(d.group()) - 1] + 1), a) for a in parts[::2]]
        seps = parts[1::2]
        if len(atoms) == 2 and all(re.search(r"\d", a) for a in atoms):
            atoms.sort(key=lambda a: int(re.search(r"\d+", a).group()))
        out = atoms[0]
        for sep, a in zip(seps, atoms[1:]):
            out += sep + a
        return f"({out})"
    return re.sub(r"\(([^()]*)\)", bond, label)


def mirror_data(rx: Reaction, frames, keep, ts_frame: int, orbitals: list[dict]) -> dict | None:
    """Mirror plane of a symmetric TS at every displayed frame (reaction.toml [mirror]).

    pairs: heavy atoms the plane swaps (0-based); hydrogens are matched by reflecting the TS.
    Each plane is [nx, ny, nz, d] with n·x = d on the plane, fitted to the atom pairs of that
    frame; rms_A is how far the reflected frame is from the frame itself (0 = exactly symmetric)."""
    cfg = rx.cfg.get("mirror")
    if not cfg:
        return None

    def plane(X, perm):
        diffs = [X[perm[i]] - X[i] for i in range(len(X)) if perm[i] > i]
        ref = max(diffs, key=np.linalg.norm)
        n = sum(d if d @ ref >= 0 else -d for d in diffs)
        n = n / np.linalg.norm(n)
        d = float(np.mean([(X[i] + X[perm[i]]) @ n / 2 for i in range(len(X))]))
        return n, d

    reflect = lambda X, n, d: X - 2 * ((X @ n) - d)[:, None] * n
    X = frames[ts_frame].coords
    perm = list(range(len(X)))
    for a, b in cfg["pairs"]:
        perm[a], perm[b] = b, a
    n, d = plane(X, perm)
    R = reflect(X, n, d)
    heavy = {i for pr in cfg["pairs"] for i in pr}
    for i, sym in enumerate(frames[0].symbols):
        if i in heavy:
            continue
        cand = [j for j, s in enumerate(frames[0].symbols) if s == sym and j not in heavy]
        perm[i] = min(cand, key=lambda j: np.linalg.norm(R[i] - X[j]))
    if sorted(perm) != list(range(len(X))) or any(perm[perm[i]] != i for i in range(len(X))):
        raise ValueError(f"[mirror] pairs do not give a mirror permutation: {perm}")
    planes, rms = [], []
    for k in keep:
        Xk = frames[k].coords
        n, d = plane(Xk, perm)
        planes.append([round(float(v), 6) for v in n] + [round(d, 6)])
        rms.append(float(np.sqrt(np.mean(np.sum((reflect(Xk, n, d) - Xk[perm]) ** 2, axis=1)))))
    print(f"mirror plane: largest deviation from mirror symmetry along the path {max(rms):.3f} Å (rms)")
    return {"perm": perm, "planes": planes, "max_rms_A": round(max(rms), 4),
            "labels": {o["key"]: _mirror_label(o["label"], perm) for o in orbitals if o["kind"] == "localized"},
            "note": cfg.get("note", "Mirror image of the orbital set through the TS mirror plane: the "
                    "same mechanism drawn with the arrows the other way round the ring. Built by "
                    "reflecting the computed orbital surfaces; no new calculation.")}


def build_meshes(rx: Reaction, iso: float | None = None, max_tris: int = 900) -> None:
    iso = iso if iso is not None else rx.cfg.get("isovalue", 0.05)
    tr = _read_json(rx.root / "03_orbitals" / "tracking.json")
    frames = read_xyz_frames(rx.root / "02_irc" / "irc_full_trj.xyz")
    csv = (rx.root / "02_irc" / "energies.csv").read_text().splitlines()
    hdr = csv[0].split(",")
    rows = [dict(zip(hdr, ln.split(","))) for ln in csv[1:]]
    meta = _read_json(rx.root / "02_irc" / "irc_meta.json")
    # IRC points flagged by `profile` (above the TS right next to it) are not displayed
    excluded = {e["frame"] for e in meta.get("display_excluded_frames", [])}
    keep = [k for k in range(len(frames)) if k not in excluded]
    show = [t for t in tr["tracked"] if t["changes"]]
    exp = tr.get("expected_changes")
    if exp:      # reaction.toml lists the orbitals to show: exactly those, in that order
        ids = [e["tracked_id"] for e in exp["items"] if e["found"]]
        show = [next(t for t in tr["tracked"] if t["id"] == j) for j in ids]
    # one hue per orbital; phase = dark (+) vs light (−) shade of that hue.
    # Nucleophile-side orbital first (blue), then the bond that breaks (orange).
    def atoms_of(carriers: str) -> set[int]:          # "Br1–C0" -> {1, 0}
        return {int("".join(ch for ch in a if ch.isdigit())) for a in carriers.split("–")}
    show.sort(key=lambda t: 0 if rx.key("nucleophile") in atoms_of(t["start"]) else 1)
    loc_colors = [["#2f6db5", "#9cc3ea"], ["#d9661f", "#f4c08e"], ["#2e8b57", "#a8dcbc"],
                  ["#8e44ad", "#d7b5e6"]]
    orbitals = [{"key": f"ibo{t['id']}", "kind": "localized",
                 "label": t.get("label") or f"{t['label_start']} → {t['label_end']}", "tracked_id": t["id"],
                 "colors": loc_colors[i % len(loc_colors)]}
                for i, t in enumerate(show)]
    orbitals += [{"key": "homo", "kind": "canonical", "label": "canonical HOMO",
                  "colors": ["#7b3fa0", "#d4b8e6"]},
                 {"key": "lumo", "kind": "canonical", "label": "canonical LUMO",
                  "colors": ["#2e8b57", "#a8dcbc"]}]

    # global box for quantization: atoms ± 6 Å covers every lobe at this isovalue
    allxyz = np.concatenate([f.coords for f in frames])
    lo, hi = allxyz.min(0) - 6.0, allxyz.max(0) + 6.0
    scale = (hi - lo) / (2 * QUANT)
    mid = (hi + lo) / 2

    # Overall sign is arbitrary: choose it once per orbital (from frame 0) so the larger lobe is
    # the + (dark) shade. Tracking already keeps the phase continuous from frame to frame.
    def frame_sign(o, k):
        if o["kind"] == "localized":
            return tr["per_frame"][k]["sign"][o["tracked_id"]]
        return tr["canonical"][o["key"].upper()]["sign"][k]

    flip = {}
    for o in orbitals:
        v0, _, _ = read_cube_gz(rx.root / "03_orbitals/cubes/f000" / f"{o['key']}.cube.gz")
        skew = float(np.sum((frame_sign(o, 0) * v0.astype(np.float64)) ** 3))
        flip[o["key"]] = 1.0 if skew >= 0 else -1.0

    # display-only σ/π rotation of a bent-bond pair (orbital_steps.orbital_display): the listed
    # frames of that one orbital are read from the rotated copy; everything else as computed
    disp_path = rx.root / "03_orbitals" / "display_rotation.json"
    disp = _read_json(disp_path) if disp_path.exists() else None
    disp_frames = set(disp["rotated_frames"]) if disp else set()

    def values(k, o):
        cd = rx.root / f"03_orbitals/cubes/f{k:03d}"
        if disp and o["kind"] == "localized" and o["tracked_id"] == disp["shown_tracked_id"]                 and k in disp_frames:
            z = np.load(cd / f"{o['key']}_display.npz")
            return z["values"], z["origin"], z["axes"]
        return read_cube_gz(cd / f"{o['key']}.cube.gz")

    def encode(max_tris):
        buf = io.BytesIO()
        index = {o["key"]: [] for o in orbitals}
        sizes = {o["key"]: 0 for o in orbitals}
        for k in keep:
            for o in orbitals:
                vals, origin, axes = values(k, o)
                sign = frame_sign(o, k) * flip[o["key"]]
                entry = []
                for verts, faces in lobes(sign * vals, origin, axes, iso, max_tris):
                    q = np.round((verts - mid) / scale).astype("<i2")
                    off_v = buf.tell(); buf.write(q.tobytes())
                    off_f = buf.tell(); buf.write(faces.astype("<u2").tobytes())
                    pad = (-buf.tell()) % 4
                    buf.write(b"\0" * pad)
                    entry.append({"v": [off_v, int(len(q))], "f": [off_f, int(len(faces))]})
                    sizes[o["key"]] += q.nbytes + faces.size * 2
                index[o["key"]].append({"pos": entry[0], "neg": entry[1]})
        return buf, index, sizes

    # the web page must stay under 10 MB: lower the triangle budget until meshes.bin fits
    budget = float(rx.cfg.get("mesh_budget_mb", 8.5)) * 1e6
    buf, index, sizes = encode(max_tris)
    while buf.tell() > budget and max_tris > 200:
        new = int(max_tris * budget / buf.tell() * 0.97)
        print(f"meshes.bin {buf.tell() / 1e6:.2f} MB > {budget / 1e6:.1f} MB at max_tris {max_tris}; "
              f"retrying with {new}")
        max_tris = new
        buf, index, sizes = encode(max_tris)
    od = rx.d("04_meshes")
    (od / "meshes.bin").write_bytes(buf.getvalue())
    ts = _read_json(rx.root / "01_ts" / "ts.json")
    manifest = {
        "id": rx.cfg["id"], "title": rx.cfg["title"], "equation": rx.cfg["equation"],
        "level": f"{rx.level()}, RIJCOSX, SMD({rx.cfg['solvent']})",
        "max_triangles_per_lobe": max_tris,
        "isovalue": iso, "isovalue_units": "e^1/2 bohr^-3/2 (orbital amplitude)",
        "quant": {"mid": mid.tolist(), "scale": scale.tolist()},
        "key_atoms": rx.cfg["key_atoms"],
        "key_bonds": rx.key_bonds(frames[0].symbols),
        "view_axis": rx.cfg.get("view_axis"),            # optional [a, b]: "down the bond" camera
        "view_axis_label": rx.cfg.get("view_axis_label"),
        "ts_frame": keep.index(meta["ts_frame"]),
        "irc_frames_not_shown": meta.get("display_excluded_frames", []),
        "imaginary_frequency_cm": ts["imaginary_cm"][0],
        "energy_reference": "separated reactants",
        "E_separated_products_kJ": meta["E_separated_products_kJ"],
        "frames": [{"irc_frame": k, "s": float(rows[k]["s_bohr_amu05"]),
                    "E_kJ": float(rows[k]["E_kJ_per_mol"]),
                    "xyz": np.round(frames[k].coords, 4).tolist()} for k in keep],
        "curves": (curves := _curves(rx, rows, meta, keep)),
        "default_curve": default_curve(rx, curves),
        "symbols": frames[0].symbols,
        "orbitals": orbitals,
        "mesh_index": index,
        "bytes_per_orbital": sizes,
        "mirror": mirror_data(rx, frames, keep, meta["ts_frame"], orbitals),
        "display_rotation": ({k_: disp[k_] for k_ in ("atoms", "shown_tracked_id", "criterion",
                                                     "rotated_frames", "overlap_across_join")}
                             if disp else None),
    }
    _write_json(od / "manifest.json", manifest)
    tot = (od / "meshes.bin").stat().st_size
    print(f"meshes.bin {tot / 1e6:.2f} MB for {len(keep)} frames × {len(orbitals)} orbitals "
          f"({tot / len(keep) / len(orbitals) / 1e3:.1f} kB per orbital per frame)")
