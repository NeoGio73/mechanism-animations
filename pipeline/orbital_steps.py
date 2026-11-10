"""Orbital steps: IBO single points at every IRC frame, tracking, cube files."""
from __future__ import annotations

import gzip
import re
import shutil
from pathlib import Path

import numpy as np

from . import orca
from .geometry import read_xyz_frames
from .orbitals import (ao_labels, export_json, load_npz, localized_range, overlap, plot_cubes,
                       read_cube, s_character, track)
from .steps import Reaction, _read_json, _write_json

LOC_BLOCK = "%loc\n  LocMet IAOIBO\nend"
POP_CUT = 0.15      # e; an atom "carries" an orbital when its Mulliken population exceeds this
PI_MAX_S = 0.10     # a two-centre orbital with less s-character than this is a π bond


def _frames(rx: Reaction):
    return read_xyz_frames(rx.root / "02_irc" / "irc_full_trj.xyz")


def _fdir(rx: Reaction, k: int) -> Path:
    return rx.d(f"03_orbitals/frames/f{k:03d}")


def orbital_sps(rx: Reaction) -> None:
    """Working-level single point with IBO localization at every IRC frame."""
    frames = _frames(rx)
    flog = orca.RunLog(rx.d("03_orbitals") / "run_log_frames.md")
    total = 0.0
    for k, g in enumerate(frames):
        wd = _fdir(rx, k)
        r = orca.run(wd, "f", rx.inp("", g, blocks=LOC_BLOCK), flog, note=f"IRC frame {k}")
        if not r.ok:
            raise orca.OrcaError(f"frame {k} single point failed; see {r.out}")
        if not (wd / "f_loc.npz").exists() or not (wd / "f_can.npz").exists():
            export_json(wd / "f.loc")
            export_json(wd / "f.gbw")
        for tmp in wd.glob("*.tmp*"):
            tmp.unlink()
        if r.wall_s == r.wall_s:      # not NaN (i.e. not reused)
            total += r.wall_s
    rx.log.note(f"IBO single points at {len(frames)} IRC frames (`03_orbitals/frames/`, per-frame "
                f"log in `03_orbitals/run_log_frames.md`): {total:.0f} s in total this run.")


def _carriers(pop_row: np.ndarray, symbols: list[str]) -> str:
    idx = [a for a in np.argsort(-pop_row) if pop_row[a] > POP_CUT]
    return "–".join(f"{symbols[a]}{a}" for a in idx)


def _label(pop_row: np.ndarray, symbols: list[str], s_frac: float | None = None,
           partner_s: list[float] = ()) -> str:
    """Chemical label from atomic populations: lone pair / σ, π or bent bond between two atoms.

    s_frac: s-character of this orbital; partner_s: s-character of the other orbitals on the
    same two atoms in the same frame. π: s-character below PI_MAX_S. A multiple bond whose
    orbitals all have s-character above PI_MAX_S has been localized as bent bonds (σ and π
    mixed; the two need not be equivalent)."""
    idx = [a for a in np.argsort(-pop_row) if pop_row[a] > POP_CUT]
    # organic convention: carbon first, hydrogen last, others alphabetical (σ(C–Br), σ(O–H))
    rank = lambda s: (0 if s == "C" else 2 if s == "H" else 1, s)
    names = sorted((symbols[a] for a in idx), key=rank)
    if len(idx) == 1:
        return f"lone pair on {names[0]}"
    if len(idx) == 2:
        if s_frac is not None and s_frac < PI_MAX_S:
            return f"π({names[0]}={names[1]})"
        if s_frac is not None and partner_s and min(partner_s) >= PI_MAX_S:
            return f"bent bond ({names[0]}={names[1]})"
        return f"σ({names[0]}–{names[1]})"
    return "delocalized over " + "–".join(names)


def _end_labels(F, perm_k, pops_k, syms, is_s) -> list[str]:
    """Labels for every tracked orbital in one frame (s-character decides σ / π / bent)."""
    m = len(perm_k)
    atoms = [frozenset(int(a) for a in np.where(pops_k[j] > POP_CUT)[0]) for j in range(m)]
    sfrac = [s_character(F, perm_k[j], is_s) for j in range(m)]
    return [_label(pops_k[j], syms, sfrac[j],
                   [sfrac[i] for i in range(m) if i != j and atoms[i] == atoms[j]])
            for j in range(m)]


def orbital_track(rx: Reaction) -> None:
    """Match IBOs between neighbouring frames, fix signs, pick the orbitals that change."""
    frames = _frames(rx)
    syms = frames[0].symbols
    F = [load_npz(_fdir(rx, k) / "f_loc.npz", syms) for k in range(len(frames))]
    Fc = [load_npz(_fdir(rx, k) / "f_can.npz", syms) for k in range(len(frames))]
    idx = []
    for k in range(len(frames)):
        a, b = localized_range((_fdir(rx, k) / "f.out").read_text(errors="replace"))
        idx.append(list(range(a, b + 1)))
    perm, sign, best = track(F, idx, ref=0)

    n = len(frames)
    m = len(idx[0])
    pops = np.array([[F[k].populations([perm[k][j]])[0] for j in range(m)] for k in range(n)])
    labs = ao_labels(_fdir(rx, 0) / "f.loc", rx.d("03_orbitals") / "ao_labels.json")
    is_s = np.array([lab.split()[1].endswith("s") for lab in labs])
    lab0 = _end_labels(F[0], perm[0], pops[0], syms, is_s)
    lab1 = _end_labels(F[-1], perm[-1], pops[-1], syms, is_s)
    tracked = []
    for j in range(m):
        start, end = _carriers(pops[0, j], syms), _carriers(pops[-1, j], syms)
        tracked.append({"id": j, "start": start, "end": end,
                        "label_start": lab0[j],
                        "label_end": lab1[j],
                        "changes": set(start.split("–")) != set(end.split("–")),
                        "min_overlap": float(min(best[k][j] for k in range(n))),
                        "frame_of_min_overlap": int(np.argmin([best[k][j] for k in range(n)]))})
    expected = _check_expected(rx, tracked, pops, best, syms)
    # canonical HOMO / LUMO: fixed index, sign from overlap with the previous frame
    nocc = int((Fc[0].occ > 1.0).sum())
    can = {"HOMO": nocc - 1, "LUMO": nocc}
    can_sign = {k: [1.0] for k in can}
    can_ovl = {k: [1.0] for k in can}
    for f in range(1, n):
        for name, i in can.items():
            o = float(overlap(Fc[f - 1], Fc[f], [i], [i])[0, 0]) * can_sign[name][-1]
            can_sign[name].append(1.0 if o >= 0 else -1.0)
            can_ovl[name].append(abs(o))
    out = {
        "isovalue_note": "orbitals are real; sign = phase applied after tracking",
        "population_cut_e": POP_CUT,
        "n_frames": n,
        "tracked": tracked,
        "expected_changes": expected,
        "per_frame": [{"frame": k,
                       "orbital_index": [int(v) for v in perm[k]],
                       "sign": sign[k],
                       "overlap_with_previous": best[k],
                       "populations": {j: {f"{syms[a]}{a}": round(float(pops[k, j, a]), 3)
                                           for a in range(len(syms)) if pops[k, j, a] > 0.05}
                                       for j in range(m)}}
                      for k in range(n)],
        "canonical": {name: {"orbital_index": i, "sign": can_sign[name],
                             "overlap_with_previous": can_ovl[name]} for name, i in can.items()},
    }
    _write_json(rx.d("03_orbitals") / "tracking.json", out)
    lines = ["# Orbital tracking log", "",
             f"{n} IRC frames, {m} localized valence IBOs per frame, matched by maximum |overlap| "
             "between neighbouring frames (C_k^T S C_k+1). Signs fixed so matched overlaps are "
             "positive. Frames where a best match falls below 0.7 are flagged.",
             f"Bond type from the s-character of each two-atom orbital: π below {PI_MAX_S:.2f}; when "
             "every orbital on a multiple bond is above it, the bond was localized as bent bonds "
             "(σ and π mixed; the two need not be equivalent). Atom numbers in brackets are 0-based.", "",
             "| Tracked orbital | Reactant end | Product end | Changes? | Lowest overlap (frame) |",
             "|---|---|---|---|---|"]
    for t in tracked:
        lines.append(f"| {t['id']} | {t['label_start']} ({t['start']}) | {t['label_end']} "
                     f"({t['end']}) | {'**yes**' if t['changes'] else 'no'} | "
                     f"{t['min_overlap']:.3f} ({t['frame_of_min_overlap']}) |")
    flags = [(k, j, best[k][j]) for k in range(n) for j in range(m) if best[k][j] < 0.7]
    lines += ["", f"Frames with a match below 0.7: {len(flags)}"]
    lines += [f"- frame {k}, orbital {j}: overlap {v:.3f}" for k, j, v in flags[:60]]
    if expected:
        lines += ["", "## Expected orbital changes (reaction.toml)", "",
                  "| Expected | Tracked orbital | Lowest overlap | Found? |", "|---|---|---|---|"]
        for e in expected["items"]:
            lines.append(f"| {e['label']} | {e['tracked_id'] if e['tracked_id'] is not None else '—'} | "
                         f"{e['min_overlap'] if e['min_overlap'] is not None else '—'} | "
                         f"{'yes' if e['found'] else '**NO**'} |")
        lines += ["", f"Orbitals that change but were not expected: {expected['unexpected'] or 'none'}",
                  f"Frames where an unexpected orbital changes its atoms: {expected['flagged_frames'] or 'none'}"]
    lines += ["", "Canonical orbitals (sign tracking only): lowest overlap "
              f"HOMO {min(can_ovl['HOMO']):.3f}, LUMO {min(can_ovl['LUMO']):.3f}"]
    (rx.root / "03_orbitals" / "tracking_log.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:12 + m]))


def _check_expected(rx: Reaction, tracked, pops, best, syms):
    """Compare the orbitals that change with `[[expected_orbital_changes]]` in reaction.toml
    (start/end = atom indices carrying the orbital at the reactant/product end).

    Each expected change must be ONE tracked orbital (continuous from end to end). Any frame
    where another orbital's set of atoms differs from its reactant-end set is flagged."""
    exp = rx.cfg.get("expected_orbital_changes")
    if not exp:
        return None
    n, m = pops.shape[0], pops.shape[1]
    atoms = lambda row: frozenset(int(a) for a in np.where(row > POP_CUT)[0])
    items, used = [], set()
    for e in exp:
        s, f = frozenset(e["start"]), frozenset(e["end"])
        hit = [j for j in range(m) if atoms(pops[0, j]) == s and atoms(pops[-1, j]) == f]
        j = hit[0] if len(hit) == 1 else None
        if j is not None:
            used.add(j)
            tracked[j]["label"] = e["label"]
        items.append({"label": e["label"], "start": sorted(s), "end": sorted(f), "found": j is not None,
                      "tracked_id": j, "candidates": hit,
                      "min_overlap": round(float(min(best[k][j] for k in range(n))), 3) if j is not None else None})
    unexpected = [t["id"] for t in tracked if t["changes"] and t["id"] not in used]
    flagged = [k for k in range(n)
               if any(atoms(pops[k, j]) != atoms(pops[0, j]) for j in range(m) if j not in used)]
    for e in items:
        if not e["found"]:
            print(f"WARNING expected orbital change not found as one tracked orbital: {e['label']} "
                  f"(candidates {e['candidates']})")
    if unexpected:
        print(f"WARNING orbitals change that were not expected: {unexpected}")
    if flagged:
        print(f"WARNING frames where an unexpected orbital changes its atoms: {flagged}")
    return {"items": items, "unexpected": unexpected, "flagged_frames": flagged}


def orbital_cubes(rx: Reaction, grid: int = 70) -> None:
    """Cube files for the changing IBOs and canonical HOMO/LUMO at every frame (gzipped)."""
    tr = _read_json(rx.root / "03_orbitals" / "tracking.json")
    show = [t["id"] for t in tr["tracked"] if t["changes"]]
    for k, pf in enumerate(tr["per_frame"]):
        wd = _fdir(rx, k)
        cd = rx.d(f"03_orbitals/cubes/f{k:03d}")
        want = {f"ibo{j}": pf["orbital_index"][j] for j in show}
        want.update({name.lower(): c["orbital_index"] for name, c in tr["canonical"].items()})
        if all((cd / f"{name}.cube.gz").exists() for name in want):
            continue
        for name, src in (("loc", wd / "f.loc"), ("gbw", wd / "f.gbw")):
            orbs = {n_: i for n_, i in want.items() if n_.startswith("ibo") == (name == "loc")}
            made = plot_cubes(src, sorted(set(orbs.values())), grid=grid)
            by_index = {int(re.search(r"\.mo(\d+)a\.cube$", p.name).group(1)): p for p in made}
            for n_, i in orbs.items():
                with open(by_index[i], "rb") as fi, gzip.open(cd / f"{n_}.cube.gz", "wb", 6) as fo:
                    shutil.copyfileobj(fi, fo)
            for p in made:
                p.unlink()
    print(f"cubes written for {len(tr['per_frame'])} frames: {['ibo%d' % j for j in show]} + homo, lumo")
    orbital_display(rx, grid=grid)


DISPLAY_FILE = "display_rotation.json"


def orbital_display(rx: Reaction, grid: int = 70) -> None:
    """Display-only σ/π rotation of a double bond that IBO localized as two bent bonds.

    `[display_rotation]` in reaction.toml names the two atoms (`atoms`) and a threshold angle
    (`max_angle_deg`). The two tracked orbitals that sit on exactly those atoms at the reactant
    end are a 2×2 orthonormal pair; the combination with the least s-character is the π
    orbital, the other the σ orbital. Bend angle = angle between the shown orbital and that π
    combination. From the reactant end, every frame with a bend angle above the threshold is
    shown as the π combination instead, until the angle first comes within the threshold.

    A rotation inside an occupied pair leaves the total density, and so every energy, exactly
    unchanged. tracking.json, the .loc files and the as-computed cubes are not modified: the
    rotated orbital is written next to them as `cubes/fNNN/ibo<id>_display.npz`, and the
    choice, angles and overlaps go to `03_orbitals/display_rotation.json`."""
    from .meshes import read_cube_gz
    cfg = rx.cfg.get("display_rotation")
    out = rx.root / "03_orbitals" / DISPLAY_FILE
    if not cfg:
        if out.exists():
            out.unlink()
        return
    tr = _read_json(rx.root / "03_orbitals" / "tracking.json")
    syms = _frames(rx)[0].symbols
    want = frozenset(cfg["atoms"])
    on_pair = [t["id"] for t in tr["tracked"]
               if frozenset(int(a[1:]) for a in tr["per_frame"][0]["populations"][str(t["id"])]
                            if tr["per_frame"][0]["populations"][str(t["id"])][a] > POP_CUT) == want]
    shown_ids = {e["tracked_id"] for e in (tr.get("expected_changes") or {}).get("items", [])}
    shown = [j for j in on_pair if j in shown_ids]
    if len(on_pair) != 2 or len(shown) != 1:
        raise ValueError(f"display_rotation: expected 2 tracked orbitals on atoms {sorted(want)} "
                         f"with 1 shown, found {on_pair} (shown {shown})")
    sh, pa = shown[0], next(j for j in on_pair if j != shown[0])
    labs = ao_labels(_fdir(rx, 0) / "f.loc", rx.d("03_orbitals") / "ao_labels.json")
    is_s = np.array([lab.split()[1].endswith("s") for lab in labs])
    Ds = np.diag(is_s.astype(float))
    limit = float(cfg["max_angle_deg"])

    def pair(k):
        pf = tr["per_frame"][k]
        F = load_npz(_fdir(rx, k) / "f_loc.npz", syms)
        sg = np.array([pf["sign"][sh], pf["sign"][pa]])
        C = F.C[:, [pf["orbital_index"][sh], pf["orbital_index"][pa]]] * sg
        Q = C.T @ (0.5 * (Ds @ F.S + F.S @ Ds)) @ C          # s-character matrix of the pair
        w, v = np.linalg.eigh(Q)
        u = v[:, 0] * np.sign(v[0, 0])                        # π combination, same phase as shown
        ang = float(np.degrees(np.arccos(min(1.0, abs(v[0, 0])))))
        return F, C, u, ang, Q

    rows, rotated, P = [], [], {}
    k = 0
    while True:
        F, C, u, ang, Q = pair(k)
        P[k] = (F, C, u)
        rot = ang > limit
        rows.append({"frame": k, "bend_angle_deg": round(ang, 2),
                     "s_character_shown": round(float(Q[0, 0]), 3),
                     "s_character_partner": round(float(Q[1, 1]), 3),
                     "rotated": rot,
                     "coefficients_shown_partner": [round(float(x), 6) for x in u] if rot else None})
        if not rot:
            break
        rotated.append(k)
        k += 1
    vec = lambda k: P[k][1] @ (P[k][2] if k in rotated else np.array([1.0, 0.0]))
    ov = lambda a, b: float(vec(a) @ (0.5 * (P[a][0].S + P[b][0].S)) @ vec(b))
    steps = [round(ov(k, k + 1), 5) for k in rotated]
    # cubes: rotated = a·(shown) + b·(partner) with tracking signs, stored in the shown orbital's
    # raw phase so meshes applies the same per-frame sign as before
    for k in rotated:
        pf = tr["per_frame"][k]
        cd = rx.root / f"03_orbitals/cubes/f{k:03d}"
        made = plot_cubes(_fdir(rx, k) / "f.loc", [pf["orbital_index"][pa]], grid=grid)
        vp, op_, ap = read_cube(made[0])
        for p in made:
            p.unlink()
        vs, os_, as_ = read_cube_gz(cd / f"ibo{sh}.cube.gz")
        if not (np.allclose(os_, op_) and np.allclose(as_, ap)):
            raise ValueError(f"frame {k}: partner cube is on a different grid")
        a, b = P[k][2]
        disp = pf["sign"][sh] * (a * pf["sign"][sh] * vs + b * pf["sign"][pa] * vp)
        np.savez_compressed(cd / f"ibo{sh}_display.npz", values=disp.astype(np.float32),
                            origin=os_, axes=as_)
    res = {
        "what": "display-only rotation of one bent-bond pair into σ + π; data files unchanged",
        "atoms": sorted(want), "shown_tracked_id": sh, "partner_tracked_id": pa,
        "criterion": (f"bend angle (shown orbital vs. the least-s-character combination of the "
                      f"pair) above {limit:g}°, from the reactant end until it first comes within "
                      f"{limit:g}°"),
        "max_angle_deg": limit,
        "rotated_frames": rotated,
        "first_unrotated_frame": rows[-1]["frame"],
        "lowest_overlap_between_rotated_frames": min(steps[:-1]) if len(steps) > 1 else None,
        "overlap_across_join": steps[-1] if steps else None,
        "join": [rotated[-1], rows[-1]["frame"]] if rotated else None,
        "frames": rows,
    }
    _write_json(out, res)
    print(f"display rotation: frames {rotated[0]}–{rotated[-1]} rotated (bend "
          f"{rows[0]['bend_angle_deg']}° at frame 0); join {res['join']} overlap "
          f"{res['overlap_across_join']}; lowest step inside {res['lowest_overlap_between_rotated_frames']}")
