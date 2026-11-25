"""Checks around the IBO step on the final M06-2X paths (J. Picione, 2026-11-25).

gate()            E2's path single points are complete: every DLPNO-CCSD(T) frame E2 must compute (its sampling, 507)
                  ends with ORCA TERMINATED NORMALLY, the level file is written, no E2 lane process is alive and no E2
                  job holds cores. If the E2 lane has ended without that (a failed frame, a missing file), the IBO
                  lanes hold and docs/IBO_HOLD.md says why; they never start while E2 holds cores.
check <rx>        Output check of the IBO single points before tracking: every IRC frame has its own clean folder with
                  f.out ending ORCA TERMINATED NORMALLY, the exported f_loc.npz and f_can.npz, and an SCF energy equal
                  to the working-level path energy of that frame (02_irc/energies.csv; same geometry, same method)
                  within 2e-5 Eh. Writes 03_orbitals/output_check.json. A failed check stops that reaction.
finish <rx>       check, then (only if it passes) orbital_track and orbital_cubes (which also runs orbital_display);
                  for the tertiary mirror path, then the mirror symmetry check. Idempotent: a step already done is not
                  redone. Writes 03_orbitals/status.json.
mirror <rx>       Symmetry check (J. Picione, 2026-11-23): every mirror frame m and its computed partner p (frame
                  2(c-1)-m) have their own IBOs; each shown orbital of m is compared with every shown orbital of p
                  mapped by the same C2 ⇄ C3 permutation and rotation that built the mirror frames (02_irc/
                  mirror_path.json, reaction.toml [mirror_path] perm). Overlap = ∫ψ_m ψ'_p / (‖ψ_m‖ ‖ψ'_p‖) on m's cube
                  grid, ψ'_p interpolated from p's cube (trilinear). Reports the best match per orbital and the worst
                  of those over all pairs. Writes 03_orbitals/mirror_symmetry.json/.md.
status            one line per reaction: single points done / check / tracking / cubes / mirror.

Usage: python -X utf8 -m pipeline.ibo_check gate|status|check <rx>|finish <rx>|mirror <rx>
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
E2 = "e2_ethoxide_2_bromopropane_m062x"
LEVEL = "DLPNO-CCSD(T)_ma-def2-TZVP"
HOLD = ROOT / "docs" / "IBO_HOLD.md"
ALL = ["e2_ethoxide_2_bromopropane_m062x", "nucadd_cyanide_acetone_m062x", "sn2_hydroxide_bromomethane_m062x",
       "hshift_2_3_dimethyl_2_butyl_m062x", "hshift_3_methyl_2_butyl_m062x", "da_butadiene_ethylene_m062x"]
E_TOL = 2e-5
BOHR = 0.529177210903


def _normal(out: Path) -> bool:
    return out.exists() and "ORCA TERMINATED NORMALLY" in out.read_text(errors="replace")[-4000:]


# ------------------------------------------------------------------ gate
def e2_state() -> dict:
    import psutil
    from .geometry import read_xyz_frames
    from .steps import Reaction, path_sp_frames
    rx = Reaction.load(ROOT / "reactions" / E2)
    n = len(read_xyz_frames(rx.root / "02_irc/irc_full_trj.xyz"))
    its = json.loads((rx.root / "02_irc/irc_meta.json").read_text())["ts_frame"]
    todo = path_sp_frames(rx, n, its, LEVEL)
    base = rx.root / "02_irc/path_sp" / LEVEL
    bad = [k for k in todo if not _normal(base / f"f{k:03d}" / "f.out")]
    alive = []
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        cl = " ".join(p.info["cmdline"] or [])
        if (p.info["name"] or "").lower().startswith(("bash", "python")) and E2 in cl and \
                ("run_singlepoints.sh" in cl or "-m pipeline reactions/" in cl):
            alive.append(p.info["pid"])
    held = []
    try:
        s = json.loads((ROOT / "docs/level_test_queue/allocation.json").read_text(encoding="utf-8"))
        held = [h["what"] for h in s.get("held", []) if E2 in (h.get("folder") or "") or E2 in h.get("what", "")]
    except (OSError, ValueError):
        held = ["allocation.json unreadable"]
    return {"frames_needed": len(todo), "frames_not_normal": bad, "level_file": (rx.root / "02_irc" /
            f"energies_{LEVEL}.json").exists(), "lane_processes": alive, "cores_held": held}


def gate(poll_s: int = 60) -> None:
    """Block until E2 is complete and gone. Holds (docs/IBO_HOLD.md) if the E2 lane ended without completing."""
    while True:
        st = e2_state()
        done = not st["frames_not_normal"] and st["level_file"] and not st["lane_processes"] and not st["cores_held"]
        if done:
            if HOLD.exists():
                HOLD.unlink()
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} gate open: all {st['frames_needed']} E2 frames ended "
                  f"normally, level file written, no E2 process, no E2 cores held", flush=True)
            return
        if not st["lane_processes"] and not st["cores_held"] and (st["frames_not_normal"] or not st["level_file"]):
            msg = (f"# IBO lanes on hold ({time.strftime('%Y-%m-%d %H:%M')})\n\nThe E2 single-point lane has ended "
                   f"but E2 is not complete: {len(st['frames_not_normal'])} of {st['frames_needed']} frames without "
                   f"ORCA TERMINATED NORMALLY (first: {st['frames_not_normal'][:10]}); level file written: "
                   f"{st['level_file']}. The IBO lanes do not start until this is resolved (the restart guard may "
                   f"restart the E2 lane; finished frames are reused). Report it to J. Picione.\n")
            if not HOLD.exists() or HOLD.read_text(encoding="utf-8").split("\n", 2)[-1] != msg.split("\n", 2)[-1]:
                HOLD.write_text(msg, encoding="utf-8")
                print(msg, flush=True)
        time.sleep(poll_s)


# ------------------------------------------------------------------ output check
def check(name: str) -> dict:
    from .geometry import read_xyz_frames
    from .steps import Reaction
    rx = Reaction.load(ROOT / "reactions" / name)
    n = len(read_xyz_frames(rx.root / "02_irc/irc_full_trj.xyz"))
    rows = (rx.root / "02_irc/energies.csv").read_text().splitlines()
    hdr = rows[0].split(",")
    ie = hdr.index("E_eh")
    e_path = [float(r.split(",")[ie]) for r in rows[1:]]
    problems, worst = [], 0.0
    for k in range(n):
        d = rx.root / f"03_orbitals/frames/f{k:03d}"
        out = d / "f.out"
        if not _normal(out):
            problems.append(f"frame {k}: no normal termination")
            continue
        for f in ("f_loc.npz", "f_can.npz"):
            if not (d / f).exists():
                problems.append(f"frame {k}: {f} missing")
        m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out.read_text(errors="replace"))
        if not m:
            problems.append(f"frame {k}: no energy")
            continue
        de = abs(float(m[-1]) - e_path[k])
        worst = max(worst, de)
        if de > E_TOL:
            problems.append(f"frame {k}: SCF energy differs from the path energy by {de:.2e} Eh")
    res = {"reaction": name, "frames": n, "pass": not problems, "problems": problems[:50],
           "n_problems": len(problems), "max_abs_dE_vs_path_eh": worst, "tolerance_eh": E_TOL,
           "checked": time.strftime("%Y-%m-%d %H:%M")}
    (rx.root / "03_orbitals").mkdir(exist_ok=True)
    (rx.root / "03_orbitals/output_check.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res


def sps_complete(name: str) -> tuple[int, int]:
    from .geometry import read_xyz_frames
    R = ROOT / "reactions" / name
    n = len(read_xyz_frames(R / "02_irc/irc_full_trj.xyz"))
    done = sum(_normal(R / f"03_orbitals/frames/f{k:03d}/f.out") for k in range(n))
    return done, n


# ------------------------------------------------------------------ mirror symmetry
def _interp(vals, origin, axes, pts_bohr):
    """Trilinear interpolation of a cube (axis-aligned grid) at points (bohr); 0 outside."""
    from scipy.ndimage import map_coordinates
    step = np.diag(axes)
    if np.abs(axes - np.diag(step)).max() > 1e-8:
        raise ValueError("cube axes are not axis-aligned")
    idx = ((pts_bohr - origin) / step).T
    return map_coordinates(vals.astype(np.float64), idx, order=1, mode="constant", cval=0.0)


def mirror(name: str) -> dict:
    from .geometry import read_xyz_frames
    from .orbitals import read_cube
    from .steps import Reaction, mirror_partners
    from .mirror_path import _kabsch
    import gzip
    import tempfile
    import os
    rx = Reaction.load(ROOT / "reactions" / name)
    frames = read_xyz_frames(rx.root / "02_irc/irc_full_trj.xyz")
    pairs = mirror_partners(rx, len(frames))
    if not pairs:
        raise RuntimeError(f"{name} has no mirror path")
    perm = list(rx.cfg["mirror_path"]["perm"])
    tr = json.loads((rx.root / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
    keys = [f"ibo{t['id']}" for t in tr["tracked"] if t["changes"]] + ["homo", "lumo"]

    def cube(k, key):
        p = rx.root / f"03_orbitals/cubes/f{k:03d}/{key}.cube.gz"
        with tempfile.NamedTemporaryFile(delete=False, suffix=".cube") as fh:
            fh.write(gzip.open(p, "rb").read())
            tmp = fh.name
        try:
            v, o, a = read_cube(Path(tmp))[:3]
        finally:
            os.unlink(tmp)
        return v, o, a

    rows = []
    for m, p in sorted(pairs.items()):
        xm = np.asarray(frames[m].coords)
        xp = np.asarray(frames[p].coords)
        R, c0, c1 = _kabsch(xp[perm], xm)            # (xp[perm] − c0) R + c1 ≈ xm
        fit = float(np.sqrt((((xp[perm] - c0) @ R + c1 - xm) ** 2).sum(1).mean()))
        cm = {kk: cube(m, kk) for kk in keys}
        cp = {kk: cube(p, kk) for kk in keys}
        for ka in keys:
            va, oa, aa = cm[ka]
            g = np.stack(np.meshgrid(*[np.arange(s) for s in va.shape], indexing="ij"), -1).reshape(-1, 3)
            pts = oa + g * np.diag(aa)                                   # bohr, m's frame
            back = ((pts * BOHR - c1) @ R.T + c0) / BOHR                 # same points in p's frame
            a = va.reshape(-1).astype(np.float64)
            best, arg = 0.0, None
            for kb in keys:
                vb, ob, ab = cp[kb]
                b = _interp(vb, ob, ab, back)
                ov = float(a @ b / np.sqrt((a @ a) * (b @ b))) if b.any() else 0.0
                if abs(ov) > abs(best):
                    best, arg = ov, kb
            rows.append({"frame": m, "partner": p, "orbital": ka, "partner_orbital": arg,
                         "overlap": round(abs(best), 4), "geometry_fit_A": round(fit, 5)})
    worst = min(rows, key=lambda r: r["overlap"])
    res = {"reaction": name, "pairs": len(pairs), "orbitals_compared": keys, "worst": worst,
           "median_overlap": float(np.median([r["overlap"] for r in rows])),
           "max_geometry_fit_A": max(r["geometry_fit_A"] for r in rows),
           "method": "cube overlap on the mirror frame's grid, partner cube interpolated (trilinear) after the "
                     "C2⇄C3 permutation and rotation that built the mirror frames", "rows": rows}
    (rx.root / "03_orbitals/mirror_symmetry.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    md = [f"# Mirror symmetry of the IBOs, {name}", "",
          f"{len(pairs)} mirror frames against their computed partners; orbitals {keys}.",
          f"- Worst best-match overlap: {worst['overlap']:.4f} (frame {worst['frame']}, {worst['orbital']} vs partner "
          f"frame {worst['partner']} {worst['partner_orbital']}).",
          f"- Median best-match overlap: {res['median_overlap']:.4f}.",
          f"- Largest geometry misfit of the mapping: {res['max_geometry_fit_A']:.5f} Å.",
          f"- Method: {res['method']}. Grid interpolation limits the overlap of identical orbitals to slightly below 1."]
    (rx.root / "03_orbitals/mirror_symmetry.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return res


# ------------------------------------------------------------------ finish / status
def finish(name: str) -> dict:
    from .orbital_steps import orbital_cubes, orbital_track
    from .steps import Reaction
    R = ROOT / "reactions" / name
    st_path = R / "03_orbitals/status.json"
    st = json.loads(st_path.read_text(encoding="utf-8")) if st_path.exists() else {}
    done, n = sps_complete(name)
    st.update(reaction=name, sps=f"{done} of {n}")
    if done < n:
        st["stage"] = "single points running"
    else:
        c = check(name)
        st["check"] = "pass" if c["pass"] else f"FAIL ({c['n_problems']} problems; see output_check.json)"
        if not c["pass"]:
            st["stage"] = "STOPPED: output check failed, not tracked"
        else:
            rx = Reaction.load(R)
            if not (R / "03_orbitals/tracking.json").exists():
                orbital_track(rx)
            tr = json.loads((R / "03_orbitals/tracking.json").read_text(encoding="utf-8"))
            exp = tr.get("expected_changes") or {}
            st["tracking"] = {"expected": [{k: e[k] for k in ("label", "found", "min_overlap")}
                                           for e in exp.get("items", [])],
                              "unexpected": exp.get("unexpected"), "flagged_frames": exp.get("flagged_frames")}
            show = [t["id"] for t in tr["tracked"] if t["changes"]]
            last = R / f"03_orbitals/cubes/f{n - 1:03d}"
            if not all((last / f"{k}.cube.gz").exists() for k in [f"ibo{j}" for j in show] + ["homo", "lumo"]):
                orbital_cubes(rx)
            st["cubes"] = "done"
            st["stage"] = "tracked, cubes done"
            if rx.cfg.get("mirror_path"):
                if not (R / "03_orbitals/mirror_symmetry.json").exists():
                    mirror(name)
                ms = json.loads((R / "03_orbitals/mirror_symmetry.json").read_text(encoding="utf-8"))
                st["mirror_check"] = {"worst_overlap": ms["worst"]["overlap"], "at": ms["worst"],
                                      "median": ms["median_overlap"]}
                st["stage"] = "tracked, cubes done, mirror check done"
    st["updated"] = time.strftime("%Y-%m-%d %H:%M")
    st_path.parent.mkdir(exist_ok=True)
    st_path.write_text(json.dumps(st, indent=1, ensure_ascii=False), encoding="utf-8")
    return st


def status() -> None:
    if HOLD.exists():
        print("IBO_HOLD:", HOLD.read_text(encoding="utf-8").strip().splitlines()[-1][:300])
    for name in ALL:
        done, n = sps_complete(name)
        sp = ROOT / "reactions" / name / "03_orbitals/status.json"
        st = json.loads(sp.read_text(encoding="utf-8")) if sp.exists() else {}
        print(f"{name}: single points {done} of {n}; {st.get('stage', 'not checked yet')}"
              + (f"; check {st['check']}" if "check" in st else "")
              + (f"; mirror worst overlap {st['mirror_check']['worst_overlap']}" if "mirror_check" in st else ""))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "gate":
        print(json.dumps(e2_state(), indent=1))
    elif cmd == "status":
        status()
    elif cmd in ("check", "finish", "mirror"):
        print(json.dumps({"check": check, "finish": finish, "mirror": mirror}[cmd](sys.argv[2]), indent=1,
                         ensure_ascii=False, default=str)[:4000])
    else:
        sys.exit(__doc__)
