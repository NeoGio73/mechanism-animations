"""Explicit-water (cluster-continuum) sensitivity test, cyanide + acetone (J. Picione, 2026-11-20).

Question: do the errors in the cyanide reaction free energy (+25.9 kJ/mol) and barrier (+9.7 kJ/mol, DLPNO-CCSD(T)/
ma-def2-TZVP against Schlesinger and Miller 1973) come from SMD under-solvating the alkoxide oxygen?

Model: n = 1, 2, 3 explicit waters, the same n on both sides, everything also in SMD water as in production:
    CN⁻·(H2O)n + acetone  →  TS·(H2O)n  →  alkoxide·(H2O)n
Acetone stays bare (its production structure and energies). For each species and n a few starting arrangements are
built (waters on the alkoxide/carbonyl oxygen, on the nitrile nitrogen, mixed; on cyanide at N and at C), optimized at
the working level (Opt Freq; OptTS Freq for the TS with the production TS settings), and the lowest in G is kept. A
cluster TS must have one imaginary frequency whose motion on the solute atoms is the bare TS mode (overlap ≥ 0.8 after
aligning the solute) and must form C–C while C–O lengthens. The kept structures get single points at every
[[diagnostic_levels]] entry; G(level) = E(level) + (G − E) at the working level, 1 atm → 1 M as in production.
Production numbers stay as they are; this is a sensitivity test on stationary points only.

Scheduling (his order): behind the path single points, so that E2 is not delayed. `run` waits until every running
reaction has its DLPNO-CCSD(T) path file and no single-point lane is left, then runs one test optimization (the
largest minimum, alkoxide·(H2O)3 with three waters on O), estimates the job time of the whole test from it, and stops
before the rest if that estimate is over the limit (04_explicit_water/STOP_ESTIMATE.md). The estimate is summed
job time (s at 4 cores); two jobs run at a time, so the clock time is about half of it.

J. Picione, 2026-11-26: limit 38 h (was 36 h; the 36.8 h estimate stopped the run on 2026-11-25 20:14); start at
07:00; run in order of n (every n = 1 structure, then the n = 1 single points as soon as its structures are kept, then
n = 2, then n = 3), so that the smaller clusters are complete if anything stops early. Not stopped at Friday 20:00
(progress.md then says what is left and the new estimate); hard stop Saturday 2026-11-28 08:00 whatever is left
(a running job is ended; it is not a result). Nothing starts after this test.

Usage: python -X utf8 -m pipeline.explicit_water <rx> build | estimate | run [--now] [--not-before T] [--hard-stop T] | md
  T as YYYY-MM-DDTHH:MM, local time
Output: <rx>/04_explicit_water/ (starts/, opt/<species>_w<n>/<arr>/, sp/<level>/<species>_w<n>/, explicit_water.json/.md)
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca
from .geometry import Geometry, read_xyz
from .steps import STD_STATE_KJ, Reaction, _read_json, _write_json, reference_names
from .units import kj

SUB = "04_explicit_water"
OH, HOH = 0.96, np.radians(104.5)
DIST = {"O": 2.75, "N": 2.95, "C": 3.15}           # acceptor···O(water) at the start, Å
CONE = {("O", 1): 55, ("O", 2): 65, ("O", 3): 70,  # half-angle of the cone around the acceptor's lone-pair axis, °
        ("N", 1): 0, ("N", 2): 40, ("N", 3): 50, ("C", 1): 0, ("C", 2): 40, ("C", 3): 50}
ARRANGEMENTS = {       # sites of the n waters
    "cyanide": {1: ["N", "C"], 2: ["NN", "CC", "NC"], 3: ["NNN", "NNC", "NCC", "CCC"]},
    "ts": {1: ["O", "N"], 2: ["OO", "NN", "ON"], 3: ["OOO", "OON", "ONN"]},
    "alkoxide": {1: ["O", "N"], 2: ["OO", "NN", "ON"], 3: ["OOO", "OON", "ONN"]},
}
TEST = ("alkoxide", 3, "OOO")
LIMIT_H = 38.0          # J. Picione 2026-11-26 (36.0 before)
PATH_RX = ["da_butadiene_ethylene_m062x", "e2_ethoxide_2_bromopropane_m062x",
           "hshift_2_3_dimethyl_2_butyl_m062x", "hshift_3_methyl_2_butyl_m062x"]
DLPNO = "DLPNO-CCSD(T)_ma-def2-TZVP"
MODE_MIN_OVERLAP = 0.8


# ------------------------------------------------------------------ bare structures
def bare(rx: Reaction) -> dict:
    """species -> (Geometry, charge, sites {symbol: atom index}); the solute atoms come first in every cluster."""
    frag = _read_json(rx.root / "01_ts/fragments.json")
    fx = lambda n: read_xyz(rx.root / f"01_ts/fragments/{n}/{frag[n].get('job', n)}.xyz")
    out = {"cyanide": fx("cyanide"), "ts": read_xyz(Path(f"{rx.ts_stem()}.xyz")),
           "alkoxide": fx("acetone_cyanohydrin_alkoxide")}
    res = {}
    for k, g in out.items():
        s = g.symbols
        sites = {"N": s.index("N")}
        if "O" in s:
            sites["O"] = s.index("O")
        else:                                                   # cyanide: the C end
            sites["C"] = s.index("C")
        res[k] = (g, -1, sites)
    return res


def _neighbors(g: Geometry, i: int) -> list[int]:
    return [j for j in range(len(g.symbols)) if j != i and np.linalg.norm(g.coords[i] - g.coords[j]) < 1.75]


def _perp(v: np.ndarray) -> np.ndarray:
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    p = np.cross(v, a)
    return p / np.linalg.norm(p)


def _rot(v: np.ndarray, axis: np.ndarray, ang: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    return v * np.cos(ang) + np.cross(axis, v) * np.sin(ang) + axis * np.dot(axis, v) * (1 - np.cos(ang))


def _water(acc: np.ndarray, ow: np.ndarray, others: np.ndarray, centroid: np.ndarray) -> np.ndarray:
    """O, H (donating to the acceptor), H (free, turned for the largest distance to the other atoms)."""
    e1 = (acc - ow) / np.linalg.norm(acc - ow)
    h1 = ow + OH * e1
    best, bd = None, -1.0
    for k in range(24):
        p = _rot(_perp(e1), e1, 2 * np.pi * k / 24)
        h2 = ow + OH * (np.cos(HOH) * e1 + np.sin(HOH) * p)
        d = np.min(np.linalg.norm(others - h2, axis=1)) + 0.01 * np.linalg.norm(h2 - centroid)
        if d > bd:
            best, bd = h2, d
    return np.array([ow, h1, best])


def build_one(g: Geometry, sites: dict, arr: str) -> tuple[Geometry, float]:
    """Waters on the sites in arr (one letter per water); returns the cluster and the closest contact of a water
    atom with any atom it is not hydrogen-bonded to (Å)."""
    centroid = g.coords.mean(axis=0)
    best = None
    for phase in np.radians(np.arange(0, 360, 15)):
        coords, syms = g.coords.copy(), list(g.symbols)
        for site in sorted(set(arr)):
            k = arr.count(site)
            i = sites[site]
            nb = _neighbors(g, i)
            b = g.coords[i] - g.coords[nb].mean(axis=0)
            b /= np.linalg.norm(b)
            theta = np.radians(CONE[(site, k)])
            ref = _perp(b)
            for j in range(k):
                phi = phase + 2 * np.pi * j / k
                u = np.cos(theta) * b + np.sin(theta) * _rot(ref, b, phi)
                ow = g.coords[i] + DIST[site] * u
                w = _water(g.coords[i], ow, coords, centroid)
                coords = np.vstack([coords, w])
                syms += ["O", "H", "H"]
        # closest contact: every water atom against every other atom, except its own water and its H-bond pair
        n0 = len(g.symbols)
        dmin = 9.9
        for wi in range(n0, len(syms), 3):
            own = set(range(wi, wi + 3))
            for a in own:
                for b2 in range(len(syms)):
                    if b2 in own:
                        continue
                    d = np.linalg.norm(coords[a] - coords[b2])
                    if a == wi + 1 and d < 2.4:          # the donating H and its acceptor
                        continue
                    if a == wi and d < DIST["C"] + 0.05 and syms[b2] in "ONC" and b2 < n0:   # O(w)···acceptor
                        continue
                    dmin = min(dmin, d)
        if best is None or dmin > best[1]:
            best = (Geometry(syms, coords), dmin)
    return best


def build(rx: Reaction) -> dict:
    B = bare(rx)
    man = {}
    d = rx.d(f"{SUB}/starts")
    for sp, byn in ARRANGEMENTS.items():
        g, q, sites = B[sp]
        for n, arrs in byn.items():
            for arr in arrs:
                cg, dmin = build_one(g, sites, arr)
                name = f"{sp}_w{n}_{arr}"
                cg.write(d / f"{name}.xyz", comment=f"{name}: {sp} + {n} H2O on {arr}; closest other contact {dmin:.2f} Å")
                man[name] = {"species": sp, "n": n, "arr": arr, "atoms": len(cg.symbols), "closest_contact_A": round(dmin, 3)}
                print(f"{name:22s} {len(cg.symbols):3d} atoms  closest other contact {dmin:.2f} Å")
    _write_json(rx.root / SUB / "starts.json", man)
    return man


# ------------------------------------------------------------------ jobs
def _opt_input(rx: Reaction, sp: str, g: Geometry) -> str:
    kw = f'{rx.cfg["method"]} CPCM({rx.cfg["solvent"]}) ' + ("OptTS Freq" if sp == "ts" else "Opt Freq")
    geom = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 150 end" if sp == "ts" else "%geom MaxIter 200 end"
    return orca.make_input(kw, g, -1, 1, blocks=rx.solvent_block + "\n" + geom, nprocs=rx.cfg["nprocs"],
                           maxcore_mb=rx.cfg["maxcore_mb"])


def _bare_mode(rx: Reaction):
    h = orca.read_hess(Path(f"{rx.ts_stem()}.hess"))
    k = int(np.argmin(h["freqs"]))
    return read_xyz(Path(f"{rx.ts_stem()}.xyz")), h["modes"][:, k].reshape(-1, 3)


def _kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotation R with (P - mean) @ R ≈ (Q - mean)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    dd = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, dd]) @ Vt


def ts_mode_check(rx: Reaction, wd: Path, name: str) -> dict:
    g = read_xyz(wd / f"{name}.xyz")
    h = orca.read_hess(wd / f"{name}.hess")
    k = int(np.argmin(h["freqs"]))
    v = h["modes"][:, k].reshape(-1, 3)
    g0, v0 = _bare_mode(rx)
    n0 = len(g0.symbols)
    R = _kabsch(g.coords[:n0], g0.coords)              # cluster solute frame → bare frame
    vs = v[:n0] @ R
    ov = float(abs(np.sum(vs * v0)) / (np.linalg.norm(vs) * np.linalg.norm(v0)))
    share_w = float((v[n0:] ** 2).sum() / (v ** 2).sum())
    rates = {}
    for lab, (a, b) in {"C(N)–C": (4, 0), "C–O": (0, 1)}.items():
        u = g.coords[a] - g.coords[b]
        u /= np.linalg.norm(u)
        rates[lab] = float(np.dot(v[a] - v[b], u))
    n_imag = int(sum(1 for f in h["freqs"] if f < 0))
    ok = n_imag == 1 and ov >= MODE_MIN_OVERLAP and np.sign(rates["C(N)–C"]) != np.sign(rates["C–O"])
    return {"n_imaginary": n_imag, "freq_cm": float(h["freqs"][k]), "overlap_with_bare_mode": round(ov, 3),
            "water_share_of_mode": round(share_w, 3), "stretch_rates": rates, "accepted": bool(ok)}


def opt_one(rx: Reaction, name: str, meta: dict, watch=None) -> dict:
    sp = meta["species"]
    wd = rx.d(f"{SUB}/opt/{sp}_w{meta['n']}/{meta['arr']}")
    g = read_xyz(rx.root / SUB / "starts" / f"{name}.xyz")
    t0 = time.time()
    r = orca.run(wd, name, _opt_input(rx, sp, g), rx.log, note=f"explicit-water test (J. Picione 2026-11-20): {name}",
                 watch=watch)
    row = {"name": name, **meta, "job": name, "wall_s": r.wall_s,
           "ok": bool(r.ok and r.opt_converged is not False and not r.stopped)}
    if row["ok"] and sp != "ts" and r.imaginary:
        # a minimum that stopped on a saddle point (a turned water): push along the mode and reoptimize once
        h = orca.read_hess(wd / f"{name}.hess")
        v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
        gp = read_xyz(wd / f"{name}.xyz")
        gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
        prev, name2 = r.imaginary, f"{name}_min"
        r = orca.run(wd, name2, _opt_input(rx, sp, gp), rx.log,
                     note=f"rerun: {name} stopped on a saddle point ({[round(f, 1) for f in prev]} cm-1); "
                          "displaced 0.1 Å along that mode, same keywords", watch=watch)
        row.update(job=name2, first_imaginary_cm=prev, wall_s=row["wall_s"] + r.wall_s,
                   ok=bool(r.ok and r.opt_converged is not False and not r.stopped))
    if r.stopped:
        row["stopped"] = r.stopped
    if row["ok"]:
        row.update(E_eh=r.energy_eh, G_eh=r.gibbs_eh, imaginary_cm=r.imaginary)
        if sp == "ts":
            row["mode"] = ts_mode_check(rx, wd, row["job"])
            row["usable"] = row["mode"]["accepted"]
        else:
            row["usable"] = not r.imaginary
    else:
        row["usable"] = False
    row["dir"] = str(wd.relative_to(rx.root)).replace("\\", "/")
    row["elapsed_s"] = round(time.time() - t0)
    return row


def sp_one(rx: Reaction, lv: dict, key: str, row: dict, watch=None) -> float:
    g = read_xyz(rx.root / row["dir"] / f"{row['job']}.xyz")
    wd = rx.d(f"{SUB}/sp/{lv['label']}/{key}")
    inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', g, -1, 1, blocks=rx.solvent_block,
                          nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
    r = orca.run(wd, key, inp, rx.log, note=f"explicit-water test: {lv['label']} single point on {row['job']}",
                 watch=watch)
    if not r.ok or r.stopped:
        raise orca.OrcaError(f"{key} {lv['label']} failed; see {r.out}")
    return r.energy_eh, r.wall_s


# ------------------------------------------------------------------ estimate
def _bare_times(rx: Reaction) -> dict:
    """Job seconds × cores of the production jobs, from the run log table."""
    t = (rx.root / "run_log.md").read_text(encoding="utf-8")
    import re
    out = {}
    for path, key in [("01_stationary/ts/ts`", "ts_opt"),
                      ("01_ts/fragments/acetone_cyanohydrin_alkoxide/acetone_cyanohydrin_alkoxide`", "alk_opt"),
                      (f"01_ts/diagnostics/{DLPNO}/ts/ts`", "ts_dlpno"),
                      ("01_ts/diagnostics/M06-2X_ma-def2-TZVP/ts/ts`", "ts_tzvp"),
                      ("01_ts/diagnostics/B3LYP_ma-def2-TZVP/ts/ts`", "ts_b3lyp")]:
        m = [x for x in re.finditer(r"\| `" + re.escape(path) + r" \| `[^`]*` \| (\d+) s", t)]
        if m:
            out[key] = int(m[-1].group(1))
    return out


def estimate(rx: Reaction, test_row: dict | None) -> dict:
    """Job time of the whole test (sum over jobs, s at 4 cores). Optimizations scale from the test optimization by
    (atoms / atoms_test)^2.5, a TS optimization costs 2.5× a minimum of the same size (Hessians every 5 cycles);
    single points scale from the bare TS single points by (atoms / 12)^3. Without the test row, the bare alkoxide
    optimization (1 core) stands in for it (a pre-estimate)."""
    man = _read_json(rx.root / SUB / "starts.json")
    bt = _bare_times(rx)
    if test_row:
        t_ref, n_ref, src = test_row["wall_s"], man["_".join([TEST[0], f"w{TEST[1]}", TEST[2]])]["atoms"], "test optimization"
    else:
        t_ref, n_ref, src = bt["alk_opt"] / 4 * 1.5, 12, "bare alkoxide Opt Freq (1 core, /4 cores, ×1.5 for loose waters)"
    opt = 0.0
    for name, m in man.items():
        f = 2.5 if m["species"] == "ts" else 1.0
        opt += t_ref * (m["atoms"] / n_ref) ** 2.5 * f
    sps = 0.0
    for sp in ARRANGEMENTS:
        for n in (1, 2, 3):
            na = {"cyanide": 2, "ts": 12, "alkoxide": 12}[sp] + 3 * n
            for k in ("ts_tzvp", "ts_b3lyp", "ts_dlpno"):
                sps += bt.get(k, 600) * (na / 12) ** 3
    return {"source": src, "ref_s": round(t_ref), "opt_h": round(opt / 3600, 1), "sp_h": round(sps / 3600, 1),
            "total_h": round((opt + sps) / 3600, 1), "limit_h": LIMIT_H, "jobs_opt": len(man), "jobs_sp": 27}


# ------------------------------------------------------------------ wait behind the path single points
def path_sp_done(root: Path) -> bool:
    import psutil
    files = all((root / "reactions" / r / f"02_irc/energies_{DLPNO}.csv").exists() for r in PATH_RX)
    lanes = any("run_singlepoints.sh" in " ".join(p.info["cmdline"] or [])
                for p in psutil.process_iter(["cmdline"]) if p.info["cmdline"])
    return files and not lanes


# ------------------------------------------------------------------ main run
def run(rx: Reaction, now: bool, not_before: str | None = None, hard_stop: str | None = None) -> None:
    root = rx.root.parents[1]
    if not_before:
        t0 = time.mktime(time.strptime(not_before, "%Y-%m-%dT%H:%M"))
        print(f"{time.strftime('%Y-%m-%d %H:%M')} waiting until {not_before} (J. Picione 2026-11-26); "
              f"limit {LIMIT_H:g} h, hard stop {hard_stop}", flush=True)
        while time.time() < t0:
            time.sleep(30)
        print(f"{time.strftime('%Y-%m-%d %H:%M')} starting", flush=True)
    if not now:
        print(f"{time.strftime('%Y-%m-%d %H:%M')} waiting behind the path single points", flush=True)
        while not path_sp_done(root):
            time.sleep(600)
        print(f"{time.strftime('%Y-%m-%d %H:%M')} path single points done; starting", flush=True)
    man = _read_json(rx.root / SUB / "starts.json")
    out_json = rx.root / SUB / "explicit_water.json"
    res = _read_json(out_json) if out_json.exists() else {"opt": {}}
    test = "_".join([TEST[0], f"w{TEST[1]}", TEST[2]])
    if test not in res["opt"]:
        res["opt"][test] = opt_one(rx, test, man[test])
        _write_json(out_json, res)
    est = estimate(rx, res["opt"][test])
    res["estimate"] = est
    _write_json(out_json, res)
    print(f"{time.strftime('%Y-%m-%d %H:%M')} test optimization {test}: {res['opt'][test]['wall_s']:.0f} s; "
          f"estimate {est}", flush=True)
    stop = rx.root / SUB / "STOP_ESTIMATE.md"
    if est["total_h"] > LIMIT_H and not (rx.root / SUB / "GO").exists():
        stop.write_text(f"# Explicit-water test stopped before the rest (J. Picione's limit {LIMIT_H:g} h)\n\n"
                        f"Test optimization {test}: {res['opt'][test]['wall_s']:.0f} s. Estimated job time of the whole "
                        f"test: {est['total_h']} h ({est['opt_h']} h optimizations, {est['sp_h']} h single points).\n"
                        "Waiting for his decision; write 04_explicit_water/GO to run it anyway.\n", encoding="utf-8")
        print("STOP: estimate over the limit", flush=True)
        return

    if stop.exists():                                   # the 2026-11-25 stop note moves aside (kept)
        stop.rename(stop.with_name("STOP_ESTIMATE_2026-11-25.md"))
    hard = time.mktime(time.strptime(hard_stop, "%Y-%m-%dT%H:%M")) if hard_stop else None
    past = lambda: hard is not None and time.time() >= hard
    watch = (lambda: f"hard stop {hard_stop} (J. Picione 2026-11-26)" if past() else None) if hard else None
    res.setdefault("kept", {})
    res.setdefault("sp", {})
    res.setdefault("not_kept", {})
    res["order"] = "by n (J. Picione 2026-11-26)"

    def opt_job(name):
        return None if past() else opt_one(rx, name, man[name], watch)

    def sp_job(j):
        if past():
            return None
        try:
            return sp_one(rx, j[0], j[1], res["opt"][res["kept"][j[1]]], watch)
        except orca.OrcaError:
            if past():
                return None
            raise

    for n in (1, 2, 3):
        todo = [k for k, m in man.items() if m["n"] == n and k not in res["opt"]]
        print(f"{time.strftime('%Y-%m-%d %H:%M')} n = {n}: {len(todo)} optimizations", flush=True)
        with ThreadPoolExecutor(max_workers=2) as ex:      # two 4-core jobs at a time, as the path lanes
            for row in ex.map(opt_job, todo):
                if row is None:
                    continue
                if row.get("stopped"):
                    res.setdefault("stopped_jobs", []).append(row["name"])
                    progress(rx, res, man)
                    continue
                res["opt"][row["name"]] = row
                _write_json(out_json, res)
                progress(rx, res, man)
                print(f"{time.strftime('%H:%M')} {row['name']}: {row['wall_s']:.0f} s, usable {row['usable']}", flush=True)
        if past():
            break
        # lowest G per species among the usable structures of this n
        for sp in ARRANGEMENTS:
            rows = [r for r in res["opt"].values() if r["species"] == sp and r["n"] == n and r.get("usable")]
            if rows:
                res["kept"][f"{sp}_w{n}"] = min(rows, key=lambda r: r["G_eh"])["name"]
            else:
                res["not_kept"][f"{sp}_w{n}"] = "no usable structure"
                print(f"{time.strftime('%H:%M')} no usable {sp}·(H2O){n}; n = {n} has no single points", flush=True)
        _write_json(out_json, res)
        keys = [f"{sp}_w{n}" for sp in ARRANGEMENTS]
        if any(k not in res["kept"] for k in keys):
            continue
        jobs = [(lv, key) for lv in rx.cfg["diagnostic_levels"] for key in keys
                if key not in res["sp"].get(lv["label"], {})]
        print(f"{time.strftime('%Y-%m-%d %H:%M')} n = {n}: {len(jobs)} single points", flush=True)
        with ThreadPoolExecutor(max_workers=2) as ex:
            for (lv, key), out in zip(jobs, ex.map(sp_job, jobs)):
                if out is None:
                    continue
                e, w = out
                res["sp"].setdefault(lv["label"], {})[key] = {"E_eh": e, "wall_s": w}
                _write_json(out_json, res)
                progress(rx, res, man)
        if past():
            break
        analyse(rx, res)
        _write_json(out_json, res)
        write_md(rx, res)
        progress(rx, res, man)
        print(f"{time.strftime('%Y-%m-%d %H:%M')} n = {n} complete", flush=True)
    progress(rx, res, man, final=True)
    print(f"{time.strftime('%Y-%m-%d %H:%M')} " + ("stopped at the hard stop" if past() else "done")
          + "; nothing else starts", flush=True)


def _model_s(rx: Reaction, res: dict, man: dict) -> tuple:
    """Model job seconds (as estimate()) per optimization and per single point (level, species_wn)."""
    test = "_".join([TEST[0], f"w{TEST[1]}", TEST[2]])
    t_ref, n_ref = res["opt"][test]["wall_s"], man[test]["atoms"]
    bt = _bare_times(rx)
    mo = {k: t_ref * (m["atoms"] / n_ref) ** 2.5 * (2.5 if m["species"] == "ts" else 1.0) for k, m in man.items()}
    lvk = {"M06-2X_ma-def2-TZVP": "ts_tzvp", "B3LYP_ma-def2-TZVP": "ts_b3lyp", DLPNO: "ts_dlpno"}
    ms = {}
    for lv in rx.cfg["diagnostic_levels"]:
        for sp in ARRANGEMENTS:
            for n in (1, 2, 3):
                na = {"cyanide": 2, "ts": 12, "alkoxide": 12}[sp] + 3 * n
                ms[(lv["label"], f"{sp}_w{n}")] = bt.get(lvk.get(lv["label"], ""), 600) * (na / 12) ** 3
    return mo, ms


def progress(rx: Reaction, res: dict, man: dict, final: bool = False) -> None:
    """04_explicit_water/progress.md/.json: what is done, what is left, and the remaining time: the model scaled by
    the measured / model ratio of the jobs already done (optimizations and single points separately)."""
    test = "_".join([TEST[0], f"w{TEST[1]}", TEST[2]])
    mo, ms = _model_s(rx, res, man)
    done_o = [k for k in res["opt"] if k != test]
    left_o = [k for k in man if k not in res["opt"]]
    ro = (sum(res["opt"][k]["wall_s"] for k in done_o) / sum(mo[k] for k in done_o)) if done_o else 1.0
    done_s = [(lv, k) for lv, d in res.get("sp", {}).items() for k in d]
    all_s = [(lv["label"], f"{sp}_w{n}") for n in (1, 2, 3) for lv in rx.cfg["diagnostic_levels"] for sp in ARRANGEMENTS]
    left_s = [j for j in all_s if j not in done_s]
    rs_ = (sum(res["sp"][lv][k]["wall_s"] for lv, k in done_s) / sum(ms[j] for j in done_s)) if done_s else 1.0
    job_h = (sum(mo[k] for k in left_o) * ro + sum(ms[j] for j in left_s) * rs_) / 3600
    per_n = {n: {"opt_done": sum(1 for k in res["opt"] if man[k]["n"] == n),
                 "opt_all": sum(1 for m in man.values() if m["n"] == n),
                 "sp_done": sum(1 for lv, k in done_s if k.endswith(f"_w{n}")),
                 "sp_all": 3 * len(rx.cfg["diagnostic_levels"])} for n in (1, 2, 3)}
    p = {"time": time.strftime("%Y-%m-%d %H:%M"), "final": final, "per_n": per_n,
         "left_opt": left_o, "left_sp": [f"{lv} {k}" for lv, k in left_s],
         "measured_over_model": {"opt": round(ro, 2), "sp": round(rs_, 2)},
         "left_job_h": round(job_h, 1), "left_clock_h_two_lanes": round(job_h / 2, 1),
         "stopped_jobs": res.get("stopped_jobs", []), "not_kept": res.get("not_kept", {})}
    _write_json(rx.root / SUB / "progress.json", p)
    L = [f"# Explicit-water test: progress at {p['time']}" + (" (final)" if final else ""), "",
         "| n | optimizations done | single points done |", "|---|---|---|"]
    L += [f"| {n} | {v['opt_done']} of {v['opt_all']} | {v['sp_done']} of {v['sp_all']} |" for n, v in per_n.items()]
    L += ["", f"Left: {len(left_o)} optimizations, {len(left_s)} single points. New estimate of what is left: "
          f"{p['left_job_h']} h job time, about {p['left_clock_h_two_lanes']} h on the clock (two jobs at a time); "
          f"jobs so far took {ro:.2f}× (optimizations) and {rs_:.2f}× (single points) the model time."]
    if p["stopped_jobs"]:
        L += ["", "Ended at the hard stop (not results): " + ", ".join(p["stopped_jobs"])]
    if p["not_kept"]:
        L += ["", "No usable structure: " + ", ".join(p["not_kept"])]
    (rx.root / SUB / "progress.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# ------------------------------------------------------------------ analysis
def analyse(rx: Reaction, res: dict) -> None:
    frag = _read_json(rx.root / "01_ts/fragments.json")
    diag = _read_json(rx.root / "01_ts/diagnostics.json")
    reac, _ = reference_names(rx, frag)
    ac = next(n for n in reac if n != "cyanide")
    wl = rx.level().replace("/", "_") + " (working)"
    levels = [wl] + [lv["label"] for lv in rx.cfg["diagnostic_levels"]]
    rows = []
    for lv in levels:
        base = next(r for r in diag["rows"] if r["level"] == lv)
        row = {"level": lv, "n0": {"dG_TS_1M_kJ": base["dG_TS_vs_separated_1M_kJ"], "dG_rxn_1M_kJ": base["dG_reaction_1M_kJ"]}}
        Eac = frag[ac]["E_eh"] if lv == wl else diag["energies_eh"][lv][ac]
        Gac = Eac + frag[ac]["G_eh"] - frag[ac]["E_eh"]
        for n in (1, 2, 3):
            keys = [f"{sp}_w{n}" for sp in ARRANGEMENTS]
            if not all(k in res.get("kept", {}) and (lv == wl or k in res["sp"].get(lv, {})) for k in keys):
                row[f"n{n}"] = None                     # not complete (order by n; J. Picione 2026-11-26)
                continue
            G = {}
            for sp in ARRANGEMENTS:
                o = res["opt"][res["kept"][f"{sp}_w{n}"]]
                E = o["E_eh"] if lv == wl else res["sp"][lv][f"{sp}_w{n}"]["E_eh"]
                G[sp] = E + o["G_eh"] - o["E_eh"]
            row[f"n{n}"] = {"dG_TS_1M_kJ": kj(G["ts"] - G["cyanide"] - Gac) - STD_STATE_KJ,
                            "dG_rxn_1M_kJ": kj(G["alkoxide"] - G["cyanide"] - Gac) - STD_STATE_KJ}
        rows.append(row)
    res["rows"] = rows
    res["measured"] = {"dG_TS_1M_kJ": 69.6, "dG_rxn_1M_kJ": 10.5}
    done_n = [n for n in (1, 2, 3) if rows and rows[-1].get(f"n{n}")]
    res["contacts"] = {k: contacts(rx.root / res["opt"][v]["dir"] / f"{res['opt'][v]['job']}.xyz", k.split('_w')[0])
                       for k, v in res["kept"].items() if int(k.split("_w")[1]) in done_n}


def contacts(xyz: Path, sp: str) -> list[dict]:
    """For each water: the solute or water atom its H points at (H···X ≤ 2.6 Å), H···X and O(w)···X."""
    g = read_xyz(xyz)
    s = g.symbols
    n0 = {"cyanide": 2, "ts": 12, "alkoxide": 12}[sp]
    # acceptors: the solute O and N (and the C end of free cyanide), and the O of the other waters
    lab = {i: {"O": "O (alkoxide/carbonyl)", "N": "N (nitrile)", "C": "C (cyanide)"}[s[i]] for i in range(n0)
           if s[i] in "ON" or (sp == "cyanide" and s[i] == "C")}
    out = []
    for w, o in enumerate(range(n0, len(s), 3)):
        for h in (o + 1, o + 2):
            cands = list(lab) + [i for i in range(n0, len(s), 3) if i != o]
            j = min(cands, key=lambda i: np.linalg.norm(g.coords[h] - g.coords[i]))
            d = float(np.linalg.norm(g.coords[h] - g.coords[j]))
            if d <= 2.6:
                out.append({"water": w + 1, "to": lab.get(j, f"O of water {(j - n0) // 3 + 1}"),
                            "H_X_A": round(d, 3), "O_X_A": round(float(np.linalg.norm(g.coords[o] - g.coords[j])), 3)})
    return out


def _lv(label: str) -> str:
    return "M06-2X/ma-def2-SVP (working level)" if label.endswith("(working)") else label.replace("_", "/")


def write_md(rx: Reaction, res: dict) -> None:
    f = lambda x: f"{x:+.1f} ({x / 4.184:+.1f})" if x is not None else "–"
    g = lambda x: f"{x:.1f} ({x / 4.184:.1f})" if x is not None else "–"
    L = ["Sensitivity test only: every other number in this summary stays on SMD water without explicit waters. "
         "Cluster-continuum model: n = 1, 2, 3 explicit waters, the same n on both sides (CN⁻·(H₂O)n + acetone → "
         "TS·(H₂O)n → alkoxide·(H₂O)n), all in SMD water; lowest-G arrangement for each species and n (working level); "
         "single points on those at the other levels with working-level G − E; 298.15 K, 1 M. kJ/mol (kcal/mol).", ""]
    rt = res.get("retry")
    old = {r["level"]: r for r in rt["before"]["rows"]} if rt else {}
    old_ref = rt["before"]["kept"]["cyanide_w3"].split("_")[-1] if rt else ""
    changed = bool(rt) and rt.get("kept_cyanide_w3") != rt["before"]["kept"]["cyanide_w3"]
    for key, title, meas in [("dG_TS_1M_kJ", "ΔG‡ from separated reactants", 69.6), ("dG_rxn_1M_kJ", "ΔrG", 10.5)]:
        fm = g if key == "dG_TS_1M_kJ" else f
        hd = "| Level | n = 0 (production) | n = 1 | n = 2 | n = 3 |" + (f" n = 3, earlier reference ({old_ref}) |"
                                                                          if changed else "")
        L += [f"**{title}:**", "", hd, "|---|---|---|---|---|" + ("---|" if changed else "")]
        for r in res["rows"]:
            L.append(f"| {_lv(r['level'])} | " + " | ".join(fm(r[k][key] if r.get(k) else None)
                                                          for k in ("n0", "n1", "n2", "n3")) + " |"
                     + (f" {fm(old[r['level']]['n3'][key])} |" if changed else ""))
        L += [f"| Measured | {fm(meas)} | | | |" + (" |" if changed else ""), ""]
    if rt:
        o_new = {k: res["opt"][k] for k in RETRY}
        L += ["**Retry of two cyanide·(H₂O)3 arrangements** (" + rt["rules"] + ")", ""]
        for k, o in o_new.items():
            b = rt["before"]["opt"][k]
            st = ("a true minimum, usable" if o.get("usable") else
                  f"still not a minimum ({', '.join(f'{-x:.0f}i' for x in o.get('imaginary_cm') or [])} cm⁻¹), not used")
            L.append(f"- {k}: was {abs(b['imaginary_cm'][0]):.0f}i cm⁻¹; after the retry {st}.")
        if changed:
            L += ["", f"The n = 3 reactant is now cyanide·(H₂O)3 {rt['kept_cyanide_w3'].split('_')[-1]}, the lowest "
                  f"working-level free energy of all usable arrangements; it replaces {old_ref}, which was kept "
                  "only because the two arrangements with water on both ends of cyanide still had an imaginary "
                  "frequency. The last column keeps the earlier n = 3 values."]
        elif rt.get("neither_minimum"):
            L += ["", "Neither became a minimum: the n = 3 values are lower bounds (the two unusable arrangements lie "
                  "5.7–7.2 kJ/mol below the kept one in working-level free energy)."]
        L.append("")
    L += ["**Lowest structures and their hydrogen bonds** (H···X and O(water)···X, Å):", "",
          "| Structure | Kept arrangement | Water → acceptor: H···X / O···X |", "|---|---|---|"]
    for k, v in res["kept"].items():
        if k not in res["contacts"]:
            continue
        o = res["opt"][v]
        c = "; ".join(f"{x['water']} → {x['to']}: {x['H_X_A']:.2f} / {x['O_X_A']:.2f}" for x in res["contacts"][k])
        extra = ""
        if o["species"] == "ts":
            m = o["mode"]
            extra = f" ({m['freq_cm']:.0f}i cm⁻¹, overlap with the bare TS mode {m['overlap_with_bare_mode']:.2f})"
        L.append(f"| {k.replace('_w', '·(H₂O)')} | {o['arr']}{extra} | {c} |")
    if res.get("lowest_freqs"):
        L += ["", "**Lowest frequency of every structure** (read by script from each final .hess; negative = "
              "imaginary), cm⁻¹:", "", "| Structure | Job | Lowest frequency | Imaginary |", "|---|---|---|---|"]
        L += [f"| {k} | {v['job']} | {v['lowest_cm']:.1f} | {v['n_imaginary']} |"
              for k, v in sorted(res["lowest_freqs"].items())]
    L += ["", f"Arrangements tried: {len(res['opt'])}; usable: {sum(1 for r in res['opt'].values() if r.get('usable'))}. "
          f"Script `pipeline/explicit_water.py`; jobs in `{SUB}/`."]
    fd = rx.root / SUB / "findings.md"
    if fd.exists():
        L += ["", fd.read_text(encoding="utf-8").strip()]
    (rx.root / SUB / "explicit_water.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# ------------------------------------------------------------------ one retry (J. Picione, 2026-11-27)
RETRY = ["cyanide_w3_NNC", "cyanide_w3_NCC"]
RETRY_RULES = ("J. Picione 2026-11-27, fixed before the run: one attempt each, started 0.1 Å along the remaining "
               "imaginary mode, same level and keywords (Opt Freq); no second attempt, no change of settings. A true "
               "minimum joins the usable set and the n = 3 reactant is again the lowest working-level G of all usable "
               "arrangements; single points on any new kept structure and the n = 3 rows rebuilt. If neither is a "
               "minimum, CCC stays and the n = 3 values are lower bounds; no energy from a structure with an imaginary "
               "mode is used.")


def lowest_freqs(rx: Reaction, res: dict) -> dict:
    """Lowest frequency of every optimized structure, read from its .hess (negative = imaginary)."""
    out = {}
    for k, o in res["opt"].items():
        h = orca.read_hess(rx.root / o["dir"] / f"{o['job']}.hess")
        f = sorted(float(x) for x in h["freqs"] if abs(x) > 1e-6)   # translations/rotations are written as 0
        out[k] = {"job": o["job"], "lowest_cm": round(f[0], 1), "n_imaginary": sum(1 for x in f if x < 0)}
    return out


def retry(rx: Reaction) -> None:
    import copy
    out_json = rx.root / SUB / "explicit_water.json"
    res = _read_json(out_json)
    if "retry" in res:
        raise SystemExit("retry already run (one attempt only, J. Picione 2026-11-27)")
    res["retry"] = {"rules": RETRY_RULES, "started": time.strftime("%Y-%m-%d %H:%M"),
                    "before": {"opt": {k: copy.deepcopy(res["opt"][k]) for k in RETRY},
                               "kept": copy.deepcopy(res["kept"]), "rows": copy.deepcopy(res["rows"]),
                               "sp_cyanide_w3": {lv: copy.deepcopy(d["cyanide_w3"]) for lv, d in res["sp"].items()}}}
    _write_json(out_json, res)
    print(f"{time.strftime('%Y-%m-%d %H:%M')} retry of {', '.join(RETRY)}: {RETRY_RULES}", flush=True)

    def one(name):
        o = res["opt"][name]
        wd = rx.root / o["dir"]
        h = orca.read_hess(wd / f"{o['job']}.hess")
        k = int(np.argmin(h["freqs"]))
        v = h["modes"][:, k].reshape(-1, 3)
        g = read_xyz(wd / f"{o['job']}.xyz")
        g = Geometry(g.symbols, g.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
        name2 = f"{name}_retry"
        r = orca.run(wd, name2, _opt_input(rx, o["species"], g), rx.log,
                     note=f"retry (J. Picione 2026-11-27): {o['job']} kept {h['freqs'][k]:.1f} cm-1; displaced 0.1 Å "
                          "along that mode, same keywords; one attempt")
        ok = bool(r.ok and r.opt_converged is not False and not r.stopped)
        row = dict(o, job=name2, wall_s=o["wall_s"] + r.wall_s, ok=ok, retry_wall_s=r.wall_s)
        if ok:
            row.update(E_eh=r.energy_eh, G_eh=r.gibbs_eh, imaginary_cm=r.imaginary, usable=not r.imaginary)
        else:
            row.update(usable=False)
            for x in ("E_eh", "G_eh", "imaginary_cm"):
                row.pop(x, None)
        return row

    with ThreadPoolExecutor(max_workers=2) as ex:
        for row in ex.map(one, RETRY):
            res["opt"][row["name"]] = row
            _write_json(out_json, res)
            print(f"{time.strftime('%H:%M')} {row['job']}: {row['retry_wall_s']:.0f} s, ok {row['ok']}, imaginary "
                  f"{row.get('imaginary_cm')}, usable {row['usable']}", flush=True)
    rows = [r for r in res["opt"].values() if r["species"] == "cyanide" and r["n"] == 3 and r.get("usable")]
    best = min(rows, key=lambda r: r["G_eh"])["name"]
    res["retry"]["kept_cyanide_w3"] = best
    if best != res["kept"]["cyanide_w3"]:
        res["kept"]["cyanide_w3"] = best
        o = res["opt"][best]
        for lv in rx.cfg["diagnostic_levels"]:      # new folder per arrangement: the CCC single points stay as they are
            g = read_xyz(rx.root / o["dir"] / f"{o['job']}.xyz")
            key = f"cyanide_w3_{o['arr']}"
            wd = rx.d(f"{SUB}/sp/{lv['label']}/{key}")
            inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', g, -1, 1, blocks=rx.solvent_block,
                                  nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
            r = orca.run(wd, key, inp, rx.log, note=f"retry (J. Picione 2026-11-27): {lv['label']} single point on {o['job']}")
            if not r.ok or r.stopped:
                raise orca.OrcaError(f"{key} {lv['label']} failed; see {r.out}")
            res["sp"][lv["label"]]["cyanide_w3"] = {"E_eh": r.energy_eh, "wall_s": r.wall_s, "job": f"sp/{lv['label']}/{key}"}
            _write_json(out_json, res)
            print(f"{time.strftime('%H:%M')} {lv['label']} {key}: {r.wall_s:.0f} s", flush=True)
    res["retry"]["neither_minimum"] = not any(res["opt"][k]["usable"] for k in RETRY)   # rule 3: lower bounds
    analyse(rx, res)
    res["lowest_freqs"] = lowest_freqs(rx, res)
    res["retry"]["ended"] = time.strftime("%Y-%m-%d %H:%M")
    _write_json(out_json, res)
    print(f"{time.strftime('%Y-%m-%d %H:%M')} retry done; kept cyanide_w3 = {best}; nothing else starts", flush=True)


def main() -> None:
    rx = Reaction.load(Path(sys.argv[1]))
    cmd = sys.argv[2]
    if cmd == "build":
        build(rx)
    elif cmd == "estimate":
        print(json.dumps(estimate(rx, None), indent=1))
    elif cmd == "run":
        a = sys.argv[3:]
        opt = lambda f: a[a.index(f) + 1] if f in a else None
        run(rx, "--now" in a, opt("--not-before"), opt("--hard-stop"))
    elif cmd == "md":
        write_md(rx, _read_json(rx.root / SUB / "explicit_water.json"))
    elif cmd == "retry":
        retry(rx)
    elif cmd == "freqs":
        print(json.dumps(lowest_freqs(rx, _read_json(rx.root / SUB / "explicit_water.json")), indent=1))


if __name__ == "__main__":
    main()
