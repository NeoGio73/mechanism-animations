"""Secondary hydride shift: the C2–C3 rotation TSs from the anti minimum at the working level (J. Picione, 2026-11-22).

Nothing in the production path folders changes; every job is in a new subfolder of 01_stationary/. Started after
the tertiary shift's last path single point, in the ordinary core queue (J. Picione 2026-11-22: "a 1-1.5 h delay to E2
is fine"); one job at a time. Order: the two torsion scans of step 2 first (whether anti and eclip interconvert at all
decides how Section 3.6 reports the barrier), then step 1, then any anti-eclip TS search.

1. ts_anti_rotation: OptTS Freq at the working level (method of reaction.toml, SMD dichloromethane, Calc_Hess true
   Recalc_Hess 5 MaxIter 100, the production TS settings), started from the gas-phase M06-2X anti TS
   (reactions/hshift_3_methyl_2_butyl/01_check/rotation/M06-2X/anti/ts.xyz). Check: one imaginary frequency. Then the
   TS is displaced both ways along that mode (largest atom displacement 0.10 Å) and each side optimized (Opt Freq):
   one side must reach secondary_anti, the other the tertiary cation.
2. ts_anti_eclip: the TS that interconverts anti and eclip directly (no hydride shift). The two minima differ in the
   C1–C2–C3–H torsion (anti −1.2°, eclip −161.9°), so the C3–H bond has to pass the empty p orbital of C2 (±90°)
   whichever way the C2–C3 bond turns. Relaxed torsion scans from secondary_anti, both directions (to −161.9° and to
   +161.9°, the mirror image of eclip), 18 points each. An arm on which the H stays on C3 and the energy has a maximum
   gives a TS guess: OptTS Freq, one imaginary frequency, displace both ways and optimize: one side must reach
   secondary_anti, the other eclip (either mirror image). An arm on which the H moves to C2 has no direct
   interconversion: that is reported, with the scan profile.
3. Single points on every TS found at M06-2X/ma-def2-TZVP, B3LYP/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP (same
   inputs as the path single points).
4. ΔG‡ from secondary_anti at four levels: E(level) difference plus the working-level thermal correction difference
   (ORCA quasi-RRHO G − E, 298.15 K; one molecule, no standard-state term), next to the production TS
   (01_ts/ts.json, the eclip-side rotation TS; DLPNO-CCSD(T) 14.0 kJ/mol).

A minimum matches a reference when the core atoms (C1, C2, its H, C3, the two C3 methyl carbons, the shifting H)
superpose within 0.10 Å RMSD (mirror image and the two C3 methyls swapped allowed) and the working-level energies
agree within 1.0 kJ/mol. A failed or unconverged job stops its branch and is reported; nothing is rerun with other
settings.

Usage: python -X utf8 reactions/hshift_3_methyl_2_butyl_m062x/00_inputs/rotation_ts_m062x.py run [after_e2|after_tertiary]
       ... check     (no ORCA job: inputs, references and the displacement on the gas-phase TS mode)
       ... report    (rewrite 01_stationary/rotation_ts.json/.md from finished jobs)
Start gate: after_e2 (default) waits until the E2 and tertiary path single points are written
(02_irc/energies_DLPNO-CCSD(T)_ma-def2-TZVP.json); after_tertiary waits for the tertiary lane only.
"""
import itertools
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
ROOT = RX.parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import orca                                   # noqa: E402
from pipeline.geometry import Geometry, read_xyz            # noqa: E402
from pipeline.steps import Reaction                         # noqa: E402
from pipeline.units import kj                               # noqa: E402

C1, C2, HC2, C3, H = 0, 1, 2, 3, 9                         # 0-based, as in reaction.toml [key_atoms]
ME = (4, 5)                                                 # the two methyl carbons on C3
CORE = [C1, C2, HC2, C3, ME[0], ME[1], H]
GAS_TS = ROOT / "reactions/hshift_3_methyl_2_butyl/01_check/rotation/M06-2X/anti/ts.xyz"
GAS_HESS = GAS_TS.with_suffix(".hess")
ST = RX / "01_stationary"
OUT_JSON, OUT_MD = ST / "rotation_ts.json", ST / "rotation_ts.md"
TS_BLOCK = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end"
ECLIP_DIH = -161.9
N_SCAN = 18
STEP_A = 0.10
RMSD_OK, DE_OK = 0.10, 1.0
GATES = {"after_e2": ["e2_ethoxide_2_bromopropane_m062x", "hshift_2_3_dimethyl_2_butyl_m062x"],
         "after_tertiary": ["hshift_2_3_dimethyl_2_butyl_m062x"]}


# ------------------------------------------------------------------ geometry
def d(g, i, j):
    return float(np.linalg.norm(g.coords[i] - g.coords[j]))


def tors(g, i, j, k, l):
    p = g.coords
    b0, b1, b2 = p[i] - p[j], p[k] - p[j], p[l] - p[k]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def describe(g):
    return {"C3_H_A": round(d(g, C3, H), 3), "C2_H_A": round(d(g, C2, H), 3),
            "C1_C2_C3_H_deg": round(tors(g, C1, C2, C3, H), 1),
            "H_on": "C3 (secondary)" if d(g, C3, H) < d(g, C2, H) else "C2 (tertiary)",
            # methyl positions (added 2026-11-23: the H-only check missed the methyl shift on scan_pos)
            **{f"C{c}_Me{m}_A": round(d(g, a, m), 3) for m in ME for c, a in ((2, C2), (3, C3))},
            "Me_moved_to_C2": [m for m in ME if d(g, C2, m) < d(g, C3, m)]}


def core_rmsd(a, b):
    """Core-atom RMSD after the best superposition; mirror image and the two C3 methyls swapped allowed."""
    best = 1e9
    for perm in (CORE, [C1, C2, HC2, C3, ME[1], ME[0], H]):
        x, y = a.coords[CORE], b.coords[perm]
        x, y = x - x.mean(0), y - y.mean(0)
        u, _, vt = np.linalg.svd(x.T @ y)
        r = u @ vt                                           # no determinant correction: reflections allowed
        best = min(best, float(np.sqrt(((x @ r - y) ** 2).sum(1).mean())))
    return best


# C2/C3-swapped labelling (degenerate 1,2-methyl shift of methyl 4): the cation centre is then C3, carrying H9 and
# methyl 5; C2 carries C1, H2 and methyl 4. Candidate atoms in the order of CORE of the unswapped cation:
SWAP = [5, C3, H, C2, C1, 4, HC2]


def swapped_rmsd(g, ref):
    """Core RMSD of candidate g, relabelled as the C2/C3-swapped cation, to ref; mirror and methyl swap allowed."""
    best = 1e9
    for perm in (SWAP, SWAP[:4] + [SWAP[5], SWAP[4]] + SWAP[6:]):
        x, y = ref.coords[CORE], g.coords[perm]
        x, y = x - x.mean(0), y - y.mean(0)
        u, _, vt = np.linalg.svd(x.T @ y)
        best = min(best, float(np.sqrt(((x @ (u @ vt) - y) ** 2).sum(1).mean())))
    return best


def displaced(ts_dir: Path, name: str, sign: int) -> Geometry:
    h = orca.read_hess(ts_dir / f"{name}.hess")
    k = int(np.argmin(h["freqs"]))
    v = h["modes"][:, k].reshape(-1, 3)
    g = read_xyz(ts_dir / f"{name}.xyz")
    return Geometry(g.symbols, g.coords + sign * STEP_A * v / np.linalg.norm(v, axis=1).max())


# ------------------------------------------------------------------ references
def references(rx):
    ep = json.loads((RX / "02_irc/endpoints.json").read_text(encoding="utf-8"))
    anti = json.loads((ST / "secondary_anti.json").read_text(encoding="utf-8"))
    refs = {"secondary_anti": {"xyz": RX / f"{anti['job']}.xyz", "E_eh": anti["E_eh"], "G_eh": anti["G_eh"],
                               "sp": anti["sp"]},
            "secondary_eclip_irc": {"xyz": RX / ep["reactant_complex"]["xyz"], "E_eh": ep["reactant_complex"]["E_eh"]},
            "tertiary_irc": {"xyz": RX / ep["product_complex"]["xyz"], "E_eh": ep["product_complex"]["E_eh"]}}
    for v in refs.values():
        v["geom"] = read_xyz(v["xyz"])
    refs["secondary_anti_swapped"] = {**refs["secondary_anti"], "swapped": True}
    return refs


def classify(g, e_eh, refs):
    out = {"geometry": describe(g), "matches": {}}
    for n, r in refs.items():
        rm = swapped_rmsd(g, r["geom"]) if r.get("swapped") else core_rmsd(g, r["geom"])
        out["matches"][n] = {"core_rmsd_A": round(rm, 3),
                             "dE_kJ": round(kj(e_eh - r["E_eh"]), 2)}
    moved = 4 in out["geometry"]["Me_moved_to_C2"]
    # anti and its C2/C3-swapped copy are only ~0.07 Å apart (the anti methyl is nearly bridging): the side is decided
    # by the carbon the methyl is bonded to, not by RMSD alone (2026-11-23)
    ok = lambda n: (n != "secondary_anti_swapped" or moved) and (n != "secondary_anti" or not moved)
    hit = [n for n, m in out["matches"].items()
           if m["core_rmsd_A"] <= RMSD_OK and abs(m["dE_kJ"]) <= DE_OK and ok(n)]
    out["reached"] = hit[0] if hit else None
    return out


# ------------------------------------------------------------------ jobs
def job(rx, folder, name, kw, g, note, blocks=""):
    wd = rx.d(f"01_stationary/{folder}")
    r = orca.run(wd, name, rx.inp(kw, g, blocks=blocks), rx.log, note=note)
    return wd, r


def sp_levels(rx, folder, g):
    res = {}
    for lv in rx.cfg["diagnostic_levels"]:
        wd = rx.d(f"01_stationary/{folder}/sp/{lv['label']}")
        inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', g, rx.cfg["charge"], rx.cfg["mult"],
                              blocks=rx.solvent_block, nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
        s = orca.run(wd, "sp", inp, rx.log, note=f"{folder}, single point {lv['label']}")
        if not s.ok:
            raise orca.OrcaError(f"single point {lv['label']} on {folder} failed; see {s.out}")
        res[lv["label"]] = s.energy_eh
    return res


def ts_branch(rx, folder, g0, refs, want, note):
    """OptTS Freq from g0, then displace both ways and optimize. want = the two minima the sides must reach."""
    rec = {"folder": f"01_stationary/{folder}", "start": note, "want": want}
    wd, r = job(rx, f"{folder}/ts", "ts", "OptTS Freq", g0, f"{folder}: OptTS Freq ({note})", TS_BLOCK)
    rec["ts_job"] = {"ok": r.ok, "opt_converged": r.opt_converged, "imaginary_cm": r.imaginary}
    if not r.ok or r.opt_converged is False:
        rec["result"] = "OptTS did not finish or converge: branch stopped (not rerun with other settings)"
        return rec
    gts = read_xyz(wd / "ts.xyz")
    rec.update(E_eh=r.energy_eh, G_eh=r.gibbs_eh, geometry=describe(gts),
               core_rmsd_to_production_ts_A=round(core_rmsd(gts, read_xyz(Path(f"{rx.ts_stem()}.xyz"))), 3))
    if len(r.imaginary) != 1:
        rec["result"] = f"{len(r.imaginary)} imaginary frequencies, not a TS: branch stopped"
        return rec
    rec["sides"] = {}
    for sign, side in ((+1, "plus"), (-1, "minus")):
        wdi, ri = job(rx, f"{folder}/{side}", side, "Opt Freq", displaced(wd, "ts", sign),
                      f"{folder}: TS displaced {'+' if sign > 0 else '−'}{STEP_A} Å along the imaginary mode, Opt Freq")
        if not ri.ok or ri.opt_converged is False:
            rec["sides"][side] = {"result": "optimization did not finish or converge"}
            continue
        c = classify(read_xyz(wdi / f"{side}.xyz"), ri.energy_eh, refs)
        c.update(E_eh=ri.energy_eh, imaginary_cm=ri.imaginary)
        rec["sides"][side] = c
    got = sorted(str(s.get("reached")) for s in rec["sides"].values())
    rec["connects"] = got
    rec["confirmed"] = got == sorted(want)
    rec["sp"] = sp_levels(rx, folder, gts)
    return rec


def scan_arm(rx, arm, target, refs):
    g = refs["secondary_anti"]["geom"]
    start = tors(g, C1, C2, C3, H)
    blk = f"%geom Scan D {C1} {C2} {C3} {H} = {start:.1f}, {target:.1f}, {N_SCAN} end end"
    wd, r = job(rx, f"ts_anti_eclip/scan_{arm}", "scan", "Opt", g,
                f"ts_anti_eclip: relaxed C1–C2–C3–H scan from secondary_anti ({start:.1f}° to {target:.1f}°)", blk)
    pts = []
    for i in range(1, N_SCAN + 1):
        p = wd / f"scan.{i:03d}.xyz"
        if not p.exists():
            break
        gi = read_xyz(p)
        m = re.search(r"E\s+(-?\d+\.\d+)", p.read_text().splitlines()[1])
        pts.append({"point": i, "E_eh": float(m.group(1)) if m else None, **describe(gi)})
    rec = {"arm": arm, "from_deg": round(start, 1), "to_deg": target, "job_ok": r.ok, "points": pts}
    shifted = [p["point"] for p in pts if p["H_on"].startswith("C2")]
    e = [p["E_eh"] for p in pts]
    rec["hydride_shifted_at_point"] = shifted[0] if shifted else None
    kmax = int(np.argmax(e)) if e and None not in e else None
    rec["max_point"] = pts[kmax]["point"] if kmax is not None else None
    rec["interior_max"] = kmax is not None and 0 < kmax < len(pts) - 1
    rec["guess"] = (str(wd / f"scan.{pts[kmax]['point']:03d}.xyz")
                  if r.ok and rec["interior_max"] and not shifted else None)
    return rec


# ------------------------------------------------------------------ report
def dg_rows(res, refs):
    anti = refs["secondary_anti"]
    tsj = json.loads((RX / "01_ts/ts.json").read_text(encoding="utf-8"))
    diag = json.loads((RX / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["energies_eh"]
    prod = {"name": "production TS (eclip side, 01_ts)", "E_eh": tsj["E_eh"], "G_eh": tsj["G_eh"],
            "sp": {lb: diag[lb]["ts"] for lb in diag}}
    labels = [None] + [lv for lv in anti["sp"]]
    rows = []
    for t in [prod] + [{"name": k, **v} for k, v in res.get("ts", {}).items() if v.get("sp")]:
        corr = (t["G_eh"] - t["E_eh"]) - (anti["G_eh"] - anti["E_eh"])
        row = {"ts": t["name"]}
        for lb in labels:
            de = (t["E_eh"] - anti["E_eh"]) if lb is None else (t["sp"][lb] - anti["sp"][lb])
            row[lb or "M06-2X_ma-def2-SVP (working)"] = {"dE_kJ": round(kj(de), 2), "dG_kJ": round(kj(de + corr), 2)}
        rows.append(row)
    return rows


def write_report(res, refs):
    res["barriers_from_secondary_anti"] = dg_rows(res, refs)
    OUT_JSON.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    f = lambda x: f"{x:.1f} ({x / 4.184:.1f})"
    L = ["# Rotation TSs from the anti minimum, secondary hydride shift (J. Picione, 2026-11-22)", "",
         "Working level M06-2X-D3(0)/ma-def2-SVP, SMD dichloromethane; G at 298.15 K (ORCA quasi-RRHO, working-level "
         "thermal correction); zero = secondary_anti. kJ/mol (kcal/mol). Written by "
         "`00_inputs/rotation_ts_m062x.py` from finished jobs.", ""]
    for k, v in res.get("ts", {}).items():
        L += [f"## {k}", f"- Start: {v.get('start')}", f"- TS job: {v.get('ts_job')}"]
        if v.get("geometry"):
            L.append(f"- TS geometry: {v['geometry']}; core RMSD to the production TS "
                     f"{v.get('core_rmsd_to_production_ts_A')} Å")
        for s, c in v.get("sides", {}).items():
            L.append(f"- {s} side: reached {c.get('reached')}; {c.get('geometry')}; "
                     + "; ".join(f"{n} RMSD {m['core_rmsd_A']} Å, ΔE {m['dE_kJ']:+.2f} kJ/mol"
                                 for n, m in c.get("matches", {}).items()))
        if "confirmed" in v:
            L.append(f"- Connects {v['connects']}; wanted {sorted(v['want'])}: "
                     f"**{'confirmed' if v['confirmed'] else 'NOT confirmed'}**")
        if v.get("result"):
            L.append(f"- **{v['result']}**")
        L.append("")
    for a in res.get("scans", []):
        L += [f"## Relaxed scan, arm {a['arm']} ({a['from_deg']}° → {a['to_deg']}°)",
              "| Point | C1–C2–C3–H (°) | E rel. secondary_anti (kJ/mol) | C3–H (Å) | C2–H (Å) | H on |", "|---|---|---|---|---|---|"]
        for p in a["points"]:
            e = kj(p["E_eh"] - refs["secondary_anti"]["E_eh"]) if p["E_eh"] is not None else float("nan")
            L.append(f"| {p['point']} | {p['C1_C2_C3_H_deg']} | {e:+.1f} | {p['C3_H_A']} | {p['C2_H_A']} | {p['H_on']} |")
        L += ["", f"- Hydride moved to C2 at point: {a['hydride_shifted_at_point']}; highest point {a['max_point']} "
              f"(interior: {a['interior_max']}); TS guess: {'yes' if a['guess'] else 'no'}", ""]
    L += ["## ΔG‡ and ΔE‡ from secondary_anti", "",
          "| TS | Level | ΔE‡ | ΔG‡ |", "|---|---|---|---|"]
    for row in res["barriers_from_secondary_anti"]:
        for lb, v in row.items():
            if lb != "ts":
                L.append(f"| {row['ts']} | {lb} | {f(v['dE_kJ'])} | {f(v['dG_kJ'])} |")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")


# ------------------------------------------------------------------ main
def gate(mode):
    need = GATES[mode]
    print(f"{time.strftime('%Y-%m-%d %H:%M')} waiting for the path single points of {', '.join(need)} ({mode})", flush=True)
    while not all((ROOT / "reactions" / n / "02_irc/energies_DLPNO-CCSD(T)_ma-def2-TZVP.json").exists() for n in need):
        time.sleep(60)
    print(f"{time.strftime('%Y-%m-%d %H:%M')} gate open; jobs now go to the core queue at low priority", flush=True)


def run(mode):
    # after_e2: low priority (takes only free cores). after_tertiary (J. Picione 2026-11-22 later: "a 1-1.5 h delay to
    # E2 is fine"): the ordinary first-come queue, so the jobs take turns with E2 instead of waiting for a gap.
    if mode == "after_e2":
        os.environ["CORES_LOW"] = "1"                        # in this process only
    gate(mode)
    rx = Reaction.load(RX)
    refs = references(rx)
    res = {"request": f"J. Picione 2026-11-22: rotation TSs from the anti minimum ({mode})", "ts": {}, "scans": []}
    # order (J. Picione 2026-11-22): the two torsion scans first, since whether anti and eclip interconvert at all
    # decides how Section 3.6 reports the barrier; then the anti-side rotation TS; then any anti-eclip TS search
    arms = []
    for arm, target in (("neg", ECLIP_DIH), ("pos", -ECLIP_DIH)):
        a = scan_arm(rx, arm, target, refs)
        res["scans"].append(a)
        arms.append(a)
        write_report(res, refs)
    res["ts"]["ts_anti_rotation"] = ts_branch(rx, "ts_anti_rotation", read_xyz(GAS_TS), refs,
                                              ["secondary_anti", "tertiary_irc"], "gas-phase M06-2X anti TS")
    write_report(res, refs)
    for a in arms:
        arm = a["arm"]
        if a["guess"]:
            res["ts"][f"ts_anti_eclip_{arm}"] = ts_branch(rx, f"ts_anti_eclip/ts_{arm}", read_xyz(Path(a["guess"])),
                                                          refs, ["secondary_anti", "secondary_eclip_irc"],
                                                          f"highest point {a['max_point']} of scan arm {arm}")
            write_report(res, refs)
    rx.log.note("Rotation TSs from the anti minimum (J. Picione, 2026-11-22): done; see `01_stationary/rotation_ts.md`.")
    print(OUT_MD.read_text(encoding="utf-8"))


def run_methyl():
    """Methyl-shift TS (J. Picione 2026-11-23): OptTS from scan_pos point 2, where the bridging methyl crosses between
    C3 and C2 (C2–CH₃ 1.763, C3–CH₃ 1.723 Å); sides must reach secondary_anti and its C2/C3-swapped copy. Ordinary
    core queue (he accepted the delay to E2). Added to rotation_ts.json/.md."""
    rx = Reaction.load(RX)
    refs = references(rx)
    res = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    res["ts"]["ts_methyl_shift"] = ts_branch(rx, "ts_methyl_shift", read_xyz(ST / "ts_anti_eclip/scan_pos/scan.002.xyz"),
                                             refs, ["secondary_anti", "secondary_anti_swapped"],
                                             "scan_pos point 2 (bridging methyl between C3 and C2)")
    write_report(res, refs)
    rx.log.note("Methyl-shift TS (J. Picione, 2026-11-23): done; see `01_stationary/rotation_ts.md`.")
    print(OUT_MD.read_text(encoding="utf-8"))


def check():
    rx = Reaction.load(RX)
    refs = references(rx)
    for n, r in refs.items():
        print(n, describe(r["geom"]), {m: round(core_rmsd(r["geom"], s["geom"]), 3) for m, s in refs.items()})
    g = read_xyz(GAS_TS)
    print("gas TS", describe(g), {m: round(core_rmsd(g, s["geom"]), 3) for m, s in refs.items()})
    h = orca.read_hess(GAS_HESS)
    print("gas TS lowest freq", round(float(np.min(h["freqs"])), 1))
    v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
    for sign in (+1, -1):
        gd = Geometry(g.symbols, g.coords + sign * STEP_A * v / np.linalg.norm(v, axis=1).max())
        print(f"displaced {sign:+d}", describe(gd))
    print(rx.inp("OptTS Freq", g, blocks=TS_BLOCK).split("* xyz")[0])
    a = refs["secondary_anti"]["geom"]
    print("scan block neg:", f"%geom Scan D {C1} {C2} {C3} {H} = {tors(a, C1, C2, C3, H):.1f}, {ECLIP_DIH:.1f}, {N_SCAN} end end")
    print("production TS check:", dg_rows({"ts": {}}, refs))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    if cmd == "run":
        run(sys.argv[2] if len(sys.argv) > 2 else "after_e2")
    elif cmd == "methyl":
        run_methyl()
    elif cmd == "check":
        check()
    else:
        sys.exit(__doc__)
