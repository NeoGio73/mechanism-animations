"""Step-size convergence test (J. Picione, 2026-11-11), SN2 rerun TS, B3LYP-D3(BJ)/ma-def2-SVP, SMD water.

Non-adaptive IRC (Adapt_Scale_Displ false, Do_SD_Corr unchanged), first step by length, Scale_Displ_SD
calibrated so the median step is about the target:
  s010  first step 0.10, median about 0.10 bohr·amu½ (the production setting)
  s005  first step 0.05, median about 0.05 bohr·amu½
Calibration from the 0.10 test (variant A): Scale_Displ_SD 0.02116 gave a median step of 0.0897;
steps scale about linearly with Scale_Displ_SD (0.05 gave about 0.21–0.24), so
Scale_Displ_SD = 0.02116 × target / 0.0897.
Both paths run far enough to cover the whole bond-changing region (on the adaptive path s = −10.7 to
+6.3 bohr·amu½; the 0.10 test reached only −3.8 to +3.6).

Comparison (criterion: within 0.2 kJ/mol and 0.005 Å throughout the bond-changing region):
  reactant side (forward, C–O long): at matching C–O distance, energy and C–Br
  product side (backward, C–Br long): at matching C–Br distance, energy and C–O
Reported both over the span of the first 80 steps of s005 (the span asked for) and over the whole
bond-changing region. Also per path: first step, median, shortest and longest step in the region,
largest change of a key bond in one step.
Run: .venv/Scripts/python.exe -X utf8 <this file> run|report|early
  early: s005 as far as it has got (first 80 steps) against the 0.10 test path (variant A).
"""
import json
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))
from pipeline import orca                                       # noqa: E402
from pipeline.geometry import read_xyz, read_xyz_frames         # noqa: E402
from pipeline.steps import Reaction, irc_init_scale             # noqa: E402

SD_PER_S = 0.02116 / 0.0897
RUNS = {"s010": {"target": 0.10, "maxiter": {"forward": 100, "backward": 160}},
        "s005": {"target": 0.05, "maxiter": {"forward": 190, "backward": 300}},
        # third path, only if the whole-region 0.05/0.10 comparison fails (J. Picione, 2026-11-11)
        "s0025": {"target": 0.025, "maxiter": {"forward": 380, "backward": 620}}}   # 2 cores per job, as 0.05 and 0.10 (J. Picione)
OUT = ROOT / "02_irc/attempts/step_convergence"
MASS = {"C": 12.011, "H": 1.008, "O": 15.999, "Br": 79.904}
CO, CBR = (0, 2), (0, 1)
H = 2625.4996
TOL_E, TOL_B = 0.2, 0.005
E_TS = -2689.43410109009   # Eh, the TS energy every IRC job here computes first (same input, same TS)


def jd(name, d):
    """Folder of one direction's job. From 19:45 each direction runs in its own folder
    (<name>/<direction>/): ORCA writes the dispersion output to a fixed file, redirect.tmp, in the job
    folder, and two IRC jobs side by side in one folder overwrote it (the likely cause of the crashes).
    Jobs that finished earlier stay in <name>/."""
    sub = OUT / name / d
    return sub if (sub / f"irc_{d[0]}.inp").exists() else OUT / name


def trj(name, d):
    return jd(name, d) / f"irc_{d[0]}_IRC_{d[0].upper()}_trj.xyz"


def job(rx, name, d):
    r = RUNS[name]
    wd = OUT / name / d
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "01_ts/ts/ts.hess", wd / "ts.hess")
    sc = irc_init_scale(ROOT / "01_ts/ts/ts.hess", r["target"])
    sd = SD_PER_S * r["target"]
    blk = (f'%irc\n  MaxIter {r["maxiter"][d]}\n  Direction {d}\n  InitHess read\n  Hess_Filename "ts.hess"\n'
           f'  Scale_Displ_SD {sd:.5f}\n  Adapt_Scale_Displ false\n  Init_Displ length\n'
           f'  Scale_Init_Displ {sc:.5f}\n  PrintLevel 1\nend')
    inp = orca.make_input(rx.keywords("IRC"), read_xyz(ROOT / "01_ts/ts/ts.xyz"), rx.cfg["charge"],
                          rx.cfg["mult"], blocks=rx.solvent_block + "\n" + blk, nprocs=r.get("nprocs", 2),
                          maxcore_mb=rx.cfg["maxcore_mb"])
    orca.run(wd, f"irc_{d[0]}", inp, rx.log,
             note=f"step-size convergence test {name}, {d}, Adapt_Scale_Displ false, Scale_Displ_SD {sd:.5f}, "
                  f"first step {r['target']} bohr·amu½, MaxIter {r['maxiter'][d]}")


def dist(c, ij):
    return float(np.linalg.norm(c[ij[0]] - c[ij[1]]))


def path(fwd, bwd, ts):
    """Frames TS-ordered: list of dict(s, step, E kJ/mol rel. TS, CO, CBr, side). Partial files are fine."""
    m = np.array([MASS[s] for s in ts.symbols])
    x0 = np.asarray(ts.coords)
    out = [{"s": 0.0, "step": 0.0, "E": 0.0, "CO": dist(x0, CO), "CBr": dist(x0, CBR), "side": "TS"}]
    Ets = None
    for sign, trj, side in ((1, fwd, "reactant"), (-1, bwd, "product")):
        if not Path(trj).exists():
            continue
        prev, s = x0 / 0.52917721, 0.0
        for f in read_xyz_frames(trj):
            c = np.asarray(f.coords)
            x = c / 0.52917721
            st = float(np.sqrt((m[:, None] * (x - prev) ** 2).sum()))
            if st < 1e-6:
                continue
            s += st
            prev = x
            out.append({"s": sign * s, "step": st, "E": float(f.comment.split()[-1]), "CO": dist(c, CO),
                        "CBr": dist(c, CBR), "side": side})
    return sorted(out, key=lambda r: r["s"])


def with_energy(p, E_ts):
    for r in p:
        if r["side"] != "TS":
            r["E"] = (r["E"] - E_ts) * H
    return p


def region_flags(p, ends):
    """In the bond-changing region: either key bond has covered 5–95 % of its change between the
    path ends (taken from the adaptive IRC of 02_irc, the same minima)."""
    (co0, cb0), (co1, cb1) = ends
    for r in p:
        pc = [(r["CO"] - co0) / (co1 - co0), (r["CBr"] - cb0) / (cb1 - cb0)]
        r["region"] = any(0.05 < q < 0.95 for q in pc)
    return p


def adaptive_ends():
    full = read_xyz_frames(ROOT / "02_irc/irc_IRC_Full_trj.xyz")
    a, b = np.asarray(full[0].coords), np.asarray(full[-1].coords)
    return (dist(a, CO), dist(a, CBR)), (dist(b, CO), dist(b, CBR))


def stats(p):
    """Steps in the region (a step counts if either frame is in it), and largest key-bond change per step."""
    st, dmax = [], (0.0, None)
    for a, b in zip(p[:-1], p[1:]):
        if a["region"] or b["region"]:
            step = abs(b["s"] - a["s"])
            st.append(step)
            for k in ("CO", "CBr"):
                if abs(b[k] - a[k]) > dmax[0]:
                    dmax = (abs(b[k] - a[k]), f"{'C–O' if k == 'CO' else 'C–Br'} at s {a['s']:+.2f}→{b['s']:+.2f}")
    first = [r["step"] for r in p if r["side"] != "TS" and abs(abs(r["s"]) - r["step"]) < 1e-9]
    return {"first_steps": first, "median_region": float(np.median(st)), "min_region": float(min(st)),
            "max_region": float(max(st)), "n_region_steps": len(st), "max_key_bond_change_A": dmax[0],
            "where": dmax[1], "s_range": [p[0]["s"], p[-1]["s"]]}


def compare(ref, fine, span=None):
    """ref frames vs the fine path interpolated at the same matching distance."""
    res = {}
    for side, key, other in (("reactant", "CO", "CBr"), ("product", "CBr", "CO")):
        F = [r for r in fine if r["side"] in (side, "TS")]
        F.sort(key=lambda r: r[key])
        fx = np.array([r[key] for r in F])
        if np.any(np.diff(fx) <= 0):
            res[side] = {"error": f"{key} not monotonic on the fine path"}
            continue
        rows = []
        for r in ref:
            if r["side"] != side or not r["region"]:
                continue
            if span and not (span[0] <= r["s"] <= span[1]):
                continue
            if not (fx[0] <= r[key] <= fx[-1]):
                continue
            dE = r["E"] - np.interp(r[key], fx, [q["E"] for q in F])
            dB = r[other] - np.interp(r[key], fx, [q[other] for q in F])
            rows.append((r["s"], r[key], float(dE), float(dB)))
        if not rows:
            res[side] = {"n": 0}
            continue
        iE = max(range(len(rows)), key=lambda i: abs(rows[i][2]))
        iB = max(range(len(rows)), key=lambda i: abs(rows[i][3]))
        res[side] = {"match_by": "C–O" if key == "CO" else "C–Br", "other": "C–Br" if other == "CBr" else "C–O",
                     "n": len(rows), "s_covered": [min(r[0] for r in rows), max(r[0] for r in rows)],
                     "max_dE_kJ": rows[iE][2], "at_s": rows[iE][0], "at_dist_A": rows[iE][1],
                     "max_dbond_A": rows[iB][3], "bond_at_s": rows[iB][0],
                     "pass": abs(rows[iE][2]) <= TOL_E and abs(rows[iB][3]) <= TOL_B,
                     "rows": rows}
    return res


def report(early=False):
    ts = read_xyz(ROOT / "01_ts/ts/ts.xyz")
    ends = adaptive_ends()
    f05 = path(trj("s005", "forward"), trj("s005", "backward"), ts)
    E_ts = E_TS
    f05 = region_flags(with_energy(f05, E_ts), ends)
    n80 = [r["s"] for r in f05 if r["side"] != "TS"]
    fw = sorted(s for s in n80 if s > 0)[:80]
    bw = sorted((s for s in n80 if s < 0), reverse=True)[:80]
    span80 = [bw[-1] if bw else 0.0, fw[-1] if fw else 0.0]
    if early:
        A = ROOT / "02_irc/attempts/constant_step/A"
        ref = region_flags(with_energy(path(A / "irc_f_IRC_F_trj.xyz", A / "irc_b_IRC_B_trj.xyz", ts), E_ts), ends)
        out = {"ref": "0.10 test path (variant A, Scale_Displ_SD 0.02116)", "span_first_80_steps_of_s005": span80,
               "compare": compare(ref, f05, span80), "stats_ref": stats(ref), "stats_s005_so_far": stats(f05)}
        name = "report_early.json"
    else:
        ref = region_flags(with_energy(path(trj("s010", "forward"), trj("s010", "backward"), ts),
                                       E_ts), ends)
        out = {"ref": "s010", "span_first_80_steps_of_s005": span80,
               "compare_span80": compare(ref, f05, span80), "compare_region": compare(ref, f05),
               "stats_s010": stats(ref), "stats_s005": stats(f05),
               "region_on_adaptive_path_s": [-10.71, 6.33]}
        name = "report.json"
    (OUT / name).write_text(json.dumps(out, indent=1), encoding="utf-8")
    def show(d, pre=""):
        for k, v in d.items():
            if isinstance(v, dict):
                print(pre + k)
                show(v, pre + "  ")
            elif k != "rows":
                print(pre + f"{k}: {np.round(v, 4).tolist() if isinstance(v, (list, float)) else v}")
    show(out)


def load(name, ts, ends):
    return region_flags(with_energy(path(trj(name, "forward"), trj(name, "backward"), ts), E_TS), ends)


def covers_region(p):
    """The path ran past the bond-changing region at both ends (its end frames are outside it)."""
    return (not p[0]["region"]) and (not p[-1]["region"])


def finished(name, d):
    o = jd(name, d) / f"irc_{d[0]}.out"
    return o.exists() and "ORCA TERMINATED NORMALLY" in o.read_text(errors="replace")


def verdict():
    """Whole-region pass/fail of s010 against s005 (both sides pass, both paths cover the region).
    Returns True (pass), False (fail) or None (incomplete: a job did not finish normally or a path
    does not cover the region; not a step-size result, so no 0.025 run on it)."""
    ts = read_xyz(ROOT / "01_ts/ts/ts.xyz")
    ends = adaptive_ends()
    p10, p05 = load("s010", ts, ends), load("s005", ts, ends)
    c = compare(p10, p05)
    cov = {"s010": covers_region(p10), "s005": covers_region(p05)}
    done = {f"{n}/{d}": finished(n, d) for n in ("s010", "s005") for d in ("forward", "backward")}
    ok = None if not (all(done.values()) and all(cov.values())) else all(v.get("pass") for v in c.values())
    out = {"pass": ok, "jobs_finished_normally": done, "covers_region": cov,
           "sides": {k: {a: b for a, b in v.items() if a != "rows"} for k, v in c.items()}}
    (OUT / "verdict.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return ok


def report3():
    """0.10, 0.05, 0.025 side by side; 0.025 is the reference. For each coarser path: largest energy
    and other-bond difference per side, where, and whether within 0.2 kJ/mol and 0.005 Å."""
    ts = read_xyz(ROOT / "01_ts/ts/ts.xyz")
    ends = adaptive_ends()
    P = {n: load(n, ts, ends) for n in RUNS}
    out = {"stats": {n: stats(p) for n, p in P.items()},
           "covers_region": {n: covers_region(p) for n, p in P.items()},
           "vs_0.025": {}, "0.10_vs_0.05": {}}
    for n in ("s010", "s005"):
        c = compare(P[n], P["s0025"])
        out["vs_0.025"][n] = {k: {a: b for a, b in v.items() if a != "rows"} for k, v in c.items()}
        # where the difference exceeds the criterion: s ranges of the failing frames
        for k, v in c.items():
            bad = [r[0] for r in v.get("rows", []) if abs(r[2]) > TOL_E or abs(r[3]) > TOL_B]
            out["vs_0.025"][n][k]["failing_s"] = [min(bad), max(bad)] if bad else None
            out["vs_0.025"][n][k]["n_failing"] = len(bad)
    c = compare(P["s010"], P["s005"])
    out["0.10_vs_0.05"] = {k: {a: b for a, b in v.items() if a != "rows"} for k, v in c.items()}
    (OUT / "report3.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))


def after(rx):
    """Verdict, then the standing instruction (J. Picione, 2026-11-11): on a whole-region fail, run
    0.025 overnight and write the three-way report. Exactly one marker (PASS, FAIL, INCOMPLETE) is
    left; the 0.025 run starts at most once (two lanes may reach this point). STEP_TEST_COMPLETE marks
    the end of the step test (the E2 arbitration lane waits for it).

    Gate (J. Picione, 2026-11-11, check 1): no verdict until every finished 0.10 and 0.05 job has
    been cleared of folder sharing: check1.json (pipeline/shared_folder_check.py) must list all four
    as "clear". Until then: HELD, and any marker written by older copies of this script is removed."""
    import os
    gate = OUT / "check1.json"
    need = [f"{n}/{d}" for n in ("s010", "s005") for d in ("forward", "backward")]
    st = json.loads(gate.read_text(encoding="utf-8")).get("status", {}) if gate.exists() else {}
    if not all(st.get(k) == "clear" for k in need):
        for m in ("PASS", "FAIL", "INCOMPLETE", "STEP_TEST_DONE", "STEP_TEST_COMPLETE"):
            (OUT / m).unlink(missing_ok=True)
        (OUT / "HELD").write_text("no verdict yet: check 1 (shared folders) has not cleared all four jobs: "
                                  + json.dumps({k: st.get(k, "not checked") for k in need}) + "\n")
        print("HELD: check 1 not cleared", flush=True)
        return
    (OUT / "HELD").unlink(missing_ok=True)
    v = verdict()
    for m in ("PASS", "FAIL", "INCOMPLETE"):
        (OUT / m).unlink(missing_ok=True)
    if v is None:
        (OUT / "INCOMPLETE").write_text("a job did not finish normally or a path does not cover the "
                                        "region: no step-size verdict, no 0.025 run\n")
    elif v:
        (OUT / "PASS").write_text("whole-region comparison passed\n")
        (OUT / "VERDICT_GIVEN").write_text("PASS\n")
        (OUT / "STEP_TEST_COMPLETE").write_text("whole-region comparison passed\n")
    else:
        (OUT / "FAIL").write_text("whole-region comparison failed; running 0.025\n")
        (OUT / "VERDICT_GIVEN").write_text("FAIL\n")      # the E2 lane starts on this (J. Picione)
        try:   # s0025.started is a block for copies of this script loaded before 19:55
            os.close(os.open(OUT / "s0025.claimed", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            print("0.025 already started by another lane", flush=True)
            return
        # J. Picione, 2026-11-11: hold the 0.025 path until the Unity window is over
        import time
        win = ROOT.parents[1] / "docs/level_test_queue/window.json"
        while win.exists() and json.loads(win.read_text(encoding="utf-8")).get("status") != "released":
            time.sleep(60)
        with ThreadPoolExecutor(2) as ex:
            list(ex.map(lambda d: job(rx, "s0025", d), ("forward", "backward")))
        report3()
        (OUT / "STEP_TEST_COMPLETE").write_text("0.025 finished; report3.json written\n")


if __name__ == "__main__":
    rx = Reaction.load(ROOT)
    if sys.argv[1] == "after":
        after(rx)
    elif sys.argv[1] == "rerun_one":
        # one job that ended without "ORCA TERMINATED NORMALLY": outputs moved aside, same settings
        # rerun, then the verdict (J. Picione, 2026-11-11: a crashed job is never a result)
        import time
        n, d = sys.argv[2], sys.argv[3]
        if not finished(n, d):
            t_crash = (jd(n, d) / f"irc_{d[0]}.out").stat().st_mtime
            dest = OUT / n / f"crashed_{time.strftime('%H%M', time.localtime(t_crash))}_{d}"
            dest.mkdir(exist_ok=True)
            for f in jd(n, d).glob(f"irc_{d[0]}*"):
                if f.is_file():
                    shutil.move(str(f), dest / f.name)
            print("rerunning", n, d, flush=True)
            job(rx, n, d)
        after(rx)
    elif sys.argv[1] == "rerun_new":
        # J. Picione, 2026-11-11: s010 backward shared its folder (15:44–16:07); rerun alone in its own
        # folder (<name>/<direction>/) with identical settings; the old run stays in <name>/ for the
        # old-vs-new comparison. jd() then points the verdict at the rerun.
        n, d = sys.argv[2], sys.argv[3]
        print("rerunning in its own folder", n, d, flush=True)
        job(rx, n, d)
    elif sys.argv[1] == "rerun_alone":
        # a finished job that failed check 1 (it shared its folder with another running job) is not a
        # result: outputs moved to shared_folder_<time>_<direction>/, rerun alone in its own folder
        import time
        n, d = sys.argv[2], sys.argv[3]
        src = jd(n, d)
        dest = OUT / n / f"shared_folder_{time.strftime('%H%M')}_{d}"
        dest.mkdir(exist_ok=True)
        for f in src.glob(f"irc_{d[0]}*"):
            if f.is_file():
                shutil.move(str(f), dest / f.name)
        if src != OUT / n:
            (src / f"irc_{d[0]}.inp").unlink(missing_ok=True)
        print("rerunning alone", n, d, flush=True)
        job(rx, n, d)
    elif sys.argv[1] == "rerun_crashed":
        # jobs that ended without "ORCA TERMINATED NORMALLY": outputs moved aside, same settings rerun
        import time
        todo = [(n, d) for n in ("s010", "s005") for d in ("forward", "backward") if not finished(n, d)]
        for n, d in todo:
            dest = OUT / n / f"crashed_{time.strftime('%H%M')}_{d}"
            dest.mkdir(exist_ok=True)
            for f in (OUT / n).glob(f"irc_{d[0]}*"):
                if f.is_file():
                    shutil.move(str(f), dest / f.name)
        print("rerunning", todo, flush=True)
        with ThreadPoolExecutor(max(len(todo), 1)) as ex:
            list(ex.map(lambda a: job(rx, *a), todo))
        report()
    elif sys.argv[1] == "report3":
        report3()
    elif sys.argv[1] == "run":
        with ThreadPoolExecutor(4) as ex:
            list(ex.map(lambda a: job(rx, *a), [(n, d) for n in RUNS for d in ("forward", "backward")]))
        report()
    elif sys.argv[1] == "early":
        report(early=True)
    else:
        report()
