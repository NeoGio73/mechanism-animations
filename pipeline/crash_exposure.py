"""Crash checks from existing logs only (J. Picione, 2026-11-11), no new jobs, no files of running jobs.

1. compare: every finished IRC job today and before, with its input, cores, path length and steps, to
   list what differs between the step-size test IRC jobs (which crashed) and the others (which did not).
2. exposure: SCF calls that printed a D3 block, per group, and crashes per 1,000 calls.
   - Single points, optimizations, frequencies: the D3 blocks printed in the .out.
   - IRC jobs: one SCF per energy/gradient evaluation, each printing its D3 block into <base>.lastscf
     (overwritten each time, so not kept). Counted exactly where the per-evaluation property file
     survives (crashed jobs: one $VdW_Correction record per evaluation); for jobs that ended normally
     ORCA deletes that file, so the count is the IRC iterations × evaluations per iteration measured on
     the crashed jobs of the same kind (stated as an estimate with its range).
Usage: python -m pipeline.crash_exposure
"""
import json
import os
import re
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions"
RUNNING = ("s010/backward/", "s005/backward/", "bridged/gas/scan", "s005/forward/")   # not read
ITER = re.compile(r"^\s+\d+\s+-\d+\.\d+\s+-?\d+\.\d+\s+\d\.\d+\s+\d\.\d+\s*$", re.M)
CRASHED = {"s010/crashed_1706_forward", "s005/crashed_1706_forward", "s005/crashed_1706_backward",
           "s005/crashed_1722_forward", "s005/crashed_1941_backward", "s005/crashed_2121_forward"}


def outputs():
    for o in RX.rglob("*.out"):
        s = str(o.relative_to(RX)).replace("\\", "/")
        if any(r in s for r in RUNNING) or s.endswith(".smd.out"):
            continue
        try:
            t = o.read_text(errors="replace")
        except OSError:
            continue
        if "O   R   C   A" not in t[:20000]:
            continue
        yield o, s, t


def inp_lines(t):
    return [x.rstrip() for x in re.findall(r"^\|\s+\d+>\s?(.*)$", t, re.M)]


def group(s, kw):
    if "step_convergence" in s:
        return "step-size test IRC"
    if re.search(r"\bIRC\b", kw):
        return "other IRC"
    return "non-IRC"


RT = re.compile(r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec")
# jobs running now: interval from file metadata only (creation time of the .out), cores from the ledger
RUNNING_NOW = {"sn2_hydroxide_bromomethane_masvp/02_irc/attempts/step_convergence/s005/backward/irc_b.out": 2,
               "sn2_hydroxide_bromomethane_masvp/02_irc/attempts/step_convergence/s010/backward/irc_b.out": 2,
               "hshift_2_3_dimethyl_2_butyl/01_ts/level_test/M06-2X/bridged/gas/scan/scan.out": 4}


def load(rows):
    """Cores in use (sum of nprocs of all ORCA jobs running) while each IRC job ran today, and crash
    rates per 1,000 D3 SCF calls with the machine at or above 8 cores vs below. Calls are spread evenly
    over each job's run time. Intervals: end = last write of the .out, start = end − ORCA's total run
    time, else the creation time of the .out."""
    from .folder_overlap import run_log_walls
    walls = run_log_walls()
    iv = []
    for r in rows:
        p = RX / r["out"]
        end = p.stat().st_mtime
        t = p.read_text(errors="replace")[-3000:]
        m = RT.search(t)
        if m:
            start = end - (int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]))
        else:
            rx_, rel = r["out"].split("/", 1)
            rel = re.sub(r"/(crashed|stopped)_[^/]+", "", rel)[:-4]
            near = [w for e, w in walls.get((rx_, rel), []) if abs(e - end) < 120]
            ct = p.stat().st_ctime
            start = end - near[0] if near else (ct if end - ct < 6 * 3600 else None)
        r["_iv"] = (start, end)
        if start is not None:
            iv.append((start, end, r["nprocs"]))
    for s, n in RUNNING_NOW.items():
        p = RX / s
        if p.exists():
            iv.append((p.stat().st_ctime, time.time(), n))
    t0 = time.mktime(time.strptime("2026-11-11 00:00", "%Y-%m-%d %H:%M"))

    def cores_at(t):
        return sum(n for a, b, n in iv if a <= t < b)

    # "hi": more cores busy than the 8 physical ones (> 8); "lo": 8 or fewer
    res = {"step-size test IRC": {"hi": [0.0, 0], "lo": [0.0, 0]}, "other IRC": {"hi": [0.0, 0], "lo": [0.0, 0]}}
    per_job = []
    for r in rows:
        if not r["irc_iterations"] or r["group"] == "non-IRC" or r["_iv"][1] < t0 or r["_iv"][0] is None:
            continue
        a, b = r["_iv"]
        ts = [a + k * 10 for k in range(max(1, int((b - a) / 10)))]
        cs = [cores_at(t) for t in ts]
        hi = sum(c > 8 for c in cs) / len(cs)
        calls = r["d3_calls"][1]
        g = res[r["group"]]
        g["hi"][0] += calls * hi
        g["lo"][0] += calls * (1 - hi)
        if r["crashed"]:
            g["hi" if cores_at(b - 1) > 8 else "lo"][1] += 1
        per_job.append({"out": r["out"], "crashed": r["crashed"], "share_of_time_above_8_cores": round(hi, 2),
                        "max_cores": max(cs), "cores_at_end": cores_at(b - 1)})
    for g in res.values():
        for k in ("hi", "lo"):
            c, n = g[k]
            g[k] = {"d3_calls": round(c), "crashes": n, "per_1000": round(1000 * n / c, 2) if c else None}
    return {"by_machine_load": res, "per_job": per_job, "by_folder_sharing": sharing(rows)}


def sharing(rows):
    """D3 SCF calls made while another ORCA job ran in the same folder (folder at run time: outputs
    moved to crashed_*/stopped_* ran in the folder above) vs calls made alone; crashes counted by
    whether another job ran in the folder at the crash (within 15 s)."""
    def folder(r):
        return re.sub(r"/(crashed|stopped)_[^/]+", "", r["out"].rsplit("/", 1)[0])

    jobs = [r for r in rows if r.get("_iv") and r["_iv"][0] is not None]
    acc = {"shared": [0.0, 0, 0], "alone": [0.0, 0, 0]}    # calls, crashes, runs touched
    per = []
    for r in jobs:
        if not r["irc_iterations"] or not r["has_d3"]:
            continue
        a, b = r["_iv"]
        others = [o["_iv"] for o in jobs if o is not r and folder(o) == folder(r)]
        ts = [a + k * 10 for k in range(max(1, int((b - a) / 10)))]
        sh = sum(any(x <= t < y for x, y in others) for t in ts) / len(ts)
        calls = r["d3_calls"][1]
        acc["shared"][0] += calls * sh
        acc["alone"][0] += calls * (1 - sh)
        if r["crashed"]:
            k = "shared" if any(x <= b and y >= b - 15 for x, y in others) else "alone"
            acc[k][1] += 1
        per.append({"out": r["out"], "crashed": r["crashed"], "share_of_time_sharing_folder": round(sh, 3),
                    "ever_shared": sh > 0})
    out = {k: {"d3_calls": round(v[0]), "crashes": v[1],
               "per_1000": round(1000 * v[1] / v[0], 2) if v[0] else None} for k, v in acc.items()}
    ever = [p for p in per if p["ever_shared"]]
    out["runs_that_ever_shared_a_folder"] = {"runs": len(ever), "crashed": sum(p["crashed"] for p in ever)}
    never = [p for p in per if not p["ever_shared"]]
    out["runs_that_never_shared"] = {"runs": len(never), "crashed": sum(p["crashed"] for p in never)}
    out["per_run"] = per
    return out


def main():
    rows = []
    for o, s, t in outputs():
        il = inp_lines(t)
        kw = next((x for x in il if x.strip().startswith("!")), "")
        g = group(s, kw)
        ok = "ORCA TERMINATED NORMALLY" in t[-3000:]
        crashed = any(c in s for c in CRASHED)
        n_it = len(ITER.findall(t)) if "IRC" in kw else 0
        d3_out = t.count("DFT DISPERSION CORRECTION")
        prop = next(o.parent.glob(f"{o.stem}_IRC_*.property.txt"), None) if "IRC" in kw else None
        evals = len(re.findall(r"^\$VdW_Correction", prop.read_text(errors="replace"), re.M)) if prop else None
        np_ = re.search(r"nprocs\s+(\d+)", "\n".join(il))
        rows.append({"out": s, "group": g, "ok": ok, "crashed": crashed, "keywords": kw.strip(),
                     "irc_block": [x.strip() for x in il if re.match(r"\s*(MaxIter|Direction|Scale_Displ_SD|"
                                   r"Adapt_Scale_Displ|Init_Displ|Scale_Init_Displ|PrintLevel|InitHess|Hess_Filename|Do_SD_Corr)", x)],
                     "nprocs": int(np_.group(1)) if np_ else 1, "path_len": len(str(o.resolve())),
                     "folder_len": len(str(o.parent.resolve())), "irc_iterations": n_it, "d3_blocks_in_out": d3_out,
                     "evals_recorded": evals, "has_d3": bool(re.search(r"D3BJ|D3ZERO|D3\b|-D3", kw)),
                     "end": time.strftime("%Y-%m-%d %H:%M", time.localtime(o.stat().st_mtime))})
    # evaluations per IRC iteration, from the jobs that kept their per-evaluation record
    ratios = [r["evals_recorded"] / r["irc_iterations"] for r in rows      # jobs with >= 20 steps
              if r["evals_recorded"] and r["irc_iterations"] >= 20 and r["group"] == "step-size test IRC"]
    rmin, rmax, rmean = min(ratios), max(ratios), sum(ratios) / len(ratios)
    for r in rows:
        if r["irc_iterations"]:
            if r["evals_recorded"]:
                r["d3_calls"] = (r["evals_recorded"],) * 3
            else:
                n = r["irc_iterations"] + 1           # + the first single point at the TS
                r["d3_calls"] = (n * rmin, n * rmean, n * rmax)
        else:
            r["d3_calls"] = (r["d3_blocks_in_out"],) * 3
        if not r["has_d3"]:
            r["d3_calls"] = (0, 0, 0)
    out = {"evals_per_irc_iteration": {"min": rmin, "mean": rmean, "max": rmax, "n_jobs": len(ratios)}, "groups": {}}
    for g in ("step-size test IRC", "other IRC", "non-IRC"):
        R = [r for r in rows if r["group"] == g]
        today = [r for r in R if r["end"].startswith("2026-11-11")]
        for label, S in (("all days", R), ("today", today)):
            calls = [sum(r["d3_calls"][k] for r in S) for k in range(3)]
            cr = sum(r["crashed"] for r in S)
            out["groups"][f"{g}, {label}"] = {
                "jobs": len(S), "crashes": cr, "d3_scf_calls": [round(c) for c in calls],
                "crashes_per_1000": [round(1000 * cr / c, 3) if c else None for c in (calls[2], calls[1], calls[0])]}
    (ROOT / "docs/level_test_queue/crash_exposure.json").write_text(
        json.dumps({"summary": out, "jobs": rows}, indent=1, default=str), encoding="utf-8")
    out["load"] = load(rows)
    (ROOT / "docs/level_test_queue/crash_exposure.json").write_text(
        json.dumps({"summary": out, "jobs": rows}, indent=1, default=str), encoding="utf-8")
    print(json.dumps(out, indent=1))
    print("\nIRC jobs:")
    for r in sorted((r for r in rows if r["irc_iterations"]), key=lambda r: r["end"]):
        print(f"{r['end']} {'CRASH' if r['crashed'] else ('ok   ' if r['ok'] else 'stop ')} np={r['nprocs']} "
              f"it={r['irc_iterations']:4d} calls~{r['d3_calls'][1]:6.0f} pathlen={r['path_len']:3d} {r['out']}")
        print(f"      {r['keywords'][:100]} | {' ; '.join(r['irc_block'])[:200]}")


if __name__ == "__main__":
    main()
