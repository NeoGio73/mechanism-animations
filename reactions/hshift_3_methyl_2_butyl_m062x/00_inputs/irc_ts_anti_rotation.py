"""IRC from ts_anti_rotation (J. Picione, 2026-11-27): settle which minima the anti form's rotation TS connects.

Until now its two sides were followed by displaced optimizations (0.10 Å along the imaginary mode, Opt Freq;
rotation_ts_m062x.py): plus side → secondary_anti, minus side → tertiary_irc. This run follows both sides by the
production IRC method instead:
- pysisyphus 1.0.0 Gonzalez–Schlegel, mass-weighted, step 0.05 bohr·amu½, exact TS Hessian from the existing
  01_stationary/ts_anti_rotation/ts/ts.hess (162.4i cm⁻¹), then Bofill updates (pipeline/true_irc_prod.run_side,
  unchanged);
- ORCA at the working level (reaction.toml method, SMD dichloromethane), 4 cores, %maxcore 2000 (J. Picione: another
  project's ORCA jobs share this computer);
- the production stopping rule of 2026-11-15 ((a) 20 flat steps below 0.05 kJ/mol after the last bond event, counted
  1 kJ/mol or more below the TS; (b) end-minimum test on 40 frames; 3000-step backstop). The [irc_end_min] entry of
  reaction.toml names the minimum of the animated path's backward side and is not used here: which minimum each side
  reaches is what this run is to find.
- the two sides one after the other (4 cores in all), forward first;
- one attempt; at 14 h of job time from the start everything stops, whatever the state (PAUSE files for the IRC sides,
  then this process's own ORCA jobs are ended after 5 minutes' grace); a stopped side is incomplete, never a result.
Then Opt Freq (working level, same 4 cores / 2000 MB) from the last frame of each side that has frames, and the match
of rotation_ts_m062x.py (core atoms within 0.10 Å RMSD and energies within 1.0 kJ/mol of secondary_anti,
secondary_eclip_irc or tertiary_irc; mirror image and the two C3 methyls swapped allowed; anti's C2/C3-swapped copy
told apart by the methyl's carbon). No single points at other levels; nothing on the card or tables changes.

Fixed before the run (J. Picione, 2026-11-27): ends secondary_anti + tertiary_irc → the present description stands;
secondary_anti + secondary_eclip_irc → the manuscript says the anti form reaches the tertiary cation through the
eclipsed form (anti barrier 18.6 kJ/mol either way); anything else → reported, nothing changed.

Output: 02_irc/ts_anti_rotation_irc/{forward,backward}/ (path.json, path.xyz, g####), .../end_<side>/ (opt), and
02_irc/ts_anti_rotation_irc/result.json/.md.
Usage: .venv-pysis/Scripts/python.exe -X utf8 -u reactions/hshift_3_methyl_2_butyl_m062x/00_inputs/irc_ts_anti_rotation.py run|report
"""
import importlib.util
import json
import sys
import threading
import time
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
ROOT = RX.parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import orca, steps, true_irc_prod                # noqa: E402
from pipeline.geometry import Geometry, read_xyz_frames         # noqa: E402
from pipeline.units import kj                                   # noqa: E402

TS_DIR = RX / "01_stationary/ts_anti_rotation/ts"
SUB = "ts_anti_rotation_irc"
OUT = RX / "02_irc" / SUB
NPROCS, MAXCORE = 4, 2000
CAP_S = 14 * 3600
GRACE_S = 300
T0 = time.time()

spec = importlib.util.spec_from_file_location("rot", RX / "00_inputs/rotation_ts_m062x.py")
rot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rot)

_load = steps.Reaction.load.__func__


def _load_capped(cls, path):
    rx = _load(cls, path)
    rx.cfg["nprocs"], rx.cfg["maxcore_mb"] = NPROCS, MAXCORE
    rx.cfg.pop("irc_end_min", None)
    return rx


steps.Reaction.load = classmethod(_load_capped)
true_irc_prod.Reaction = steps.Reaction
true_irc_prod.ts_paths = lambda rx: (TS_DIR / "ts.xyz", TS_DIR / "ts.hess", TS_DIR / "ts.out")
STOP = {"hit": False}


def log(msg):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "run.log", "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}\n")
    print(msg, flush=True)


def watchdog():
    """At 14 h: PAUSE every side folder (the IRC stops after its current step), then after GRACE_S end this process's
    own child processes (its ORCA jobs) if any are still running. Never touches any other process."""
    import psutil
    while time.time() - T0 < CAP_S:
        time.sleep(30)
    STOP["hit"] = True
    log("14 h of job time reached: stopping (PAUSE files written)")
    for d in ("forward", "backward"):
        if (OUT / d).exists():
            (OUT / d / "PAUSE").write_text("14 h limit of job time reached (J. Picione, 2026-11-27)", encoding="utf-8")
    time.sleep(GRACE_S)
    kids = psutil.Process().children(recursive=True)
    for p in kids:
        try:
            p.terminate()
        except psutil.Error:
            pass
    if kids:
        log(f"ended {len(kids)} of this run's own child processes after the grace period")


def side(direction):
    if STOP["hit"]:
        log(f"{direction}: not started (14 h limit)")
        return None
    log(f"{direction}: start (pysisyphus GS, step 0.05, {NPROCS} cores, %maxcore {MAXCORE})")
    try:
        res = true_irc_prod.run_side(RX, direction, check_dir=SUB)
    except Exception as e:                                   # a crash is recorded, never a result
        log(f"{direction}: FAILED ({type(e).__name__}: {e})")
        return {"ended": f"FAILED: {e}"}
    log(f"{direction}: ended '{res['ended']}', {res['n_points']} points, {res['n_engrad']} gradients, "
        f"{res['wall_total_s'] / 3600:.2f} h")
    return res


def end_opt(rx, direction, refs):
    pj = OUT / direction / "path.json"
    if not pj.exists():
        return None
    js = json.loads(pj.read_text(encoding="utf-8"))
    fr = read_xyz_frames(OUT / direction / "path.xyz")
    if not fr:
        return None
    rec = {"side_ended": js["ended"], "n_points": js["n_points"], "last_frame_E_minus_TS_kJ":
           kj(js["energies"][-1] - js["E_TS_eh"]), "last_frame": rot.classify(fr[-1], js["energies"][-1], refs)}
    if STOP["hit"]:
        rec["opt"] = "not run (14 h limit)"
        return rec
    wd = rx.d(f"02_irc/{SUB}/end_{direction}")
    g = Geometry(fr[-1].symbols, np.asarray(fr[-1].coords))
    r = orca.run(wd, "opt", rx.inp("Opt Freq", g), rx.log,
                 note=f"{SUB}: Opt Freq from the last frame of the {direction} side (J. Picione 2026-11-27)")
    if not r.ok:
        rec["opt"] = f"did not end normally ({r.out}): not a result"
        return rec
    imag = [f for f in orca.read_hess(wd / "opt.hess")["freqs"] if f < 0]
    go = read_xyz_frames(wd / "opt.xyz")[-1]
    rec["opt"] = {"E_eh": r.energy_eh, "imaginary_cm": [round(float(f), 1) for f in imag],
                  **rot.classify(go, r.energy_eh, refs)}
    return rec


def verdict(ends):
    got = sorted(e["opt"]["reached"] or "none" for e in ends.values()
                 if e and isinstance(e.get("opt"), dict) and not e["opt"]["imaginary_cm"])
    if len(got) < 2:
        return "incomplete: fewer than two optimized ends without an imaginary frequency"
    if got == ["secondary_anti", "tertiary_irc"]:
        return "secondary_anti + tertiary_irc: the present description stands"
    if got == ["secondary_anti", "secondary_eclip_irc"]:
        return ("secondary_anti + secondary_eclip_irc: the anti form reaches the tertiary cation through the eclipsed "
                "form (manuscript wording changes; anti barrier 18.6 kJ/mol either way)")
    return f"{' + '.join(got)}: neither fixed case; reported, nothing changed"


def write(res, ends):
    out = {"what": __doc__.split("\n")[0], "started": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(T0)),
           "job_time_h": (time.time() - T0) / 3600, "cap_hit": STOP["hit"], "nprocs": NPROCS, "maxcore_mb": MAXCORE,
           "sides": {d: ({k: v for k, v in r.items() if k not in ("energies", "gradients", "calls")} if r else None)
                     for d, r in res.items()}, "ends": ends, "verdict": verdict(ends)}
    (OUT / "result.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    L = ["# IRC from ts_anti_rotation (J. Picione, 2026-11-27)", "",
         f"Started {out['started']}; job time {out['job_time_h']:.2f} h; 14 h limit {'reached' if STOP['hit'] else 'not reached'}. "
         f"Working level, SMD dichloromethane, {NPROCS} cores, %maxcore {MAXCORE}; pysisyphus Gonzalez–Schlegel, step 0.05 "
         "bohr·amu½, exact TS Hessian from ts.hess; production stopping rule (no [irc_end_min]).", "",
         "| Side | Ended | Points | Last frame − TS | Opt Freq: imaginary | Reached | Core RMSD / ΔE to each minimum |",
         "|---|---|---|---|---|---|---|"]
    for d, e in ends.items():
        if not e:
            L.append(f"| {d} | no frames | | | | | |")
            continue
        o = e["opt"]
        if isinstance(o, dict):
            m = "; ".join(f"{n} {v['core_rmsd_A']:.3f} Å / {v['dE_kJ']:+.2f} kJ/mol" for n, v in o["matches"].items())
            L.append(f"| {d} | {e['side_ended']} | {e['n_points']} | {e['last_frame_E_minus_TS_kJ']:+.2f} kJ/mol | "
                     f"{o['imaginary_cm'] or 'none'} | {o['reached'] or 'none'} | {m} |")
        else:
            L.append(f"| {d} | {e['side_ended']} | {e['n_points']} | {e['last_frame_E_minus_TS_kJ']:+.2f} kJ/mol | {o} | | |")
    L += ["", f"**Reading (rule fixed before the run):** {out['verdict']}.", ""]
    (OUT / "result.md").write_text("\n".join(L), encoding="utf-8")
    log("result written: " + out["verdict"])


def main():
    rx = steps.Reaction.load(RX)
    refs = rot.references(rx)
    threading.Thread(target=watchdog, daemon=True).start()
    log(f"start; TS {TS_DIR.relative_to(RX).as_posix()}; limit 14 h of job time")
    res = {d: side(d) for d in ("forward", "backward")}
    ends = {}
    for d in ("forward", "backward"):
        try:
            ends[d] = end_opt(rx, d, refs)
        except Exception as e:
            ends[d] = {"side_ended": (res.get(d) or {}).get("ended"), "n_points": None,
                       "last_frame_E_minus_TS_kJ": float("nan"), "opt": f"FAILED ({type(e).__name__}: {e})"}
    write(res, ends)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "run":
        main()
    elif cmd == "check":                       # no ORCA job: settings and the first input
        rx = steps.Reaction.load(RX)
        print(true_irc_prod.ts_paths(rx), rx.cfg["nprocs"], rx.cfg["maxcore_mb"], rx.cfg.get("irc_end_min"))
        print(rx.inp("EnGrad NoPop", read_xyz_frames(TS_DIR / "ts.xyz")[-1]).split("* xyz")[0])
        print({k: str(v["xyz"]) for k, v in rot.references(rx).items()})
    else:
        sys.exit(__doc__)
