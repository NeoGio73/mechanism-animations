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
        "s005": {"target": 0.05, "maxiter": {"forward": 190, "backward": 300}}}
OUT = ROOT / "02_irc/attempts/step_convergence"
MASS = {"C": 12.011, "H": 1.008, "O": 15.999, "Br": 79.904}
CO, CBR = (0, 2), (0, 1)
H = 2625.4996
TOL_E, TOL_B = 0.2, 0.005
E_TS = -2689.43410109009   # Eh, the TS energy every IRC job here computes first (same input, same TS)


def job(rx, name, d):
    r = RUNS[name]
    wd = OUT / name
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "01_ts/ts/ts.hess", wd / "ts.hess")
    sc = irc_init_scale(ROOT / "01_ts/ts/ts.hess", r["target"])
    sd = SD_PER_S * r["target"]
    blk = (f'%irc\n  MaxIter {r["maxiter"][d]}\n  Direction {d}\n  InitHess read\n  Hess_Filename "ts.hess"\n'
           f'  Scale_Displ_SD {sd:.5f}\n  Adapt_Scale_Displ false\n  Init_Displ length\n'
           f'  Scale_Init_Displ {sc:.5f}\n  PrintLevel 1\nend')
    inp = orca.make_input(rx.keywords("IRC"), read_xyz(ROOT / "01_ts/ts/ts.xyz"), rx.cfg["charge"],
                          rx.cfg["mult"], blocks=rx.solvent_block + "\n" + blk, nprocs=2,
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
    f05 = path(OUT / "s005/irc_f_IRC_F_trj.xyz", OUT / "s005/irc_b_IRC_B_trj.xyz", ts)
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
        ref = region_flags(with_energy(path(OUT / "s010/irc_f_IRC_F_trj.xyz", OUT / "s010/irc_b_IRC_B_trj.xyz", ts),
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


if __name__ == "__main__":
    rx = Reaction.load(ROOT)
    if sys.argv[1] == "run":
        with ThreadPoolExecutor(4) as ex:
            list(ex.map(lambda a: job(rx, *a), [(n, d) for n in RUNS for d in ("forward", "backward")]))
        report()
    elif sys.argv[1] == "early":
        report(early=True)
    else:
        report()
