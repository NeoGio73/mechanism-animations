"""Constant-step IRC test (J. Picione, 2026-11-11), on the SN2 rerun TS (B3LYP-D3(BJ)/ma-def2-SVP).

Target step 0.10 bohr·amu½ for every step including the first (J. Picione, 2026-11-11). First a
calibration (5 forward steps per variant at Scale_Displ_SD 0.05) to find the Scale_Displ_SD that
gives 0.10, then 40 steps forward and 40 backward per variant, first step 0.10 (Init_Displ length):
  A  Adapt_Scale_Displ false
  B  Adapt_Scale_Displ false, SD_ParabolicFit false (no line search, so no step-length choice)
Report: step lengths in s (are they equal, including the first?), and agreement with the
adaptive path in 02_irc (energy at matching s within 0.2 kJ/mol, key bonds within 0.01 Å).
Run: .venv/Scripts/python.exe -X utf8 reactions/sn2_hydroxide_bromomethane_masvp/00_inputs/constant_step_test.py [run|report]
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

TARGET = 0.10      # bohr·amu½, every step including the first (J. Picione, 2026-11-11)
SD0 = 0.05
VARIANTS = {"A": "  Adapt_Scale_Displ false\n",
            "B": "  Adapt_Scale_Displ false\n  SD_ParabolicFit false\n"}
MASS = {"C": 12.011, "H": 1.008, "O": 15.999, "Br": 79.904}
BOND = {"C–O": (0, 2), "C–Br": (0, 1)}


def run_variant(rx, v, calib=False):
    wd = ROOT / "02_irc/attempts/constant_step" / (f"calib_{v}" if calib else v)
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "01_ts/ts/ts.hess", wd / "ts.hess")
    sc = irc_init_scale(ROOT / "01_ts/ts/ts.hess", TARGET)
    sd = SD0 if calib else json.loads((ROOT / "02_irc/attempts/constant_step/calibration.json").read_text())[v]["sd"]
    for d in (("forward",) if calib else ("forward", "backward")):
        blk = (f'%irc\n  MaxIter {5 if calib else 40}\n  Direction {d}\n  InitHess read\n  Hess_Filename "ts.hess"\n'
               f'  Scale_Displ_SD {sd:.5f}\n{VARIANTS[v]}  Init_Displ length\n  Scale_Init_Displ {sc:.5f}\n'
               '  PrintLevel 1\nend')
        inp = orca.make_input(rx.keywords("IRC"), read_xyz(ROOT / "01_ts/ts/ts.xyz"), rx.cfg["charge"],
                              rx.cfg["mult"], blocks=rx.solvent_block + "\n" + blk, nprocs=2,
                              maxcore_mb=rx.cfg["maxcore_mb"])
        orca.run(wd, f"irc_{d[0]}", inp, rx.log,
                 note=f"constant-step IRC test {v} ({VARIANTS[v].strip()}), {d}, "
                      f"{'calibration, 5 steps' if calib else '40 steps'}, Scale_Displ_SD {sd:.5f}, first step {TARGET}")


def calibrate():
    ts = read_xyz(ROOT / "01_ts/ts/ts.xyz")
    out = {}
    for v in VARIANTS:
        fr = frames_s([(1, ROOT / f"02_irc/attempts/constant_step/calib_{v}/irc_f_IRC_F_trj.xyz")], ts)
        st = [t[1] for t in fr][1:]                     # steps after the first (displacement) step
        med = float(np.median(st))
        out[v] = {"sd0": SD0, "steps_s": st, "median_s": med, "sd": SD0 * TARGET / med}
        print(v, "steps", np.round(st, 4).tolist(), "→ Scale_Displ_SD", round(out[v]["sd"], 5), flush=True)
    (ROOT / "02_irc/attempts/constant_step/calibration.json").write_text(json.dumps(out, indent=1))


def frames_s(trjs, ts):
    """Frames ordered backward end → TS → forward end, with s (bohr·amu½) and E (kJ/mol rel. TS)."""
    m = np.array([MASS[s] for s in ts.symbols])
    x0 = np.asarray(ts.coords) / 0.52917721
    out = []
    for sign, trj in trjs:
        prev, s = x0, 0.0
        for f in read_xyz_frames(trj):          # first frame = the displaced point, not the TS
            x = np.asarray(f.coords) / 0.52917721
            step = float(np.sqrt((m[:, None] * (x - prev) ** 2).sum()))
            if step < 1e-6:
                continue
            s += step
            prev = x
            E = float(f.comment.split()[-1])   # 'Coordinates from ORCA-job ... E -2689.43...'
            out.append((sign * s, step, E, f))
    return sorted(out, key=lambda t: t[0])


def report():
    rx = Reaction.load(ROOT)
    ts = read_xyz(ROOT / "01_ts/ts/ts.xyz")
    # adaptive path: the joined IRC (reactant → TS → product), s accumulated frame to frame, s = 0 at the TS
    full = read_xyz_frames(ROOT / "02_irc/irc_IRC_Full_trj.xyz")
    m = np.array([MASS[s] for s in ts.symbols])
    X = [np.asarray(f.coords) / 0.52917721 for f in full]
    cum = np.concatenate([[0.0], np.cumsum([np.sqrt((m[:, None] * (b - a) ** 2).sum()) for a, b in zip(X[:-1], X[1:])])])
    its = int(np.argmin([np.abs(np.asarray(f.coords) - np.asarray(ts.coords)).max() for f in full]))
    rs = cum - cum[its]
    rE = np.array([float(f.comment.split()[-1]) for f in full])
    rd = {b: np.array([np.linalg.norm(np.asarray(f.coords)[i] - np.asarray(f.coords)[j]) for f in full])
          for b, (i, j) in BOND.items()}
    Ets = float(rE[its])
    res = {}
    for v in VARIANTS:
        wd = ROOT / "02_irc/attempts/constant_step" / v
        fr = frames_s([(1, wd / "irc_f_IRC_F_trj.xyz"), (-1, wd / "irc_b_IRC_B_trj.xyz")], ts)
        steps_f = [t[1] for t in fr if t[0] > 0]
        steps_b = [t[1] for t in fr if t[0] < 0][::-1]
        inside = [t for t in fr if rs.min() <= t[0] <= rs.max()]
        dE = [abs((t[2] - Ets) * 2625.4996 - np.interp(t[0], rs, (rE - Ets) * 2625.4996)) for t in inside]
        dB = [max(abs(np.linalg.norm(np.asarray(t[3].coords)[i] - np.asarray(t[3].coords)[j]) - np.interp(t[0], rs, rd[b]))
                  for b, (i, j) in BOND.items()) for t in inside]
        res[v] = {"forward_steps_s": steps_f, "backward_steps_s": steps_b,
                  "first_step_s": [steps_f[0], steps_b[0]],
                  "step_cv_excluding_first": float(np.std(steps_f[1:] + steps_b[1:]) / np.mean(steps_f[1:] + steps_b[1:])),
                  "step_min_max_s": [float(min(steps_f + steps_b)), float(max(steps_f + steps_b))],
                  "max_dE_vs_adaptive_kJ": float(np.nanmax(dE)), "max_dbond_vs_adaptive_A": float(max(dB)),
                  "n_compared": len(inside)}
    (ROOT / "02_irc/attempts/constant_step/report.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    for v, r in res.items():
        print(v, {k: (np.round(x, 3).tolist() if isinstance(x, list) else round(x, 4) if isinstance(x, float) else x)
                  for k, x in r.items() if not k.endswith("steps_s")})
        print("   forward steps", np.round(r["forward_steps_s"], 3).tolist())
        print("   backward steps", np.round(r["backward_steps_s"], 3).tolist())


if __name__ == "__main__":
    rx = Reaction.load(ROOT)
    if sys.argv[1] in ("all", "calibrate"):
        with ThreadPoolExecutor(2) as ex:
            list(ex.map(lambda v: run_variant(rx, v, calib=True), VARIANTS))
        calibrate()
    if sys.argv[1] in ("all", "run"):
        with ThreadPoolExecutor(2) as ex:
            list(ex.map(lambda v: run_variant(rx, v), VARIANTS))
    if sys.argv[1] in ("all", "run", "report"):
        report()
