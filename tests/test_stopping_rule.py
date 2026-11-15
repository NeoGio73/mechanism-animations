"""Stopping rule as amended by J. Picione on 2026-11-15 (pipeline.steps.irc_stop_prod).
1. Replays the rule, point by point as the running check does, on every finished production side: each
   side must stop at the frame recorded under the old rule; the secondary hydride shift sides (stopped at
   frame 20 on its flat TS by the old rule) must not stop within their 20 computed points.
2. The safeguard: the end-minimum test must hold on 40 consecutive frames (39 is not enough).
Run: .venv/Scripts/python.exe -X utf8 tests/test_stopping_rule.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.geometry import Geometry as G  # noqa: E402
from pipeline.path_angles import _xyz_frames  # noqa: E402
from pipeline.steps import END_TEST_STEPS, Reaction, irc_stop_prod  # noqa: E402

EXPECT = {("sn2_hydroxide_bromomethane_m062x", "backward"): 153, ("sn2_hydroxide_bromomethane_m062x", "forward"): 188,
          ("nucadd_cyanide_acetone_m062x", "backward"): 158, ("nucadd_cyanide_acetone_m062x", "forward"): 94,
          ("da_butadiene_ethylene_m062x", "backward"): 179,
          ("hshift_2_3_dimethyl_2_butyl_m062x", "backward"): 113, ("hshift_2_3_dimethyl_2_butyl_m062x", "forward"): 129,
          ("hshift_3_methyl_2_butyl_m062x", "backward"): None, ("hshift_3_methyl_2_butyl_m062x", "forward"): None}
FIRST_RUN = {"hshift_3_methyl_2_butyl_m062x": "02_irc/attempts/first_run_early_stop_flat_ts/true_irc"}


def frames_of(rxn, d):
    R = ROOT / "reactions" / rxn
    base = R / FIRST_RUN.get(rxn, "02_irc/true_irc")
    if not base.exists():
        base = R / "02_irc/true_irc"
    sym, xts = _xyz_frames(R / "01_stationary/ts/ts.xyz")[0]
    js = json.loads((base / d / "path.json").read_text(encoding="utf-8"))
    fr = _xyz_frames(base / d / "path.xyz")
    frames = [G(sym, xts, f"TS E {js['E_TS_eh']:.12f}")] + [G(s, x, f"x E {e:.12f}") for (s, x), e in zip(fr, js["energies"])]
    return Reaction.load(R), frames, js["gradients"]


def replay(rx, frames, grads):
    for n in range(2, len(frames) + 1):
        h = irc_stop_prod(rx, frames[:n], grads[:n - 1])
        if h:
            return h["stop_frame"]
    return None


def main():
    for (rxn, d), want in EXPECT.items():
        got = replay(*frames_of(rxn, d))
        assert got == want, f"{rxn} {d}: stop {got}, expected {want}"
        print(f"ok  {rxn} {d}: stop {got}")
    # safeguard: a side resting at a minimum less than 1 kJ/mol below the TS (energy flat, gradient tiny)
    rx, frames, grads = frames_of("hshift_2_3_dimethyl_2_butyl_m062x", "forward")
    ts = frames[0]
    e = float(ts.comment.split(" E ")[-1]) - 0.3 / 2625.499639
    # frame 1 lies 0.3 kJ/mol below the TS (not a test frame); frames 2 .. n+1 are test frames (n in a row)
    for n, want in ((END_TEST_STEPS - 1, None), (END_TEST_STEPS, END_TEST_STEPS + 1)):
        fr = [ts] + [G(ts.symbols, ts.coords, f"x E {e:.12f}") for _ in range(n + 1)]
        h = irc_stop_prod(rx, fr, [[1e-6] * (3 * len(ts.symbols))] * (n + 1))
        assert (h and h["stop_frame"]) == want or (h is None and want is None), (n, h)
        print(f"ok  safeguard with {n} test frames in a row: {'stops at ' + str(h['stop_frame']) if h else 'no stop'}")
    print("all passed")


if __name__ == "__main__":
    main()
