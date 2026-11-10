"""Builds ts_guess.xyz for cyanide adding to acetone. Run once: python make_ts_guess.py

Atom order: 0 carbonyl C, 1 O, 2 and 3 methyl C, 4 cyanide C, 5 N, 6–11 H.
Cyanide carbon 2.10 Å from the carbonyl carbon, approaching at 107° to the C=O bond (a typical
Bürgi–Dunitz angle) in the plane perpendicular to the acetone plane; C=O stretched to 1.25 Å;
methyls and O bent slightly away from the incoming nucleophile."""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent


def unit(v):
    return v / np.linalg.norm(v)


C0 = np.zeros(3)
O = np.array([1.25, 0.0, -0.15])
M1 = np.array([-0.76, 1.28, -0.25])
M2 = np.array([-0.76, -1.28, -0.25])
a = np.radians(107.0)
d = np.array([np.cos(a), 0.0, np.sin(a)])          # 107° from the C→O direction (+x)
CN_C = 2.10 * d
CN_N = CN_C + 1.18 * d


def methyl_h(m, c):
    """Three H on methyl carbon m bonded to c, staggered, 1.09 Å, tetrahedral."""
    ax = unit(m - c)
    p = unit(np.cross(ax, [0, 0, 1.0]) if abs(ax[2]) < 0.9 else np.cross(ax, [1.0, 0, 0]))
    q = np.cross(ax, p)
    th = np.radians(180 - 109.47)
    return [m + 1.09 * (np.cos(th) * ax + np.sin(th) * (np.cos(phi) * p + np.sin(phi) * q))
            for phi in np.radians([90.0, 210.0, 330.0])]


H = methyl_h(M1, C0) + methyl_h(M2, C0)
xyz = [C0, O, M1, M2, CN_C, CN_N] + H
sym = ["C", "O", "C", "C", "C", "N"] + ["H"] * 6
lines = [str(len(sym)), "CN- + acetone TS guess: C(N)...C 2.10 A, Nu...C=O 107 deg"]
lines += [f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(sym, xyz)]
(HERE / "ts_guess.xyz").write_text("\n".join(lines) + "\n")
