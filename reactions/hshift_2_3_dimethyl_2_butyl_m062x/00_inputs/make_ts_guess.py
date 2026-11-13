"""Builds ts_guess.xyz for the degenerate 1,2-hydride shift in the 2,3-dimethyl-2-butyl cation,
and the classical cation start structure. Run once: python make_ts_guess.py

Atom order: 0 C2, 1 C3 (the two central carbons), 2 H (the migrating hydrogen), 3 and 4 methyl C
on C2, 5 and 6 methyl C on C3, then the methyl H atoms. TS guess: the hydrogen bridges the C2–C3
bond symmetrically (both C···H 1.30 Å, C2–C3 1.45 Å), the four methyls splayed away from it.
Cation: the same skeleton with H on C3 (C3–H 1.10 Å, sp3 C3) and C2 planar."""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent


def unit(v):
    return v / np.linalg.norm(v)


def methyl_h(c, m, phase):
    """Three staggered H on methyl carbon m bonded to c."""
    u = unit(m - c)
    a = unit(np.cross(u, [0.0, 0.0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0.0, 0.0]))
    b = np.cross(u, a)
    return [m + 1.09 * (0.3338 * u + 0.9426 * (np.cos(p) * a + np.sin(p) * b))
            for p in np.radians([phase, phase + 120, phase + 240])]


def build(cc, h_pos, tilt3):
    C2, C3 = np.array([-cc / 2, 0, 0]), np.array([cc / 2, 0, 0])
    # methyls ~120° from the C–C bond in the xy plane, tilted away from the hydrogen (−z);
    # tilt3 is the extra tilt at C3 (sp3 when the H sits on C3)
    ms = []
    for c, sx, tilt in ((C2, -1, 0.25), (C3, 1, tilt3)):
        for sy in (1, -1):
            ms.append(c + 1.50 * unit(np.array([0.5 * sx, 0.866 * sy, -tilt])))
    sym = ["C", "C", "H", "C", "C", "C", "C"]
    xyz = [C2, C3, h_pos, *ms]
    for k, m in enumerate(ms):
        c = C2 if k < 2 else C3
        for h in methyl_h(c, m, 60.0):
            sym.append("H")
            xyz.append(h)
    return sym, np.array(xyz)


def write(path, sym, xyz, title):
    lines = [str(len(sym)), title] + [f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(sym, xyz)]
    path.write_text("\n".join(lines) + "\n")


cc, ch = 1.45, 1.30
sym, xyz = build(cc, np.array([0.0, 0.0, np.sqrt(ch ** 2 - (cc / 2) ** 2)]), 0.25)
write(HERE / "ts_guess.xyz", sym, xyz, "2,3-dimethyl-2-butyl cation, symmetric H-bridged TS guess (C-H 1.30, C-C 1.45 A)")
cc = 1.50
C3 = np.array([cc / 2, 0, 0])
sym, xyz = build(cc, C3 + 1.10 * unit(np.array([0.35, 0, 1.0])), 0.6)
write(HERE / "cation_start.xyz", sym, xyz, "2,3-dimethyl-2-butyl cation, classical start (H on C3)")
d = lambda i, j: np.linalg.norm(xyz[i] - xyz[j])
print("cation: C2-H", round(d(0, 2), 2), "C3-H", round(d(1, 2), 2), "atoms", len(sym))
