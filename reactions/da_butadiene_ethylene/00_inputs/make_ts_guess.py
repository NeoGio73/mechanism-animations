"""Builds ts_guess.xyz (synchronous Diels–Alder TS guess, mirror plane x = 0) and the two
butadiene starting conformers. Run once: python make_ts_guess.py

Atom order: C1 C2 C3 C4 (diene), C5 C6 (ethylene), then H. New σ bonds: C1–C6 and C4–C5.
Ethylene lies along x at z = 0; the diene termini sit 2.27 Å away (both forming bonds equal),
the diene tilted away from ethylene (envelope, as in the known concerted TS)."""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
CH = 1.09


def unit(v):
    return v / np.linalg.norm(v)


def sp2_pair(c, nb, away, normal):
    """Two H on carbon c: bisector opposite to neighbour nb (tilted along `away`), ±60° in the
    plane perpendicular to `normal`."""
    b = unit(-unit(nb - c) + 0.35 * away)
    q = unit(np.cross(b, normal))
    return [c + CH * (0.5 * b + np.sqrt(3) / 2 * q), c + CH * (0.5 * b - np.sqrt(3) / 2 * q)]


def write(path, sym, xyz, title):
    lines = [str(len(sym)), title] + [f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(sym, xyz)]
    path.write_text("\n".join(lines) + "\n")


# ---- TS guess
d_form, dx = 2.27, 0.70
C5, C6 = np.array([0.69, 0, 0]), np.array([-0.69, 0, 0])
zt = np.sqrt(d_form ** 2 - dx ** 2)
C1, C4 = np.array([-1.39, 0, zt]), np.array([1.39, 0, zt])
C2, C3 = np.array([-0.70, 1.10, zt + 0.46]), np.array([0.70, 1.10, zt + 0.46])
n_diene = unit(np.cross(C2 - C1, C3 - C2))
H = []
H += sp2_pair(C1, C2, unit(C1 - C6), n_diene)
H += [C2 + CH * unit(-(unit(C1 - C2) + unit(C3 - C2)))]
H += [C3 + CH * unit(-(unit(C4 - C3) + unit(C2 - C3)))]
H += sp2_pair(C4, C3, unit(C4 - C5), n_diene)
H += sp2_pair(C5, C6, unit(C5 - C4), np.array([0, 0, 1.0]))
H += sp2_pair(C6, C5, unit(C6 - C1), np.array([0, 0, 1.0]))
xyz = np.array([C1, C2, C3, C4, C5, C6] + H)
write(HERE / "ts_guess.xyz", ["C"] * 6 + ["H"] * 10, xyz,
      "Diels-Alder TS guess: C1-C6 = C4-C5 = 2.27 A, synchronous (mirror plane x = 0)")


# ---- butadiene conformers: C1=C2-C3=C4 with the C1=C2-C3=C4 dihedral set
def butadiene(dihedral_deg):
    c2, c3 = np.array([0, 0, 0.0]), np.array([1.46, 0, 0])
    a = np.radians(180 - 124.0)
    c1 = c2 + 1.34 * np.array([-np.cos(a), np.sin(a), 0])
    t = np.radians(dihedral_deg)
    # C4 from C3: same bond angle, rotated by the dihedral about the C2–C3 axis
    v = 1.34 * np.array([np.cos(a), np.sin(a), 0])
    R = np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    c4 = c3 + R @ v
    n12 = unit(np.cross(c1 - c2, c3 - c2))
    n34 = unit(np.cross(c2 - c3, c4 - c3))
    h = sp2_pair(c1, c2, np.zeros(3), n12)
    h += [c2 + CH * unit(-(unit(c1 - c2) + unit(c3 - c2)))]
    h += [c3 + CH * unit(-(unit(c4 - c3) + unit(c2 - c3)))]
    h += sp2_pair(c4, c3, np.zeros(3), n34)
    return np.array([c1, c2, c3, c4] + h)


# s-cis start slightly off planar (15°) so it relaxes to the cisoid minimum, not the planar saddle
write(HERE / "butadiene_s-cis_start.xyz", ["C"] * 4 + ["H"] * 6, butadiene(15.0),
      "1,3-butadiene, s-cis start (C1=C2-C3=C4 dihedral 15 deg)")
write(HERE / "butadiene_s-trans_start.xyz", ["C"] * 4 + ["H"] * 6, butadiene(180.0),
      "1,3-butadiene, s-trans start (dihedral 180 deg)")
