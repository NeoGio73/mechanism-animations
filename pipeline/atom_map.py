"""Atom mapping between a structure and a reference minimum whose atoms are numbered differently (2026-11-16).

Needed for the tertiary hydride shift (mirror path): its open-side end and the optimized open cation are the same
cation with C2 and C3 (and their methyls) numbered the other way, so a plain RMSD matches the end to the bridged
structure instead (end_check_open_side.md, 2026-11-15, done by hand there). Here: every relabelling of the heavy atoms
that keeps elements and the heavy-atom bond graph is tried, the one with the lowest heavy-atom RMSD (after alignment)
is kept, and each hydrogen is then matched within the heavy atom it is bonded to (closest after alignment).
best_map returns P with x_ref[P] numbered as x.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import linear_sum_assignment

from .path_minimum import bonds


def _kabsch(a: np.ndarray, b: np.ndarray):
    """Rotation and centres that put b onto a."""
    ca, cb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((b - cb).T @ (a - ca))
    d = np.sign(np.linalg.det(u @ vt))
    R = u @ np.diag([1, 1, d]) @ vt
    return R, ca, cb


def best_map(sym, x, x_ref, max_heavy: int = 9) -> tuple[np.ndarray, float]:
    x, x_ref = np.asarray(x, float), np.asarray(x_ref, float)
    heavy = [i for i, s in enumerate(sym) if s != "H"]
    if len(heavy) > max_heavy:
        raise ValueError(f"{len(heavy)} heavy atoms: too many for the permutation search")
    hg = {p for p in bonds(sym, x) if p[0] in heavy and p[1] in heavy}
    hg_ref = {p for p in bonds(sym, x_ref) if p[0] in heavy and p[1] in heavy}
    best, best_r = None, 9e9
    for perm in itertools.permutations(heavy):
        m = dict(zip(heavy, perm))
        if any(sym[i] != sym[m[i]] for i in heavy):
            continue
        if {tuple(sorted((m[i], m[j]))) for i, j in hg} != hg_ref:
            continue
        a, b = x[heavy], x_ref[[m[i] for i in heavy]]
        R, ca, cb = _kabsch(a, b)
        r = float(np.sqrt((((b - cb) @ R + ca - a) ** 2).sum(1).mean()))
        if r < best_r:
            best, best_r = m, r
    if best is None:
        raise ValueError("no relabelling keeps the heavy-atom bond graph")
    R, ca, cb = _kabsch(x[heavy], x_ref[[best[i] for i in heavy]])
    y = (x_ref - cb) @ R + ca                                   # reference aligned onto x
    hyd = [i for i, s in enumerate(sym) if s == "H"]
    owner = lambda xx, h: min(heavy, key=lambda k: np.linalg.norm(xx[h] - xx[k]))
    P = np.arange(len(sym))
    for i in heavy:
        P[i] = best[i]
    for k in heavy:
        hx = [h for h in hyd if owner(x, h) == k]
        hr = [h for h in hyd if owner(x_ref, h) == best[k]]
        if len(hx) != len(hr):
            raise ValueError(f"hydrogens on atom {k + 1}: {len(hx)} here, {len(hr)} on its partner in the reference")
        if hx:
            cost = np.linalg.norm(x[hx][:, None, :] - y[hr][None, :, :], axis=2)
            ri, ci = linear_sum_assignment(cost)
            for a_, b_ in zip(ri, ci):
                P[hx[a_]] = hr[b_]
    return P, best_r
