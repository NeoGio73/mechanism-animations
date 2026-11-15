"""Second half of a degenerate rearrangement as the mirror image of the computed half (J. Picione, plan of
2026-11-12, confirmed 2026-11-14; tertiary hydride shift, 2,3-dimethyl-2-butyl cation in DCM).

The computed path runs open cation (H on C3) → TS → bridged structure. The bridged structure is C₂-symmetric
(the C2 ⇄ C3 swap maps it onto itself to 0.0002 Å), so the second half, bridged → TS′ → open cation (H on C2),
is the image of the computed half under that swap. Each mirror frame is a computed frame with its atoms
relabelled by the swap permutation and rotated by the one rotation that maps the relabelled bridged structure
onto itself; its energy is the computed frame's energy. Mirror frames are labelled as such in irc_path.json
("frame_source") and from there in the viewer data and on the card. They are not IRC points.

The side that leads to the bridged structure ends at its frame closest (RMSD) to the optimized bridged
structure (reaction.toml [mirror_path] reference). If the computed side carries on past that point toward the
other open cation (C3–H − C2–H above drift_A), the frames beyond it are not used as IRC points (they are kept in
true_irc/<side>/path.xyz). The energy span of that side, from the TS to the end point, is reported against the
SCF noise floor (0.034 kJ/mol, docs/protocol.md), because this side is noise-limited.
"""
from __future__ import annotations

import json

import numpy as np

from .geometry import Geometry as PGeom
from .path_angles import HARTREE_KJ, _minima, rmsd

NOISE_FLOOR_KJ = 0.034


def _kabsch(a, b):
    """Rotation R and centres such that (a − ca) R + cb ≈ b."""
    ca, cb = a.mean(0), b.mean(0)
    U, S, Vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, d]) @ Vt, ca, cb


def build(rx, fb, eb, ff, ef, ts, e_ts) -> dict:
    """fb/ff: backward/forward frames from the TS outward (TS excluded), eb/ef their energies (Eh)."""
    mp = rx.cfg["mirror_path"]
    perm = list(mp["perm"])
    ref = next(m for m in _minima(rx.root) if m["name"] == mp["reference"])
    xb = ref["xyz"]
    kb = [b for b in rx.key_bonds(ts.symbols) if b["role"] in ("forms", "breaks")]
    asym = lambda x: (np.linalg.norm(x[kb[0]["atoms"][0]] - x[kb[0]["atoms"][1]])
                      - np.linalg.norm(x[kb[1]["atoms"][0]] - x[kb[1]["atoms"][1]]))
    # the symmetry of the reference under the swap, and the rotation that realizes it
    R, c0, c1 = _kabsch(xb[perm], xb)
    sym_rmsd = rmsd(xb[perm], xb)
    if sym_rmsd > 0.01:
        raise RuntimeError(f"reference {ref['name']} is not symmetric under the swap (RMSD {sym_rmsd:.4f} Å)")
    mirror = lambda x: (x[perm] - c0) @ R + c1

    sides = {"B": (fb, eb), "F": (ff, ef)}
    # The side that leads to the bridged structure is the one that comes closest to it at any frame. (Fixed
    # 2026-11-15: this used the last frame of each side, which fails when the bridged side, on the noise-limited
    # plateau, slides back down to the open cation it started from; both sides then end at the same open cation.)
    d_min = {t: min(rmsd(np.asarray(g.coords), xb) for g in f) for t, (f, _) in sides.items()}
    d_end = {t: rmsd(np.asarray(f[-1].coords), xb) for t, (f, _) in sides.items()}
    tb = min(d_min, key=d_min.get)
    to = "F" if tb == "B" else "B"
    fr, en = sides[tb]
    X = [np.asarray(f.coords) for f in fr]
    d = [rmsd(x, xb) for x in X]
    k = int(np.argmin(d))
    a = [asym(x) for x in X]
    a_ts = asym(np.asarray(ts.coords))
    sign = np.sign(a_ts) or 1.0                          # the TS side of the asymmetry; past the bridge it flips
    drift = [i for i in range(k + 1, len(X)) if -sign * a[i] > mp.get("drift_A", 0.05)]
    # frames past the end point that go back toward the open cation on the TS side (more asymmetric than the TS)
    back = [i for i in range(k + 1, len(X)) if sign * a[i] > sign * a_ts]
    span = (max([e_ts] + list(en[:k + 1])) - min([e_ts] + list(en[:k + 1]))) * HARTREE_KJ
    info = {"reference": ref["name"], "reference_swap_rmsd_A": sym_rmsd, "bridged_side": tb, "open_side": to,
            "bridged_side_frames_computed": len(X), "end_frame_index": k, "end_rmsd_to_reference_A": d[k],
            "end_asymmetry_A": a[k], "end_energy_minus_reference_kJ": (en[k] - ref["E_eh"]) * HARTREE_KJ,
            "drifted_past_reference": bool(drift), "frames_beyond_end_not_used": len(X) - 1 - k,
            "largest_drift_A": float(max((-sign * a[i] for i in drift), default=0.0)),
            "closest_approach_to_reference_A": {t: float(v) for t, v in d_min.items()},
            "last_frame_rmsd_to_reference_A": {t: float(v) for t, v in d_end.items()},
            "returned_toward_start_after_end": bool(back), "first_frame_back_past_TS_asymmetry": (back[0] if back else None),
            "ts_asymmetry_A": float(a_ts),
            "energy_span_TS_to_end_kJ": span, "noise_floor_kJ": NOISE_FLOOR_KJ,
            "span_over_noise": span / NOISE_FLOOR_KJ}
    keep = fr[:k + 1]
    ekeep = list(en[:k + 1])
    other, eother = sides[to]
    # open (computed side, reversed) + TS + computed half to the bridge + mirror back out
    frames = [f.coords for f in other[::-1]] + [ts.coords] + [f.coords for f in keep]
    E = list(eother[::-1]) + [e_ts] + ekeep
    src = ([f"computed IRC ({'backward' if to == 'B' else 'forward'} side)"] * len(other) + ["computed TS"]
           + [f"computed IRC ({'backward' if tb == 'B' else 'forward'} side)"] * len(keep))
    n_comp = len(frames)
    # mirror: the computed half in reverse order, without the bridged end point itself
    comp = [np.asarray(c) for c in frames]
    for j in range(n_comp - 2, -1, -1):
        frames.append(mirror(comp[j]))
        E.append(E[j])
        src.append(f"mirror image (C2⇄C3 swap) of computed frame {j}")
    info.update(n_frames=len(frames), n_computed=n_comp, n_mirror=len(frames) - n_comp, ts_index=len(other),
                mirror_ts_index=len(frames) - 1 - len(other), bridged_index=n_comp - 1)
    geoms = [PGeom(ts.symbols, np.asarray(c), "") for c in frames]
    return {"frames": geoms, "energies": E, "frame_source": src, "info": info}
