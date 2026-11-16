"""Continuation test (J. Picione, 2026-11-16): continue a side downhill from its last frame, or rerun it from the TS?

Test case: the tertiary hydride shift's open (backward) side, which stopped 4.1 kJ/mol above the open cation.
  rerun        02_irc/true_irc/<side>                    from the TS, rule (c) (flat steps within 0.5 kJ/mol of the
                                                           named end minimum), everything else as in production
  continuation 02_irc/attempts/continuation_test/<side>  the earlier side's frames, then downhill from its last
                                                           frame (exact Hessian there, same step, patches, rule (c))
The new frames of the continuation are compared with the rerun at equal distance along the path: arc length s from
the TS (mass-weighted, bohr·amu½, as in the angle check), the rerun's coordinates and energy interpolated linearly
between its two frames that bracket each s. Heavy-atom RMSD after alignment, and the energy difference.
The earlier frames are compared too, frame by frame (same settings from the TS, so they should match to SCF noise).

Agreement (set by Claude on 2026-11-16 before either run had finished, stated in the 2026-11-16 report):
  - every compared new frame: heavy-atom RMSD <= 0.05 Å and |ΔE| <= 0.1 kJ/mol;
  - the two side ends: each within 0.5 kJ/mol of the end minimum (rule (c)), and the end structures the same
    minimum by the end-check criteria where they apply to two path frames: |ΔE| <= 0.1 kJ/mol, heavy-atom RMSD <= 0.2 Å.
Usage: python -m pipeline.continuation_compare <reaction> <side> [rerun folder] [continuation folder]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from .path_angles import HARTREE_KJ, ROOT, _xyz_frames, rmsd, side_angles

RMSD_A, DE_KJ, END_DE_KJ, END_RMSD_A = 0.05, 0.1, 0.1, 0.2


def _side(base: Path, x_ts, sym):
    js = json.loads((base / "path.json").read_text(encoding="utf-8"))
    fr = _xyz_frames(base / "path.xyz")
    keep = [k for k, (_, x) in enumerate(fr) if np.max(np.abs(x - x_ts)) > 1e-6]
    X = [fr[k][1] for k in keep]
    E = [js["energies"][k] for k in keep]
    G = [np.asarray(js["gradients"][k]) for k in keep]
    s = np.array([r["s_abs"] for r in side_angles(X, G, x_ts, sym)] + [np.nan])
    # side_angles gives s for frames 1 .. n-1; the last frame's s from the last step
    sq = None
    if len(X) >= 2:
        from .path_angles import MASS, BOHR
        sq = np.sqrt(np.repeat([MASS[a] for a in sym], 3))
        s[-1] = s[-2] + float(np.linalg.norm(((X[-1] - X[-2]) / BOHR).ravel() * sq))
    return X, np.array(E), s, js


def compare(rx: str, d: str, rerun_dir: str = "true_irc", cont_dir: str = "attempts/continuation_test") -> dict:
    R = ROOT / "reactions" / rx
    sym, x_ts = _xyz_frames(R / "01_stationary/ts/ts.xyz")[0]
    heavy = [i for i, a in enumerate(sym) if a != "H"]
    Xr, Er, sr, jr = _side(R / "02_irc" / rerun_dir / d, x_ts, sym)
    Xc, Ec, sc, jc = _side(R / "02_irc" / cont_dir / d, x_ts, sym)
    n_old = jc["continuation"]["earlier_points"]          # earlier points incl. the start frame
    em = jc.get("end_min") or jr.get("end_min")
    e_min = em["E_eh"] if em else None
    # earlier frames, frame by frame
    early = [{"frame": k + 1, "heavy_rmsd_A": rmsd(Xc[k][heavy], Xr[k][heavy]), "dE_kJ": (Ec[k] - Er[k]) * HARTREE_KJ}
             for k in range(min(n_old - 1, len(Xr)))]
    rows = []
    for k in range(n_old - 1, len(Xc)):
        s = sc[k]
        j = int(np.searchsorted(sr, s))
        if j == 0 or j >= len(sr):
            rows.append({"frame": k + 1, "s": float(s), "beyond_rerun": True})
            continue
        w = (s - sr[j - 1]) / (sr[j] - sr[j - 1])
        xr = (1 - w) * Xr[j - 1] + w * Xr[j]
        er = (1 - w) * Er[j - 1] + w * Er[j]
        rows.append({"frame": k + 1, "s": float(s), "heavy_rmsd_A": rmsd(Xc[k][heavy], xr[heavy]),
                     "dE_kJ": (Ec[k] - er) * HARTREE_KJ, "E_above_min_kJ": (Ec[k] - e_min) * HARTREE_KJ if e_min else None})
    cmp_rows = [r for r in rows if not r.get("beyond_rerun")]
    ends = {"rerun": {"frames": len(Xr), "s_end": float(sr[-1]), "E_above_min_kJ": (Er[-1] - e_min) * HARTREE_KJ if e_min else None,
                      "ended": jr["ended"], "by": (jr.get("own_stop") or {}).get("by")},
            "continuation": {"frames": len(Xc), "new_frames": len(Xc) - n_old + 1, "s_end": float(sc[-1]),
                             "E_above_min_kJ": (Ec[-1] - e_min) * HARTREE_KJ if e_min else None,
                             "ended": jc["ended"], "by": (jc.get("own_stop") or {}).get("by")},
            "end_to_end": {"heavy_rmsd_A": rmsd(Xc[-1][heavy], Xr[-1][heavy]), "dE_kJ": (Ec[-1] - Er[-1]) * HARTREE_KJ}}
    agree_frames = bool(cmp_rows) and all(r["heavy_rmsd_A"] <= RMSD_A and abs(r["dE_kJ"]) <= DE_KJ for r in cmp_rows)
    agree_ends = (e_min is not None and ends["rerun"]["E_above_min_kJ"] <= 0.5 and ends["continuation"]["E_above_min_kJ"] <= 0.5
                  and abs(ends["end_to_end"]["dE_kJ"]) <= END_DE_KJ and ends["end_to_end"]["heavy_rmsd_A"] <= END_RMSD_A)
    out = {"reaction": rx, "side": d, "rerun_dir": rerun_dir, "continuation_dir": cont_dir, "end_minimum": em,
           "criteria": {"frame_heavy_rmsd_A": RMSD_A, "frame_dE_kJ": DE_KJ, "end_dE_kJ": END_DE_KJ, "end_heavy_rmsd_A": END_RMSD_A},
           "continuation": jc["continuation"], "earlier_frames": {
               "n": len(early), "max_heavy_rmsd_A": max((r["heavy_rmsd_A"] for r in early), default=None),
               "max_abs_dE_kJ": max((abs(r["dE_kJ"]) for r in early), default=None)},
           "new_frames_compared": len(cmp_rows), "new_frames_beyond_rerun": len(rows) - len(cmp_rows),
           "max_heavy_rmsd_A": max((r["heavy_rmsd_A"] for r in cmp_rows), default=None),
           "median_heavy_rmsd_A": float(np.median([r["heavy_rmsd_A"] for r in cmp_rows])) if cmp_rows else None,
           "max_abs_dE_kJ": max((abs(r["dE_kJ"]) for r in cmp_rows), default=None),
           "ends": ends, "agree_frames": agree_frames, "agree_ends": agree_ends, "agree": agree_frames and agree_ends,
           "rows": rows}
    base = R / "02_irc" / cont_dir
    (base / f"compare_{d}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    f = lambda v, p=3: "—" if v is None else f"{v:.{p}f}"
    L = [f"# Continuation test, {rx}, {d} side", "",
         f"Continuation from the last frame of `{jc['continuation']['from']}` vs rerun from the TS (`{rerun_dir}/{d}`), "
         f"compared at equal arc length s (rerun interpolated). Criteria: each new frame heavy-atom RMSD ≤ {RMSD_A} Å and "
         f"|ΔE| ≤ {DE_KJ} kJ/mol; both ends within 0.5 kJ/mol of {em['minimum'] if em else '—'}, end to end |ΔE| ≤ {END_DE_KJ} "
         f"kJ/mol and heavy-atom RMSD ≤ {END_RMSD_A} Å. **Agree: {'yes' if out['agree'] else 'no'}**", "",
         f"- Earlier frames ({len(early)}, same settings from the TS): max heavy-atom RMSD {f(out['earlier_frames']['max_heavy_rmsd_A'], 4)} Å, "
         f"max |ΔE| {f(out['earlier_frames']['max_abs_dE_kJ'], 4)} kJ/mol.",
         f"- Start frame recomputed: {jc['continuation']['start_frame_energy_recomputed_minus_earlier_kJ']:+.4f} kJ/mol.",
         f"- New frames compared: {len(cmp_rows)} (beyond the rerun's end: {out['new_frames_beyond_rerun']}); heavy-atom RMSD "
         f"max {f(out['max_heavy_rmsd_A'], 4)} Å, median {f(out['median_heavy_rmsd_A'], 4)} Å; max |ΔE| {f(out['max_abs_dE_kJ'], 4)} kJ/mol.",
         f"- Ends: rerun {ends['rerun']['frames']} frames, s = {ends['rerun']['s_end']:.2f}, {f(ends['rerun']['E_above_min_kJ'])} kJ/mol above "
         f"the minimum; continuation {ends['continuation']['frames']} frames ({ends['continuation']['new_frames']} new), s = "
         f"{ends['continuation']['s_end']:.2f}, {f(ends['continuation']['E_above_min_kJ'])} kJ/mol above; end to end heavy-atom RMSD "
         f"{ends['end_to_end']['heavy_rmsd_A']:.3f} Å, ΔE {ends['end_to_end']['dE_kJ']:+.3f} kJ/mol.", "",
         "| frame | s (bohr·amu½) | heavy-atom RMSD (Å) | ΔE continuation − rerun (kJ/mol) | E − minimum (kJ/mol) |", "|---|---|---|---|---|"]
    L += [f"| {r['frame']} | {r['s']:.3f} | " + ("beyond the rerun's end | — | — |" if r.get("beyond_rerun") else
          f"{r['heavy_rmsd_A']:.4f} | {r['dE_kJ']:+.4f} | {f(r['E_above_min_kJ'])} |") for r in rows]
    (base / f"compare_{d}.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    o = compare(a[0], a[1], *(a[2:4]))
    print(json.dumps({k: v for k, v in o.items() if k != "rows"}, indent=1))
