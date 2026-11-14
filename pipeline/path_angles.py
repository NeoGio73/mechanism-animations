"""Angle check on every production path (J. Picione, 2026-11-13).

For each side of a production true IRC (02_irc/true_irc/<direction>/path.json and path.xyz, written
by pipeline/true_irc_prod.py), at every interior frame: the angle between the step through the frame,
M^1/2 Δx (central chord), and the mass-weighted downhill direction −M^−1/2 g. Reported as median and
90th percentile within |s| ≤ 1 bohr·amu½ and beyond, per side. The rule (the tail rule of 2026-11-13):
a zone fails if its median is 5° or more or its 90th percentile 15° or more; a path with any failing
zone is flagged.

Also, for each end: the nearest optimized M06-2X minimum of the reaction (docs/stationary_m062x.json)
by RMSD after alignment, and the energy of the end above it, so a path that ends in another conformer
(secondary hydride shift) is reported plainly.

Writes 02_irc/true_irc/angles.json and angles.md. Called by true_irc_prod.join; also
Usage: python -m pipeline.path_angles <reaction folder>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BOHR = 0.529177210903
HARTREE_KJ = 2625.499639
MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "Br": 79.904}   # as pysisyphus 1.0.0
MEDIAN_MAX, P90_MAX = 5.0, 15.0
S_NEAR = 1.0
# frames at the end minimum are left out of the check (J. Picione, 2026-11-14): the trailing frames of a
# side whose energy changes by less than END_DE_KJ from the previous frame and whose RMS Cartesian
# gradient is below END_RMS_G; there the downhill direction is not defined
END_DE_KJ, END_RMS_G = 0.005, 5e-5


def _xyz_frames(path: Path) -> list[tuple[list[str], np.ndarray]]:
    L = path.read_text().splitlines()
    out, i = [], 0
    while i < len(L) and L[i].strip():
        n = int(L[i])
        rows = [ln.split() for ln in L[i + 2:i + 2 + n]]
        out.append(([r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])))
        i += n + 2
    return out


def _ang(a, b) -> float:
    c = float(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, c)))))


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD (Å) after centring and optimal rotation (Kabsch), same atom order."""
    A, B = a - a.mean(0), b - b.mean(0)
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(np.mean(np.sum((A @ R - B) ** 2, axis=1))))


def side_angles(X: list[np.ndarray], G: list[np.ndarray], x_ts: np.ndarray, sym: list[str]) -> list[dict]:
    """X: frames from the TS outward (TS excluded), Å; G: gradients (Eh/bohr) at the same frames."""
    sq = np.sqrt(np.repeat([MASS[a] for a in sym], 3))
    P = [x_ts] + X
    rows, s = [], 0.0
    for i in range(1, len(P) - 1):
        a = ((P[i] - P[i - 1]) / BOHR).ravel() * sq
        b = ((P[i + 1] - P[i]) / BOHR).ravel() * sq
        s += float(np.linalg.norm(a))
        g = np.asarray(G[i - 1]).ravel()
        rows.append({"frame_from_ts": i, "s_abs": s, "mw_deg": _ang(a + b, -g / sq)})
    return rows


def _stats(v: list[float]) -> dict | None:
    if not v:
        return None
    a = np.array(v)
    st = {"n": len(v), "median_deg": float(np.median(a)), "p90_deg": float(np.percentile(a, 90)),
          "max_deg": float(a.max())}
    st["passes"] = st["median_deg"] < MEDIAN_MAX and st["p90_deg"] < P90_MAX
    return st


def _minima(rx_root: Path) -> list[dict]:
    js = json.loads((ROOT / "docs/stationary_m062x.json").read_text(encoding="utf-8"))
    out = []
    for v in js.values():
        if v.get("kind") == "min" and v["folder"].split("/")[0] == rx_root.name and v.get("accepted"):
            f = ROOT / "reactions" / v["folder"] / f"{v['job']}.xyz"
            out.append({"name": v["name"], "E_eh": v["E_eh"], "xyz": _xyz_frames(f)[-1][1]})
    return out


def check(rx_root: Path, minima: list[dict] | None = None) -> dict:
    """minima: override for tests (name, E_eh, xyz)."""
    wd = rx_root / "02_irc/true_irc"
    ts_sym, x_ts = _xyz_frames(rx_root / "01_stationary/ts/ts.xyz")[0]
    mins = _minima(rx_root) if minima is None else minima
    out = {"reaction": rx_root.name, "rule": f"each zone: median < {MEDIAN_MAX:g}°, 90th percentile < {P90_MAX:g}°",
           "excluded": f"trailing frames at the end minimum: |dE| < {END_DE_KJ:g} kJ/mol from the previous frame and "
                       f"RMS gradient < {END_RMS_G:g} Eh/bohr (J. Picione, 2026-11-14)",
           "zones": f"abs(s) <= {S_NEAR:g} bohr*amu^0.5 and beyond", "sides": {}}
    for d in ("backward", "forward"):
        js = json.loads((wd / d / "path.json").read_text(encoding="utf-8"))
        fr = _xyz_frames(wd / d / "path.xyz")
        keep = [k for k, (_, x) in enumerate(fr) if np.max(np.abs(x - x_ts)) > 1e-6]
        X = [fr[k][1] for k in keep]
        G = [np.asarray(js["gradients"][k]) for k in keep]
        E = [js["energies"][k] for k in keep]
        rows = side_angles(X, G, x_ts, ts_sym)
        # trailing frames at the end minimum (index k in X; E before X[0] is the TS energy)
        Eprev = [js["E_TS_eh"]] + E[:-1] if "E_TS_eh" in js else [None] + E[:-1]
        at_min = [Eprev[k] is not None and abs(E[k] - Eprev[k]) * HARTREE_KJ < END_DE_KJ
                  and float(np.sqrt(np.mean(np.asarray(G[k]) ** 2))) < END_RMS_G for k in range(len(X))]
        k = len(X) - 1
        while k >= 0 and at_min[k]:
            k -= 1
        end_block = set(range(k + 1, len(X)))
        for r in rows:
            r["at_end_minimum"] = (r["frame_from_ts"] - 1) in end_block
        used = [r for r in rows if not r["at_end_minimum"]]
        near = _stats([r["mw_deg"] for r in used if r["s_abs"] <= S_NEAR])
        far = _stats([r["mw_deg"] for r in used if r["s_abs"] > S_NEAR])
        n_excl = sum(r["at_end_minimum"] for r in rows)
        end = None
        if mins and X:
            cand = sorted(({"minimum": m["name"], "rmsd_A": rmsd(X[-1], m["xyz"]),
                            "end_above_min_kJ": (E[-1] - m["E_eh"]) * HARTREE_KJ} for m in mins),
                          key=lambda c: c["rmsd_A"])
            end = {"nearest": cand[0], "all": cand}
        out["sides"][d] = {"ended": js["ended"], "n_frames": len(X), "n_excluded_at_end_minimum": n_excl,
                           "angles_near": near, "angles_beyond": far,
                           "passes": all(z is None or z["passes"] for z in (near, far)), "end": end,
                           "frames": rows}
    out["passes"] = all(v["passes"] for v in out["sides"].values())
    (wd / "angles.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    z = lambda st: "—" if st is None else f"{st['median_deg']:.2f}° / {st['p90_deg']:.2f}° ({st['n']})"
    L = [f"# Angle check, production path: {rx_root.name}", "",
         f"Mass-weighted angle, M^1/2 Δx (central chord) vs −M^−1/2 g; median / 90th percentile (frames). "
         f"Rule: {out['rule']}. Left out: {out['excluded']}. Path: **{'PASSES' if out['passes'] else 'FAILS (flagged)'}**", "",
         "| Side | ended | frames | left out at the end minimum | abs(s) ≤ 1 | beyond | passes | nearest optimized minimum at the end | RMSD (Å) | end above it (kJ/mol) |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for d, v in out["sides"].items():
        n = v["end"]["nearest"] if v["end"] else None
        tail = f"{n['minimum']} | {n['rmsd_A']:.3f} | {n['end_above_min_kJ']:+.3f} |" if n else "— | — | — |"
        L.append(f"| {d} | {v['ended']} | {v['n_frames']} | {v['n_excluded_at_end_minimum']} | {z(v['angles_near'])} | {z(v['angles_beyond'])} | "
                 f"{'yes' if v['passes'] else '**NO**'} | " + tail)
    (wd / "angles.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return out


if __name__ == "__main__":
    p = Path(sys.argv[1])
    check(p if p.is_absolute() else ROOT / "reactions" / p)
