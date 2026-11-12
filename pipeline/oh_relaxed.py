"""Gradient-check summary near and far from the TS, and the relaxed O–H length at each matched point
(J. Picione, 2026-11-12). Relaxed O–H = O–H + F/k: F the force stretching O–H at the point (from the
gradient check), k the O–H stretching force constant from the TS Hessian (B3LYP-D3(BJ)/ma-def2-SVP,
SMD water), every other atom held where it is. No new jobs.
Result: <step_convergence>/angle_oh_summary.json
Usage: python -m pipeline.oh_relaxed
"""
from __future__ import annotations

import json

import numpy as np

from .geometry import read_xyz
from .irc_angle import BOHR, st
from .mass_weight_evidence import read_hess_matrix


def main() -> dict:
    a = json.loads((st.OUT / "angle_check.json").read_text(encoding="utf-8"))["rows"]
    oh = json.loads((st.OUT / "oh_check.json").read_text(encoding="utf-8"))["points"]
    Hm = read_hess_matrix(st.ROOT / "01_ts/ts/ts.hess")
    x = np.asarray(read_xyz(st.ROOT / "01_ts/ts/ts.xyz").coords)
    u = (x[3] - x[2]) / np.linalg.norm(x[3] - x[2])
    B = np.zeros(Hm.shape[0])
    B[9:12], B[6:9] = u, -u                          # stretch O(2)–H(3)
    k = float(B @ Hm @ B) / 4                        # Eh/bohr²; B moves H and O by 1 each, so O–H changes by 2
    summ = {}
    for n in st.R4:
        R = [r for r in a if r["path"] == n]
        d = {}
        for zone, sel in (("near_TS_|s|<=1", [r for r in R if abs(r["s"]) <= 1.0]),
                          ("far_|s|>1", [r for r in R if abs(r["s"]) > 1.0]), ("all", R)):
            d[zone] = {"n": len(sel), **{f"{m}_{f}": float(getattr(np, f)([r[m] for r in sel]))
                                          for m in ("cart_deg", "mw_deg") for f in ("median", "max")}}
        summ[n] = d
    pts = []
    for r in oh:
        rel = r["OH_A"] + r["F_OH"] / k * BOHR
        pts.append({"path": r["path"], "side": r["side"], "target_A": r["target_A"], "OH_A": r["OH_A"],
                    "F_OH": r["F_OH"], "OH_relaxed_A": rel, "OH_minus_relaxed_A": r["OH_A"] - rel})
    out = {"k_OH_Eh_per_bohr2": k, "k_OH_mdyn_per_A": k * 15.569141, "angles": summ, "oh_points": pts}
    (st.OUT / "angle_oh_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    o = main()
    print("k_OH", round(o["k_OH_Eh_per_bohr2"], 4), "Eh/bohr2 =", round(o["k_OH_mdyn_per_A"], 2), "mdyn/A")
    for n, d in o["angles"].items():
        print(n, " | ".join(f"{z}: n={v['n']} cart {v['cart_deg_median']:.1f}/{v['cart_deg_max']:.1f} mw {v['mw_deg_median']:.1f}/{v['mw_deg_max']:.1f}" for z, v in d.items()))
    for side in ("product", "reactant"):
        for t in sorted({p["target_A"] for p in o["oh_points"] if p["side"] == side}):
            row = [p for p in o["oh_points"] if p["side"] == side and p["target_A"] == t]
            print(side, t, " ".join(f"{p['path']}:{p['OH_A']:.4f}->{p['OH_relaxed_A']:.4f}" for p in row))
