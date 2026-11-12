"""Is ORCA 6.1.1's IRC path mass-weighted? Evidence from the outputs themselves (J. Picione, 2026-11-12).

1. Initial displacement: angle between the first displacement from the TS (each step-test path, both
   directions) and (a) the lowest eigenvector of the plain Cartesian Hessian, (b) the lowest eigenvector
   of the mass-weighted Hessian M^-1/2 H M^-1/2, taken back to Cartesian displacements (M^-1/2 q),
   (c) the imaginary mode as stored in ts.hess ($normal_modes).
2. Step length: each step's length in plain Cartesian bohr and in bohr·amu½, against the SD step length
   ORCA prints ("Length of steepest descent step ... bohr").
The angle of each path to the plain and to the mass-weighted gradient is in pipeline/irc_angle.py.
Result: <step_convergence>/mass_weight_evidence.json
Usage: python -m pipeline.mass_weight_evidence
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, st, ang
from .orca import read_hess


def read_hess_matrix(path: Path) -> np.ndarray:
    lines = Path(path).read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.strip() == "$hessian")
    n = int(lines[i + 1])
    H = np.zeros((n, n))
    j = i + 2
    while j < len(lines) and not lines[j].startswith("$") and lines[j].strip():
        cols = [int(c) for c in lines[j].split()]
        for r in range(n):
            vals = lines[j + 1 + r].split()[1:]
            for c, v in zip(cols, vals):
                H[r, c] = float(v)
        j += 1 + n
    return H


def main() -> dict:
    hp = st.ROOT / "01_ts/ts/ts.hess"
    h = read_hess(hp)
    H = read_hess_matrix(hp)
    m = np.repeat(h["masses"], 3)
    wc, vc = np.linalg.eigh(H)
    v_cart = vc[:, 0]
    Hm = H / np.sqrt(np.outer(m, m))
    wm, vm = np.linalg.eigh(Hm)
    v_mw = vm[:, 0] / np.sqrt(m)
    v_mw /= np.linalg.norm(v_mw)
    v_file = h["modes"][:, int(np.argmin(h["freqs"]))]
    v_file = v_file / np.linalg.norm(v_file)
    sa = lambda a, b: min(ang(a, b), ang(a, -b))           # sign of a mode is arbitrary
    out = {"hess_modes": {"plain_vs_massweighted_deg": sa(v_cart, v_mw), "file_vs_plain_deg": sa(v_file, v_cart),
                          "file_vs_massweighted_deg": sa(v_file, v_mw)},
           "first_displacement": {}, "steps": {}}
    x0 = np.asarray(read_xyz(st.ROOT / "01_ts/ts/ts.xyz").coords).ravel() / BOHR
    for n in st.R4:
        for d in ("forward", "backward"):
            X = [np.asarray(f.coords).ravel() / BOHR for f in read_xyz_frames(st.trj(n, d))]
            X = [x0] + [x for x in X if np.max(np.abs(x - x0)) > 1e-6]
            X = [X[0]] + [b for a, b in zip(X, X[1:]) if np.max(np.abs(b - a)) > 1e-6]
            d1 = X[1] - X[0]
            out["first_displacement"][f"{n}/{d}"] = {
                "vs_plain_mode_deg": sa(d1, v_cart), "vs_massweighted_mode_deg": sa(d1, v_mw),
                "length_bohr": float(np.linalg.norm(d1)), "length_bohr_amu05": float(np.sqrt((m * d1 ** 2).sum()))}
            dx = [b - a for a, b in zip(X[1:], X[2:])]
            cart = [float(np.linalg.norm(v)) for v in dx]
            mw = [float(np.sqrt((m * v ** 2).sum())) for v in dx]
            o = (st.jd(n, d) / f"irc_{d[0]}.out").read_text(errors="replace")
            printed = float(re.search(r"Length of steepest descent step\s+\.+\s+([\d.]+) bohr", o).group(1))
            out["steps"][f"{n}/{d}"] = {"printed_SD_step_bohr": printed, "median_step_bohr": float(np.median(cart)),
                                        "median_step_bohr_amu05": float(np.median(mw)), "n": len(cart)}
    # what the IRC output prints per iteration and in its summary (units)
    o = (st.jd("s010", "forward") / "irc_f.out").read_text(errors="replace")
    out["printed_columns"] = re.findall(r"^Iteration .*$", o, re.M)[:1] + re.findall(r"^\s*Step\s+E\(Eh\).*$", o, re.M)[:1]
    (st.OUT / "mass_weight_evidence.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
