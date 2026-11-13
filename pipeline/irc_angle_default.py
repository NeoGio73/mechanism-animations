"""Fairness check, ORCA default IRC settings (J. Picione, 2026-11-12, item 1).

The three-angle check of pipeline/true_irc_report.py on ORCA IRCs run with ORCA's default IRC
settings (Scale_Displ_SD 0.15, Adapt_Scale_Displ true, Init_Displ energy 2 mEh):
  sn2_default   a new run on the SN2 TS, B3LYP-D3(BJ)/ma-def2-SVP, SMD water (the step-test level), every
                IRC keyword at its default except the TS Hessian read and MaxIter
  sn2_pilot     the SN2 pilot (B3LYP-D3(BJ)/def2-SVP, SMD water): defaults except Scale_Displ_SD 0.05
  mb_gas_eclip, mb_gas_anti   3-methyl-2-butyl cation rotation IRCs (B3LYP-D3(BJ)/ma-def2-SVP, gas),
                all IRC keywords at their defaults
Angles at interior frames (up to 15 per side, evenly spaced), central chord Δx through the frame,
gradient from one EnGrad job at the IRC's own level:
  plain angle(Δx, −g); mass-weighted angle(M^1/2 Δx, −M^−1/2 g); Cartesian form angle(Δx, −M^−1 g).
Result: <step_convergence>/default_irc/angles.json and angles.md
Usage: python -m pipeline.irc_angle_default
"""
from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca
from .geometry import read_xyz, read_xyz_frames
from .irc_angle import BOHR, MASS, ang, st

R = st.ROOT.parent
OUT = st.OUT / "default_irc"
MB = R / "hshift_3_methyl_2_butyl/01_check/rotation/gas"
B3_SVP = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF"
SMD = '%cpcm smd true SMDsolvent "WATER" end\n'
PATHS = {
    "sn2_default": dict(dir=OUT / "sn2_default", ts=st.ROOT / "01_ts/ts/ts.xyz",
                        level=B3_SVP + " CPCM(water)", blocks=SMD, q=-1,
                        what="new run, all IRC keywords default"),
    "sn2_pilot": dict(dir=R / "sn2_hydroxide_bromomethane/02_irc", ts=R / "sn2_hydroxide_bromomethane/01_ts/ts/ts.xyz",
                      level="B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF CPCM(water)", blocks=SMD, q=-1,
                      what="SN2 pilot, defaults except Scale_Displ_SD 0.05"),
    "mb_gas_eclip": dict(dir=MB / "eclip", ts=MB / "eclip/ts.xyz", level=B3_SVP, blocks="", q=1,
                         what="3-methyl-2-butyl rotation, all IRC keywords default"),
    "mb_gas_anti": dict(dir=MB / "anti", ts=MB / "anti/ts.xyz", level=B3_SVP, blocks="", q=1,
                        what="3-methyl-2-butyl rotation, all IRC keywords default"),
}
MASS = {**MASS}


def run_default_sn2():
    """ORCA IRC with default settings on the SN2 TS (input = the s010 input without the step keywords)."""
    p = PATHS["sn2_default"]["dir"]
    src = (st.jd("s010", "forward") / "irc_f.inp").read_text()
    keep = [ln for ln in src.splitlines() if not re.match(
        r"\s*(Scale_Displ_SD|Adapt_Scale_Displ|Init_Displ|Scale_Init_Displ|Direction)\b", ln)]
    inp = "\n".join(keep).replace("MaxIter 100", "MaxIter 300\n  Direction both") + "\n"
    p.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy(st.jd("s010", "forward") / "ts.hess", p / "ts.hess")
    r = orca.run(p, "irc", inp, orca.RunLog(OUT / "run_log_default.md"),
                 note="fairness check: ORCA IRC, default IRC settings, SN2 TS, step-test level")
    if not r.ok:
        raise orca.OrcaError("default IRC did not end normally")


def frames(trj: Path, x0):
    fr = [x0] + [np.asarray(f.coords) for f in read_xyz_frames(trj)]
    return [fr[0]] + [c for a, c in zip(fr, fr[1:]) if np.max(np.abs(c - a)) > 1e-7]


def engrad(name, p, side, k, sym, x):
    wd = OUT / "angles" / name / f"{side}_{k:04d}"
    pf, of = wd / "e.property.txt", wd / "e.out"
    if not (pf.exists() and of.exists() and "ORCA TERMINATED NORMALLY" in of.read_text(errors="replace")):
        wd.mkdir(parents=True, exist_ok=True)
        geo = "\n".join(f"{s:2s} {a:16.10f} {b:16.10f} {c:16.10f}" for s, (a, b, c) in zip(sym, x))
        inp = f"! {p['level']} EnGrad NoPop\n%pal nprocs 2 end\n%maxcore 3000\n{p['blocks']}\n* xyz {p['q']} 1\n{geo}\n*\n"
        r = orca.run(wd, "e", inp, orca.RunLog(OUT / "run_log_default.md"), note=f"angle check (defaults): {name} {side} frame {k}")
        if not r.ok:
            raise orca.OrcaError(f"{wd} did not end normally")
        for t in wd.glob("*.tmp*"):
            t.unlink()
    t = pf.read_text()
    gb = re.search(r"\$SCF_Nuc_Gradient(.*?)\$End", t, re.S).group(1)
    return np.array([float(v) for v in re.findall(r"^\s*\d+\s+(-?\d\.\d+e[-+]\d+)\s*$", gb, re.M)])


def main() -> dict:
    run_default_sn2()
    todo = []
    for name, p in PATHS.items():
        ts = read_xyz(p["ts"])
        sym = list(ts.symbols)
        m3 = np.repeat([MASS[a] for a in sym], 3)
        x0 = np.asarray(ts.coords)
        for side, f in (("forward", "irc_IRC_F_trj.xyz"), ("backward", "irc_IRC_B_trj.xyz")):
            X = frames(p["dir"] / f, x0)
            s = np.concatenate([[0.0], np.cumsum([np.linalg.norm(((X[i + 1] - X[i]) / BOHR).ravel() * np.sqrt(m3))
                                                 for i in range(len(X) - 1)])])
            ks = list(range(1, len(X) - 1))
            if len(ks) > 15:
                ks = sorted({ks[int(round(j))] for j in np.linspace(0, len(ks) - 1, 15)})
            for k in ks:
                todo.append(dict(name=name, side=side, k=k, sym=sym, m3=m3, s=float(s[k]), x=X[k],
                                 ch=X[k + 1] - X[k - 1], first_step=float(s[1])))
    with ThreadPoolExecutor(2) as ex:
        G = list(ex.map(lambda d: engrad(d["name"], PATHS[d["name"]], d["side"], d["k"], d["sym"], d["x"]), todo))
    rows = []
    for d, g in zip(todo, G):
        dx = (d["ch"] / BOHR).ravel()
        m3 = d["m3"]
        rows.append({"path": d["name"], "side": d["side"], "frame": d["k"], "s": d["s"],
                     "rms_grad": float(np.sqrt(np.mean(g ** 2))),
                     "plain": ang(dx, -g), "mw": ang(dx * np.sqrt(m3), -g / np.sqrt(m3)), "cart_Minv": ang(dx, -g / m3)})
    summ = {}
    for name in PATHS:
        for zone, f in (("|s|<=1", lambda r: r["s"] <= 1.0), ("|s|>1", lambda r: r["s"] > 1.0), ("all", lambda r: True)):
            rr = [r for r in rows if r["path"] == name and f(r)]
            if rr:
                summ.setdefault(name, {})[zone] = {"n": len(rr), **{k: [float(np.median([r[k] for r in rr])),
                                                                        float(max(r[k] for r in rr))]
                                                                    for k in ("plain", "mw", "cart_Minv")}}
        summ[name]["first_step_s"] = sorted({d["first_step"] for d in todo if d["name"] == name})
    out = {"paths": {k: v["what"] for k, v in PATHS.items()}, "summary": summ, "rows": rows}
    (OUT / "angles.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = ["# Three-angle check on ORCA IRCs with default settings", "",
         "Degrees, median / largest. plain: Δx vs −g; mass-weighted: M^½Δx vs −M^−½g; Cartesian form: Δx vs −M⁻¹g.", "",
         "| Path | Settings | Zone | n | plain | mass-weighted | Cartesian form |", "|---|---|---|---|---|---|---|"]
    for name, z in summ.items():
        for zone in ("|s|<=1", "|s|>1"):
            if zone in z:
                v = z[zone]
                L.append(f"| {name} | {PATHS[name]['what']} | {zone} | {v['n']} | " + " | ".join(
                    f"{v[k][0]:.2f} / {v[k][1]:.2f}" for k in ("plain", "mw", "cart_Minv")) + " |")
    (OUT / "angles.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return out


if __name__ == "__main__":
    main()
