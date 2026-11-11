"""Report on the C2–C3 rotation transition states of the 3-methyl-2-butyl cation (J. Picione,
2026-11-11), from finished jobs only (no ORCA run):

  - imaginary frequency and what the mode does (change of the C2–C3 torsions, the C2–C3–CH₃ angle
    and the key distances for a small step along the mode);
  - barrier from the turned-away minimum the TS was searched from (ΔE and ΔG at 298.15 K, ORCA's
    quasi-RRHO G, 1 M not needed: unimolecular), and from the IRC backward end;
  - the two IRC ends (optimized) with C3–H, C2···H and the stretched C3–CH₃ bond;
  - whether the "anti" and "eclip" starts give the same TS (RMSD after best superposition,
    mirror images and the two C3 methyls swapped allowed);
  - gas phase: the same at M06-2X/ma-def2-SVP, side by side with B3LYP.

Run: .venv/Scripts/python.exe -X utf8 reactions/hshift_3_methyl_2_butyl/00_inputs/rotation_report.py dcm|gas
Writes 01_check/rotation_report_<phase>.md and .json.
"""
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pipeline import orca                          # noqa: E402
from pipeline.geometry import read_xyz              # noqa: E402
from pipeline.units import kj                       # noqa: E402
from followups import H, describe                   # noqa: E402

CHECK = ROOT / "01_check"
MIN = {"B3LYP": {ph: {s: CHECK / ph / s / s for s in ("anti", "eclip")} for ph in ("gas", "dcm")},
       "M06-2X": {"gas": {"anti": CHECK / "method/M06-2X/anti/anti",
                          "eclip": CHECK / "method/M06-2X/eclip/eclip_min"}}}
TS = {"B3LYP": lambda ph, s: CHECK / "rotation" / ph / s / "ts",
      "M06-2X": lambda ph, s: CHECK / "rotation" / "M06-2X" / s / "ts"}


def eg(stem):
    t = Path(f"{stem}.out").read_text(errors="replace")
    E = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
    G = re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", t)
    return E, (float(G[-1]) if G else None)


def coords(stem):
    return np.asarray(read_xyz(Path(f"{stem}.xyz")).coords)


def tors(x, a, b, c, d):
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def ang(x, a, b, c):
    u, v = x[a] - x[b], x[c] - x[b]
    return float(np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))))


def internals(x, me):
    d = lambda i, j: float(np.linalg.norm(x[i] - x[j]))
    return {"C1–C2–C3–H torsion (°)": tors(x, 0, 1, 3, H),
            "C1–C2–C3–CH₃(a) torsion (°)": tors(x, 0, 1, 3, 4),
            "C1–C2–C3–CH₃(b) torsion (°)": tors(x, 0, 1, 3, 5),
            "H(C2)–C2–C3–H torsion (°)": tors(x, 2, 1, 3, H),
            "C2–C3–CH₃(stretched) angle (°)": ang(x, 1, 3, me),
            "C3–CH₃(stretched) (Å)": d(3, me), "C3–H (Å)": d(3, H), "C2···H (Å)": d(1, H),
            "C2–C3 (Å)": d(1, 3)}


def mode_analysis(stem, me):
    """Internal-coordinate change for a 0.05 step along the normalized imaginary mode (both signs,
    half the difference), and the share of the mode's (Cartesian) amplitude on each group."""
    h = orca.read_hess(Path(f"{stem}.hess"))
    i = int(np.argmin(h["freqs"]))
    v = h["modes"][:, i].reshape(-1, 3)
    v = v / np.linalg.norm(v)
    x = coords(stem)
    p, m = internals(x + 0.05 * v, me), internals(x - 0.05 * v, me)
    dlt = {k: (p[k] - m[k]) / 2 for k in p}
    for k in dlt:                                      # torsion wrap-around
        if "torsion" in k and abs(dlt[k]) > 90:
            dlt[k] -= np.sign(dlt[k]) * 180
    amp = (v ** 2).sum(1)
    groups = {"C1 methyl (C1 + its H)": [0, 6, 7, 8], "C2 + its H": [1, 2], "C3 + its H": [3, H],
              "CH₃(a) on C3": [4, 10, 11, 12], "CH₃(b) on C3": [5, 13, 14, 15]}
    share = {g: float(amp[ix].sum()) for g, ix in groups.items()}
    return {"freq_cm": float(h["freqs"][i]), "delta_per_step": dlt, "amplitude_share": share}


def kabsch_rmsd(a, b):
    a, b = a - a.mean(0), b - b.mean(0)
    best = np.inf
    for mirror in (1, -1):
        bb = b * np.array([1, 1, mirror])
        u, _, vt = np.linalg.svd(bb.T @ a)
        dd = np.sign(np.linalg.det(u @ vt))
        r = u @ np.diag([1, 1, dd]) @ vt
        best = min(best, float(np.sqrt(((bb @ r - a) ** 2).sum(1).mean())))
    return best


def same_ts(sa, sb):
    """RMSD over C1, C2, H(C2), C3, both C3 methyl carbons and H(C3); C3 methyls may swap."""
    xa, xb = coords(sa), coords(sb)
    core = [0, 1, 2, 3, 4, 5, H]
    out = []
    for perm in ([0, 1, 2, 3, 4, 5, H], [0, 1, 2, 3, 5, 4, H]):
        out.append(kabsch_rmsd(xa[core], xb[perm]))
    return min(out)


def kjk(x):
    return f"{x:.1f} kJ/mol ({x / 4.184:.1f} kcal/mol)"


def main(phase):
    rot = json.loads((CHECK / "rotation_ts.json").read_text(encoding="utf-8"))
    rep = {"phase": phase, "ts": {}}
    L = [f"# C2–C3 rotation transition states, {'SMD dichloromethane' if phase == 'dcm' else 'gas phase'}", "",
         "Working level B3LYP-D3(BJ)/ma-def2-SVP; G at 298.15 K (ORCA, quasi-RRHO). Atoms: C1 (methyl on C2), "
         "C2 (cationic), C3, the two methyls on C3, H = the hydrogen on C3 that shifts.", ""]
    levels = ["B3LYP"] + (["M06-2X"] if phase == "gas" else [])
    for start in ("anti", "eclip"):
        row = rot.get(f"{phase}/{start}", {})
        d = {"start": start, "levels": {}}
        L.append(f"## Start: turned-away minimum \"{start}\"")
        for lev in levels:
            ts = TS[lev](phase, start)
            if not Path(f"{ts}.out").exists() or not Path(f"{ts}.hess").exists():
                L += [f"- {lev}: TS job not finished.", ""]
                continue
            imag = orca.read_hess(Path(f"{ts}.hess"))["freqs"]
            nimag = int((imag < 0).sum())
            Ets, Gts = eg(ts)
            Em, Gm = eg(MIN[lev][phase][start])
            gts = describe(read_xyz(Path(f"{ts}.xyz")))
            me = gts["stretched_methyl"]
            ma = mode_analysis(ts, me)
            gm = describe(read_xyz(Path(f"{MIN[lev][phase][start]}.xyz")))
            e = {"n_imag": nimag, "imag_cm": [float(f) for f in imag[imag < 0]], "mode": ma,
                 "dE_kJ": kj(Ets - Em), "dG_kJ": kj(Gts - Gm) if Gm and Gts else None,
                 "ts_geom": gts, "min_geom": gm}
            L.append(f"### {lev}")
            L.append(f"- Imaginary modes: {nimag} ({', '.join(f'{abs(f):.0f}i' for f in e['imag_cm'])} cm⁻¹).")
            dl = ma["delta_per_step"]
            L.append("- Mode, change for a small step (0.05, normalized Cartesian): "
                     + "; ".join(f"{k} {v:+.2f}" for k, v in dl.items()) + ".")
            L.append("- Share of the mode amplitude: " + "; ".join(f"{k} {v:.0%}" for k, v in ma["amplitude_share"].items()) + ".")
            L.append(f"- Barrier from the turned-away minimum \"{start}\": ΔE‡ = {kjk(e['dE_kJ'])}, "
                     + (f"ΔG‡(298 K) = {kjk(e['dG_kJ'])}." if e["dG_kJ"] is not None else "ΔG‡ not available."))
            L.append(f"- TS: C1–C2–C3–H {gts['C1_C2_C3_H_deg']:.1f}°, C3–H {gts['C3_H_A']:.3f} Å, C2···H {gts['C2_H_A']:.3f} Å, "
                     f"C3–CH₃ (stretched) {gts['C3_CH3_A']:.3f} Å, C2–C3–CH₃ {gts['C2_C3_CH3_deg']:.1f}°.")
            L.append(f"- Minimum \"{start}\": C1–C2–C3–H {gm['C1_C2_C3_H_deg']:.1f}°, C3–H {gm['C3_H_A']:.3f} Å, "
                     f"C2···H {gm['C2_H_A']:.3f} Å, C3–CH₃ (stretched) {gm['C3_CH3_A']:.3f} Å, C2–C3–CH₃ {gm['C2_C3_CH3_deg']:.1f}°.")
            if lev == "B3LYP":
                ends = {}
                for tag, end in (row.get("irc_ends") or {}).items():
                    st = CHECK / "rotation" / phase / start / f"end_{tag}" / f"end_{tag}"
                    if not Path(f"{st}.out").exists():
                        continue
                    Ee, Ge = eg(st)
                    ge = describe(read_xyz(Path(f"{st}.xyz")))
                    ends[tag] = {"dE_vs_ts_kJ": kj(Ee - Ets), "dE_vs_min_kJ": kj(Ee - Em),
                                 "dG_vs_min_kJ": kj(Ge - Gm) if Ge and Gm else None,
                                 "n_imag": len(end.get("imaginary_cm") or []), "geom": ge,
                                 "irc_end_geom": end.get("irc_end"),
                                 "rmsd_to_start_min_A": same_ts(st, MIN[lev][phase][start])}
                    o = ends[tag]
                    L.append(f"- IRC {'forward' if tag == 'F' else 'backward'} end, optimized: {ge['kind']}; "
                             f"C3–H {ge['C3_H_A']:.3f} Å, C2···H {ge['C2_H_A']:.3f} Å, C3–CH₃ (stretched) {ge['C3_CH3_A']:.3f} Å, "
                             f"C2–C3–CH₃ {ge['C2_C3_CH3_deg']:.1f}°, C1–C2–C3–H {ge['C1_C2_C3_H_deg']:.1f}°; "
                             f"E − TS = {kjk(o['dE_vs_ts_kJ'])}; relative to minimum \"{start}\": ΔE {kjk(o['dE_vs_min_kJ'])}"
                             + (f", ΔG {kjk(o['dG_vs_min_kJ'])}" if o["dG_vs_min_kJ"] is not None else "")
                             + f"; {o['n_imag']} imaginary; RMSD to minimum \"{start}\" {o['rmsd_to_start_min_A']:.3f} Å "
                             f"(core atoms). Unoptimized IRC end: C3–H {o['irc_end_geom']['C3_H_A']:.3f} Å, "
                             f"C2···H {o['irc_end_geom']['C2_H_A']:.3f} Å, C1–C2–C3–H {o['irc_end_geom']['C1_C2_C3_H_deg']:.1f}°.")
                e["irc_ends"] = ends
            d["levels"][lev] = e
            L.append("")
        rep["ts"][start] = d
    L.append("## Same transition state from both starts?")
    for lev in levels:
        a, b = TS[lev](phase, "anti"), TS[lev](phase, "eclip")
        if Path(f"{a}.xyz").exists() and Path(f"{b}.xyz").exists() and Path(f"{a}.out").exists() and Path(f"{b}.out").exists():
            r = same_ts(a, b)
            dE = kj(eg(b)[0] - eg(a)[0])
            rep.setdefault("same_ts", {})[lev] = {"rmsd_A": r, "dE_eclip_minus_anti_kJ": dE}
            L.append(f"- {lev}: RMSD between the two TSs {r:.3f} Å (core atoms, mirror image and methyl swap allowed); "
                     f"E(eclip TS) − E(anti TS) = {kjk(dE)}.")
    out = CHECK / f"rotation_report_{phase}"
    Path(f"{out}.json").write_text(json.dumps(rep, indent=1, default=float), encoding="utf-8")
    Path(f"{out}.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1])
