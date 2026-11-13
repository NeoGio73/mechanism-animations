"""E2 arbitration report (J. Picione, 2026-11-11; parts a, b and c of pipeline.level_test.e2_arbitration).

a. The 7.4 kJ/mol between the DLPNO-CCSD(T) barriers (ΔG‡, 298.15 K, 1 M) on B3LYP and M06-2X
   geometries, split into the electronic part (CC ΔE‡ from separated reactants on each geometry set)
   and the thermal part (G − E of that geometry level's own frequencies).
b. ωB97X-D3(BJ) and RI-MP2 TS (each from the B3LYP and the M06-2X TS) and reactants: key bonds,
   imaginary frequency, acceptance, CC ΔE‡, thermal part, CC ΔG‡.
c. CC single points near the TS: on the B3LYP production path (02_irc, existing) and on the short
   M06-2X IRC (0.10 step); the highest frame and the parabola through it and its two neighbours,
   energies relative to the TS geometry at the same level; s from the TS in bohr·amu½.
No ORCA jobs. Result: 01_ts/level_test/arbitration/report.json and report.md
Usage: python -m pipeline.e2_arbitration_report
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .level_test import RX, STD_STATE_KJ, eg, kj

R = RX / "e2_ethoxide_2_bromopropane"
A = R / "01_ts/level_test/arbitration"
LT = R / "01_ts/level_test/M06-2X"
REAC = ["ethoxide", "2-bromopropane"]
KCAL = 4.184


def both(x):
    return f"{x:.1f} kJ/mol ({x / KCAL:.1f} kcal/mol)"


def peak(s, e, bonds):
    """s sorted ascending, + toward the reactants; bonds: {label: array} at the same frames."""
    i = int(np.argmax(e))
    out = {"frame_s": float(s[i]), "frame_E_kJ": float(e[i]), "parabola": None,
           "frame_bonds_A": {k: float(v[i]) for k, v in bonds.items()}}
    if 0 < i < len(s) - 1:
        c = np.polyfit(s[i - 1:i + 2], e[i - 1:i + 2], 2)
        if c[0] < 0:
            sm = -c[1] / (2 * c[0])
            out["parabola"] = {"s": float(sm), "E_kJ": float(np.polyval(c, sm)),
                               "bonds_A": {k: float(np.interp(sm, s[i - 1:i + 2], v[i - 1:i + 2])) for k, v in bonds.items()},
                               "neighbour_s": [float(s[i - 1]), float(s[i + 1])]}
    return out


def toward_reactants(s, oh):
    """Sign so that s grows toward the reactants (O–H, which forms, is longer there)."""
    return s if np.polyfit(s, oh, 1)[0] > 0 else -s


def main() -> dict:
    res = json.loads((A / "results.json").read_text(encoding="utf-8"))
    sens = json.loads((Path(__file__).resolve().parents[1] / "docs/sensitivity_M06-2X.json")
                      .read_text(encoding="utf-8"))["e2_ethoxide_2_bromopropane"]
    # ---- a
    diag = next(r for r in json.loads((R / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["rows"]
                if r["level"].startswith("DLPNO"))
    b3 = {"dE": diag["dE_TS_vs_separated_kJ"], "dG": diag["dG_TS_vs_separated_1M_kJ"]}
    b3["thermal"] = b3["dG"] - b3["dE"] + STD_STATE_KJ
    S = json.loads((LT / "summary.json").read_text(encoding="utf-8"))
    E = S["dlpno_smd_eh"]
    st = {k: R / v["stem"].replace("\\", "/") for k, v in S["structures"].items()}
    gc = {k: (lambda t: t[1] - t[0])(eg(v)) for k, v in st.items()}
    m6 = {"dE": kj(E["ts"] - sum(E[n] for n in REAC)), "thermal": kj(gc["ts"] - sum(gc[n] for n in REAC))}
    m6["dG"] = m6["dE"] + m6["thermal"] - STD_STATE_KJ
    a = {"B3LYP_geom": b3, "M06-2X_geom": m6, "check_vs_sensitivity_table": sens["CC_dG298_kJ"],
         "difference_B3LYP_minus_M06-2X": {k: b3[k] - m6[k] for k in ("dE", "thermal", "dG")}}
    # ---- c
    rows = list(csv.DictReader(open(R / "02_irc/energies.csv", encoding="utf-8")))
    cc = list(csv.DictReader(open(R / "02_irc/energies_DLPNO-CCSD(T)_ma-def2-TZVP.csv", encoding="utf-8")))
    ts_k = json.loads((R / "02_irc/irc_meta.json").read_text(encoding="utf-8"))["ts_frame"]
    col = {lab: next(k for k in rows[0] if k.startswith("d_" + lab.replace("–", "-") + "_")) for lab in ("O–H", "Cβ–H", "Cα–Br")}
    bB = {lab: np.array([float(r[c]) for r in rows]) for lab, c in col.items()}
    sB = toward_reactants(np.array([float(r["s_bohr_amu05"]) for r in rows]), bB["O–H"])
    eB = np.array([float(r["E_eh"]) for r in cc])
    eBk = kj(eB - eB[ts_k])
    sel = np.array([r["computed"] == "1" for r in cc]) & (np.abs(sB) <= 1.2)
    o = np.argsort(sB[sel])
    pB = peak(sB[sel][o], eBk[sel][o], {k: v[sel][o] for k, v in bB.items()})
    pB["n_points"] = int(sel.sum())
    pB["frames"] = [{"s": float(s), "E_kJ": float(e)} for s, e in zip(sB[sel][o], eBk[sel][o])]
    fr = res["M06-2X_short_IRC"]["frames"]
    bM = {lab: np.array([f["bonds_A"][lab] for f in fr]) for lab in ("O–H", "Cβ–H", "Cα–Br")}
    sM = toward_reactants(np.array([f["s"] for f in fr]), bM["O–H"])
    eM = np.array([f["E_cc_eh"] for f in fr])
    eMk = kj(eM - eM[int(np.argmin(np.abs(sM)))])
    sel = np.abs(sM) <= 1.2
    o = np.argsort(sM[sel])
    pM = peak(sM[sel][o], eMk[sel][o], {k: v[sel][o] for k, v in bM.items()})
    pM["n_points"] = int(sel.sum())
    oa = np.argsort(sM)
    pM["frames"] = [{"s": float(sM[i]), "E_kJ": float(eMk[i]), **{k: float(bM[k][i]) for k in bM}} for i in oa]
    c = {"B3LYP_path": pB, "M06-2X_short_IRC": pM, "M06-2X_step_settings": {"Scale_Displ_SD": res["M06-2X_short_IRC"]["Scale_Displ_SD"]}}
    out = {"a": a, "b": {k: res[k] for k in ("wB97X-D3(BJ)", "RI-MP2")}, "c": c}
    (A / "report.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

    L = ["# E2 arbitration (EtO⁻ + 2-bromopropane, SMD ethanol): parts a, b and c", "",
         "DLPNO-CCSD(T) TightPNO/ma-def2-TZVP SMD ethanol single points throughout. ΔG‡ at 298.15 K, 1 M, from separated reactants.", "",
         "## a. The 7.4 kJ/mol between the CC barriers on B3LYP and M06-2X geometries", "",
         "| Geometry set | CC ΔE‡ (electronic) | thermal part (own frequencies) | 1 atm → 1 M | CC ΔG‡ |", "|---|---|---|---|---|"]
    for k, v in (("B3LYP", b3), ("M06-2X", m6)):
        L.append(f"| {k} | {both(v['dE'])} | {both(v['thermal'])} | −{STD_STATE_KJ:.2f} kJ/mol | {both(v['dG'])} |")
    d = a["difference_B3LYP_minus_M06-2X"]
    L.append(f"| B3LYP − M06-2X | {both(d['dE'])} | {both(d['thermal'])} | 0 | {both(d['dG'])} |")
    L += ["", "## b. ωB97X-D3(BJ) and RI-MP2 TS (ma-def2-SVP, SMD ethanol), CC barrier on each", "",
          "| Level | TS from | O–H (Å) | Cβ–H (Å) | Cα–Br (Å) | imaginary (cm⁻¹) | accepted | CC ΔE‡ | thermal | CC ΔG‡ |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    L.append(f"| B3LYP-D3(BJ) | (production) | {sens['bonds_A']['O–H']['B3LYP']:.3f} | {sens['bonds_A']['Cβ–H']['B3LYP']:.3f} | "
             f"{sens['bonds_A']['Cα–Br']['B3LYP']:.3f} | {abs(sens['imag_cm']['B3LYP'][0]):.0f}i | yes | {both(b3['dE'])} | {both(b3['thermal'])} | {both(b3['dG'])} |")
    L.append(f"| M06-2X | (level test) | {sens['bonds_A']['O–H']['M06-2X']:.3f} | {sens['bonds_A']['Cβ–H']['M06-2X']:.3f} | "
             f"{sens['bonds_A']['Cα–Br']['M06-2X']:.3f} | {abs(sens['imag_cm']['M06-2X'][0]):.0f}i | yes | {both(m6['dE'])} | {both(m6['thermal'])} | {both(m6['dG'])} |")
    for lev in ("wB97X-D3(BJ)", "RI-MP2"):
        v = res[lev]
        if "INCOMPLETE" in v:
            L.append(f"| {lev} | — | INCOMPLETE: {v['INCOMPLETE']} |")
            continue
        for src, t in v["ts"].items():
            b = t["bonds_A"]
            im = ", ".join(f"{abs(x):.0f}i" for x in t["imag_cm"]) or "none"
            L.append(f"| {lev} | {src.replace('from_', '')} | {b['O–H']:.3f} | {b['Cβ–H']:.3f} | {b['Cα–Br']:.3f} | {im} | "
                     f"{'yes' if t['accepted'] else 'NO: ' + '; '.join(t['reasons'])} | {both(t['CC_dE_kJ'])} | {both(t['thermal_kJ'])} | {both(t['CC_dG298_kJ'])} |")
    L += ["", "Reactants at each level: " + "; ".join(
        f"{lev}: " + ", ".join(f"{n} imaginary {len(r['imag_cm'])}" for n, r in res[lev]["reactants"].items())
        for lev in ("wB97X-D3(BJ)", "RI-MP2") if "reactants" in res[lev]) + ".", "",
        "## c. Where CC peaks near the TS (energies relative to the TS geometry at the same level)", "",
        "s in bohr·amu½ from the TS, + toward the reactants on both paths. Peak from the parabola through the highest computed frame and its two neighbours.", "",
        "| Path | points, abs(s) ≤ 1.2 | neighbours used (s) | peak s | peak above TS geometry | O–H at peak (Å) | Cβ–H (Å) | Cα–Br (Å) | at the TS: O–H / Cβ–H / Cα–Br (Å) |",
        "|---|---|---|---|---|---|---|---|---|"]
    tsb = {"B3LYP production path": {k: float(v[ts_k]) for k, v in bB.items()},
           "M06-2X short IRC, 0.10 step": {k: float(v[int(np.argmin(np.abs(sM)))]) for k, v in bM.items()}}
    for k, p in (("B3LYP production path", pB), ("M06-2X short IRC, 0.10 step", pM)):
        pa = p["parabola"]
        if pa is None:
            L.append(f"| {k} | {p['n_points']} | — | {p['frame_s']:+.3f} (frame) | {p['frame_E_kJ']:+.3f} kJ/mol | — | — | — | — |")
            continue
        b, t0 = pa["bonds_A"], tsb[k]
        L.append(f"| {k} | {p['n_points']} | {pa['neighbour_s'][0]:+.3f}, {pa['neighbour_s'][1]:+.3f} | {pa['s']:+.3f} | "
                 f"{pa['E_kJ']:+.2f} kJ/mol ({pa['E_kJ'] / KCAL:+.2f} kcal/mol) | {b['O–H']:.3f} | {b['Cβ–H']:.3f} | {b['Cα–Br']:.3f} | "
                 f"{t0['O–H']:.3f} / {t0['Cβ–H']:.3f} / {t0['Cα–Br']:.3f} |")
    L += ["", "M06-2X short IRC frames (CC relative to the TS geometry):", "",
          "| s | E (kJ/mol) | O–H (Å) | Cβ–H (Å) | Cα–Br (Å) |", "|---|---|---|---|---|"]
    for f in pM["frames"]:
        L.append(f"| {f['s']:+.3f} | {f['E_kJ']:+.3f} | {f['O–H']:.3f} | {f['Cβ–H']:.3f} | {f['Cα–Br']:.3f} |")
    (A / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return out


if __name__ == "__main__":
    main()
