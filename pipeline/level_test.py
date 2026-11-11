"""M06-2X geometry work decided by J. Picione, 2026-11-11 (no orbitals, no cards):

  hshift_fallback   2,3-dimethyl-2-butyl (tertiary → tertiary): TS and lowest conformer at
                    M06-2X/ma-def2-SVP, SMD dichloromethane, Opt Freq; single points at
                    B3LYP/ma-def2-TZVP, M06-2X/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP on the
                    B3LYP geometries (DCM and gas) and DLPNO-CCSD(T) on the M06-2X geometries.
  methylbutyl       3-methyl-2-butyl, eclip-side rotation path at the new working level for this
                    reaction, M06-2X/ma-def2-SVP: minima ("eclip", "anti") and TS in DCM and gas,
                    IRC both ways in DCM with a small first step, end-point optimizations;
                    DLPNO-CCSD(T)/ma-def2-TZVP single points on minima and TS at the B3LYP and
                    M06-2X geometries, gas and DCM.
  sensitivity <id>  Functional sensitivity test for a finished reaction: TS and separated
                    reactants re-optimized at M06-2X/ma-def2-SVP (same SMD solvent, Opt Freq),
                    one imaginary mode with the same motion (overlap with the B3LYP mode),
                    DLPNO-CCSD(T)/ma-def2-TZVP single points on the M06-2X geometries.

Integration grid for M06-2X: DefGrid3 (ORCA 6 default is DefGrid2; meta-GGAs such as M06-2X are
grid-sensitive, and DefGrid3 is ORCA's finer grid recommended for them). All M06-2X geometry
work here uses it, including gas-phase re-optimizations of the 3-methyl-2-butyl structures that
were first optimized on DefGrid2 (method check, rotation TS), so every compared geometry is on
the same grid.

Run as a queue: .venv/Scripts/python.exe -X utf8 -m pipeline.level_test queue
A second worker started the same way takes the next task not yet taken (lock files in
docs/level_test_queue/). Each task is restartable: finished ORCA jobs are reused.
"""
from __future__ import annotations

import json
import math
import os
import sys
import traceback
from pathlib import Path

import numpy as np

from . import orca, thermo
from .geometry import Geometry, read_xyz, read_xyz_frames
from .steps import Reaction, _read_json, irc_init_scale
from .units import kj

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions"
M062X = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX TightSCF DefGrid3"
LEVELS_SP = {"B3LYP_ma-def2-TZVP": "B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF",
             "M06-2X_ma-def2-TZVP": "M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF",   # as on the cards
             "DLPNO-CCSD(T)_ma-def2-TZVP": "DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX TightSCF"}
DLPNO = "DLPNO-CCSD(T)_ma-def2-TZVP"
TSB = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 150 end"
STD_STATE_KJ = 7.92            # 1 atm → 1 M, per species (as on the cards)
R_KJ = 8.314462618e-3
QUEUE = ROOT / "docs" / "level_test_queue"


def make(rx, kw, g, q=None, m=None, solv=True, blocks=""):
    if solv:
        kw += f" CPCM({rx.cfg['solvent']})"
        blocks = rx.solvent_block + ("\n" + blocks if blocks else "")
    return orca.make_input(kw, g, rx.cfg["charge"] if q is None else q, rx.cfg["mult"] if m is None else m,
                           blocks=blocks, nprocs=4, maxcore_mb=rx.cfg["maxcore_mb"])


def nimag(stem: Path) -> list[float]:
    h = Path(f"{stem}.hess")
    if not h.exists():
        return []
    f = orca.read_hess(h)["freqs"]
    return [float(x) for x in f[f < 0]]


def opt_min(rx, wd: Path, name: str, g, q=None, m=None, solv=True, note=""):
    """Opt Freq; a minimum left with an imaginary mode is displaced 0.1 Å (largest atom) along it
    and re-optimized once with the same keywords (as in the earlier method check)."""
    r = orca.run(wd, name, make(rx, M062X + " Opt Freq", g, q, m, solv), rx.log, note=note)
    stem = wd / name
    im = nimag(stem)
    if r.ok and im:
        h = orca.read_hess(Path(f"{stem}.hess"))
        v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
        v *= 0.1 / np.linalg.norm(v, axis=1).max()
        g1 = read_xyz(Path(f"{stem}.xyz"))
        g2 = Geometry(g1.symbols, np.asarray(g1.coords) + v)
        r = orca.run(wd, name + "_min", make(rx, M062X + " Opt Freq", g2, q, m, solv), rx.log,
                     note=note + f"; displaced 0.1 Å along the {abs(im[0]):.0f}i mode and re-optimized")
        stem = wd / (name + "_min")
    return stem, r


def opt_ts(rx, wd: Path, name: str, g, solv=True, note=""):
    r = orca.run(wd, name, make(rx, M062X + " OptTS Freq", g, solv=solv, blocks=TSB), rx.log, note=note)
    return wd / name, r


def sp(rx, wd: Path, name: str, level: str, xyz: Path, q=None, m=None, solv=True, note=""):
    r = orca.run(wd, name, make(rx, LEVELS_SP[level], read_xyz(xyz), q, m, solv), rx.log,
                 note=note or f"single point {level}")
    if not r.ok:
        raise orca.OrcaError(f"single point failed: {wd / name}")
    return r.energy_eh


def eg(stem: Path):
    import re
    t = Path(f"{stem}.out").read_text(errors="replace")
    E = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
    G = re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", t)
    return E, (float(G[-1]) if G else None)


def mode_overlap(stem_a: Path, stem_b: Path) -> float:
    """|overlap| of the two normalized Cartesian imaginary modes after superposing geometry b on a
    (same atom order; rotation only)."""
    ha, hb = orca.read_hess(Path(f"{stem_a}.hess")), orca.read_hess(Path(f"{stem_b}.hess"))
    va = ha["modes"][:, int(np.argmin(ha["freqs"]))].reshape(-1, 3)
    vb = hb["modes"][:, int(np.argmin(hb["freqs"]))].reshape(-1, 3)
    xa = np.asarray(read_xyz(Path(f"{stem_a}.xyz")).coords)
    xb = np.asarray(read_xyz(Path(f"{stem_b}.xyz")).coords)
    xa, xb = xa - xa.mean(0), xb - xb.mean(0)
    u, _, vt = np.linalg.svd(xb.T @ xa)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    vb = vb @ rot
    va, vb = va / np.linalg.norm(va), vb / np.linalg.norm(vb)
    return float(abs((va * vb).sum()))


def h_minus_e(stem: Path, T: float) -> float:
    """H(T) − E, kJ/mol, harmonic (as pipeline/gas_phase.thermal_H_kj)."""
    sigma = thermo.symmetry_number(Path(f"{stem}.out"))
    return thermo.thermo(Path(f"{stem}.hess"), T, sigma)["H_minus_E_kJ"]


def g_minus_e(stem: Path, T: float, sigma: int | None = None) -> float:
    sigma = sigma or thermo.symmetry_number(Path(f"{stem}.out"))
    return thermo.thermo(Path(f"{stem}.hess"), T, sigma)["G_minus_E_kJ"]


def write(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=float), encoding="utf-8")


# ------------------------------------------------------------------ tertiary → tertiary shift
def hshift_fallback() -> None:
    rx = Reaction.load(RX / "hshift_2_3_dimethyl_2_butyl")
    base = rx.d("01_ts/level_test")
    B = {"dcm": {"ts": rx.root / "01_ts/ts/ts",
                 "min": rx.root / "02_irc/endpoints/reactant_complex/reactant_complex"},
         "gas": {"ts": rx.root / "01_ts/gas_phase_opt/ts/ts", "min": rx.root / "01_ts/gas_phase_opt/cation/cation"}}
    M = {}
    M["ts"], r = opt_ts(rx, base / "M06-2X/dcm/ts", "ts", read_xyz(Path(f"{B['dcm']['ts']}.xyz")),
                        note="M06-2X/ma-def2-SVP (DefGrid3) TS, SMD DCM, from the B3LYP TS")
    M["min"], r2 = opt_min(rx, base / "M06-2X/dcm/min", "min", read_xyz(Path(f"{B['dcm']['min']}.xyz")),
                           note="M06-2X/ma-def2-SVP (DefGrid3) lowest cation conformer, SMD DCM, from B3LYP")
    out = {"M06-2X_dcm": {k: {"stem": str(v.relative_to(rx.root)), "imag_cm": nimag(v)} for k, v in M.items()},
           "sp": {}}
    for ph, solv in (("dcm", True), ("gas", False)):
        for k, stem in B[ph].items():
            for lev in LEVELS_SP:
                out["sp"].setdefault(f"B3LYP_geom/{ph}", {}).setdefault(lev, {})[k] = sp(
                    rx, base / f"sp/B3LYP_geom/{ph}/{lev}", k, lev, Path(f"{stem}.xyz"), solv=solv)
    for k, stem in M.items():
        for lev in (DLPNO, "M06-2X_ma-def2-TZVP"):
            out["sp"].setdefault("M06-2X_geom/dcm", {}).setdefault(lev, {})[k] = sp(
                rx, base / f"sp/M06-2X_geom/dcm/{lev}", k, lev, Path(f"{stem}.xyz"))
    # 135 K and 298 K barriers, each geometry with its own thermal corrections
    tab = []
    geo = {"B3LYP": B["dcm"], "M06-2X": M}
    for gname, st in geo.items():
        own_E = {k: eg(v)[0] for k, v in st.items()}
        for T in (135.0, 298.15):
            corr = g_minus_e(st["ts"], T, 1) - g_minus_e(st["min"], T, 1)    # σ = 1 for both (low_T.py)
            rows = {f"{gname}/ma-def2-SVP (own level)": own_E}
            for lev, E in out["sp"][f"{gname}_geom/dcm"].items():
                rows[lev] = E
            for lev, E in rows.items():
                dE = kj(E["ts"] - E["min"])
                tab.append({"geometry": gname, "phase": "dcm", "level": lev, "T_K": T, "dE_kJ": dE,
                            "dG_kJ": dE + corr})
    out["barriers"] = tab
    out["geometry"] = {g: _ts_bonds_hshift(Path(f"{st['ts']}.xyz")) for g, st in geo.items()}
    write(base / "summary.json", out)
    # levels for low_T.py (B3LYP geometries, DCM and gas)
    lv = {}
    for lev in LEVELS_SP:
        lv[lev] = {ph: out["sp"][f"B3LYP_geom/{ph}"][lev] for ph in ("dcm", "gas")}
    write(rx.root / "01_ts/low_T_levels.json", lv)


def _ts_bonds_hshift(xyz: Path) -> dict:
    """C–H distances of the migrating H (the H closest to both cationic carbons)."""
    g = read_xyz(xyz)
    x = np.asarray(g.coords)
    C = [i for i, s in enumerate(g.symbols) if s == "C"]
    H = [i for i, s in enumerate(g.symbols) if s == "H"]
    best = min(((h, a, b) for h in H for a in C for b in C if a < b),
               key=lambda t: np.linalg.norm(x[t[0]] - x[t[1]]) + np.linalg.norm(x[t[0]] - x[t[2]]))
    h, a, b = best
    return {"C_H_A": sorted([round(float(np.linalg.norm(x[h] - x[a])), 3), round(float(np.linalg.norm(x[h] - x[b])), 3)]),
            "C_C_A": round(float(np.linalg.norm(x[a] - x[b])), 3)}


# ------------------------------------------------------------------ 3-methyl-2-butyl, eclip side
def _median_step(trj: Path) -> float:
    fr = read_xyz_frames(trj)
    m = np.array([{"C": 12.011, "H": 1.008}[s] for s in fr[0].symbols])
    X = [np.asarray(f.coords) / 0.52917721 for f in fr]
    st = [math.sqrt((m[:, None] * (b - a) ** 2).sum()) for a, b in zip(X[:-1], X[1:])]
    return float(np.median(st))


def methylbutyl() -> None:
    rx = Reaction.load(RX / "hshift_3_methyl_2_butyl")
    C = rx.root / "01_check"
    base = rx.d("01_check/m062x")
    start = {"eclip": C / "method/M06-2X/eclip/eclip_min", "anti": C / "method/M06-2X/anti/anti",
             "ts": C / "rotation/M06-2X/eclip/ts"}
    st = {}
    for ph, solv in (("dcm", True), ("gas", False)):
        for k in ("eclip", "anti"):
            st[(ph, k)], _ = opt_min(rx, base / ph / k, k, read_xyz(Path(f"{start[k]}.xyz")), solv=solv,
                                     note=f"M06-2X/ma-def2-SVP (DefGrid3) {k} minimum, {ph}")
        st[(ph, "ts")], r = opt_ts(rx, base / ph / "ts", "ts", read_xyz(Path(f"{start['ts']}.xyz")), solv=solv,
                                   note=f"M06-2X/ma-def2-SVP (DefGrid3) eclip-side rotation TS, {ph}")
    out = {"structures": {f"{ph}/{k}": {"stem": str(v.relative_to(rx.root)), "imag_cm": nimag(v)}
                          for (ph, k), v in st.items()}}
    # IRC both ways in DCM, first step about the median step of the B3LYP DCM eclip IRC
    ts = st[("dcm", "ts")]
    if len(nimag(ts)) == 1:
        wd = base / "dcm/irc"
        wd.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copy(f"{ts}.hess", wd / "irc_ts.hess")
        s1 = _median_step(C / "rotation/dcm/eclip/irc_IRC_Full_trj.xyz")
        sc = irc_init_scale(Path(f"{ts}.hess"), s1)
        blk = ('%irc\n  MaxIter 300\n  Direction both\n  InitHess read\n  Hess_Filename "irc_ts.hess"\n'
               f'  Init_Displ length\n  Scale_Init_Displ {sc:.5f}\n  PrintLevel 1\nend')
        ri = orca.run(wd, "irc", make(rx, M062X + " IRC", read_xyz(Path(f"{ts}.xyz")), blocks=blk), rx.log,
                      note=f"IRC both ways, M06-2X DCM, first step {s1:.3f} bohr·amu½ (median step of the B3LYP DCM path)")
        out["irc"] = {"ok": ri.ok, "first_step_target_s": s1, "Scale_Init_Displ": sc}
        for tag in ("F", "B"):
            p = wd / f"irc_IRC_{tag}.xyz"
            if p.exists():
                e, _ = opt_min(rx, wd / f"end_{tag}", f"end_{tag}", read_xyz(p),
                               note=f"M06-2X DCM, optimization of the IRC {tag} end (steps appended to the path, marked as optimization steps)")
                out["irc"][f"end_{tag}"] = {"stem": str(e.relative_to(rx.root)), "imag_cm": nimag(e)}
    # DLPNO-CCSD(T) single points: B3LYP and M06-2X geometries, gas and DCM
    B = {"dcm": {"eclip": C / "dcm/eclip/eclip", "anti": C / "dcm/anti/anti", "ts": C / "rotation/dcm/eclip/ts"},
         "gas": {"eclip": C / "gas/eclip/eclip", "anti": C / "gas/anti/anti", "ts": C / "rotation/gas/eclip/ts"}}
    out["sp"] = {}
    for ph, solv in (("dcm", True), ("gas", False)):
        for gname, geo in (("B3LYP", B[ph]), ("M06-2X", {k: st[(ph, k)] for k in ("eclip", "anti", "ts")})):
            for k, stem in geo.items():
                out["sp"].setdefault(f"{gname}/{ph}", {})[k] = sp(
                    rx, base / f"sp/{gname}_geom/{ph}", k, DLPNO, Path(f"{stem}.xyz"), solv=solv)
    rows = []
    for ph in ("dcm", "gas"):
        for gname, geo in (("B3LYP", B[ph]), ("M06-2X", {k: st[(ph, k)] for k in ("eclip", "anti", "ts")})):
            own = {k: eg(v) for k, v in geo.items()}
            for ref in ("eclip", "anti"):
                corr = kj(own["ts"][1] - own["ts"][0]) - kj(own[ref][1] - own[ref][0])
                for lev, E in ((f"{gname}/ma-def2-SVP (own level)", {k: v[0] for k, v in own.items()}),
                               ("DLPNO-CCSD(T)/ma-def2-TZVP", out["sp"][f"{gname}/{ph}"])):
                    dE = kj(E["ts"] - E[ref])
                    rows.append({"phase": ph, "geometry": gname, "from_minimum": ref, "level": lev,
                                 "dE_kJ": dE, "dG298_kJ": dE + corr})
    out["barriers"] = rows
    write(base / "summary.json", out)


# ------------------------------------------------------------------ sensitivity test
MEASURED = {"sn2_hydroxide_bromomethane": ("ΔG‡", 95.0, "water, 25 °C; Mabey and Mill 1978"),
            "e2_ethoxide_2_bromopropane": ("ΔG‡", 105.2, "ethanol, 25 °C, E2 share; Dhar, Hughes and Ingold 1948"),
            "nucadd_cyanide_acetone": ("ΔG‡", 69.6, "water, 25 °C; Schlesinger and Miller 1973"),
            "da_butadiene_ethylene": ("Ea (840 K, gas)", 115.0, "gas phase, 760–921 K")}


def sensitivity(folder: str) -> None:
    rx = Reaction.load(RX / folder)
    base = rx.d("01_ts/level_test/M06-2X")
    frag = _read_json(rx.root / "01_ts/fragments.json")
    from .geometry import formal_charge
    pts = {s["name"]: (rx.root / f"01_ts/fragments/{s['name']}/{frag[s['name']].get('job', s['name'])}.xyz",
                       formal_charge(s["smiles"]), s.get("mult", 1)) for _, s in rx.species() if s["name"] in frag}
    names = [n for n, v in frag.items() if v["role"] in ("reactants", "references")]
    tsb3 = rx.root / "01_ts/ts/ts"
    ts, _ = opt_ts(rx, base / "ts", "ts", read_xyz(Path(f"{tsb3}.xyz")),
                   note="sensitivity test: M06-2X/ma-def2-SVP (DefGrid3) TS from the B3LYP TS")
    st = {"ts": ts}
    for n in names:
        xyz, q, m = pts[n]
        st[n], _ = opt_min(rx, base / n, n, read_xyz(xyz), q, m,
                           note=f"sensitivity test: M06-2X/ma-def2-SVP (DefGrid3) {n}")
    E = {k: sp(rx, base / "sp_dlpno_smd", k, DLPNO, Path(f"{v}.xyz"),
               *((pts[k][1], pts[k][2]) if k != "ts" else (None, None))) for k, v in st.items()}
    out = {"structures": {k: {"stem": str(v.relative_to(rx.root)), "imag_cm": nimag(v)} for k, v in st.items()},
           "dlpno_smd_eh": E, "ts_mode_overlap_with_B3LYP": mode_overlap(tsb3, ts)}
    from .steps import accept_ts_stem
    acc = accept_ts_stem(rx, ts)
    out["ts_acceptance"] = {k: acc[k] for k in ("accepted", "classification", "reasons", "imaginary_cm")}
    if not acc["accepted"]:
        out["FLAG"] = f"ts: not a valid TS ({'; '.join(acc['reasons']) or acc['classification']})"
    if folder.startswith("da_"):
        Eg = {k: sp(rx, base / "sp_dlpno_gas", k, DLPNO, Path(f"{st[k]}.xyz"),
                    *((pts[k][1], pts[k][2]) if k != "ts" else (None, None)), solv=False)
              for k in ("ts", "butadiene_s-trans", "ethylene")}
        out["dlpno_gas_eh"] = Eg
    write(base / "summary.json", out)


# ------------------------------------------------------------------ TS acceptance for queue summaries
def check_ts_acceptance() -> dict:
    """J. Picione, 2026-11-11: every job meant to be a TS must pass the pipeline's own acceptance test
    (steps.accept_ts_stem: exactly one imaginary mode, key bonds moving the expected way, the
    reaction's [ts_mode] rules); summaries flag any that do not. Applied after the fact to every
    level-test summary, since queue worker 1 runs code loaded before this check existed."""
    from .steps import accept_ts_stem
    jobs = {"hshift_2_3_dimethyl_2_butyl": ["01_ts/level_test/M06-2X/dcm/ts/ts"],
            "nucadd_cyanide_acetone": ["01_ts/level_test/M06-2X/ts/ts"],
            "da_butadiene_ethylene": ["01_ts/level_test/M06-2X/ts/ts"],
            "e2_ethoxide_2_bromopropane": ["01_ts/level_test/M06-2X/ts/ts"],
            "sn2_hydroxide_bromomethane_masvp": ["01_ts/level_test/M06-2X/ts/ts"]}
    res = {}
    for rid, stems in jobs.items():
        rx = Reaction.load(RX / rid)
        for s in stems:
            stem = rx.root / s
            if not Path(f"{stem}.hess").exists():
                res[f"{rid}/{s}"] = {"status": "not finished"}
                continue
            a = accept_ts_stem(rx, stem)
            res[f"{rid}/{s}"] = {"accepted": a["accepted"], "classification": a["classification"],
                                 "reasons": a["reasons"], "imaginary_cm": a["imaginary_cm"]}
            summ = stem.parent.parent / "summary.json" if stem.parent.name == "ts" else None
            for cand in (summ, stem.parents[2] / "summary.json"):
                if cand and cand.exists():
                    d = json.loads(cand.read_text(encoding="utf-8"))
                    d.setdefault("ts_acceptance", {})[s] = res[f"{rid}/{s}"]
                    if not a["accepted"]:
                        d["FLAG"] = f"{s}: not a valid TS ({'; '.join(a['reasons']) or a['classification']})"
                    write(cand, d)
                    break
    write(QUEUE / "ts_acceptance.json", res)
    return res


# ------------------------------------------------------------------ bridged ion (tertiary shift)
def _sym_swap(xyz: Path) -> dict:
    """Does the structure map onto itself with C2 ↔ C3 (atoms 0 ↔ 1), H 2 fixed? Atom matching by
    distances to (0, 1, 2), then the best rotation or reflection (Kabsch); RMSD in Å."""
    from scipy.optimize import linear_sum_assignment
    g = read_xyz(xyz)
    x = np.asarray(g.coords)
    sig = np.array([[np.linalg.norm(x[i] - x[a]) for a in (0, 1, 2)] for i in range(len(x))])
    img = sig[:, [1, 0, 2]]
    cost = np.abs(sig[:, None, :] - img[None, :, :]).sum(-1)
    for i in range(len(x)):
        for j in range(len(x)):
            if g.symbols[i] != g.symbols[j]:
                cost[i, j] = 1e6
    _, perm = linear_sum_assignment(cost)
    a, b = x - x.mean(0), x[perm] - x[perm].mean(0)
    out = {}
    u, _, vt = np.linalg.svd(b.T @ a)
    s = np.sign(np.linalg.det(u @ vt))
    for kind, det in (("rotation", 1), ("reflection", -1)):
        rot = u @ np.diag([1, 1, det * s]) @ vt          # det(rot) = det
        out[kind] = float(np.sqrt(((b @ rot - a) ** 2).sum(1).mean()))
    return {"perm": [int(p) for p in perm], "rmsd_A": out}


def _mode_overlaps(stem_ts: Path, stem_min: Path) -> list[tuple[float, float]]:
    """Overlap of the TS imaginary mode with every vibration of the minimum (after superposing the
    minimum on the TS); returns (overlap, frequency cm⁻¹) sorted by overlap."""
    ha, hb = orca.read_hess(Path(f"{stem_ts}.hess")), orca.read_hess(Path(f"{stem_min}.hess"))
    va = ha["modes"][:, int(np.argmin(ha["freqs"]))].reshape(-1, 3)
    xa = np.asarray(read_xyz(Path(f"{stem_ts}.xyz")).coords)
    xb = np.asarray(read_xyz(Path(f"{stem_min}.xyz")).coords)
    xa, xb = xa - xa.mean(0), xb - xb.mean(0)
    u, _, vt = np.linalg.svd(xb.T @ xa)
    rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    va = va / np.linalg.norm(va)
    res = []
    for k, f in enumerate(hb["freqs"]):
        if abs(f) < 1.0:
            continue
        vb = hb["modes"][:, k].reshape(-1, 3) @ rot
        res.append((float(abs((va * vb).sum() / np.linalg.norm(vb))), float(f)))
    return sorted(res, reverse=True)


def _kind(xyz: Path) -> dict:
    g = read_xyz(xyz)
    d0, d1 = g.distance(0, 2), g.distance(1, 2)
    kind = "bridged" if abs(d0 - d1) < 0.10 else ("open, H on C3" if d1 < d0 else "open, H on C2")
    return {"C2–H_A": round(d0, 3), "C3–H_A": round(d1, 3), "C2–C3_A": round(g.distance(0, 1), 3), "kind": kind}


def bridged() -> None:
    """J. Picione, 2026-11-11, items 1–3: M06-2X TS between the open cation and the bridged ion (DCM and
    gas, DefGrid3, Freq), mirror-image check, depth of the bridged well; RI-MP2/ma-def2-SVP gas
    optimization + numerical frequencies of the bridged structure from the M06-2X geometry."""
    from concurrent.futures import ThreadPoolExecutor
    from .steps import accept_ts_stem
    rx = Reaction.load(RX / "hshift_2_3_dimethyl_2_butyl")
    base = rx.d("01_ts/level_test/M06-2X/bridged")
    lt = rx.root / "01_ts/level_test/M06-2X"
    bridged_dcm = lt / "dcm/ts/ts"                 # the OptTS that converged to a minimum

    def mp2():
        wd = base / "rimp2_gas"
        r = orca.run(wd, "bridged", make(rx, "RI-MP2 ma-def2-SVP AutoAux TightSCF Opt NumFreq",
                                         read_xyz(Path(f"{bridged_dcm}.xyz")), solv=False), rx.log,
                     note="RI-MP2/ma-def2-SVP gas, Opt NumFreq of the bridged ion from the M06-2X (DCM) geometry")
        st = wd / "bridged"
        return {"ok": r.ok, "opt_converged": r.opt_converged, "imag_cm": nimag(st), **_kind(Path(f"{st}.xyz")),
                "symmetry": _sym_swap(Path(f"{st}.xyz")), "E_eh": eg(st)[0]}

    def phase(ph):
        solv = ph == "dcm"
        out = {}
        if solv:
            opn, brg = lt / "dcm/min/min", bridged_dcm
        else:
            opn, _ = opt_min(rx, base / "gas/open", "open", read_xyz(rx.root / "01_ts/gas_phase_opt/cation/cation.xyz"),
                             solv=False, note="M06-2X gas, open cation (lowest conformer) from B3LYP gas")
            brg, _ = opt_min(rx, base / "gas/bridged", "bridged", read_xyz(Path(f"{bridged_dcm}.xyz")), solv=False,
                             note="M06-2X gas, bridged ion from the DCM M06-2X bridged minimum")
        out["open"] = {"stem": str(opn.relative_to(rx.root)), "imag_cm": nimag(opn), **_kind(Path(f"{opn}.xyz"))}
        out["bridged"] = {"stem": str(brg.relative_to(rx.root)), "imag_cm": nimag(brg), **_kind(Path(f"{brg}.xyz")),
                          "symmetry": _sym_swap(Path(f"{brg}.xyz"))}
        if not solv:
            out["bridged"]["first_optimization_imag_cm"] = nimag(base / "gas/bridged/bridged")
        if out["bridged"]["kind"] != "bridged":
            out["note"] = "bridged structure did not stay bridged on optimization; no well, no scan"
            return out
        # relaxed scan of C2–H from the bridged ion toward the open cation (H stays on C3)
        r0 = read_xyz(Path(f"{brg}.xyz")).distance(0, 2)
        wd = base / ph / "scan"
        blk = f"%geom\n  Scan B 0 2 = {r0 + 0.02:.3f}, 1.95, 16 end\nend"
        r = orca.run(wd, "scan", make(rx, M062X + " Opt", read_xyz(Path(f"{brg}.xyz")), solv=solv, blocks=blk),
                     rx.log, note=f"M06-2X {ph}, relaxed scan of C2–H {r0 + 0.02:.2f}→1.95 Å from the bridged ion")
        dat = [ln.split() for ln in (wd / "scan.relaxscanact.dat").read_text().splitlines() if ln.strip()]
        sc = [(float(a), float(b)) for a, b in dat]
        out["scan"] = [{"C2–H_A": a, "dE_vs_bridged_kJ": kj(b - eg(brg)[0])} for a, b in sc]
        k = int(np.argmax([b for _, b in sc]))
        out["scan_max_point"] = k + 1
        if k == 0 or k == len(sc) - 1:
            out["note"] = f"scan maximum at the end point {k + 1}; no interior maximum"
            return out
        gmax = read_xyz(wd / f"scan.{k + 1:03d}.xyz")
        ts, r = opt_ts(rx, base / ph / "ts", "ts", gmax, solv=solv,
                       note=f"M06-2X {ph}, TS open ↔ bridged from scan point {k + 1}")
        acc = accept_ts_stem(rx, ts)
        out["ts"] = {"stem": str(ts.relative_to(rx.root)), "imag_cm": nimag(ts), **_kind(Path(f"{ts}.xyz")),
                     "accepted": acc["accepted"], "reasons": acc["reasons"],
                     "symmetry_of_ts": _sym_swap(Path(f"{ts}.xyz"))}
        # which minima does the TS connect? ±0.1 Å along the imaginary mode, then optimize
        if len(nimag(ts)) == 1:
            h = orca.read_hess(Path(f"{ts}.hess"))
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            v *= 0.1 / np.linalg.norm(v, axis=1).max()
            g0 = read_xyz(Path(f"{ts}.xyz"))
            ends = {}
            for sgn, tag in ((1, "plus"), (-1, "minus")):
                g1 = Geometry(g0.symbols, np.asarray(g0.coords) + sgn * v)
                e, _ = opt_min(rx, base / ph / f"ts_{tag}", tag, g1, solv=solv,
                               note=f"M06-2X {ph}, TS displaced {tag} 0.1 Å along its mode, optimized")
                ends[tag] = {**_kind(Path(f"{e}.xyz")), "imag_cm": nimag(e), "dE_vs_open_kJ": kj(eg(e)[0] - eg(opn)[0])}
            out["ts_connects"] = ends
        # energies, ZPE, G(135 K); σ = 1 throughout (as low_T.py)
        T = 135.0
        def th(stem):
            t = thermo.thermo(Path(f"{stem}.hess"), T, 1)
            return eg(stem)[0], t["ZPE_kJ"], t["G_minus_E_kJ"]
        E = {k: th(s) for k, s in (("open", opn), ("bridged", brg), ("ts", ts))}
        def diff(a, b):
            dE = kj(E[a][0] - E[b][0])
            return {"dE_kJ": dE, "dE_ZPE_kJ": dE + E[a][1] - E[b][1], "dG135_kJ": dE + E[a][2] - E[b][2]}
        out["ts_minus_bridged"] = diff("ts", "bridged")
        out["ts_minus_open"] = diff("ts", "open")
        out["bridged_minus_open"] = diff("bridged", "open")
        ov = _mode_overlaps(ts, brg)
        out["transfer_mode_in_bridged"] = {"freq_cm": ov[0][1], "overlap_with_ts_mode": ov[0][0],
                                           "next": ov[1], "ZPE_kJ": 0.5 * ov[0][1] * 0.0119627}
        return out

    with ThreadPoolExecutor(2) as ex:
        fm = ex.submit(mp2)
        res = {"dcm": phase("dcm"), "gas": phase("gas")}
        res["rimp2_gas_bridged"] = fm.result()
    write(base / "summary.json", res)
    print(json.dumps(res, indent=1, default=float, ensure_ascii=False))


# ------------------------------------------------------------------ queue
TASKS = [("hshift_fallback", hshift_fallback, ()),
         ("methylbutyl", methylbutyl, ()),
         ("sens_cyanide", sensitivity, ("nucadd_cyanide_acetone",)),
         ("sens_da", sensitivity, ("da_butadiene_ethylene",)),
         ("sens_e2", sensitivity, ("e2_ethoxide_2_bromopropane",)),
         ("sens_sn2", sensitivity, ("sn2_hydroxide_bromomethane_masvp",))]


def queue() -> None:
    QUEUE.mkdir(parents=True, exist_ok=True)
    for name, fn, args in TASKS:
        lock = QUEUE / f"{name}.lock"
        if (QUEUE / f"{name}.done").exists():
            continue
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
        print(f"=== task {name}", flush=True)
        try:
            fn(*args)
            (QUEUE / f"{name}.done").write_text("ok")
            print(f"=== task {name} done", flush=True)
        except Exception:
            traceback.print_exc()
            (QUEUE / f"{name}.failed").write_text(traceback.format_exc())
            print(f"=== task {name} FAILED", flush=True)
        finally:
            lock.unlink(missing_ok=True)


if __name__ == "__main__":
    if sys.argv[1] == "sensitivity":
        sensitivity(sys.argv[2])
    else:
        {"queue": queue, "hshift_fallback": hshift_fallback, "methylbutyl": methylbutyl,
         "bridged": bridged, "check_ts": check_ts_acceptance}[sys.argv[1]]()
