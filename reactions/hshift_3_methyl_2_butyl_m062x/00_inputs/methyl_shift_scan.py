"""Relaxed scan of the degenerate 1,2-methyl shift in the 3-methyl-2-butyl cation (J. Picione, 2026-11-23).

Coordinate q = d(C3–CH₃) − d(C2–CH₃) for the bridging methyl (carbon 4, hydrogens 10, 11, 12). At the anti minimum
q = 1.674 − 1.804 = −0.130 Å; at the C2/C3-swapped anti q = +0.130 Å. Scan from anti, q = −0.13 → +0.13 Å in 0.02 Å steps
(14 points). Each point starts from the previous point's optimized geometry, with the methyl moved rigidly along the
C3→C2 direction until q hits the target, and is then optimized with q held fixed (pysisyphus 1.0.0 RFOptimizer,
redundant internals, DIST_FUNC constraint, thresh "gau") on M06-2X-D3(0)/ma-def2-SVP DefGrid3 VeryTightSCF, SMD
dichloromethane: the working level, one ORCA gradient job per optimization step, each in its own folder.

Methyl rotamer (his rule): the H10–C4–C2–C3 dihedral is logged at every point. A change of more than 40° from the
previous point is flagged as a jump to the other rotamer: the point is rerun from the previous point with a smaller
trust radius; if it still jumps, it is rerun after the pass from the next point's geometry. A flagged point is never
read as part of the shift profile.

At the highest accepted point: Freq (working level) and a DLPNO-CCSD(T)/ma-def2-TZVP single point; the imaginary mode,
if any, is checked for whether it moves the methyl between C3 and C2.

Cores: one 4-core slot is taken through the ordinary queue (pipeline/cores.py) and held for the whole scan, because
each optimization step is a short ORCA job (about 30 s); queuing each behind 20-minute E2 jobs would take days.
Run with the pysisyphus environment: .venv-pysis/Scripts/python.exe -X utf8 reactions/hshift_3_methyl_2_butyl_m062x/
00_inputs/methyl_shift_scan.py run   (check: no ORCA job; tests the constraint with a dummy calculator)
"""
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

RX = Path(__file__).resolve().parents[1]
ROOT = RX.parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import cores, orca                                          # noqa: E402
from pipeline.geometry import read_xyz                                    # noqa: E402
from pysisyphus.calculators.Calculator import Calculator                  # noqa: E402
from pysisyphus.Geometry import Geometry                                  # noqa: E402
from pysisyphus.optimizers.RFOptimizer import RFOptimizer                 # noqa: E402

C1, C2, C3, ME, MEH, HDIH = 0, 1, 3, 4, (10, 11, 12), 10
GROUP = [ME, *MEH]
BOHR = 0.529177210903
EH_KJ = 2625.499639
LEVEL = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX VeryTightSCF DefGrid3 CPCM(dichloromethane)"
DLPNO = "DLPNO-CCSD(T) TightPNO ma-def2-TZVP AutoAux RIJCOSX VeryTightSCF CPCM(dichloromethane)"
SMD = '%cpcm smd true SMDsolvent "DICHLOROMETHANE" end'
NPROCS, MAXCORE = 4, 3000
OUT = RX / "01_stationary/methyl_shift_scan"
ANTI = RX / "01_stationary/secondary_anti/secondary_anti.xyz"
E_ANTI = -196.679490940788
TARGETS = [round(-0.13 + 0.02 * i, 3) for i in range(14)]
PRIM = ["DIST_FUNC", C3, ME, C2, ME, -1]                                  # d(C3–C4) − d(C2–C4)
JUMP_DEG = 40.0


def d(x, i, j):
    return float(np.linalg.norm(x[i] - x[j]))


def q_of(x):
    return d(x, C3, ME) - d(x, C2, ME)


def tors(x, i, j, k, l):
    b0, b1, b2 = x[i] - x[j], x[k] - x[j], x[l] - x[k]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def dih(x):
    return tors(x, HDIH, ME, C2, C3)


def wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def shifted(x, target):
    """Move the methyl (C4 + its H) rigidly along C3→C2 until q = target (Å)."""
    u = x[C2] - x[C3]
    u = u / np.linalg.norm(u)
    lo, hi = -0.6, 0.6
    for _ in range(80):
        mid = (lo + hi) / 2
        y = x.copy()
        y[GROUP] += mid * u
        if q_of(y) < target:            # moving toward C2 lengthens C3–C4 and shortens C2–C4: q increases
            lo = mid
        else:
            hi = mid
    y = x.copy()
    y[GROUP] += (lo + hi) / 2 * u
    return y


def log(msg):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "scan_log.txt", "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
    print(msg, flush=True)


class HeldSlot:
    """ORCA jobs on one 4-core slot held for the whole scan (taken once through the ordinary queue)."""

    def __init__(self):
        self.rid = None
        self.env = os.environ.copy()
        self.env["PATH"] = f"{orca.MSMPI_BIN};{orca.ORCA_DIR};" + self.env.get("PATH", "")

    def take(self):
        self.rid = cores.acquire(NPROCS, NPROCS * MAXCORE, what="hshift_3_methyl_2_butyl_m062x/01_stationary/"
                                 "methyl_shift_scan (whole scan, held)", folder=str(OUT.resolve()))
        log(f"core slot taken ({NPROCS} cores), id {self.rid}")

    def give_back(self):
        if self.rid:
            cores.release(self.rid)
            log("core slot released")
            self.rid = None

    def run(self, wd: Path, name: str, inp: str) -> str:
        wd.mkdir(parents=True, exist_ok=False)                    # one job per clean folder
        (wd / f"{name}.inp").write_text(inp)
        with open(wd / f"{name}.out", "w") as fh:
            subprocess.run([str(orca.ORCA_DIR / "orca.exe"), f"{name}.inp"], cwd=wd, stdout=fh,
                           stderr=subprocess.STDOUT, env=self.env)
        for t in wd.glob("*.tmp*"):
            t.unlink()
        out = (wd / f"{name}.out").read_text(errors="replace")
        if "ORCA TERMINATED NORMALLY" not in out:
            raise orca.OrcaError(f"{wd / name}.out did not end normally (not a result)")
        return out


def xyz_block(atoms, x):
    return "\n".join(f"{a:2s} {p[0]:20.12f} {p[1]:20.12f} {p[2]:20.12f}" for a, p in zip(atoms, x))


class OrcaHeld(Calculator):
    def __init__(self, slot: HeldSlot, wd: Path, **kw):
        super().__init__(charge=1, mult=1, **kw)
        self.slot, self.wd, self.n = slot, wd, 0

    def _engrad(self, atoms, coords):
        wd = self.wd / f"g{self.n:04d}"
        self.n += 1
        c = np.asarray(coords).reshape(-1, 3)
        inp = (f"! {LEVEL} EnGrad Bohrs NoPop\n%pal nprocs {NPROCS} end\n%maxcore {MAXCORE}\n{SMD}\n\n"
               f"* xyz 1 1\n{xyz_block(atoms, c)}\n*\n")
        self.slot.run(wd, "g", inp)
        vals = [ln.strip() for ln in (wd / "g.engrad").read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        nat = int(vals[0])
        return float(vals[1]), np.array([float(v) for v in vals[2:2 + 3 * nat]])

    def get_energy(self, atoms, coords, **kw):
        return {"energy": self._engrad(atoms, coords)[0]}

    def get_forces(self, atoms, coords, **kw):
        E, g = self._engrad(atoms, coords)
        return {"energy": E, "forces": -g}


def optimize(atoms, x_ang, calc, out_dir: Path, trust_max=0.3):
    geom = Geometry(atoms, (x_ang / BOHR).ravel(), coord_type="redund",
                    coord_kwargs={"define_prims": [PRIM], "constrain_prims": [PRIM]})
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau", max_cycles=150, trust_max=trust_max, dump=False, out_dir=str(out_dir))
    opt.run()
    x = np.asarray(geom.coords3d) * BOHR
    return x, float(geom.energy), bool(opt.is_converged), int(opt.cur_cycle) + 1


def point(slot, atoms, start, target, tag, trust_max=0.3):
    wd = OUT / tag
    wd.mkdir(parents=True, exist_ok=False)
    x0 = shifted(start, target)
    calc = OrcaHeld(slot, wd)
    x, E, conv, cyc = optimize(atoms, x0, calc, wd, trust_max)
    (wd / "opt.xyz").write_text(f"{len(atoms)}\nq {q_of(x):+.4f} E {E:.12f}\n{xyz_block(atoms, x)}\n")
    r = {"tag": tag, "target_q": target, "q": round(q_of(x), 4), "C3_CH3": round(d(x, C3, ME), 4),
         "C2_CH3": round(d(x, C2, ME), 4), "dihedral_H10_C4_C2_C3": round(dih(x), 1), "E_eh": E,
         "dE_kJ": round((E - E_ANTI) * EH_KJ, 3), "converged": conv, "cycles": cyc, "gradient_jobs": calc.n,
         "trust_max": trust_max}
    log(f"{tag}: q {r['q']:+.3f} Å, C3–CH₃ {r['C3_CH3']:.3f}, C2–CH₃ {r['C2_CH3']:.3f}, dihedral "
        f"{r['dihedral_H10_C4_C2_C3']:.1f}°, ΔE {r['dE_kJ']:+.2f} kJ/mol, converged {conv} in {cyc} cycles")
    return r, x


def mode_check(out_text, hess: Path, x):
    h = orca.read_hess(hess)
    f = np.asarray(h["freqs"])
    imag = [float(v) for v in f if v < -1]
    res = {"imaginary_cm": imag, "lowest_real_cm": [round(float(v), 1) for v in sorted(v for v in f if v > 1)[:4]]}
    if imag:
        k = int(np.argmin(f))
        v = h["modes"][:, k].reshape(-1, 3)
        amp = (v ** 2).sum(1)
        share = amp / amp.sum()
        y = x + 0.05 * v / np.sqrt(amp.max())
        res.update(share_methyl_H=round(float(share[list(MEH)].sum()), 3), share_methyl_C=round(float(share[ME]), 3),
                   dq_for_0p05_step=round(q_of(y) - q_of(x), 4),
                   moves_methyl_between_C3_and_C2=bool(abs(q_of(y) - q_of(x)) > 0.02))
    return res


def run():
    OUT.mkdir(parents=True, exist_ok=False)
    g = read_xyz(ANTI)
    atoms, x = list(g.symbols), np.asarray(g.coords, dtype=float)
    log(f"start: anti q {q_of(x):+.4f} Å, dihedral {dih(x):.1f}°; targets {TARGETS}")
    slot = HeldSlot()
    slot.take()
    res = {"request": "J. Picione 2026-11-23: relaxed scan of q = d(C3–CH3) − d(C2–CH3), anti → swapped anti",
           "level": LEVEL, "smd": SMD, "points": [], "flagged": []}
    try:
        prev = x
        prev_dih = dih(x)
        for i, t in enumerate(TARGETS):
            r, xi = point(slot, atoms, prev, t, f"p{i:02d}")
            if abs(wrap(r["dihedral_H10_C4_C2_C3"] - prev_dih)) > JUMP_DEG:
                log(f"p{i:02d}: rotamer jump ({prev_dih:.1f}° → {r['dihedral_H10_C4_C2_C3']:.1f}°), rerun from the "
                    f"previous point with trust radius 0.1")
                r["flag"] = "rotamer jump"
                res["flagged"].append(r)
                r, xi = point(slot, atoms, prev, t, f"p{i:02d}_retry", trust_max=0.1)
                if abs(wrap(r["dihedral_H10_C4_C2_C3"] - prev_dih)) > JUMP_DEG:
                    r["flag"] = "rotamer jump after retry"
                    res["flagged"].append(r)
                    log(f"p{i:02d}: still jumps; rerun from the next point after the pass")
            res["points"].append(r)
            prev, prev_dih = xi, r["dihedral_H10_C4_C2_C3"]
            (OUT / "scan.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
        # flagged points that still jumped: rerun from the next point (other neighbour)
        for i, r in enumerate(res["points"]):
            if r.get("flag") == "rotamer jump after retry" and i + 1 < len(res["points"]) and \
                    not res["points"][i + 1].get("flag"):
                nxt = read_xyz(OUT / res["points"][i + 1]["tag"] / "opt.xyz")
                r2, _ = point(slot, atoms, np.asarray(nxt.coords, dtype=float), r["target_q"], f"p{i:02d}_from_next")
                ok = abs(wrap(r2["dihedral_H10_C4_C2_C3"] - res["points"][i + 1]["dihedral_H10_C4_C2_C3"])) <= JUMP_DEG
                r2["flag"] = None if ok else "rotamer jump from both neighbours"
                res["points"][i] = r2
        good = [r for r in res["points"] if not r.get("flag")]
        top = max(good, key=lambda r: r["E_eh"])
        res["highest"] = top["tag"]
        xt = np.asarray(read_xyz(OUT / top["tag"] / "opt.xyz").coords, dtype=float)
        fq = OUT / "highest_freq"
        out = slot.run(fq, "freq", f"! {LEVEL} Freq\n%pal nprocs {NPROCS} end\n%maxcore {MAXCORE}\n{SMD}\n\n"
                                   f"* xyz 1 1\n{xyz_block(atoms, xt)}\n*\n")
        res["highest_freq"] = mode_check(out, fq / "freq.hess", xt)
        sp = OUT / "highest_dlpno"
        out = slot.run(sp, "sp", f"! {DLPNO}\n%pal nprocs {NPROCS} end\n%maxcore {MAXCORE}\n{SMD}\n\n"
                                 f"* xyz 1 1\n{xyz_block(atoms, xt)}\n*\n")
        e_sp = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)[-1])
        t1 = re.findall(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", out)
        anti_sp = json.loads((RX / "01_stationary/secondary_anti.json").read_text(encoding="utf-8"))["sp"]
        res["highest_dlpno"] = {"E_eh": e_sp, "dE_kJ_vs_anti": round((e_sp - anti_sp["DLPNO-CCSD(T)_ma-def2-TZVP"]) * EH_KJ, 3),
                                "T1": float(t1[-1]) if t1 else None}
    finally:
        slot.give_back()
        (OUT / "scan.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    write_md(res)


def write_md(res):
    L = ["# Relaxed scan of the degenerate methyl shift (J. Picione, 2026-11-23)", "",
         "q = d(C3–CH₃) − d(C2–CH₃); anti q = −0.130 Å, swapped anti q = +0.130 Å. Working level M06-2X-D3(0)/ma-def2-SVP, "
         "SMD dichloromethane; ΔE electronic, zero = secondary_anti. Each point optimized with q fixed, started from the "
         "previous point.", "",
         "| Point | q (Å) | C3–CH₃ (Å) | C2–CH₃ (Å) | H10–C4–C2–C3 (°) | ΔE kJ/mol (kcal/mol) | Converged | Flag |",
         "|---|---|---|---|---|---|---|---|"]
    for r in res["points"]:
        L.append(f"| {r['tag']} | {r['q']:+.3f} | {r['C3_CH3']:.3f} | {r['C2_CH3']:.3f} | {r['dihedral_H10_C4_C2_C3']:.1f} | "
                 f"{r['dE_kJ']:+.2f} ({r['dE_kJ'] / 4.184:+.2f}) | {r['converged']} | {r.get('flag') or ''} |")
    if res.get("flagged"):
        L += ["", "Flagged (rotamer jumps, not part of the profile): " +
              "; ".join(f"{r['tag']} dihedral {r['dihedral_H10_C4_C2_C3']}°, ΔE {r['dE_kJ']:+.2f} kJ/mol" for r in res["flagged"])]
    if res.get("highest"):
        L += ["", f"Highest accepted point: {res['highest']}.", f"- Freq: {res.get('highest_freq')}",
              f"- DLPNO-CCSD(T)/ma-def2-TZVP: {res.get('highest_dlpno')}"]
    (OUT / "scan.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def check():
    """No ORCA: the DIST_FUNC constraint holds q while a dummy harmonic calculator relaxes the rest."""
    g = read_xyz(ANTI)
    atoms, x = list(g.symbols), np.asarray(g.coords, dtype=float)
    print("anti q", round(q_of(x), 4), "dihedral", round(dih(x), 1))
    ref = x / BOHR

    class Dummy(Calculator):
        def get_forces(self, atoms, coords, **kw):
            dx = np.asarray(coords) - ref.ravel()
            return {"energy": 0.5 * float(dx @ dx), "forces": -dx}

        def get_energy(self, atoms, coords, **kw):
            return {"energy": self.get_forces(atoms, coords)["energy"]}

    import tempfile
    for t in (-0.09, 0.0, 0.07):
        x0 = shifted(x, t)
        geom = Geometry(atoms, (x0 / BOHR).ravel(), coord_type="redund",
                        coord_kwargs={"define_prims": [PRIM], "constrain_prims": [PRIM]})
        geom.set_calculator(Dummy(charge=1, mult=1))
        with tempfile.TemporaryDirectory() as td:
            opt = RFOptimizer(geom, thresh="gau", max_cycles=100, dump=False, out_dir=td)
            opt.run()
        xs = np.asarray(geom.coords3d) * BOHR
        print(f"target {t:+.3f}: start q {q_of(x0):+.4f}, after opt q {q_of(xs):+.4f}, converged {opt.is_converged}, "
              f"max move from anti {np.abs(xs - x).max():.3f} Å")


if __name__ == "__main__":
    (run if (sys.argv[1:] or ["check"])[0] == "run" else check)()
