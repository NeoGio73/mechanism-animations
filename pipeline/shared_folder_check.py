"""Check 1 (J. Picione, 2026-11-11): did any job that shared its folder with another running ORCA job
use a number that was not its own?

ORCA writes text output of sub-steps (the D3 printout among them) to a fixed file, redirect.tmp, in the
job folder (orca_tools/qcmsg.cpp). Two jobs in one folder share that file. The check recomputes, alone
in a clean folder, what each job used, and compares:

Part A, jobs whose per-evaluation property file survives (crashed jobs; a snapshot of the running
s005 forward rerun). That file holds, for every energy/gradient evaluation, the exact geometry (Bohr,
12 decimals), the D3 energy used, the final energy and, for gradient evaluations, the gradient.
  A1  D3 alone at every evaluation's geometry (B3LYP-D3(BJ)/STO-3G; D3 depends only on geometry and
      functional) vs the D3 value the job used.                       tolerance 1e-9 Eh
  A2  energy and gradient at the job's own level, at every gradient evaluation that overlapped another
      job in the folder, vs the values the job used.                  1e-7 Eh, 1e-5 Eh/bohr
Part B, finished jobs (ORCA deletes the per-evaluation file at normal termination; only the
trajectory is left, coordinates to 6 decimals in Å):
  B   energy and gradient at the job's own level at every trajectory frame vs the trajectory energy.
      tolerance per frame 1e-7 Eh + sum|g_i| × 0.5e-6 Å (the coordinate rounding)
      Also D3 alone at every frame, to give the size an exchanged D3 value would have had.

Every check job runs alone in its own folder (one ORCA job per folder) through the core ledger.
Results: <step_convergence>/check1/results.json; gate for the verdict: <step_convergence>/check1.json.
Usage: python -m pipeline.shared_folder_check partA | partB <job> | gate
"""
from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reactions/sn2_hydroxide_bromomethane_masvp/02_irc/attempts/step_convergence"
CK = OUT / "check1"
LEVEL = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF CPCM(water)"
SMD = '%cpcm smd true SMDsolvent "WATER" end'
CHARGE, MULT = -1, 1
BOHR = 0.529177210903
TOL_D3, TOL_E, TOL_G = 1e-9, 1e-7, 1e-5
ROUND_A = 0.5e-6                      # Å, half the last printed digit of the trajectory coordinates
LOG = orca.RunLog(CK / "run_log_check1.md") if CK.exists() else None

# Part A sources: (tag, property file, which gradient evaluations overlapped another job in the folder)
# "all": the job shared the folder for its whole run; ("last", n) / ("first", n): only at its end/start.
PART_A = [
    ("s010_forward_crashed1607", OUT / "s010/crashed_1706_forward/irc_f_IRC_F.property.txt", "all"),
    ("s005_forward_crashed1607", OUT / "s005/crashed_1706_forward/irc_f_IRC_F.property.txt", "all"),
    ("s005_backward_crashed1607", OUT / "s005/crashed_1706_backward/irc_b_IRC_B.property.txt", "all"),
    ("s005_forward_crashed1722", OUT / "s005/crashed_1722_forward/irc_f_IRC_F.property.txt", "all"),
    ("s005_backward_crashed1941", OUT / "s005/crashed_1941_backward/irc_b_IRC_B.property.txt", ("last", 6)),
    ("s005_forward_rerun_first_steps", CK / "snapshots/s005_forward_rerun_IRC_F.property.txt", ("first", 6)),
]
# Part B: finished jobs that shared a folder (check 1a), with the trajectory and the partner's record
PART_B = {
    "s010/backward": {"trj": OUT / "s010/irc_b_IRC_B_trj.xyz",
                      "partner": "s010_forward_crashed1607", "partner_slice": [0, None], "own_slice": [0, None]},
    "s005/forward": {"trj": OUT / "s005/irc_f_IRC_F_trj.xyz",
                     "partner": "s005_backward_crashed1941", "partner_slice": [-12, None], "own_slice": [0, 6]},
}


def read_property(p: Path) -> list[dict]:
    """Per GeometryIndex: symbols, coords (Bohr), vdw, final energy, gradient (21 values or None)."""
    t = p.read_text(errors="replace")
    ev: dict[int, dict] = {}
    for blk in re.finditer(r"^\$(\w+)\n(.*?)^\$End", t, re.S | re.M):
        name, body = blk.group(1), blk.group(2)
        gi = re.search(r"&GeometryIndex (\d+)", body)
        if not gi:
            continue
        e = ev.setdefault(int(gi.group(1)), {"idx": int(gi.group(1)), "grad": None})
        if name == "Geometry":
            rows = re.findall(r"^\s+([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", body, re.M)
            e["symbols"] = [r[0] for r in rows]
            e["bohr"] = np.array([[float(x) for x in r[1:]] for r in rows])
        elif name == "VdW_Correction":
            e["vdw"] = float(re.search(r"&vdW \[&Type \"Double\"\]\s+(\S+)", body).group(1))
        elif name == "Single_Point_Data":
            m = re.search(r"&FinalEnergy \[&Type \"Double\"\]\s+(\S+)", body)
            if m:
                e["final"] = float(m.group(1))
        elif name == "SCF_Nuc_Gradient":
            vals = re.findall(r"^\s*\d+\s+(-?\d\.\d+e[-+]\d+)\s*$", body, re.M)
            e["grad"] = np.array([float(v) for v in vals])
    return [ev[k] for k in sorted(ev) if "bohr" in ev[k] and "vdw" in ev[k]]


def read_trj(p: Path) -> list[dict]:
    lines = p.read_text().splitlines()
    out, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        e = float(re.search(r"E\s+(-?\d+\.\d+)", lines[i + 1]).group(1))
        rows = [ln.split() for ln in lines[i + 2:i + 2 + n]]
        out.append({"frame": len(out), "symbols": [r[0] for r in rows],
                    "bohr": np.array([[float(x) for x in r[1:4]] for r in rows]) / BOHR, "E": e})
        i += 2 + n
    return out


def xyz_bohr(symbols, bohr) -> str:
    return "\n".join(f"{s:2s} {x:20.12f} {y:20.12f} {z:20.12f}" for s, (x, y, z) in zip(symbols, bohr))


def d3_batch(tag: str, geoms: list[dict]) -> list[float]:
    """D3 alone at each geometry: one ORCA run with one sub-job per geometry ($new_job), alone in
    its own folder; values from each sub-job's property file (16 digits)."""
    wd = CK / "d3" / tag
    parts = []
    for k, g in enumerate(geoms):
        parts.append(f"! B3LYP D3BJ STO-3G Bohrs NoPop\n%maxcore 1000\n%base \"g{k:04d}\"\n%scf MaxIter 300 end\n"
                     f"* xyz {CHARGE} {MULT}\n{xyz_bohr(g['symbols'], g['bohr'])}\n*\n")
    inp = "\n$new_job\n".join(parts)
    r = orca.run(wd, "d3", inp, LOG, note=f"check 1: D3 alone at {len(geoms)} geometries of {tag}")
    if not r.ok:
        raise orca.OrcaError(f"D3 batch {tag} did not end normally")
    vals = []
    for k in range(len(geoms)):
        t = (wd / f"g{k:04d}.property.txt").read_text()
        vals.append(float(re.search(r"&vdW \[&Type \"Double\"\]\s+(\S+)", t).group(1)))
    return vals


def engrad(tag: str, k: int, g: dict) -> dict:
    """Energy and gradient at the job's own level, alone in its own folder."""
    wd = CK / "engrad" / tag / f"e{k:04d}"
    inp = (f"! {LEVEL} EnGrad Bohrs NoPop\n%pal nprocs 2 end\n%maxcore 3000\n{SMD}\n\n"
           f"* xyz {CHARGE} {MULT}\n{xyz_bohr(g['symbols'], g['bohr'])}\n*\n")
    r = orca.run(wd, "e", inp, LOG, note=f"check 1: energy and gradient, {tag} #{k}")
    if not r.ok:
        raise orca.OrcaError(f"{tag} #{k} did not end normally")
    t = (wd / "e.property.txt").read_text()
    fe = float(re.search(r"&FinalEnergy \[&Type \"Double\"\]\s+(\S+)", t).group(1))
    gb = re.search(r"\$SCF_Nuc_Gradient(.*?)\$End", t, re.S).group(1)
    gr = np.array([float(v) for v in re.findall(r"^\s*\d+\s+(-?\d\.\d+e[-+]\d+)\s*$", gb, re.M)])
    for tmp in wd.glob("*.tmp*"):
        tmp.unlink()
    return {"E": fe, "grad": gr}


def part_a() -> dict:
    res = {}
    for tag, prop, which in PART_A:
        ev = read_property(prop)
        d3 = d3_batch(tag, ev)
        dd = [abs(a - e["vdw"]) for a, e in zip(d3, ev)]
        gev = [e for e in ev if e["grad"] is not None and "final" in e]
        if which != "all":
            gev = gev[-which[1]:] if which[0] == "last" else gev[:which[1]]
        with ThreadPoolExecutor(2) as ex:
            eg = list(ex.map(lambda a: engrad(tag, a[0], a[1]), enumerate(gev)))
        dE = [abs(x["E"] - e["final"]) for x, e in zip(eg, gev)]
        dG = [float(np.max(np.abs(x["grad"] - e["grad"]))) for x, e in zip(eg, gev)]
        res[tag] = {
            "property_file": str(prop.relative_to(OUT)), "n_evaluations": len(ev),
            "d3_vs_used": {"max_abs_Eh": max(dd), "at_index": ev[int(np.argmax(dd))]["idx"],
                           "n_over_tol": sum(x > TOL_D3 for x in dd), "tol_Eh": TOL_D3},
            "energy_gradient_vs_used": {"evaluations": [e["idx"] for e in gev], "which": str(which),
                                        "max_abs_dE_Eh": max(dE) if dE else None,
                                        "max_abs_dgrad_Eh_per_bohr": max(dG) if dG else None,
                                        "n_over_tol": sum(a > TOL_E or b > TOL_G for a, b in zip(dE, dG)),
                                        "tol": [TOL_E, TOL_G]},
            "d3_used_range_Eh": [min(e["vdw"] for e in ev), max(e["vdw"] for e in ev)],
            "d3_used": [e["vdw"] for e in ev],
        }
        print(tag, json.dumps({k: v for k, v in res[tag].items() if k != "d3_used"}), flush=True)
    return res


def part_b(job: str) -> dict:
    spec = PART_B[job]
    tag = job.replace("/", "_")
    fr = read_trj(spec["trj"])
    d3 = d3_batch(tag + "_frames", fr)
    with ThreadPoolExecutor(2) as ex:
        eg = list(ex.map(lambda a: engrad(tag, a[0], a[1]), enumerate(fr)))
    rows = []
    for f, x, d in zip(fr, eg, d3):
        tol = TOL_E + float(np.sum(np.abs(x["grad"]))) * ROUND_A / BOHR
        rows.append({"frame": f["frame"], "E_trj": f["E"], "E_recomputed": x["E"], "dE": x["E"] - f["E"],
                     "tol": tol, "d3_own": d})
    worst = max(rows, key=lambda r: abs(r["dE"]) / r["tol"])
    out = {"trajectory": str(spec["trj"].relative_to(OUT)), "n_frames": len(rows),
           "max_abs_dE_Eh": max(abs(r["dE"]) for r in rows),
           "worst_relative_to_tol": {k: worst[k] for k in ("frame", "dE", "tol")},
           "n_over_tol": sum(abs(r["dE"]) > r["tol"] for r in rows), "rows": rows}
    # how large an exchanged D3 value would have been: own D3 vs every D3 the partner job used
    pa = CK / "results.json"
    if pa.exists():
        partner = json.loads(pa.read_text(encoding="utf-8"))["partA"].get(spec["partner"])
        if partner:
            ps, os_ = spec["partner_slice"], spec["own_slice"]
            diffs = [abs(r["d3_own"] - p) for r in rows[os_[0]:os_[1]]
                     for p in partner["d3_used"][ps[0]:ps[1]]]
            out["exchange_size_Eh"] = {"smallest": min(diffs), "largest": max(diffs),
                                       "note": "|own D3 − a D3 value the job sharing the folder used|"}
    print(job, json.dumps({k: v for k, v in out.items() if k != "rows"}), flush=True)
    return out


def save(key: str, val) -> None:
    p = CK / "results.json"
    d = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    d.setdefault(key if key == "partA" else "partB", {})
    if key == "partA":
        d["partA"] = val
    else:
        d["partB"][key] = val
    p.write_text(json.dumps(d, indent=1, default=float), encoding="utf-8")


def gate() -> dict:
    """check1.json for the verdict. s010 forward and s005 backward ran alone (check 1a: clear);
    s010 backward and s005 forward are clear only if part B finds no mismatch."""
    d = json.loads((CK / "results.json").read_text(encoding="utf-8")) if (CK / "results.json").exists() else {}
    pb = d.get("partB", {})
    st = {"s010/forward": "clear", "s005/backward": "clear"}
    for j in PART_B:
        st[j] = "not checked" if j not in pb else ("clear" if pb[j]["n_over_tol"] == 0 else "mismatch")
    (OUT / "check1.json").write_text(json.dumps({"status": st}, indent=1), encoding="utf-8")
    print(st)
    return st


def _wait_for(pattern: str) -> None:
    import psutil
    while any(pattern in " ".join(p.info["cmdline"] or []) for p in psutil.process_iter(["cmdline"])
              if p.info["cmdline"] and p.pid != __import__("os").getpid() and "python" in (p.info["cmdline"][0] or "").lower()):
        __import__("time").sleep(60)


def lane() -> None:
    """The whole check-1 lane, in order:
    part A; part B for s010 backward; part B for s005 forward once it has finished; reruns alone of
    any job with a mismatch; the step-test verdict (after the s005 backward rerun has finished); then
    the redirect repro (point 3b, 4 cores, when the queue has them)."""
    import subprocess
    py = [sys.executable, "-X", "utf8"]
    step = str(OUT.parents[2] / "00_inputs/step_convergence_test.py")
    d = json.loads((CK / "results.json").read_text(encoding="utf-8")) if (CK / "results.json").exists() else {}
    if "partA" not in d:
        save("partA", part_a())
    if "s010/backward" not in d.get("partB", {}):
        save("s010/backward", part_b("s010/backward"))
    gate()
    _wait_for("rerun_one s005 forward")
    if "s005/forward" not in d.get("partB", {}):
        save("s005/forward", part_b("s005/forward"))
    st = gate()
    for j, s in st.items():
        if s == "mismatch":
            n, dr = j.split("/")
            subprocess.run(py + [step, "rerun_alone", n, dr], cwd=ROOT, check=True)
            st[j] = "clear"
            st[j + " note"] = "mismatch in check 1; rerun alone in its own folder"
            (OUT / "check1.json").write_text(json.dumps({"status": st}, indent=1), encoding="utf-8")
    _wait_for("rerun_one s005 backward")
    subprocess.run(py + [step, "after"], cwd=ROOT)
    # point 3b, when cores are free: the queue gives the repro its 4 cores in turn
    if not (OUT.parent / "redirect_repro/results.json").exists():
        subprocess.run(py + ["-m", "pipeline.redirect_repro", "run"], cwd=ROOT)


if __name__ == "__main__":
    CK.mkdir(parents=True, exist_ok=True)
    LOG = orca.RunLog(CK / "run_log_check1.md")
    if sys.argv[1] == "lane":
        lane()
    if sys.argv[1] == "partA":
        save("partA", part_a())
    elif sys.argv[1] == "partB":
        save(sys.argv[2], part_b(sys.argv[2]))
        gate()
    elif sys.argv[1] == "gate":
        gate()
