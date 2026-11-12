"""SCF noise floor (J. Picione, 2026-11-12): how much do energy and gradient at one geometry depend on
the starting orbitals, at TightSCF and VeryTightSCF, and do the solvent model or the grid contribute?

Geometries: the six control evaluations of check 1 (s005 forward rerun, B3LYP step-test path, gradient
evaluations 13, 14, 17, 18, 21, 22; none overlapped another job). For each setting and each geometry:
  scratch  energy + gradient from ORCA's default guess, alone in its own new folder
  moread   the same, starting from the orbitals of the IRC's previous evaluation (12, 13, 16, 17, 20, 21),
           as the IRC itself does; that previous evaluation is also run from scratch first
Scatter at one setting = |moread - scratch| (energy, Eh; largest gradient component, Eh/bohr). For the
B3LYP SMD settings also |scratch - value the IRC used|. A repeated from-scratch run is deterministic
(checked on one geometry per setting), so the starting orbitals are the only thing that differs.
Settings: B3LYP step-test level (Tight, VeryTight; gas phase Tight and VeryTight; DefGrid3 Tight;
Tight with no incremental Fock build) and the M06-2X working level (Tight, VeryTight).
Every job: 2 cores, alone in its own new folder, through the core ledger.
Result: <step_convergence>/check1/scf_noise.json
Usage: python -m pipeline.scf_noise
"""
from __future__ import annotations

import json
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import orca
from . import shared_folder_check as sfc

W = sfc.CK / "scf_noise"
LOG = orca.RunLog(sfc.CK / "run_log_check1.md")
SMD = sfc.SMD
B3 = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX"
M6 = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX DefGrid3"
SET = {  # name: (keywords, extra blocks, solvent)
    "B3LYP_Tight": (f"{B3} TightSCF CPCM(water)", "", True),
    "B3LYP_VeryTight": (f"{B3} VeryTightSCF CPCM(water)", "", True),
    "B3LYP_Tight_gas": (f"{B3} TightSCF", "", False),
    "B3LYP_VeryTight_gas": (f"{B3} VeryTightSCF", "", False),
    "B3LYP_Tight_DefGrid3": (f"{B3} TightSCF DefGrid3 CPCM(water)", "", True),
    "B3LYP_Tight_noIncFock": (f"{B3} TightSCF CPCM(water)", "%scf DirectResetFreq 1 end\n", True),
    "M062X_Tight": (f"{M6} TightSCF CPCM(water)", "", True),
    "M062X_VeryTight": (f"{M6} VeryTightSCF CPCM(water)", "", True),
}
TARGET = [13, 14, 17, 18, 21, 22]
PREV = {13: 12, 14: 13, 17: 16, 18: 17, 21: 20, 22: 21}


def inp(setting: str, g: dict, moread: bool) -> str:
    kw, extra, solv = SET[setting]
    head = f"! {kw} EnGrad Bohrs NoPop" + (" MORead" if moread else "") + "\n"
    return (head + "%pal nprocs 2 end\n%maxcore 3000\n" + (SMD + "\n" if solv else "") + extra
            + ('%moinp "prev.gbw"\n' if moread else "")
            + f"\n* xyz {sfc.CHARGE} {sfc.MULT}\n{sfc.xyz_bohr(g['symbols'], g['bohr'])}\n*\n")


def read(wd: Path) -> dict:
    t = (wd / "e.property.txt").read_text()
    fe = float(re.search(r"&FinalEnergy \[&Type \"Double\"\]\s+(\S+)", t).group(1))
    gb = re.search(r"\$SCF_Nuc_Gradient(.*?)\$End", t, re.S).group(1)
    gr = np.array([float(v) for v in re.findall(r"^\s*\d+\s+(-?\d\.\d+e[-+]\d+)\s*$", gb, re.M)])
    out = (wd / "e.out").read_text(errors="replace")
    cyc = [int(x) for x in re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", out)]
    scf_s = re.search(r"SCF iterations\s+\.+\s+([\d.]+) sec", out)
    m = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", out)
    wall = (int(m[1]) * 86400 + int(m[2]) * 3600 + int(m[3]) * 60 + int(m[4]) + int(m[5]) / 1000) if m else None
    return {"E": fe, "grad": gr, "cycles": cyc[-1] if cyc else None,
            "scf_iter_s": float(scf_s.group(1)) if scf_s else None, "wall_s": wall}


def job(setting: str, idx: int, g: dict, moread: bool, tag: str = "") -> dict:
    wd = W / setting / f"{'moread' if moread else 'scratch'}{tag}_{idx:03d}"
    if (wd / "e.property.txt").exists() and "ORCA TERMINATED NORMALLY" in (wd / "e.out").read_text(errors="replace"):
        return read(wd)                                   # finished earlier in this check
    if wd.exists():
        shutil.rmtree(wd)                                 # an unfinished attempt of this check job only
    wd.mkdir(parents=True)
    if moread:
        shutil.copy(W / setting / f"scratch_{PREV[idx]:03d}" / "e.gbw", wd / "prev.gbw")
    r = orca.run(wd, "e", inp(setting, g, moread), LOG,
                 note=f"SCF noise: {setting}, evaluation {idx}, {'from the previous evaluation' if moread else 'from scratch'}")
    if not r.ok:
        raise orca.OrcaError(f"{setting} {idx} did not end normally")
    d = read(wd)
    for f in list(wd.glob("*.tmp*")) + list(wd.glob("prev.gbw")):
        f.unlink()
    return d


def main() -> dict:
    ev = {e["idx"]: e for e in sfc.read_property(sfc.CK / "snapshots/s005_forward_rerun_IRC_F.property.txt")}
    geoms = sorted(set(TARGET) | set(PREV.values()))
    res = {}
    for setting in SET:
        with ThreadPoolExecutor(2) as ex:
            sc = dict(zip(geoms, ex.map(lambda i: job(setting, i, ev[i], False), geoms)))
            mo = dict(zip(TARGET, ex.map(lambda i: job(setting, i, ev[i], True), TARGET)))
        again = job(setting, TARGET[0], ev[TARGET[0]], False, tag="_repeat")
        dE = [mo[i]["E"] - sc[i]["E"] for i in TARGET]
        dG = [float(np.max(np.abs(mo[i]["grad"] - sc[i]["grad"]))) for i in TARGET]
        r = {"target_evaluations": TARGET,
             "moread_minus_scratch_Eh": dE, "moread_vs_scratch_max_grad_Eh_per_bohr": dG,
             "max_abs_dE_Eh": max(abs(x) for x in dE), "max_dgrad_Eh_per_bohr": max(dG),
             "repeat_scratch_dE_Eh": again["E"] - sc[TARGET[0]]["E"],
             "repeat_scratch_dgrad": float(np.max(np.abs(again["grad"] - sc[TARGET[0]]["grad"]))),
             "cycles_scratch": [sc[i]["cycles"] for i in TARGET], "cycles_moread": [mo[i]["cycles"] for i in TARGET],
             "wall_scratch_s": [round(sc[i].get("wall_s", 0), 1) for i in TARGET],
             "wall_moread_s": [round(mo[i].get("wall_s", 0), 1) for i in TARGET],
             "scf_iter_s_scratch": [sc[i]["scf_iter_s"] for i in TARGET],
             "scf_iter_s_moread": [mo[i]["scf_iter_s"] for i in TARGET],
             "E_scratch": [sc[i]["E"] for i in TARGET], "E_moread": [mo[i]["E"] for i in TARGET]}
        if setting.startswith("B3LYP") and SET[setting][2] and "DefGrid3" not in setting:
            used = [ev[i]["final"] for i in TARGET]
            r["scratch_minus_irc_used_Eh"] = [sc[i]["E"] - u for i, u in zip(TARGET, used)]
            r["moread_minus_irc_used_Eh"] = [mo[i]["E"] - u for i, u in zip(TARGET, used)]
            r["scratch_vs_irc_used_max_grad"] = max(float(np.max(np.abs(sc[i]["grad"] - ev[i]["grad"]))) for i in TARGET)
            r["moread_vs_irc_used_max_grad"] = max(float(np.max(np.abs(mo[i]["grad"] - ev[i]["grad"]))) for i in TARGET)
        res[setting] = r
        print(setting, json.dumps({k: v for k, v in r.items() if not k.startswith("E_")}), flush=True)
        (sfc.CK / "scf_noise.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res


def inrun() -> dict:
    """No new jobs: the same geometry evaluated twice inside one IRC run (identical to 1e-11 bohr in the
    per-evaluation property files of check 1), energy and gradient differences at TightSCF.
    Result: <step_convergence>/check1/scf_noise_inrun.json"""
    allE, allG, per = [], [], {}
    for tag, prop, _ in sfc.PART_A:
        ev = sfc.read_property(prop)
        pairs = []
        for i in range(len(ev)):
            for j in range(i + 1, min(i + 12, len(ev))):
                if "final" in ev[i] and "final" in ev[j] and np.max(np.abs(ev[i]["bohr"] - ev[j]["bohr"])) < 1e-11:
                    dG = (float(np.max(np.abs(ev[j]["grad"] - ev[i]["grad"])))
                          if ev[i]["grad"] is not None and ev[j]["grad"] is not None else None)
                    pairs.append({"a": ev[i]["idx"], "b": ev[j]["idx"], "dE_Eh": ev[j]["final"] - ev[i]["final"], "dgrad": dG})
        per[tag] = pairs
        allE += [abs(p["dE_Eh"]) for p in pairs]
        allG += [p["dgrad"] for p in pairs if p["dgrad"] is not None]
    a, g = np.array(allE), np.array(allG)
    out = {"n_pairs": len(a), "dE_median_Eh": float(np.median(a)), "dE_p95_Eh": float(np.percentile(a, 95)),
           "dE_max_Eh": float(a.max()), "n_grad_pairs": len(g), "dgrad_median": float(np.median(g)),
           "dgrad_p95": float(np.percentile(g, 95)), "dgrad_max": float(g.max()), "pairs": per}
    (sfc.CK / "scf_noise_inrun.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "pairs"}, indent=1))
    return out


if __name__ == "__main__":
    import sys
    inrun() if sys.argv[1:] == ["inrun"] else main()
