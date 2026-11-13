"""Path-quality rule (decided by J. Picione, 2026-11-11): within the bond-changing region, flag any
IRC step longer than 2.0 times (first set at 1.5) the median step of that region (not of the whole path, whose long
flat approach lowers the median); the card states the flagged steps compactly and
02_irc/path_quality.json keeps the full list.

Steps are differences of the mass-weighted path length s between consecutive IRC frames
(02_irc/energies.csv). The bond-changing region is the set of frames at which at least one
forming or breaking key bond has covered more than 5 % and less than 95 % of its total change
from the first to the last frame; a step counts if either of its frames is in the region.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

FACTOR = 2.0       # J. Picione, 2026-11-11 (was 1.5); with Adapt_Scale_Displ false steps vary ±50 %
LO, HI = 0.05, 0.95


def region(rx) -> tuple[list[int], np.ndarray, np.ndarray]:
    """IRC frames, their s, and whether each is in the bond-changing region (the definition above);
    also used for the coupled-cluster frame choice (J. Picione, 2026-11-13)."""
    rows, s, frames, inreg, _ = _region_rows(rx)
    return frames, s, inreg


def _region_rows(rx):
    rows = list(csv.DictReader(open(rx.root / "02_irc/energies.csv", encoding="utf-8")))
    s = np.array([float(r["s_bohr_amu05"]) for r in rows])
    frames = [int(r["frame"]) for r in rows]
    if "key_bonds" in rx.cfg:
        kb = [b for b in rx.key_bonds() if b["role"] in ("forms", "breaks")]
        pairs = [tuple(sorted(b["atoms"])) for b in kb]
    else:
        ka = rx.cfg["key_atoms"]
        pairs = [tuple(sorted((ka["center"], ka["nucleophile"]))), tuple(sorted((ka["center"], ka["leaving_group"])))]
    cols = []
    for c in rows[0]:
        if c.startswith("d_") and c.endswith("_A"):
            parts = c[:-2].split("_")
            if tuple(sorted((int(parts[-2]), int(parts[-1])))) in pairs:
                cols.append(c)
    D = np.array([[float(r[c]) for c in cols] for r in rows])
    prog = (D - D[0]) / np.where(np.abs(D[-1] - D[0]) > 1e-9, D[-1] - D[0], 1.0)
    inreg = ((prog > LO) & (prog < HI)).any(axis=1)
    return rows, s, frames, inreg, cols


def quality(rx) -> dict:
    rows, s, frames, inreg, cols = _region_rows(rx)
    ds = np.diff(s)
    step_in = np.array([inreg[k] or inreg[k + 1] for k in range(len(ds))])
    med = float(np.median(ds[step_in])) if step_in.any() else float(np.median(ds))
    flagged = []
    for k, d in enumerate(ds):
        if step_in[k] and d > FACTOR * med:
            flagged.append({"from_frame": frames[k], "to_frame": frames[k + 1], "step_s": float(d),
                            "ratio_to_median": float(d / med), "s_from": float(s[k]), "s_to": float(s[k + 1])})
    out = {"median_step_s": med, "median_of": "steps in the bond-changing region",
           "median_step_whole_path_s": float(np.median(ds)), "factor": FACTOR, "region_frames": [frames[i] for i in np.where(inreg)[0][[0, -1]]]
           if inreg.any() else None, "bonds": cols, "flagged": flagged,
           "first_steps_at_TS": None}
    ts = int(np.argmin(np.abs(s)))
    if 0 < ts < len(s) - 1:
        out["first_steps_at_TS"] = {"backward_s": float(s[ts] - s[ts - 1]), "forward_s": float(s[ts + 1] - s[ts])}
    wd = rx.root / "02_irc"
    outs = [o for o in (wd / "irc_f.out", wd / "irc_b.out", wd / "irc.out") if o.exists()]
    out["completeness"] = path_completeness(outs, "02_irc") if outs else None
    (rx.root / "02_irc/path_quality.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def card_text(q: dict) -> str:
    """Card sentence: count, the largest step, the two steps at the TS; full list in path_quality.json."""
    base = (f"bond-changing region IRC frames {q['region_frames'][0]}–{q['region_frames'][1]}, median step there "
            f"{q['median_step_s']:.3f} bohr·amu½; steps longer than {q['factor']:g} × that median: ")
    c = q.get("completeness")
    if c and c["verdict"] != "complete":     # J. Picione, 2026-11-12: never shown as finished
        why = "; ".join(f"{s}: {', '.join(v['why_incomplete'])}" for s, v in c["sides"].items() if v["why_incomplete"])
        base = (f"**IRC INCOMPLETE, not a result** ({why or 'missing side: ' + ', '.join(c['missing_sides'])}); "
                + base)
    fl = q["flagged"]
    if not fl:
        return base + "none"
    one = lambda f: (f"frames {f['from_frame']}→{f['to_frame']} (s {f['s_from']:+.2f} → {f['s_to']:+.2f}), "
                     f"{f['step_s']:.3f} bohr·amu½ = {f['ratio_to_median']:.1f} × median")
    if len(fl) <= 4:
        return base + "; ".join(one(f) for f in fl)
    big = max(fl, key=lambda f: f["step_s"])
    at_ts = [f for f in fl if abs(f["s_from"]) < 1e-6 or abs(f["s_to"]) < 1e-6]
    txt = base + f"{len(fl)} steps, longest {one(big)}"
    rest = [f for f in at_ts if f is not big]
    if rest:
        txt += "; at the TS also " + "; ".join(one(f) for f in rest)
    return txt + " (full list in `02_irc/path_quality.json`)"


# ---------------------------------------------------------------------------------------------------
# Completeness rule (J. Picione, 2026-11-12): an IRC is a result only if it reached its end point by the
# stopping rule, on both sides. A side has reached its end point when
#   (a) ORCA ended it with "THE IRC HAS CONVERGED" (gradient criterion), the job terminated normally,
#       and it is not a false stop next to the TS: at least MIN_STEPS steps and an end at least
#       MIN_DROP_KJ below the TS (the SN2 0.025 runs of 2026-11-11/12 stopped after one step that way); or
#   (b) the project's IRC stopping rule (docs/protocol.md; pipeline.steps.irc_stop_index: energy change
#       < 0.05 kJ/mol per step for 20 steps after the last bond change) is met within the side, and the
#       job did not crash.
# Everything else (MaxIter before the rule is met, a crash, a stop without the rule) is incomplete,
# never finished. A job still running is listed as running and not judged.

import re
import time

MIN_STEPS = 10
MIN_DROP_KJ = 1.0
EH_KJ = 2625.4996
_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")
_HEAD = re.compile(r"^\s*\*\s+(FORWARD|BACKWARD) IRC\s+\*\s*$")


def _rule_hit(rx, trj: Path):
    from .geometry import read_xyz_frames
    from .steps import irc_stop_index
    if rx is None or not trj.exists():
        return None
    try:
        return irc_stop_index(rx, read_xyz_frames(trj))
    except (ValueError, IndexError, KeyError):
        return None


def irc_completeness(out: Path, rx=None) -> dict:
    """Per side of one ORCA IRC output: how it ended, steps, energy drop from the TS, complete or not
    (and why). rx (a Reaction) is needed for the stopping rule; without it only (a) is tested."""
    out = Path(out)
    text = out.read_text(errors="replace")
    lines = text.splitlines()
    normal = "ORCA TERMINATED NORMALLY" in text
    crashed = "error termination" in text or "aborting the run" in text
    running = not normal and not crashed and time.time() - out.stat().st_mtime < 900
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    e_ts = float(e[0]) if e else None
    base = out.stem
    stopped_txt = out.with_suffix(".stopped.txt")
    stop_json = out.parent / "irc_stop.json"
    decided = {}
    if stop_json.exists():
        decided = json.loads(stop_json.read_text(encoding="utf-8")).get("directions", {})
    heads = [(i, m.group(1).lower()) for i, ln in enumerate(lines) if (m := _HEAD.match(ln))]
    sides = {}
    for k, (i, side) in enumerate(heads):
        j_end = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        block = lines[i:j_end]
        cut = next((n for n, ln in enumerate(block) if "IRC PATH SUMMARY" in ln), None)
        if cut is not None:
            block = block[:cut]
        rows = [m for m in (_ROW.match(ln) for ln in block) if m]
        btxt = "\n".join(block)
        steps = len(rows)
        e_end = float(rows[-1].group(2)) if rows else None
        drop = (e_ts - e_end) * EH_KJ if (e_ts is not None and e_end is not None) else None
        tag = side[0].upper()
        hit = _rule_hit(rx, out.parent / f"{base}_IRC_{tag}_trj.xyz")
        if "THE IRC HAS CONVERGED" in btxt:
            how = "ORCA converged"
        elif "MAXIMUM NUMBER OF ITERATIONS REACHED" in btxt:
            how = "MaxIter reached"
        elif stopped_txt.exists() and "stopping rule" in stopped_txt.read_text(encoding="utf-8"):
            how = "stopped by the IRC stopping rule"
        elif tag in decided:
            how = f"stopped by decision ({decided[tag].get('decided_by', '?')}, {decided[tag].get('date', '?')})"
        elif crashed:
            how = "crashed"
        elif running:
            how = "running"
        else:
            how = "ended without an IRC status"
        why = []
        if how == "running":
            verdict = "running (not judged)"
        else:
            if crashed:
                why.append("job crashed (a crashed job is never a result)")
            elif how == "ORCA converged":
                # a both-direction job stopped on purpose later (rule or recorded decision) keeps this side
                if not normal and not (stopped_txt.exists() or decided):
                    why.append("job did not terminate normally")
                if steps < MIN_STEPS:
                    why.append(f"only {steps} steps (< {MIN_STEPS})")
                if drop is None or drop < MIN_DROP_KJ:
                    why.append("end only " + (f"{drop:.3f}" if drop is not None else "?")
                               + f" kJ/mol below the TS (< {MIN_DROP_KJ}): gradient criterion met next to the TS")
            elif hit is None:
                why.append(f"{how}; the IRC stopping rule is not met within the side")
            verdict = "complete" if not why else "INCOMPLETE"
        sides[side] = {"ended": how, "steps": steps, "drop_from_TS_kJ": drop,
                       "stopping_rule_met_at_step": hit[1] if hit else None,
                       "verdict": verdict, "why_incomplete": why, "file": str(out)}
    return sides


def _rx_for(path: Path):
    from .steps import Reaction
    for q in [path, *path.parents]:
        if (q / "00_inputs/reaction.toml").exists():
            return Reaction.load(q)
    return None


def path_completeness(outs: list, label: str) -> dict:
    """One path from one or more IRC outputs (irc_f + irc_b, one irc.out with both, or split folders)."""
    sides = {}
    for o in outs:
        sides.update(irc_completeness(o, _rx_for(Path(o).parent)))
    missing = [s for s in ("forward", "backward") if s not in sides]
    vs = [v["verdict"] for v in sides.values()]
    if any(v.startswith("running") for v in vs):
        verdict = "running (not judged)"
    else:
        verdict = "complete" if not missing and all(v == "complete" for v in vs) else "INCOMPLETE"
    return {"path": label, "verdict": verdict, "missing_sides": missing, "sides": sides}


def scan_all(root: Path) -> list:
    """Every IRC output under reactions/ (old and new). A folder holding both sides is one path; the
    step-size test's split folders (<run>/forward, <run>/backward) are paired with their run folder,
    and a side rerun in <run>/<side>/ replaces the old side in <run>/."""
    root = Path(root)
    outs = sorted(p for p in root.glob("reactions/**/irc*.out")
                  if re.fullmatch(r"irc(_[fb])?\.out", p.name) and "redirect_repro" not in p.parts)
    groups: dict = {}
    for o in outs:
        key = o.parent.parent if o.parent.name in ("forward", "backward") else o.parent
        groups.setdefault(key, []).append(o)
    res = []
    for key, os_ in groups.items():
        use = {}
        for o in os_:
            if o.name == "irc.out":
                use["both"] = o
                continue
            s = "forward" if o.name == "irc_f.out" else "backward"
            if s not in use or o.parent.name == s:      # the split rerun wins over the old side
                use[s] = o
        res.append(path_completeness(sorted(use.values()), str(key.relative_to(root / "reactions")).replace("\\", "/")))
    return res


def _side_text(s: str, v: dict) -> str:
    d = v["drop_from_TS_kJ"]
    txt = f"{s}: {v['ended']}, {v['steps']} steps, " + (f"end {d:.2f} kJ/mol below TS" if d is not None else "no energy")
    if v["stopping_rule_met_at_step"] is not None:
        txt += f", stopping rule met at step {v['stopping_rule_met_at_step']}"
    if v["why_incomplete"]:
        txt += " -> " + "; ".join(v["why_incomplete"])
    return txt


def line(r: dict) -> str:
    det = "; ".join(_side_text(s, v) for s, v in r["sides"].items())
    if r["missing_sides"]:
        det += "; missing side: " + ", ".join(r["missing_sides"])
    return f"{r['verdict']:<22} {r['path']}: {det}"


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    res = scan_all(root)
    (root / "docs/irc_completeness.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    for r in res:
        print(line(r))
