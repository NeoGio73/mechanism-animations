"""Energy zero of a one-molecule reaction (J. Picione, 2026-11-19).

Rule: the zero is the lowest-free-energy conformer of the reactant (working-level G), the same rule as for separated
species (docs/protocol.md, "Reference conformers"), not automatically the form the path reaches (frame 0 and its
optimized end minimum, `reactant_complex` in 02_irc/endpoints.json). Candidates are the reactant end minimum and every
`[[zero_conformers]]` entry in reaction.toml (name, result = a JSON with E_eh, G_eh, job and sp = {level label: E_eh}).

`zero(rx)` returns the chosen conformer and its energy at every level. `python -m pipeline.unimolecular_zero <folder>`
writes 02_irc/zero.json (choice, offsets of frame 0 and of the end minimum) and rewrites the kJ/mol and kcal/mol
columns of 02_irc/energies*.csv and the kJ/mol fields of 02_irc/energies_<level>.json from their E_eh values with
that zero. E_eh is never changed and nothing is recomputed. Idempotent: run again after a path level finishes.
"""
import csv
import io
import json
import sys
from pathlib import Path

from .steps import Reaction
from .units import kj

KCAL = 4.184


def zero(rx: Reaction) -> dict:
    R = rx.root
    ep = json.loads((R / "02_irc/endpoints.json").read_text(encoding="utf-8"))["reactant_complex"]
    diag = json.loads((R / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["energies_eh"]
    labels = [lv["label"] for lv in rx.cfg.get("diagnostic_levels", [])]
    cands = {"reactant_complex": {"minimum": ep["minimum"], "E_eh": ep["E_eh"], "G_eh": ep["G_eh"],
                                  "job": ep["xyz"].rsplit(".", 1)[0],
                                  "sp": {lb: diag[lb]["reactant_complex"] for lb in labels}}}
    for c in rx.cfg.get("zero_conformers", []):
        j = json.loads((R / c["result"]).read_text(encoding="utf-8"))
        if j.get("G_eh") is None or j.get("imaginary_cm") or not j.get("still_secondary", True):
            raise RuntimeError(f"zero conformer {c['name']} is not a usable minimum ({c['result']})")
        missing = [lb for lb in labels if lb not in j["sp"]]
        if missing:
            raise RuntimeError(f"zero conformer {c['name']} has no single point at {missing}")
        cands[c["name"]] = {"minimum": c["name"], "E_eh": j["E_eh"], "G_eh": j["G_eh"], "job": j["job"], "sp": j["sp"]}
    best = min(cands, key=lambda k: cands[k]["G_eh"])
    z = cands[best]
    return {"name": best, "minimum": z["minimum"], "job": z["job"], "G_eh": z["G_eh"],
            "E_eh": {None: z["E_eh"], **z["sp"]}, "candidates": cands}


def _rewrite_csv(p: Path, e0: float) -> None:
    rows = list(csv.DictReader(io.StringIO(p.read_text(encoding="utf-8"))))
    head = p.read_text(encoding="utf-8").splitlines()[0].split(",")
    for r in rows:
        if r.get("E_eh"):
            rel = kj(float(r["E_eh"]) - e0)
            r["E_kJ_per_mol"], r["E_kcal_per_mol"] = f"{rel:.3f}", f"{rel / KCAL:.3f}"
    out = io.StringIO()
    wr = csv.DictWriter(out, fieldnames=head, lineterminator="\n")
    wr.writeheader()
    wr.writerows(rows)
    p.write_text(out.getvalue(), encoding="utf-8")


def main(folder: str) -> None:
    rx = Reaction.load(Path(folder))
    R = rx.root
    if json.loads((R / "01_ts/fragments.json").read_text(encoding="utf-8")):
        sys.exit("not a one-molecule reaction (fragments.json lists separated species)")
    z = zero(rx)
    ep = z["candidates"]["reactant_complex"]
    first = lambda p: next(csv.DictReader(io.StringIO(p.read_text(encoding="utf-8"))))
    nm = lambda k: ep["minimum"] if k == "reactant_complex" else k   # in zero.json the path end is named by its minimum
    rep = {"rule": "zero = lowest-free-energy conformer of the reactant, working-level G (J. Picione, 2026-11-19)",
           "zero": nm(z["name"]), "zero_minimum": z["minimum"], "zero_job": z["job"],
           "candidates_G_kJ_vs_zero": {nm(k): kj(v["G_eh"] - z["G_eh"]) for k, v in z["candidates"].items()},
           "path_end_minimum": ep["minimum"], "levels": {}}
    files = {None: R / "02_irc/energies.csv",
             **{lb: R / f"02_irc/energies_{lb}.csv" for lb in z["E_eh"] if lb is not None}}
    for lb, p in files.items():
        e0 = z["E_eh"][lb]
        e_end = ep["E_eh"] if lb is None else ep["sp"][lb]
        d = {"E_zero_eh": e0, "end_minimum_dE_kJ": kj(e_end - e0)}
        if p.exists():
            r0 = first(p)
            assert r0["frame"] == "0"
            d["frame0_dE_kJ"] = kj(float(r0["E_eh"]) - e0) if r0.get("E_eh") else None
            _rewrite_csv(p, e0)
            jp = p.with_suffix(".json")
            if lb is not None and jp.exists():
                j = json.loads(jp.read_text(encoding="utf-8"))
                rows = list(csv.DictReader(io.StringIO(p.read_text(encoding="utf-8"))))
                Emax = float(rows[j["max_frame"]]["E_kJ_per_mol"])
                j.update(E_max_kJ=Emax, zero=f"{z['name']} (lowest-G conformer of the reactant, J. Picione 2026-11-19)")
                j.pop("E_separated_products_kJ", None)
                jp.write_text(json.dumps(j, indent=1, ensure_ascii=False), encoding="utf-8")
            d["files"] = "rewritten"
        else:
            o = R / "02_irc/path_sp" / str(lb) / "f000/f.out"
            if lb is not None and o.exists():
                import re
                m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", o.read_text(errors="replace"))
                d["frame0_dE_kJ"] = kj(float(m[-1]) - e0) if m else None
            d["files"] = "level still running; rerun this after it finishes"
        rep["levels"]["working" if lb is None else lb] = d
    gc = lambda k: z["candidates"][k]["G_eh"] - z["candidates"][k]["E_eh"]
    for lb, d in rep["levels"].items():
        d["end_minimum_dG_kJ"] = d["end_minimum_dE_kJ"] + kj(gc("reactant_complex") - gc(z["name"]))
    (R / "02_irc/zero.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(rep, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main(sys.argv[1])
