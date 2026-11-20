"""Sensitivity check: SMD18 Coulomb radius for Br (J. Picione, 2026-11-20). Single points only, on existing geometries.

SMD18 (Engelage, Mohamed, Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652) changed the SMD Coulomb radius of
Br from 3.06 Å to 2.60 Å. Every bromine-containing structure of the energetics (TS, reactant and product complexes,
separated species; pipeline/steps._stationary_points) is recomputed at the working level and at every
[[diagnostic_levels]] entry with `radius[35] 2.60` added to the %cpcm block; everything else is the production
input. The working-level thermal corrections (G - E) are kept. Production numbers stay on standard SMD.

Standard-radius energies: the diagnostics single points (01_ts/diagnostics.json) at the extra levels; at the working
level a single point with the production input on the same geometry, run here, so that both radii are single points
(its difference from the stored optimization energy is recorded as a check).
Every new output is read back: ORCA must print "Radius for Br used is ... 2.6000 Ang." or the script stops.

Output: <rx>/03_br_radius/<level>/<name>/ (jobs), <rx>/03_br_radius/br_radius.json and br_radius.md.
Usage: python -X utf8 -m pipeline.br_radius_check <reaction folder> [--test]
  --test: the bromide jobs at the working level only (both radii), radius lines printed.
"""
import json
import re
import sys
from pathlib import Path

from . import orca
from .geometry import read_xyz
from .steps import Reaction, _read_json, _stationary_points, _write_json, reference_names
from .units import kj

BR_RADIUS_A = 2.60
SOURCE = "SMD18: Engelage, Mohamed, Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652 (Br 3.06 → 2.60 Å)"

_RAD = re.compile(r"Radius for Br used is\s+([\d.]+) Bohr \(=\s*([\d.]+) Ang\.\)")
_DIEL = re.compile(r"CPCM Dielectric\s*:\s*(-?[\d.]+) Eh")
_CDS = re.compile(r"SMD CDS \(Gcds\)\s*:\s*(-?[\d.]+) Eh")
_AREA = re.compile(r"Cavity Surface-area\s*\.\.\.\s*([\d.]+)")


def solv_info(out: Path) -> dict:
    """Cavity radius of Br, cavity area, CPCM dielectric and SMD-CDS energy as printed in an ORCA output."""
    t = out.read_text(errors="replace")
    r, d, c, a = _RAD.findall(t), _DIEL.findall(t), _CDS.findall(t), _AREA.findall(t)
    return {"br_radius_A": float(r[-1][1]) if r else None,
            "cavity_area_A2": float(a[-1]) if a else None,
            "dielectric_kJ": kj(float(d[-1])) if d else None,
            "cds_kJ": kj(float(c[-1])) if c else None}


def main() -> None:
    rx = Reaction.load(Path(sys.argv[1]))
    test = "--test" in sys.argv[2:]
    solv = rx.cfg["solvent"]
    std_block = rx.solvent_block
    new_block = std_block.replace(" end", f" radius[35] {BR_RADIUS_A:.2f} end")
    assert new_block != std_block
    pts = {n: v for n, v in _stationary_points(rx).items() if "Br" in read_xyz(v[0]).symbols}
    if test:
        pts = {n: v for n, v in pts.items() if n == "bromide"}
    wlabel = rx.level().replace("/", "_") + " (working)"
    levels = [{"label": wlabel, "dir": "working", "method": rx.cfg["method"]}]
    if not test:
        levels += [{"label": lv["label"], "dir": lv["label"], "method": lv["method"]}
                   for lv in rx.cfg["diagnostic_levels"]]

    frag = _read_json(rx.root / "01_ts/fragments.json")
    tsj = _read_json(rx.root / "01_ts/ts.json")
    ep = _read_json(rx.root / "02_irc/endpoints.json")
    stored_working = {**{n: v["E_eh"] for n, v in frag.items()}, "ts": tsj["E_eh"],
                      "reactant_complex": ep["reactant_complex"]["E_eh"],
                      "product_complex": ep["product_complex"]["E_eh"]}
    diag = _read_json(rx.root / "01_ts/diagnostics.json")

    res = {"rule": "SMD single points with the Br Coulomb radius of SMD18 (2.60 Å), everything else as in "
                   "production; working-level thermal corrections (J. Picione, 2026-11-20). Production numbers stay "
                   "on standard SMD (3.06 Å).",
           "source": SOURCE, "solvent": solv, "structures": list(pts), "levels": {}}
    out_json = rx.root / "03_br_radius" / ("br_radius_test.json" if test else "br_radius.json")
    for lv in levels:
        L = res["levels"][lv["label"]] = {}
        for name, (xyz, q, m) in pts.items():
            g = read_xyz(xyz)
            kw = f'{lv["method"]} CPCM({solv})'
            row = {}
            if lv["dir"] == "working":
                wd = rx.d(f"03_br_radius/working_std/{name}")
                r0 = orca.run(wd, name, orca.make_input(kw, g, q, m, blocks=std_block, nprocs=rx.cfg["nprocs"],
                                                        maxcore_mb=rx.cfg["maxcore_mb"]),
                              rx.log, note="Br radius check: standard SMD single point, working level", front=True)
                if not r0.ok:
                    raise orca.OrcaError(f"standard single point {name} failed; see {r0.out}")
                row["E_std_eh"] = r0.energy_eh
                row["E_stored_minus_sp_kJ"] = kj(stored_working[name] - r0.energy_eh)
                row["std"] = solv_info(wd / f"{name}.out")
            else:
                row["E_std_eh"] = diag["energies_eh"][lv["label"]][name]
                row["std"] = solv_info(rx.root / f"01_ts/diagnostics/{lv['label']}/{name}/{name}.out")
            wd = rx.d(f"03_br_radius/{lv['dir']}/{name}")
            r = orca.run(wd, name, orca.make_input(kw, g, q, m, blocks=new_block, nprocs=rx.cfg["nprocs"],
                                                   maxcore_mb=rx.cfg["maxcore_mb"]),
                         rx.log, note=f"Br radius check: SMD with Br {BR_RADIUS_A:.2f} Å, {lv['label']}", front=True)
            if not r.ok:
                raise orca.OrcaError(f"Br 2.60 single point {name} {lv['label']} failed; see {r.out}")
            row["new"] = solv_info(wd / f"{name}.out")
            if row["new"]["br_radius_A"] is None or abs(row["new"]["br_radius_A"] - BR_RADIUS_A) > 1e-4:
                raise RuntimeError(f"{wd / name}.out: Br radius {row['new']['br_radius_A']} Å, not {BR_RADIUS_A}")
            if abs(row["std"]["br_radius_A"] - 3.06) > 1e-4:
                raise RuntimeError(f"{name} {lv['label']}: standard Br radius {row['std']['br_radius_A']} Å, not 3.06")
            row["E_new_eh"] = r.energy_eh
            row["delta_kJ"] = kj(r.energy_eh - row["E_std_eh"])
            L[name] = row
            _write_json(out_json, res)
            print(f"{lv['label']:40s} {name:18s} radius {row['std']['br_radius_A']:.2f} -> "
                  f"{row['new']['br_radius_A']:.2f} Å  delta {row['delta_kJ']:+.2f} kJ/mol  "
                  f"CDS {row['std']['cds_kJ']:+.3f} -> {row['new']['cds_kJ']:+.3f}", flush=True)
    if test:
        return

    # barriers and reaction energies, both radii: the production rows plus the shift from the Br-containing species
    reac, prod = reference_names(rx, frag)
    rows = []
    for lv in levels:
        d = {n: v["delta_kJ"] for n, v in res["levels"][lv["label"]].items()}
        sh = lambda names: sum(d.get(n, 0.0) for n in names)
        base = next(r for r in diag["rows"] if r["level"] == lv["label"])
        s = {"dG_TS_vs_separated_1M_kJ": d["ts"] - sh(reac),
             "dG_TS_vs_complex_kJ": d["ts"] - d["reactant_complex"],
             "dG_reaction_1M_kJ": sh(prod) - sh(reac),
             "dG_reverse_from_product_complex_kJ": d["ts"] - d["product_complex"]}
        rows.append({"level": lv["label"], "shift_kJ": s,
                     "std": {k: base[k] for k in ("dG_TS_vs_separated_1M_kJ", "dG_TS_vs_complex_kJ",
                                                  "dG_reaction_1M_kJ")}})
        rows[-1]["new"] = {k: rows[-1]["std"][k] + s[k] for k in rows[-1]["std"]}
    res["reference_species"] = {"reactants": reac, "products": prod}
    res["rows"] = rows
    _write_json(out_json, res)
    write_md(rx, res)
    print(json.dumps(rows, indent=1, ensure_ascii=False))


def _lvname(label: str) -> str:
    if label.endswith("(working)"):
        return "M06-2X/ma-def2-SVP (working level)"
    return label.replace("_", "/")


def _f(x: float, sign: bool = False) -> str:
    s = "+" if sign else ""
    return f"{x:{s}.1f} ({x / 4.184:{s}.1f})"


_SNAME = {"ts": "TS", "reactant_complex": "Reactant complex", "product_complex": "Product complex"}


def write_md(rx: Reaction, res: dict) -> None:
    """03_br_radius/br_radius.md: the section for the results summary (tables; findings from findings.md if present)."""
    import tomllib
    notes = rx.root / "00_inputs"
    meas = {}
    for fn, sec in (("summary_notes.toml", "measured"), ("card_notes.toml", "measured_barrier")):
        p = notes / fn
        if p.exists():
            meas.update({f"{fn}:{k}": v for k, v in tomllib.loads(p.read_text(encoding="utf-8")).get(sec, {}).items()})
    m_sep = meas.get("summary_notes.toml:dG_TS_vs_separated_1M_kJ", meas.get("card_notes.toml:dG_kJ"))
    m_rxn = meas.get("summary_notes.toml:dG_reaction_1M_kJ")
    order = lambda x: (x.startswith("DLPNO"), 0)             # as in the summary tables: DLPNO-CCSD(T) last
    lv = sorted(res["levels"], key=order)
    res = dict(res, rows=sorted(res["rows"], key=lambda r: order(r["level"])))
    names = res["structures"]
    n_out = sum(len(v) for v in res["levels"].values())
    cds = max(abs(r["new"]["cds_kJ"] - r["std"]["cds_kJ"]) for v in res["levels"].values() for r in v.values())
    brm = res["levels"][lv[0]].get("bromide")
    chk = max(abs(r["E_stored_minus_sp_kJ"]) for r in res["levels"][lv[0]].values())
    L = [f"Sensitivity check only: every other number in this summary stays on standard SMD (Br radius 3.06 Å). "
         f"SMD18 (Engelage, Mohamed and Huber, Chem. Eur. J. 2018, doi:10.1002/chem.201803652) changed the SMD "
         f"Coulomb radius of Br from 3.06 to 2.60 Å. The {len(names)} bromine-containing structures ("
         + ", ".join(_SNAME.get(n, n[:1].upper() + n[1:]) for n in names) +
         ") were recomputed at all four levels with only that radius changed: single points on the existing "
         "M06-2X/ma-def2-SVP geometries, SMD " + res["solvent"] + ", working-level thermal corrections. "
         "Species without bromine do not change. kJ/mol (kcal/mol).", "",
         f"**Radius check.** In all {n_out} new outputs ORCA 6.1.1 printed \"Radius for Br used is 4.9133 Bohr "
         f"(= 2.6000 Ang.)\", against 3.0600 Å in the standard outputs. The SMD-CDS term did not change (largest "
         f"change {cds:.3f} kJ/mol): only the electrostatic cavity uses this radius."
         + (f" Bromide's cavity area went from {brm['std']['cavity_area_A2']:.1f} to {brm['new']['cavity_area_A2']:.1f} Å². "
            if brm else " ")
         + f"At the working level the standard single points reproduce the stored optimization energies within "
           f"{chk:.3f} kJ/mol.", "",
         "**Energy change with the 2.60 Å radius, E(2.60 Å) − E(3.06 Å):**", "",
         "| Structure | " + " | ".join(_lvname(x) for x in lv) + " |", "|---" * (len(lv) + 1) + "|"]
    for n in names:
        L.append(f"| {_SNAME.get(n, n[:1].upper() + n[1:])} | " + " | ".join(_f(res["levels"][x][n]["delta_kJ"], True) for x in lv) + " |")
    for key, title, m in (
            ("dG_TS_vs_separated_1M_kJ", "ΔG‡ from separated reactants (298.15 K, 1 M)", m_sep),
            ("dG_TS_vs_complex_kJ", "ΔG‡ from the reactant complex (298.15 K)", None),
            ("dG_reaction_1M_kJ", "Reaction free energy ΔrG (298.15 K, 1 M)", m_rxn)):
        L += ["", f"**{title}:**", "", "| Level | Standard SMD (Br 3.06 Å) | SMD18 Br (2.60 Å) | Change |",
              "|---|---|---|---|"]
        for r in res["rows"]:
            L.append(f"| {_lvname(r['level'])} | {_f(r['std'][key])} | {_f(r['new'][key])} | "
                     f"{_f(r['shift_kJ'][key], True)} |")
        if m is not None:
            L.append(f"| Measured | {_f(m)} | | |")
    fnd = rx.root / "03_br_radius" / "findings.md"
    if fnd.exists():
        L += ["", fnd.read_text(encoding="utf-8").strip()]
    (rx.root / "03_br_radius" / "br_radius.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if "--md" in sys.argv[2:]:
        _rx = Reaction.load(Path(sys.argv[1]))
        write_md(_rx, _read_json(_rx.root / "03_br_radius" / "br_radius.json"))
    else:
        main()
