"""Per-reaction card (card.md) filled from the computed results.

Free-text fields that need the chemist (course mapping, extra limitations) are read from
00_inputs/card_notes.toml when present, so a rerun never overwrites his wording."""
from __future__ import annotations

import tomllib

from .steps import Reaction, _read_json, level_name
from .units import KJ_TO_KCAL, fmt_kj


def card(rx: Reaction) -> None:
    r = rx.root
    ts = _read_json(r / "01_ts/ts.json")
    bar = _read_json(r / "02_irc/barriers.json")
    meta = _read_json(r / "02_irc/irc_meta.json")
    tr = _read_json(r / "03_orbitals/tracking.json")
    man = _read_json(r / "04_meshes/manifest.json")
    diag = _read_json(r / "01_ts/diagnostics.json") if (r / "01_ts/diagnostics.json").exists() else None
    notes_path = r / "00_inputs/card_notes.toml"
    notes = tomllib.loads(notes_path.read_text(encoding="utf-8")) if notes_path.exists() else {}

    tz, svp = bar["sp"], bar["working"]
    wb, sb = rx.basis(), rx.basis(sp=True)
    changing = [t for t in tr["tracked"] if t["changes"]]
    if tr.get("expected_changes"):       # the names given in reaction.toml, in that order
        ids = [e["tracked_id"] for e in tr["expected_changes"]["items"] if e["found"]]
        changing = [next(t for t in tr["tracked"] if t["id"] == j) for j in ids]
    orb_txt = "; ".join(t.get("label") or f"{t['label_start']} → {t['label_end']}" for t in changing)
    min_ovl = min(t["min_overlap"] for t in changing) if changing else float("nan")

    lines = [
        f"# Reaction card: {rx.cfg['title']}", "",
        "| | |", "|---|---|",
        f"| Reaction | {rx.cfg['title']} |",
        f"| Model system | {rx.cfg['equation']} |",
        f"| Level of theory | {_level(rx.cfg['method'])}; barrier at {_level(rx.cfg['sp_method'])} "
        f"(single points on the working-level geometries) |",
        f"| Solvent model | CPCM/SMD, {rx.cfg['solvent']} |",
        f"| Barrier ({sb}, electronic) | **{fmt_kj(tz['barrier_from_reactant_complex'])}** relative to "
        f"the reactant complex; {fmt_kj(tz['barrier_from_separated_reactants'])} relative to separated reactants |",
        f"| Barrier ({wb}, electronic, the path shown) | {fmt_kj(svp['barrier_from_reactant_complex'])} "
        f"relative to the reactant complex; {fmt_kj(svp['barrier_from_separated_reactants'])} relative to "
        f"separated reactants |",
        f"| Reaction energy ({sb}, electronic) | {fmt_kj(tz['reaction_energy'])} |",
        f"| Imaginary frequency at TS | {abs(ts['imaginary_cm'][0]):.0f}i cm⁻¹ (exactly one) |",
        f"| TS geometry | {_ts_geometry(rx, ts)} |",
        f"| IRC frames | {meta['n_frames']}; reaction coordinate s from {meta['s_min']:.2f} to "
        f"{meta['s_max']:.2f} bohr·amu^1/2 (TS at s = 0) |",
        f"| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±{man['isovalue']} "
        f"{man['isovalue_units']}; tracked frame to frame by maximum overlap with fixed phase "
        f"(lowest frame-to-frame overlap of the shown orbitals {min_ovl:.2f}); canonical HOMO/LUMO "
        f"optional in the viewer |",
        f"| Orbitals that change | {orb_txt} |",
        "| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, "
        "no explicit solvent molecules |",
        f"| Course mapping | {notes.get('course_mapping', 'McMurry chapter __ (to be filled in)')} |",
        "",
    ]
    if diag:
        lines += ["## Barrier at other levels (same geometries)", "",
                  "| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |",
                  "|---|---|---|---|---|"]
        for row in diag["rows"]:
            g = row["dG_TS_vs_separated_1M_kJ"]
            lines.append(f"| {row['level'].replace('_', '/')} | {fmt_kj(row['dE_TS_vs_separated_kJ'])} | "
                         f"{fmt_kj(row['dE_TS_vs_complex_kJ'])} | {fmt_kj(row['dE_reaction_kJ'])} | "
                         f"{fmt_kj(g) if g is not None else '—'} |")
        lines += ["", "ΔG‡ adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M "
                  "standard-state correction to each level's electronic barrier.", ""]
    if notes.get("experimental"):
        lines += [notes["experimental"], ""]
    lines += ["## Known artifacts and limitations", ""]
    lines += [f"- {x}" for x in _auto_limitations(rx, tz, svp, ts, diag, tr)]
    lines += [f"- {x}" for x in notes.get("limitations", [])]
    lines += ["", "## Files", "",
              "- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)",
              "- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product",
              "- `02_irc/energy_profile.png`: energy profile",
              "- `03_orbitals/tracking_log.md`: orbital matching along the path",
              f"- `05_render/{rx.cfg['id']}.mp4`: rendered video, energy inset at "
              f"{'the working level' if rx.cfg.get('default_curve', 'working') == 'working' else rx.cfg['default_curve'].replace('_', '/')}"
              f" (`{rx.cfg['id']}_<level>.mp4`: same video with the other energy curves)",
              f"- `web/viewer/index.html?r={rx.cfg['id']}`: interactive viewer",
              "- `run_log.md`: every ORCA job, wall time, and reruns", "",
              "Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.", ""]
    (r / "card.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


def _level(method: str) -> str:
    """'B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF' → 'B3LYP-D3(BJ)/def2-SVP, RIJCOSX'"""
    return level_name(method) + (", RIJCOSX" if "RIJCOSX" in method.split() else "")


def _ts_geometry(rx: Reaction, ts: dict) -> str:
    if "key_bonds" not in rx.cfg:                    # SN2-type: two bonds and the Nu–C–LG angle
        return (f"{_bond(rx, 'nucleophile')} {ts['d_nu_center_A']:.2f} Å, {_bond(rx, 'leaving_group')} "
                f"{ts['d_center_lg_A']:.2f} Å, angle {ts['angle_nu_center_lg_deg']:.1f}°")
    roles = {"forms": "forming", "breaks": "breaking", "shortens": "becoming double"}
    parts = [f"{b['label']} {ts['bonds_A'][b['label']]:.2f} Å ({roles[b['role']]})" for b in rx.key_bonds()]
    for name, atoms in rx.cfg.get("ts_dihedrals", {}).items():
        parts.append(f"{name} {_dihedral(rx, atoms):.1f}°")
    return ", ".join(parts)


def _dihedral(rx: Reaction, atoms) -> float:
    import numpy as np
    from .geometry import read_xyz
    p = read_xyz(rx.root / "01_ts/ts/ts.xyz").coords[list(atoms)]
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 /= np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _anions(rx: Reaction) -> str:
    from .geometry import formal_charge
    names = [s["smiles"] for _, s in rx.species() if formal_charge(s["smiles"]) < 0]
    pretty = {"[OH-]": "HO⁻", "[Br-]": "Br⁻", "CC[O-]": "EtO⁻", "[Cl-]": "Cl⁻", "[I-]": "I⁻"}
    return ", ".join(pretty.get(n, n) for n in names) or "anions"


def _bond(rx: Reaction, end: str) -> str:
    sy = _read_json(rx.root / "04_meshes/manifest.json")["symbols"]
    ka = rx.cfg["key_atoms"]
    a, b = (ka["nucleophile"], ka["center"]) if end == "nucleophile" else (ka["center"], ka["leaving_group"])
    return f"{sy[b] if end == 'nucleophile' else sy[a]}–{sy[a] if end == 'nucleophile' else sy[b]}"


def _auto_limitations(rx, tz, svp, ts, diag, tr) -> list[str]:
    """Limitations written from the numbers, so they stay true when results change."""
    out = []
    wb, sb = rx.basis(), rx.basis(sp=True)
    rows = {r["level"]: r for r in diag["rows"]} if diag else {}
    best = diag["rows"][-1] if diag else None          # last level = highest (reference)
    proto = next((r for r in rows.values() if "protocol" in r["level"]), None)
    if svp["barrier_from_separated_reactants"] < 0:
        out.append(f"On the working-level curve the TS lies {fmt_kj(-svp['barrier_from_separated_reactants'])} "
                   f"*below* separated reactants, which would suggest no barrier from separated "
                   f"reactants at all. This is mostly the small basis set describing the anion poorly "
                   f"(at {rx.basis(sp=True)} the TS is {fmt_kj(tz['barrier_from_separated_reactants'])} above). "
                   + (f"The viewer and the main video therefore show the same path with "
                      f"{rx.cfg['default_curve'].replace('_', '/')} energies; the working-level "
                      f"curve stays selectable in the viewer."
                      if rx.cfg.get("default_curve", "working") != "working" else
                      "The viewer can show the same path with energies at a higher level."))
    if best and proto and best["dG_TS_vs_separated_1M_kJ"] and proto["dG_TS_vs_separated_1M_kJ"]:
        gap = best["dG_TS_vs_separated_1M_kJ"] - proto["dG_TS_vs_separated_1M_kJ"]
        out.append(f"{rx.level().split('-')[0]} {'underestimates' if gap > 0 else 'overestimates'} this barrier: ΔG‡ = {fmt_kj(proto['dG_TS_vs_separated_1M_kJ'])} "
                   f"at the protocol level versus {fmt_kj(best['dG_TS_vs_separated_1M_kJ'])} at "
                   f"{best['level'].replace('_', '/')} ({abs(gap):.0f} kJ/mol {'lower' if gap > 0 else 'higher'}). The experimental value "
                   f"is given above; the close agreement of the coupled-cluster number is partly "
                   f"fortuitous (implicit solvent, working-level geometries and thermal corrections).")
    rc_svp, rc_tz = svp["reactant_complex"], tz["reactant_complex"]
    txt = (f"Ion–dipole wells: the reactant complex lies {fmt_kj(-rc_svp)} below separated reactants "
           f"on the working-level curve but only {fmt_kj(-rc_tz)} below at {rx.basis(sp=True)}")
    if best:
        rc_best = best["dE_TS_vs_separated_kJ"] - best["dE_TS_vs_complex_kJ"]
        txt += (f", and {fmt_kj(abs(rc_best))} {'above' if rc_best > 0 else 'below'} at "
                f"{best['level'].replace('_', '/')}")
    txt += (f". In {rx.cfg['solvent']} the deep double well of gas-phase ion–molecule reactions is "
            "largely washed out; what remains here is a shallow, partly artificial well (a teaching "
            "point, not a failure).")
    out.append(txt)
    if not wb.startswith(("ma-", "aug-")):
        out.append(f"The {wb} and {sb} basis sets have no diffuse functions; small anions "
                   f"({_anions(rx)}) are described less well than neutral molecules.")
    out.append(f"Implicit solvent only: no explicit {rx.cfg['solvent']} molecules, so no specific "
               "hydrogen bonds to the base/nucleophile or leaving group.")
    meta = _read_json(rx.root / "02_irc/irc_meta.json")
    path_E = [float(ln.split(",")[3]) for ln in
              (rx.root / "02_irc/energies.csv").read_text(encoding="utf-8").splitlines()[1:]]
    for end, e_path in (("reactant", path_E[0]), ("product", path_E[-1])):
        e_opt = svp[f"{end}_complex"]
        if abs(e_opt - e_path) >= 4.0:      # ~1 kcal/mol: a different minimum, not just loose IRC convergence
            out.append(f"The optimized {end} complex used for the energies ({e_opt:.1f} kJ/mol) is a "
                       f"different local minimum on a flat surface from the {'first' if end == 'reactant' else 'last'} "
                       f"frame of the path ({e_path:.1f} kJ/mol, the IRC end shown in the viewer and "
                       f"video), {fmt_kj(abs(e_opt - e_path))} {'lower' if e_opt < e_path else 'higher'} "
                       f"(working level, relative to separated reactants).".replace("-", "−")
                       .replace("−def2", "-def2"))
    h_transfer = str(ts.get("mode_check", {}).get("largest_atom", "")).startswith("H") \
        and not str(ts.get("mode_check", {}).get("largest_atom", "")).startswith("Hg")
    for e in meta.get("display_excluded_frames", []):
        side = "reactant" if e["s_bohr_amu05"] < 0 else "product"
        out.append(f"The IRC point {e['steps_from_TS']} step{'s' if e['steps_from_TS'] > 1 else ''} "
                   f"from the TS on the {side} side (s = {e['s_bohr_amu05']:.2f} bohr·amu½) lies "
                   f"{fmt_kj(e['dE_above_TS_kJ'])} above the TS: a straight-line first step off a "
                   f"curved {'hydrogen-transfer ' if h_transfer else ''}path. It is kept in the data "
                   f"files (02_irc/energies.csv) but left out of the viewer path and the video, so "
                   f"the displayed maximum is the TS.")
    for t in meta.get("truncated", []):
        who = ("by decision" if t.get("stopped_by") == "decision" else
               "by the IRC stopping rule (docs/protocol.md)")
        out.append(f"The {t['end']} end of the path was truncated where the energy had flattened: the "
                   f"IRC was stopped {who} at s = {t['s_bohr_amu05']:.2f} bohr·amu½ (IRC step "
                   f"{t['last_step']}), after the bonds had finished changing; beyond it only the "
                   f"separating fragments drift apart. The {t['end']} complex was optimized from that "
                   f"last point.")
    if "key_bonds" in rx.cfg:
        out.append("One conformer only: the path starts from a single hand-built TS guess; other "
                   "rotamers of the reactants were not searched.")
    hom = tr["canonical"]["HOMO"]["overlap_with_previous"]
    if min(hom) < 0.5:
        k = int(min(range(len(hom)), key=lambda i: hom[i]))
        out.append(f"The canonical HOMO changes character abruptly along the path (near frame {k}); "
                   f"canonical orbitals are delocalized and re-ordered by energy, which is why the "
                   f"localized IBOs are the default view.")
    return out
