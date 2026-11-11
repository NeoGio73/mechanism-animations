"""Per-reaction card (card.md) filled from the computed results.

Free-text fields that need the chemist (course mapping, extra limitations) are read from
00_inputs/card_notes.toml when present, so a rerun never overwrites his wording."""
from __future__ import annotations

import json
import tomllib

from .steps import Reaction, _read_json, level_name
from .units import KJ_TO_KCAL, fmt_kj


def _shown(tr: dict) -> list[dict]:
    """The tracked orbitals a card reports: the expected ones when reaction.toml names them
    (and all were found), otherwise every orbital whose atoms change."""
    ec = tr.get("expected_changes") or {}
    if ec and all(e["found"] for e in ec["items"]):
        return [next(t for t in tr["tracked"] if t["id"] == e["tracked_id"]) for e in ec["items"]]
    return [t for t in tr["tracked"] if t["changes"]]


def _overlap_text(rx: Reaction, tr: dict, man: dict, changing: list[dict], min_ovl) -> str:
    """Lowest frame-to-frame overlap of each shown orbital, where on the path it occurs, and the
    same number for the other reactions that already have a card (reported as computed)."""
    if not changing:
        return "lowest frame-to-frame overlap of the shown orbitals —"
    s_of = {f.get("irc_frame", i): f["s"] for i, f in enumerate(man["frames"])}
    if len(changing) == 1 or len({t["frame_of_min_overlap"] for t in changing}) > 1 or min_ovl >= 0.95:
        per = "; ".join(f"{t.get('label') or t['label_start'] + ' → ' + t['label_end']} {t['min_overlap']:.3f}"
                        f" (IRC frames {t['frame_of_min_overlap'] - 1}→{t['frame_of_min_overlap']})"
                        for t in changing)
        txt = f"lowest frame-to-frame overlap of each shown orbital: {per}"
    else:   # all at the same step: say where that step is
        k = changing[0]["frame_of_min_overlap"]
        per = ", ".join(f"{t.get('label') or t['label_start'] + ' → ' + t['label_end']} {t['min_overlap']:.3f}"
                        for t in changing)
        txt = (f"lowest frame-to-frame overlap of each shown orbital, all at the step from IRC frame {k - 1} "
               f"to {k} (s = {s_of.get(k - 1, float('nan')):+.2f} to {s_of.get(k, float('nan')):+.2f} "
               f"bohr·amu^1/2; TS at frame {man['ts_frame']}): {per}")
        per_frame = tr.get("per_frame") or []
        if len(per_frame) > k + 1:
            ids = [t["id"] for t in changing]
            nb = [per_frame[j]["overlap_with_previous"][i] for j in (k - 1, k + 1) for i in ids]
            txt += f"; the neighbouring steps are {min(nb):.2f}–{max(nb):.2f}"
    others, other_min = [], []
    for d in sorted(rx.root.parent.iterdir()):
        if d == rx.root or not (d / "card.md").exists() or not (d / "03_orbitals/tracking.json").exists():
            continue
        o = _shown(_read_json(d / "03_orbitals/tracking.json"))
        if o:
            title = tomllib.loads((d / "00_inputs/reaction.toml").read_text(encoding="utf-8"))["title"]
            m = min(t['min_overlap'] for t in o)
            others.append(f"{title.split(',')[0]} {m:.3f}")
            other_min.append((title.split(',')[0], m))
    if others:
        txt += "; for comparison, lowest overlap of the shown orbitals in " + ", ".join(others)
    txt += _step_test_text(rx, other_min)
    return txt


def _step_test_text(rx: Reaction, other_min: list[tuple[str, float]]) -> str:
    """Intermediate IBO single points inside a low-overlap step (03_orbitals/step_test_A_B/)."""
    out = ""
    for p in sorted((rx.root / "03_orbitals").glob("step_test_*_*/result.json")):
        a, b = (int(x) for x in p.parent.name.split("_")[-2:])
        pts = sorted(_read_json(p)["points"].values(), key=lambda q: q["fraction_of_step"])
        first = pts[0]
        to_pt = first[f"frame{a}_to_point"].values()
        lo1, hi1 = min(to_pt), max(to_pt)
        on = first[f"point_to_frame{b}"].values()
        lo2, hi2 = min(on), max(on)
        rng = lambda lo, hi, n: f"{lo:.{n}f}" + (f"–{hi:.{n}f}" if round(hi, n) != round(lo, n) else "")
        part = {0.25: "¼", 0.5: "½"}
        frac = {0.25: "first quarter", 0.5: "first half"}.get(first["fraction_of_step"],
                                                              f"first {first['fraction_of_step']:g}")
        names = " and ".join(part.get(q["fraction_of_step"], f"{q['fraction_of_step']:g}") for q in pts)
        pt = part.get(first["fraction_of_step"], f"{first['fraction_of_step']:g}")
        out += (f"; tested with intermediate IBO single points at {names} of step {a}→{b} "
                f"(Cartesian interpolation, `03_orbitals/{p.parent.name}/`): the localized orbitals "
                f"switch abruptly within its {frac}, overlap {rng(lo1, hi1, 2)} from frame {a} to the "
                f"{pt} point, then {rng(lo2, hi2, 3)} on to frame {b}, so more IRC frames would not "
                f"spread the change")
        if other_min and all(m >= 0.99 for _, m in other_min):
            nm = [n for n, _ in other_min]
            who = nm[0] if len(nm) == 1 else ", ".join(nm[:-1]) + " and " + nm[-1]
            out += (f"; {who} show no such switch: every "
                    f"frame-to-frame overlap of their shown orbitals is at least 0.99")
    return out


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
    ec = tr.get("expected_changes") or {}
    if ec and not all(e["found"] for e in ec["items"]):
        moved = [t for t in tr["tracked"] if t["changes"]]
        orb_txt = ("**Expected assignment not matched** (reaction.toml: "
                   + "; ".join(e["label"] for e in ec["items"]) + "). Tracked orbitals that change: "
                   + "; ".join(f"{t['start']} → {t['end']}" for t in moved) + " (0-based atom indices)"
                   + ". See 03_orbitals/tracking_log.md; awaiting J. Picione's decision.")
        changing = moved
    min_ovl = min(t["min_overlap"] for t in changing) if changing else None
    if notes.get("orbital_note"):
        orb_txt += ". " + notes["orbital_note"]

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
        f"({_overlap_text(rx, tr, man, changing, min_ovl)}); canonical HOMO/LUMO "
        f"optional in the viewer |",
        *([f"| Display rotation | {_display_rotation_text(rx)} |"] if _display_rotation_text(rx) else []),
        f"| Orbitals that change | {orb_txt} |",
        "| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, "
        "no explicit solvent molecules |",
        f"| Course mapping | {notes.get('course_mapping', 'McMurry chapter __ (to be filled in)')} |",
        "",
    ]
    if diag:
        f = lambda x: fmt_kj(x) if x is not None else "—"
        full = "dG_TS_vs_complex_kJ" in diag["rows"][0]    # older diagnostics.json: one ΔG column
        lines += ["## Barrier at other levels (same geometries)", ""]
        if full:
            lines += ["| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE "
                      "| ΔG‡ from separated reactants (298 K, 1 M) | ΔG‡ from reactant complex (298 K) "
                      "| Reaction free energy ΔG (298 K, 1 M) |",
                      "|---|---|---|---|---|---|---|"]
        else:
            lines += ["| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |",
                      "|---|---|---|---|---|"]
        for row in diag["rows"]:
            cells = [row["level"].replace("_", "/"), fmt_kj(row["dE_TS_vs_separated_kJ"]),
                     fmt_kj(row["dE_TS_vs_complex_kJ"]), fmt_kj(row["dE_reaction_kJ"]),
                     f(row["dG_TS_vs_separated_1M_kJ"])]
            if full:
                cells += [f(row["dG_TS_vs_complex_kJ"]), f(row["dG_reaction_1M_kJ"])]
            lines.append("| " + " | ".join(cells) + " |")
        xr = notes.get("experiment_row")          # measured values in the same table
        if xr and full:
            cells = [f"**{xr['label']}**", "—", "—", "—", f(xr.get("dG_TS_vs_separated_1M_kJ")), "—",
                     f(xr.get("dG_reaction_1M_kJ"))]
            lines.append("| " + " | ".join(cells) + " |")
        lines += ["", "ΔG adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M "
                  "standard-state correction to each level's electronic energies"
                  + ("; the reactant complex → TS step is unimolecular and needs no standard-state "
                     "correction" if full else "") + ".", ""]
        if xr and full:
            sg = lambda x: f"{x:+.1f} kJ/mol ({x / 4.184:+.1f} kcal/mol)"
            lines += ["### Computed − measured", "",
                      "| Level | ΔG‡ from separated reactants | Reaction free energy ΔG |", "|---|---|---|"]
            for row in diag["rows"]:
                a = row["dG_TS_vs_separated_1M_kJ"]
                b = row["dG_reaction_1M_kJ"]
                ea, eb = xr.get("dG_TS_vs_separated_1M_kJ"), xr.get("dG_reaction_1M_kJ")
                lines.append(f"| {row['level'].replace('_', '/')} | "
                             f"{sg(a - ea) if a is not None and ea is not None else '—'} | "
                             f"{sg(b - eb) if b is not None and eb is not None else '—'} |")
            lines += [""]
    gp_path = r / "01_ts/gas_phase.json"
    if gp_path.exists():
        lines += _gas_phase_section(rx, _read_json(gp_path))
    if notes.get("experimental"):
        lines += [notes["experimental"], ""]
        # stays on the card until the value has been checked in the primary paper
        if not notes.get("experimental_confirmed", False):
            lines += ["*Value to be confirmed against the primary paper by J. Picione.*", ""]
    lines += ["## Known artifacts and limitations", ""]
    lines += [f"- {x}" for x in _auto_limitations(rx, tz, svp, ts, diag, tr,
                                                  has_exp=bool(notes.get("experimental")),
                                                  gas=gp_path.exists(),
                                                  exp_dg=(notes.get("experiment_row") or {}).get(
                                                      "dG_TS_vs_separated_1M_kJ"))]
    lines += [f"- {x}" for x in notes.get("limitations", [])]
    lines += ["", "## Files", "",
              "- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)",
              "- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product",
              "- `02_irc/energy_profile.png`: energy profile",
              "- `03_orbitals/tracking_log.md`: orbital matching along the path",
              f"- `05_render/{rx.cfg['id']}.mp4`: rendered video, energy inset at "
              f"{'the working level' if man.get('default_curve', 'working') == 'working' else man['default_curve'].replace('_', '/')}"
              f" (`{rx.cfg['id']}_<level>.mp4`: same video with the other energy curves)",
              f"- `web/viewer/index.html?r={rx.cfg['id']}`: interactive viewer",
              "- `run_log.md`: every ORCA job, wall time, and reruns", "",
              "Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.", ""]
    (r / "card.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


def _display_rotation_text(rx: Reaction) -> str | None:
    """One sentence on the display-only σ/π rotation (orbital_steps.orbital_display), or None."""
    p = rx.root / "03_orbitals" / "display_rotation.json"
    if not p.exists():
        return None
    d = _read_json(p)
    sym = json.loads((rx.root / "04_meshes" / "manifest.json").read_text(encoding="utf-8"))["symbols"]
    a, b = d["atoms"]
    bond = f"{sym[a]}{a + 1}={sym[b]}{b + 1}"
    f0, fr = d["frames"][0], d["rotated_frames"]
    return (f"Near the reactant end the IBO localization gave the isolated {bond} double bond as two "
            f"bent bonds, not equivalent (s-character {f0['s_character_shown']:.2f} and "
            f"{f0['s_character_partner']:.2f} at IRC frame 0); for display only, the pair was rotated "
            f"within itself into the σ + π pair that spans the same two orbitals, and the π orbital is "
            f"shown, in IRC frames {fr[0]}–{fr[-1]} (where the computed orbital is more than "
            f"{d['max_angle_deg']:g}° from that split; {f0['bend_angle_deg']:.0f}° at frame 0), joining "
            f"the unrotated frames with overlap {d['overlap_across_join']:.4f}; a rotation inside an "
            f"occupied pair changes no energy or density, and the unrotated orbitals are kept in the "
            f"data files (`03_orbitals/display_rotation.json`)")


def _level(method: str) -> str:
    """'B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF' → 'B3LYP-D3(BJ)/def2-SVP, RIJCOSX'"""
    return level_name(method) + (", RIJCOSX" if "RIJCOSX" in method.split() else "")


def _ts_geometry(rx: Reaction, ts: dict) -> str:
    if "key_bonds" not in rx.cfg:                    # SN2-type: two bonds and the Nu–C–LG angle
        return (f"{_bond(rx, 'nucleophile')} {ts['d_nu_center_A']:.2f} Å, {_bond(rx, 'leaving_group')} "
                f"{ts['d_center_lg_A']:.2f} Å, angle {ts['angle_nu_center_lg_deg']:.1f}°")
    roles = {"forms": "forming", "breaks": "breaking", "shortens": "becoming double",
             "lengthens": "becoming single"}
    parts = [f"{b['label']} {ts['bonds_A'][b['label']]:.2f} Å ({roles[b['role']]})" for b in rx.key_bonds()]
    for name, atoms in rx.cfg.get("ts_dihedrals", {}).items():
        parts.append(f"{name} {_dihedral(rx, atoms):.1f}°")
    for name, val in ts.get("angles_deg", {}).items():
        parts.append(f"{name} {val:.1f}°")
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


def _gas_phase_section(rx, gp: dict) -> list[str]:
    """Card section: gas phase against solution at each level, Ea = ΔH‡(T) + 2RT vs experiment."""
    T = f"{gp['T_K']:.0f} K" + (", provisional" if gp["T_provisional"] else "")
    rng = gp.get("T_range_K")
    T_long = T + (f", the midpoint of the measured range {rng[0]}–{rng[1]} K" if rng else "")
    ref = " + ".join(n.replace("_", " ") for n in gp["reference"])
    reac = " + ".join(n.replace("_", " ") for n in gp["reactants"])
    out = ["## Comparison with the gas-phase experiment", "",
           f"Stationary points recomputed without the solvent model at every level (same geometries). "
           f"Barriers from {ref} (the measured reference) unless marked; Ea = ΔH‡(T) + 2RT for a "
           f"gas-phase bimolecular reaction, with ΔH‡ recomputed at T = {T_long}.", "",
           f"| Level | Phase | ΔE‡ from {reac} | ΔE‡ | ΔH‡ (298.15 K) | ΔG‡ (298.15 K) | Ea ({T}) |",
           "|---|---|---|---|---|---|---|"]
    for row in gp["rows"]:
        out.append(f"| {row['level'].replace('_', '/')} | {row['phase']} | "
                   f"{fmt_kj(row['dE_from_reactants_kJ'])} | {fmt_kj(row['dE_from_reference_kJ'])} | "
                   f"{fmt_kj(row['dH_298_kJ'])} | {fmt_kj(row['dG_298_kJ'])} ({row['standard_state']}) | "
                   f"**{fmt_kj(row['Ea_T_kJ'])}** |")
    if gp.get("experimental_Ea_kJ"):
        out += ["", f"Experiment (gas phase): Ea ≈ {fmt_kj(gp['experimental_Ea_kJ'])}."]
    st = [r for r in gp["rows"] if r.get("s_trans_to_s_cis_dE_kJ") is not None]
    if st:
        best = st[-1]
        out += [f"{gp['reference'][0].replace('_', ' ')} → {gp['reactants'][0].replace('_', ' ')}, ΔE: "
                f"{fmt_kj(best['s_trans_to_s_cis_dE_kJ'])} at {best['level'].replace('_', '/')} "
                f"({best['phase']})."]
    out += ["Thermal corrections: working-level harmonic frequencies (SMD), quasi-RRHO entropy for "
            "ΔG‡; the s-trans reference is why ΔE‡ here exceeds the barrier from the s-cis conformer.", ""]
    tst = [r for r in gp["rows"] if r.get("k_TST_L_mol_s")]
    if tst:
        from .gas_phase import ratio_txt, sci
        from .thermo import NA
        out += [f"### Rate constant at {T}", "",
                "Transition-state theory, gas phase: k = (k_B T/h)(RT/p°) exp(−ΔG‡/RT), with ΔG‡ at "
                f"{T} for a 1 atm (101 325 Pa) standard state, no tunnelling correction, and the "
                "reaction-path degeneracy carried by the rotational symmetry numbers. ΔG‡(T) uses each "
                "level's gas-phase ΔE‡ and the working-level thermal corrections (quasi-RRHO entropy, "
                "checked against ORCA's own values at 298.15 K).", "",
                f"| Level | ΔG‡ ({T}, 1 atm) | k(TST) (L mol⁻¹ s⁻¹) | k(TST) (cm³ molecule⁻¹ s⁻¹) | "
                "k(TST) / k(experiment) |", "|---|---|---|---|---|"]
        for r in tst:
            ratio = ratio_txt(r["k_ratio_to_exp"]) if r.get("k_ratio_to_exp") else "—"
            out.append(f"| {r['level'].replace('_', '/')} | {fmt_kj(r['dG_T_kJ'])} | "
                       f"{r['k_TST_L_mol_s']:.3g} | {sci(r['k_TST_L_mol_s'] * 1000 / NA)} | **{ratio}** |")
        if gp.get("k_exp_L_mol_s"):
            k = gp["k_exp_L_mol_s"]
            A = gp["experimental_A_cm3_molecule_s"]
            out += ["", f"Experiment at {T}, from the measured Arrhenius parameters (A = {sci(A)} "
                        f"cm³ molecule⁻¹ s⁻¹ = {sci(A * NA / 1000)} L mol⁻¹ s⁻¹, Ea = "
                        f"{gp['experimental_Ea_kJ']:g} kJ/mol): k = **{k:.3g} L mol⁻¹ s⁻¹** "
                        f"({sci(k * 1000 / NA)} cm³ molecule⁻¹ s⁻¹)."]
        out += [""]
    return out


def _auto_limitations(rx, tz, svp, ts, diag, tr, has_exp: bool = True, gas: bool = False,
                      exp_dg: float | None = None) -> list[str]:
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
                   f"{best['level'].replace('_', '/')} ({abs(gap):.0f} kJ/mol {'lower' if gap > 0 else 'higher'})."
                   + (" The comparison with experiment is the gas-phase Ea table above." if gas else
                      " The measured value is in the table above, with the differences below it."
                      if exp_dg is not None and abs(best["dG_TS_vs_separated_1M_kJ"] - exp_dg) > 4 else
                      " The experimental value is given above; the close agreement of the "
                      "coupled-cluster number is partly fortuitous (implicit solvent, working-level "
                      "geometries and thermal corrections)." if has_exp else ""))
    rc_svp, rc_tz = svp["reactant_complex"], tz["reactant_complex"]
    charged = any(s.get("smiles", "").count("-]") or s.get("smiles", "").count("+]")
                  for _, s in rx.species())
    txt = (f"{'Ion–dipole wells: the reactant complex' if charged else 'Reactant complex (a weak, dispersion-bound pair of the neutral reactants)'} lies {fmt_kj(-rc_svp)} below separated reactants "
           f"on the working-level curve but only {fmt_kj(-rc_tz)} below at {rx.basis(sp=True)}")
    if best:
        rc_best = best["dE_TS_vs_separated_kJ"] - best["dE_TS_vs_complex_kJ"]
        txt += (f", and {fmt_kj(abs(rc_best))} {'above' if rc_best > 0 else 'below'} at "
                f"{best['level'].replace('_', '/')}")
    txt += (f". In {rx.cfg['solvent']} the deep double well of gas-phase ion–molecule reactions is "
            "largely washed out; what remains here is a shallow, partly artificial well (a teaching "
            "point, not a failure)." if charged else
            "; the barriers above are also given from separated reactants.")
    out.append(txt)
    if not wb.startswith(("ma-", "aug-")):
        out.append(f"The {wb} and {sb} basis sets have no diffuse functions; small anions "
                   f"({_anions(rx)}) are described less well than neutral molecules.")
    kb = rx.cfg.get("key_bonds")      # the SN2 pilot has none: its key_atoms define Nu–C–LG
    breaks = any(b.get("role") == "breaks" for b in kb) if kb else "leaving_group" in rx.cfg.get("key_atoms", {})
    partners = "the base/nucleophile or leaving group" if breaks else "the nucleophile or the product anion"
    out.append(f"Implicit solvent only: no explicit {rx.cfg['solvent']} molecules"
               + (f", so no specific hydrogen bonds to {partners}." if charged else "."))
    meta = _read_json(rx.root / "02_irc/irc_meta.json")
    path_E = [float(ln.split(",")[3]) for ln in
              (rx.root / "02_irc/energies.csv").read_text(encoding="utf-8").splitlines()[1:]]
    for end, e_path in (("reactant", path_E[0]), ("product", path_E[-1])):
        e_opt = svp[f"{end}_complex"]
        side = [s["name"] for role, s in rx.species() if role == f"{end}s"]
        if abs(e_opt - e_path) >= 4.0 and len(side) == 1:
            # one molecule on this side: the path ends at a higher-energy conformer of it
            mol = side[0].replace("_", " ")
            num = lambda x: f"{x:.1f}".replace("-", "−")
            out.append(f"The last frame of the path on the {end} side (the IRC end shown in the viewer "
                       f"and video, {num(e_path)} kJ/mol) is a higher-energy conformer of {mol}: the "
                       f"IRC follows the steepest-descent path from the TS and stops where the "
                       f"gradient along it vanishes. The optimized {mol} used for the energies "
                       f"({num(e_opt)} kJ/mol) is {fmt_kj(abs(e_opt - e_path))} "
                       f"{'lower' if e_opt < e_path else 'higher'} (working level, relative to "
                       f"separated reactants).")
        elif abs(e_opt - e_path) >= 4.0:    # ~1 kcal/mol: a different minimum, not just loose IRC convergence
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
