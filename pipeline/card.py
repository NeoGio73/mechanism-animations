"""Per-reaction card (card.md) filled from the computed results.

Free-text fields that need the chemist (course mapping, extra limitations) are read from
00_inputs/card_notes.toml when present, so a rerun never overwrites his wording."""
from __future__ import annotations

import json
import tomllib
from pathlib import Path

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
        # production reactions only (the M06-2X folders): the earlier B3LYP folders' orbitals are not results (2026-11-26)
        if not d.name.endswith("_m062x"):
            continue
        if d == rx.root or not (d / "card.md").exists() or not (d / "03_orbitals/tracking.json").exists():
            continue
        o = _shown(_read_json(d / "03_orbitals/tracking.json"))
        if o:
            title = tomllib.loads((d / "00_inputs/reaction.toml").read_text(encoding="utf-8"))["title"]
            sn = d / "00_inputs/summary_notes.toml"         # short name as in the SI ("Secondary hydride shift")
            short = (tomllib.loads(sn.read_text(encoding="utf-8")).get("si_row") if sn.exists() else None) \
                or title.split(", ")[0]
            m = min(t['min_overlap'] for t in o)
            others.append(f"{short} {m:.3f}")
            other_min.append((short, m))
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


def _one_molecule(rx: Reaction) -> dict | None:
    """One molecule throughout (the hydride shifts): barriers measured from the zeros of pipeline/results_summary.py
    (lowest-free-energy conformer, J. Picione 2026-11-19), in the shape of diagnostics.json rows and barriers.json
    tables, so the rest of the card reads them unchanged. None for reactions of two or more molecules.

    barriers.json and diagnostics.json were written for two-molecule reactions; for the secondary shift (no species
    in reaction.toml) their "separated reactants" columns hold absolute energies, so the card never reads them
    here (2026-11-26)."""
    import contextlib
    import io
    import tempfile
    from . import results_summary as rs
    from .steps import reference_names
    frag = _read_json(rx.root / "01_ts/fragments.json")
    if frag:
        reac, prod = reference_names(rx, frag)
        if len(reac) != 1 or len(prod) != 1:
            return None
    with tempfile.TemporaryDirectory() as td, contextlib.redirect_stdout(io.StringIO()):
        rs.main(rx.root.name, str(Path(td) / "summary.md"))      # same numbers and cross-checks as docs/results
    Q = rs.LAST["Q"]
    z = lambda x: 0.0 if x is not None and abs(x) < 0.05 else x      # a degenerate shift prints 0.0, not -0.0
    Q = {lv: {k: z(v) for k, v in q.items()} for lv, q in Q.items()}
    tab = lambda q: {"barrier_from_reactant_complex": q["dE_cx"], "barrier_from_separated_reactants": q["dE_sep"],
                     "reaction_energy": q["dEr"], "reactant_complex": q["dE_cx_sep"],
                     "product_complex": q["dE_pc_sep"]}
    diag = _read_json(rx.root / "01_ts/diagnostics.json")
    key = lambda lb: (None if lb.endswith("(working)") else "M06-2X_ma-def2-TZVP" if lb.endswith("(protocol)")
                      else lb)
    rows = [{"level": d["level"], "dE_TS_vs_separated_kJ": Q[key(d["level"])]["dE_sep"],
             "dE_TS_vs_complex_kJ": Q[key(d["level"])]["dE_cx"], "dE_reaction_kJ": Q[key(d["level"])]["dEr"],
             "dG_TS_vs_separated_1M_kJ": Q[key(d["level"])]["dG_sep"],
             "dG_TS_vs_complex_kJ": Q[key(d["level"])]["dG_cx"], "dG_reaction_1M_kJ": Q[key(d["level"])]["dGr"]}
            for d in diag["rows"]]
    # a conformer with its own TS (pipeline/own_barrier.py): the barrier of the path shown is measured from the form
    # the path starts from, and that conformer's barrier through its own TS is given beside it; the TS shown measured
    # from that conformer is not the barrier of any path and is not printed (J. Picione, 2026-11-26)
    from .own_barrier import own_barriers
    own = own_barriers(rx)
    if own:
        for d, row in zip(diag["rows"], rows):
            row.update(dE_TS_vs_separated_kJ=row["dE_TS_vs_complex_kJ"], dG_TS_vs_separated_1M_kJ=row["dG_TS_vs_complex_kJ"],
                       own_dE=own["rows"][key(d["level"])]["dE"], own_dG=own["rows"][key(d["level"])]["dG"])
        for q in (Q[None], Q["M06-2X_ma-def2-TZVP"]):
            q["dE_sep"] = q["dE_cx"]
    return {"working": tab(Q[None]), "sp": tab(Q["M06-2X_ma-def2-TZVP"]), "diag": {"rows": rows},
            "zero": rs.LAST["react"][0], "own": own}


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
    one = _one_molecule(rx)
    if one:                    # card_notes.toml names the two zeros (zero_name, start_name)
        tz, svp, diag = one["sp"], one["working"], one["diag"]
        z_nm = notes.get("zero_name", f"the reactant, lowest-free-energy conformer (`{one['zero']}`)")
        s_nm = notes.get("start_name", "the start of the path (frame 0)")
        o_nm = notes.get("own_name", f"`{one['own']['name']}` through its own TS" if one.get("own") else "")
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
    orb_base = orb_txt                   # the one-page card carries its own plain-words note (orbital_card_note)
    if notes.get("orbital_note"):
        orb_txt += ". " + notes["orbital_note"]
    # J. Picione, 2026-11-26: changes found by tracking that the textbook drawing does not show have their own row,
    # described in card_notes.toml documented_note (hyperconjugation withdrawn or re-formed), never as extra arrows
    doc = [e for e in ec.get("documented", []) if e["found"]]
    doc_txt = (notes.get("documented_note") or "; ".join(f"{e['label']} ({e['note']})" for e in doc)) if doc else ""

    lines = [
        f"# Reaction card: {rx.cfg['title']}", "",
        "| | |", "|---|---|",
        f"| Reaction | {rx.cfg['title']} |",
        f"| Model system | {rx.cfg['equation']} |",
        f"| Level of theory | {_level(rx.cfg['method'])}; barrier at {_level(rx.cfg['sp_method'])} "
        f"(single points on the working-level geometries) |",
        f"| Solvent model | CPCM/SMD, {rx.cfg['solvent']} |",
        *([f"| Barrier ({sb}, electronic) | **{fmt_kj(tz['barrier_from_reactant_complex'])}** for the path shown, from "
           f"{s_nm}; {fmt_kj(one['own']['rows']['M06-2X_ma-def2-TZVP']['dE'])} for {o_nm} |",
           f"| Barrier ({wb}, electronic) | {fmt_kj(svp['barrier_from_reactant_complex'])} for the path shown, from "
           f"{s_nm}; {fmt_kj(one['own']['rows'][None]['dE'])} for {o_nm} |",
           f"| Reaction energy ({sb}, electronic) | {fmt_kj(tz['reaction_energy'])}, from {z_nm} |"]
          if one and one.get("own") else
          [f"| Barrier ({sb}, electronic) | **{fmt_kj(tz['barrier_from_separated_reactants'])}** from {z_nm}; "
           f"{fmt_kj(tz['barrier_from_reactant_complex'])} from {s_nm} |",
           f"| Barrier ({wb}, electronic, the path shown) | {fmt_kj(svp['barrier_from_separated_reactants'])} "
           f"from {z_nm}; {fmt_kj(svp['barrier_from_reactant_complex'])} from {s_nm} |",
           f"| Reaction energy ({sb}, electronic) | {fmt_kj(tz['reaction_energy'])}, from {z_nm} |"] if one else
          [f"| Barrier ({sb}, electronic) | **{fmt_kj(tz['barrier_from_reactant_complex'])}** relative to "
           f"the reactant complex; {fmt_kj(tz['barrier_from_separated_reactants'])} relative to separated reactants |",
           f"| Barrier ({wb}, electronic, the path shown) | {fmt_kj(svp['barrier_from_reactant_complex'])} "
           f"relative to the reactant complex; {fmt_kj(svp['barrier_from_separated_reactants'])} relative to "
           f"separated reactants |",
           f"| Reaction energy ({sb}, electronic) | {fmt_kj(tz['reaction_energy'])} |"]),
        f"| Imaginary frequency at TS | {abs(ts['imaginary_cm'][0]):.0f}i cm⁻¹ (exactly one) |",
        f"| TS geometry | {_ts_geometry(rx, ts)} |",
        f"| IRC frames | {meta['n_frames']}; reaction coordinate s from {meta['s_min']:.2f} to "
        f"{meta['s_max']:.2f} bohr·amu^1/2 (TS at s = 0) |",
        f"| Path quality | {_path_quality(rx)} |",
        f"| Orbitals shown | IBOs (IAO/IBO localization, occupied valence), isovalue ±{man['isovalue']} "
        f"{man['isovalue_units']}; tracked frame to frame by maximum overlap with fixed phase "
        f"({_overlap_text(rx, tr, man, changing, min_ovl)}); canonical HOMO/LUMO "
        f"optional in the viewer |",
        *([f"| Display rotation | {_display_rotation_text(rx)} |"] if _display_rotation_text(rx) else []),
        f"| Orbitals that change | {orb_txt} |",
        *([f"| Also found by tracking (not arrows) | {doc_txt} |"] if doc_txt else []),
        # card_notes.toml mechanism_note: one approved plain-words line on the computed path (J. Picione, E2 2026-11-25)
        *([f"| Along the path | {notes['mechanism_note']}"
           + (" *(Draft wording, to be approved by J. Picione.)*" if notes.get("mechanism_note_draft") else "")
           + " |"] if notes.get("mechanism_note") else []),
        "| What this is | The minimum-energy path, not a trajectory: no temperature, no dynamics, "
        "no explicit solvent molecules |",
        f"| Course mapping | {notes.get('course_mapping', 'McMurry chapter __ (to be filled in)')} |",
        "",
    ]
    if diag:
        f = lambda x: fmt_kj(x) if x is not None else "—"
        full = "dG_TS_vs_complex_kJ" in diag["rows"][0]    # older diagnostics.json: one ΔG column
        lines += ["## Barrier at other levels (same geometries)", ""]
        if one and one.get("own"):
            lines += [f"| Level | ΔE‡, path shown (from {s_nm}) | ΔG‡ (298 K), path shown | ΔE‡, {o_nm} "
                      f"| ΔG‡ (298 K), {o_nm} | Reaction energy ΔE, from {z_nm} | Reaction free energy ΔG (298 K) |",
                      "|---|---|---|---|---|---|---|"]
        elif one:
            lines += [f"| Level | ΔE‡ from {z_nm} | ΔE‡ from {s_nm} | Reaction energy, ΔE | ΔG‡ from {z_nm} (298 K) "
                      f"| ΔG‡ from {s_nm} (298 K) | Reaction free energy ΔG (298 K) |",
                      "|---|---|---|---|---|---|---|"]
        elif full:
            lines += ["| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE "
                      "| ΔG‡ from separated reactants (298 K, 1 M) | ΔG‡ from reactant complex (298 K) "
                      "| Reaction free energy ΔG (298 K, 1 M) |",
                      "|---|---|---|---|---|---|---|"]
        else:
            lines += ["| Level | TS − separated reactants, ΔE | TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |",
                      "|---|---|---|---|---|"]
        for row in diag["rows"]:
            if one and one.get("own"):
                lines.append("| " + " | ".join([row["level"].replace("_", "/"), fmt_kj(row["dE_TS_vs_complex_kJ"]),
                             fmt_kj(row["dG_TS_vs_complex_kJ"]), fmt_kj(row["own_dE"]), fmt_kj(row["own_dG"]),
                             fmt_kj(row["dE_reaction_kJ"]), fmt_kj(row["dG_reaction_1M_kJ"])]) + " |")
                continue
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
        lines += ["", ("ΔG adds the working-level thermal corrections (298.15 K) to each level's electronic energies; "
                       "there is one molecule throughout, so no standard-state correction applies." if one else
                       "ΔG adds the working-level thermal corrections (298.15 K) and the 1 atm → 1 M "
                       "standard-state correction to each level's electronic energies"
                       + ("; the reactant complex → TS step is unimolecular and needs no standard-state "
                          "correction" if full else "") + "."), ""]
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
            er = xr.get("dG_reverse_kJ")        # product → TS, first order, no standard state
            if er is not None:
                lines += ["### Reverse barrier: product → TS (ΔG‡, 298 K)", "",
                          "| Level | ΔG‡(reverse) = ΔG‡ − ΔG | Computed − measured |", "|---|---|---|"]
                for row in diag["rows"]:
                    a, b = row["dG_TS_vs_separated_1M_kJ"], row["dG_reaction_1M_kJ"]
                    if a is None or b is None:
                        continue
                    lines.append(f"| {row['level'].replace('_', '/')} | {fmt_kj(a - b)} | {sg(a - b - er)} |")
                lines += [f"| **{xr.get('reverse_label', xr['label'])}** | {fmt_kj(er)} | — |", "",
                          "The reverse step is first order, so it needs no standard-state correction; its "
                          "error equals the barrier error minus the reaction-free-energy error.", ""]
    cm = _curve_max(man, diag, notes)
    if cm:
        cm["note"] = notes.get("curve_max_note")
        cm["zero"] = z_nm if one else "separated reactants"
        lines += _curve_max_section(cm, notes.get("curve_max_note"))
    gp_path = r / "01_ts/gas_phase.json"
    gp_new = r / "03_gas_phase/gas_phase.md"      # gas-phase geometries (pipeline/gas_opt.py, 2026-11-20)
    if gp_new.exists():
        lines += ["## Comparison with the gas-phase experiment", "", gp_new.read_text(encoding="utf-8").strip(), ""]
    elif gp_path.exists():
        lines += _gas_phase_section(rx, _read_json(gp_path))
    if notes.get("experimental"):
        lines += [notes["experimental"], ""]
        # stays on the card until the value has been checked in the primary paper
        if not notes.get("experimental_confirmed", False):
            lines += ["*Value to be confirmed against the primary paper by J. Picione.*", ""]
    if notes.get("context_note"):       # no measured value for this reaction: context only, never a comparison
        lines += ["## Measured value", "", notes["context_note"], ""]
    lines += ["## Known artifacts and limitations", ""]
    auto = _auto_limitations(rx, tz, svp, ts, diag, tr, has_exp=bool(notes.get("experimental")),
                             gas=gp_path.exists() or gp_new.exists(),
                             exp_dg=(notes.get("experiment_row") or {}).get("dG_TS_vs_separated_1M_kJ"),
                             curve_max=cm, one=bool(one))
    lines += [f"- {x}" for x in auto]
    lines += [f"- {x}" for x in notes.get("limitations", [])]
    lines += ["", "## Files", "",
              "- `02_irc/energies.csv`: frame, s (bohr·amu^1/2), E (Eh, kJ/mol, kcal/mol), key distances (Å)",
              "- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product",
              "- `02_irc/energy_profile.png`: energy profile",
              "- `03_orbitals/tracking_log.md`: orbital matching along the path",
              f"- `05_render/{rx.cfg['id']}.mp4`: rendered video, energy inset at "
              f"{'the working level' if man.get('default_curve', 'working') == 'working' else man['default_curve'].replace('_', '/')}"
              f" (`{rx.cfg['id']}_<level>.mp4`: same video with the other energy curves)",
              f"- `web/viewer/index.html?r={rx.root.name}`: interactive viewer (this folder's data)",
              "- `run_log.md`: every ORCA job, wall time, and reruns", "",
              "Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.", ""]
    if notes.get("layout") == "split":      # one-page card + card_details.md (J. Picione, 2026-11-26; E2 first)
        ctx = {"one": one, **({"z_nm": z_nm, "s_nm": s_nm, "o_nm": o_nm} if one else {})}
        _write_split(rx, notes, ts, meta, man, tr, diag, changing, min_ovl, cm,
                     [a.replace(FORTUITOUS, "") for a in auto], orb_base, doc_txt, ctx)
        return
    (r / "card.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


WHAT_THIS_IS = ("The minimum-energy path, not a recording of how one molecule moves: no temperature, no vibration, "
                "no solvent molecules.")          # J. Picione, 2026-11-26


def _true_minus(text: str) -> str:
    """A hyphen before a digit that is not part of a name (M06-2X, 2-bromopropane, p. 2-38) is a minus sign."""
    import re
    return re.sub(r"(?<![\w)\]])-(?=\d)", "−", text)


def _sgn(x: float) -> str:
    """'+5.5 kJ/mol (+1.3 kcal/mol)'; a kcal/mol value below 0.05 keeps two decimals (−0.05, not −0.0)."""
    k = x * KJ_TO_KCAL
    return f"{x:+.1f} kJ/mol ({k:+.2f} kcal/mol)" if 0 < abs(k) < 0.05 else f"{x:+.1f} kJ/mol ({k:+.1f} kcal/mol)"


def _pm(v: float, e: float) -> str:
    return f"{v:.1f} ± {e:.1f} kJ/mol ({v * KJ_TO_KCAL:.1f} ± {e * KJ_TO_KCAL:.1f} kcal/mol)"


def _details_edit(text: str, notes: dict) -> str:
    """card_notes.toml [[details_edit]] old/new: removes project history from texts read from other files (the
    gas-phase report). Each `old` must occur exactly once, so an edit can never silently miss or hit twice."""
    for e in notes.get("details_edit", []):
        n = text.count(e["old"])
        if n != 1:
            raise ValueError(f"details_edit: {n} matches for {e['old'][:60]!r}")
        text = text.replace(e["old"], e["new"])
    return text


def _write_split(rx, notes, ts, meta, man, tr, diag, changing, min_ovl, cm, auto, orb_base, doc_txt, ctx) -> None:
    """The card as one page a student reads (card.md) and the rest in card_details.md (J. Picione, 2026-11-26).

    card.md: summary rows at the default curve's level (DLPNO-CCSD(T)/ma-def2-TZVP), each with the measured value and
    the difference on the same zero where there is one; the orbitals; the approved "along the path" line; at most six
    limitations in plain words (card_notes.toml card_limitations). Summary kinds (card_notes.toml `summary`):
    "two" (two or more molecules: ΔG‡ and ΔrG from separated reactants, 298 K, 1 M), "gas" (Diels–Alder: the
    measurement is in the gas phase, so the comparison row is the gas-phase Ea from the conformer mixture,
    03_gas_phase/conformers.json, and the solution ΔG‡ is computed only), "low_T" (tertiary shift: ΔG‡ at the
    temperature of the measurement, 01_ts/rotamer_zero.json, from the lowest cation conformer), "own" (secondary
    shift: one barrier per form, pipeline/own_barrier.py). No dates, names or project history on either page; that
    stays in card_notes.toml, run_log.md and docs/si_corrections.md. Every number comes from the same files as the
    full card, the results summary or the files named here; none is recomputed from ORCA output."""
    from .steps import level_name
    r = rx.root
    one = ctx.get("one")
    sol = rx.cfg["solvent"]
    work = level_name(rx.cfg["method"])
    lab = {lv["label"]: level_name(lv["method"]) for lv in rx.cfg.get("diagnostic_levels", [])}
    dflt = man.get("default_curve", "working")
    cc = next(row for row in diag["rows"] if row["level"] == dflt)
    cc_name = lab[dflt]
    kind = notes.get("summary", "own" if one and one.get("own") else "one" if one else "two")
    ms = notes.get("summary_measured") or {}
    lim = notes.get("card_limitations", [])
    if len(lim) > 6:
        raise ValueError("card_limitations: at most six on the one-page card")
    mtxt = lambda key, val, fm=fmt_kj: (f"; measured {fm(ms[key])}, {ms.get(key.replace('kJ', 'label'), ms['label'])}; "
                             f"computed − measured {_sgn(val - ms[key])}"
                             if ms.get(key) is not None else "; no measured value")
    course = notes.get("course_mapping", "")
    course = "" if "__" in course else course
    level_txt = (f"Geometries, reaction path and orbitals: {work}, with the SMD model of {sol}. Energies (the default "
                 f"curve in the viewer and video): {cc_name} single points on those geometries, SMD {sol}; thermal "
                 f"corrections from {work} frequencies.")
    if kind == "gas":
        level_txt += (f" Comparison with the gas-phase measurement: separate gas-phase geometries and frequencies at "
                      f"{work}, {cc_name} single points.")
    # ---------------------------------------------------------------- summary rows
    if kind == "two":
        zero = notes.get("summary_zero", "from separated reactants, 298 K, 1 M")
        summ = [f"| Free energy of activation, ΔG‡ ({zero}) | **{fmt_kj(cc['dG_TS_vs_separated_1M_kJ'])}** at "
                f"{cc_name}{mtxt('dG_TS_kJ', cc['dG_TS_vs_separated_1M_kJ'])} |",
                # ΔrG always with its sign (+36.4, +10.5; J. Picione, 2026-11-27)
                f"| Reaction free energy, ΔrG ({zero}) | {_sgn(cc['dG_reaction_1M_kJ'])} at {cc_name}"
                f"{mtxt('dG_reaction_kJ', cc['dG_reaction_1M_kJ'], _sgn)} |"]
    elif kind == "gas":
        conf = _read_json(r / "03_gas_phase/conformers.json")
        g = next(x for x in conf["rows"] if x["level"] == dflt)["mixture"]
        x = conf["experiment"]
        zero = notes.get("summary_zero", "from separated reactants, 298 K, 1 M")
        summ = [f"| Activation energy, Ea ({notes['gas_label']}) | **{fmt_kj(g['Ea_kJ'])}** at {cc_name}; measured "
                f"{_pm(x['Ea_kJ'], x['Ea_err_kJ'])}, {notes['gas_measured_label']}; computed − measured "
                f"{_sgn(g['Ea_minus_exp_kJ'])} |",
                f"| Free energy of activation in {sol}, ΔG‡ ({zero}) | {fmt_kj(cc['dG_TS_vs_separated_1M_kJ'])} at "
                f"{cc_name}; computed only (the measurement is in the gas phase) |",
                f"| Reaction free energy in {sol}, ΔrG ({zero}) | {fmt_kj(cc['dG_reaction_1M_kJ'])} at {cc_name}; "
                "no measured value |"]
    elif kind == "low_T":
        rz = _read_json(r / "01_ts/rotamer_zero.json")
        x = rz["experiment"]
        g = rz["by_temperature"][f"{float(x['T_K'])}"]["levels"][dflt]["dG_ts_kJ"]
        z = notes["summary_zero"]
        summ = [f"| Free energy of activation at {x['T_K']:.0f} K, ΔG‡ ({z}) | **{fmt_kj(g)}** at {cc_name}; measured "
                f"{_pm(x['dG_kJ'], x['err_kJ'])}, {ms['label']}; computed − measured {_sgn(g - x['dG_kJ'])} |",
                f"| Free energy of activation at 298 K, ΔG‡ (same zero) | {fmt_kj(cc['dG_TS_vs_separated_1M_kJ'])} at "
                f"{cc_name}; computed only |",
                f"| Reaction free energy, ΔrG | {fmt_kj(cc['dG_reaction_1M_kJ'])}: {notes['reaction_note']} |"]
    elif kind == "own":
        summ = [f"| Free energy of activation, ΔG‡ (298 K) | **{fmt_kj(cc['dG_TS_vs_complex_kJ'])}** for the path "
                f"shown, from {ctx['s_nm']}; {fmt_kj(cc['own_dG'])} for {ctx['o_nm']}; both at {cc_name} |",
                f"| Reaction free energy, ΔrG (298 K) | {fmt_kj(cc['dG_reaction_1M_kJ'])} at {cc_name}, from "
                f"{ctx['z_nm']} |"]
    else:
        raise ValueError(f"summary kind {kind!r}")
    summ += [f"| {row['label']} | {row['text']} |" for row in notes.get("summary_rows", [])]
    orb_card = (orb_base + (". Localized orbitals (IBOs), followed frame to frame" if orb_base.endswith(")") else
                            " (localized orbitals, IBOs, followed frame to frame)")
                + (". " + notes["orbital_card_note"] if notes.get("orbital_card_note") else ""))
    card = [
        f"# {rx.cfg['title']}", "",
        "| | |", "|---|---|",
        f"| Reaction | {rx.cfg['equation']}, in {sol} |",
        f"| Level of theory | {level_txt} |",
        *summ,
        f"| Transition state | one imaginary frequency, {abs(ts['imaginary_cm'][0]):.0f}i cm⁻¹; {_ts_geometry(rx, ts)} |",
        f"| Orbitals that change | {orb_card} |",
        *([f"| Also found by tracking (not arrows) | {doc_txt} |"] if doc_txt else []),
        *([f"| Along the path | {notes['mechanism_note']} |"] if notes.get("mechanism_note") else []),
        f"| What this is | {WHAT_THIS_IS} |",
        f"| Course mapping | {course} |", "",
        "## Limitations", "",
        *[f"{i}. {x}" for i, x in enumerate(lim, 1)], "",
        "Energies at every level, " + ("how the measured value was derived, " if notes.get("experimental_details")
                                       else "the context for the missing measured value, "
                                       if notes.get("context_details") else "")
        + "and technical notes: `card_details.md`.", "",
        "Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.", ""]
    # ---------------------------------------------------------------- details
    rows = [row for row in diag["rows"] if not row["level"].endswith("(protocol)")]   # same level as M06-2X/TZVP
    nm = lambda row: (f"{work} (geometries)" if row["level"].endswith("(working)") else
                      lab.get(row["level"], row["level"].replace("_", "/"))
                      + (" (default curve)" if row["level"] == dflt else ""))
    f = lambda x: fmt_kj(x) if x is not None else "—"
    # reaction energies with their sign (J. Picione, 2026-11-27); an exact 0.0 (degenerate shift) stays unsigned
    fs = lambda x: "—" if x is None else fmt_kj(x) if abs(x) < 0.05 else _sgn(x)
    det = [f"# Card details: {rx.cfg['title']}", "",
           f"## Energies at every level (same {work} geometries, SMD {sol})", ""]
    if kind == "own":
        s_nm, o_nm, z_nm = ctx["s_nm"], ctx["o_nm"], ctx["z_nm"]
        det += [f"| Level | ΔE‡, path shown (from {s_nm}) | ΔG‡ (298 K), path shown | ΔE‡, {o_nm} | ΔG‡ (298 K), "
                f"{o_nm} | Reaction energy ΔrE, from {z_nm} | Reaction free energy ΔrG (298 K) |",
                "|---|---|---|---|---|---|---|"]
        det += ["| " + " | ".join([nm(row), fmt_kj(row["dE_TS_vs_complex_kJ"]), fmt_kj(row["dG_TS_vs_complex_kJ"]),
                                   fmt_kj(row["own_dE"]), fmt_kj(row["own_dG"]), fs(row["dE_reaction_kJ"]),
                                   fs(row["dG_reaction_1M_kJ"])]) + " |" for row in rows]
    elif one:
        z_nm, s_nm = ctx["z_nm"], ctx["s_nm"]
        det += [f"| Level | ΔE‡ from {z_nm} | ΔE‡ from {s_nm} | Reaction energy ΔrE | ΔG‡ from {z_nm} (298 K) "
                f"| ΔG‡ from {s_nm} (298 K) | Reaction free energy ΔrG (298 K) |", "|---|---|---|---|---|---|---|"]
        det += ["| " + " | ".join([nm(row), fmt_kj(row["dE_TS_vs_separated_kJ"]), fmt_kj(row["dE_TS_vs_complex_kJ"]),
                                   fs(row["dE_reaction_kJ"]), f(row["dG_TS_vs_separated_1M_kJ"]),
                                   f(row["dG_TS_vs_complex_kJ"]), fs(row["dG_reaction_1M_kJ"])]) + " |" for row in rows]
    else:
        det += ["| Level | ΔE‡ from separated reactants | ΔE‡ from the reactant complex | Reaction energy ΔrE "
                "| ΔG‡ from separated reactants (298 K, 1 M) | ΔG‡ from the reactant complex (298 K) "
                "| Reaction free energy ΔrG (298 K, 1 M) |", "|---|---|---|---|---|---|---|"]
        det += ["| " + " | ".join([nm(row), fmt_kj(row["dE_TS_vs_separated_kJ"]), fmt_kj(row["dE_TS_vs_complex_kJ"]),
                                   fs(row["dE_reaction_kJ"]), f(row["dG_TS_vs_separated_1M_kJ"]),
                                   f(row["dG_TS_vs_complex_kJ"]), fs(row["dG_reaction_1M_kJ"])]) + " |" for row in rows]
        if kind == "two" and ms.get("dG_TS_kJ") is not None:
            det.append(f"| **Measured** | — | — | — | {f(ms.get('dG_TS_kJ'))} | — | {fs(ms.get('dG_reaction_kJ'))} |")
    det += ["", (f"ΔG adds the {work} thermal corrections (298.15 K, quasi-RRHO) to each level's electronic energies; "
                 "there is one molecule throughout, so no standard-state correction applies." if one else
                 f"ΔG adds the {work} thermal corrections (298.15 K, quasi-RRHO) to each level's electronic energies, "
                 "and the 1 atm → 1 M standard-state correction for steps that change the number of molecules; the "
                 "reactant complex → TS step needs none."), ""]
    if kind == "two" and ms.get("dG_TS_kJ") is not None:
        det += [f"Measured ΔG‡: {ms['label']}"
                + (f"; measured ΔrG: {ms.get('dG_reaction_label', ms['label'])}" if ms.get("dG_reaction_kJ")
                   is not None else "") + ".", "",
                "## Computed − measured", ""]
        both = ms.get("dG_reaction_kJ") is not None
        det += (["| Level | ΔG‡ from separated reactants | Reaction free energy ΔrG |", "|---|---|---|"] if both else
                ["| Level | ΔG‡ from separated reactants |", "|---|---|"])
        for row in rows:
            a, b = row["dG_TS_vs_separated_1M_kJ"], row["dG_reaction_1M_kJ"]
            det.append(f"| {nm(row)} | {_sgn(a - ms['dG_TS_kJ'])} |"
                       + (f" {_sgn(b - ms['dG_reaction_kJ'])} |" if both else ""))
        det.append("")
        er = ms.get("dG_reverse_kJ")        # product → TS, first order, no standard state
        if er is not None:
            det += ["## Reverse barrier: product → TS (ΔG‡, 298 K)", "",
                    "| Level | ΔG‡(reverse) = ΔG‡ − ΔrG | Computed − measured |", "|---|---|---|"]
            det += [f"| {nm(row)} | {fmt_kj(row['dG_TS_vs_separated_1M_kJ'] - row['dG_reaction_1M_kJ'])} | "
                    f"{_sgn(row['dG_TS_vs_separated_1M_kJ'] - row['dG_reaction_1M_kJ'] - er)} |" for row in rows]
            det += [f"| **Measured** ({ms['dG_reverse_label']}) | {fmt_kj(er)} | — |", "",
                    "The reverse step is first order, so it needs no standard-state correction; its error equals the "
                    "barrier error minus the reaction-free-energy error.", ""]
    if kind == "low_T":
        T = f"{float(x['T_K'])}"
        det += [f"## Comparison with the measurement at {x['T_K']:.0f} K", "",
                f"ΔG‡ at {x['T_K']:.0f} K, {notes['summary_zero']} (`01_ts/rotamer_zero.json`; same thermal model "
                f"as at 298.15 K). Measured: {_pm(x['dG_kJ'], x['err_kJ'])}, {ms['label']}.", "",
                f"| Level | ΔG‡ ({x['T_K']:.0f} K) | Computed − measured |", "|---|---|---|"]
        for row in rows:
            key = "working" if row["level"].endswith("(working)") else row["level"]
            gv = rz["by_temperature"][T]["levels"][key]["dG_ts_kJ"]
            det.append(f"| {nm(row)} | {fmt_kj(gv)} | {_sgn(gv - x['dG_kJ'])} |")
        det.append("")
    if cm:
        sec = _curve_max_section(dict(cm, exp_label="Measured"), None)
        sec = [ln.replace("| **Measured** |", f"| **Measured**, {ms['label']} |") if ms.get("label") else ln
               for ln in sec]
        det += [ln.replace("### ", "## ").replace(cm["level"], lab.get(dflt, cm["level"])) for ln in sec]
    gp = r / "03_gas_phase/gas_phase.md"
    if kind == "gas" and gp.exists():
        det += ["## Comparison with the gas-phase experiment", "", _details_edit(gp.read_text(encoding="utf-8").strip(),
                                                                                 notes), ""]
    if notes.get("experimental_details"):
        det += ["## The measured value: source and derivation", "", notes["experimental_details"], ""]
    if notes.get("context_details"):
        det += ["## Measured value", "", notes["context_details"], ""]
    # extra sections from files next to the data (e.g. the explicit-water test; J. Picione, 2026-11-27)
    for sec in notes.get("details_sections", []):
        det += [f"## {sec['title']}", "", (r / sec["file"]).read_text(encoding="utf-8").strip(), ""]
    rot = _display_rotation_text(rx)
    det += ["## Path and orbitals", "",
            f"- IRC frames: {meta['n_frames']}; reaction coordinate s from {meta['s_min']:.2f} to {meta['s_max']:.2f} "
            "bohr·amu½ (TS at s = 0).",
            f"- Path quality: {_path_quality(rx)}.",
            f"- Orbitals: IBOs (IAO/IBO localization, occupied valence), isovalue ±{man['isovalue']} "
            f"{man['isovalue_units']}; tracked frame to frame by maximum overlap with fixed phase "
            f"({_overlap_text(rx, tr, man, changing, min_ovl)}); canonical HOMO/LUMO optional in the viewer.",
            *([f"- Direction of the arrows: {notes['orbital_details']}"] if notes.get("orbital_details") else []),
            *([f"- Display rotation: {rot}."] if rot else []), "",
            "## Technical notes", "",
            *[f"- {x}" for x in auto], *[f"- {x}" for x in notes.get("details_notes", [])], "",
            "## Files", "",
            "- `02_irc/energies.csv`: frame, s (bohr·amu½), E (Eh, kJ/mol, kcal/mol), key distances (Å)",
            "- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product",
            "- `02_irc/energy_profile.png`: energy profile",
            "- `03_orbitals/tracking_log.md`: orbital matching along the path",
            f"- `05_render/{rx.cfg['id']}.mp4`: rendered video, energy inset at {cc_name} "
            f"(`{rx.cfg['id']}_<level>.mp4`: the same video with the other energy curves)",
            f"- `web/viewer/index.html?r={rx.root.name}`: interactive viewer (this folder's data)",
            "- `run_log.md`: every ORCA job and rerun", "",
            "Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.", ""]
    (r / "card.md").write_text(_true_minus("\n".join(card)), encoding="utf-8")
    (r / "card_details.md").write_text(_true_minus("\n".join(det)), encoding="utf-8")
    print(_true_minus("\n".join(card)))


CURVE_MAX_MIN_KJ = 1.0      # J. Picione, 2026-11-11
# full card only; the split card's limitation says why the agreement should not be over-read (J. Picione, 2026-11-26)
FORTUITOUS = (" The experimental value is given above; the close agreement of the coupled-cluster number is partly "
              "fortuitous (implicit solvent, working-level geometries and thermal corrections).")


def _curve_max(man: dict, diag: dict | None, notes: dict) -> dict | None:
    """Default-level curve maximum off the working-level TS frame by more than CURVE_MAX_MIN_KJ
    (decided by J. Picione, 2026-11-11; the practice of Malick, Petersson and Montgomery 1998,
    IRCMax). Returns the numbers for the card, or None."""
    key = man.get("default_curve", "working")
    curve = next((c for c in man.get("curves", []) if c["key"] == key), None)
    if curve is None or key == "working" or diag is None:
        return None
    E = curve["E_kJ"]
    ts = man["ts_frame"]
    if E[ts] is None:
        raise RuntimeError("default curve not computed at the TS frame")
    k = max((i for i in range(len(E)) if E[i] is not None), key=lambda i: E[i])   # computed frames only
    d = E[k] - E[ts]
    if k == ts or d <= CURVE_MAX_MIN_KJ:
        return None
    row = next((r for r in diag["rows"] if r["level"] == key), None)
    if row is None or row.get("dG_TS_vs_separated_1M_kJ") is None:
        return None
    irc = lambda i: man["frames"][i].get("irc_frame", i)
    mb = notes.get("measured_barrier") or {}
    xr = notes.get("experiment_row") or {}
    sm = notes.get("summary_measured") or {}           # the split card's measured value (E2)
    exp = mb.get("dG_kJ", xr.get("dG_TS_vs_separated_1M_kJ", sm.get("dG_TS_kJ")))
    g_ts = row["dG_TS_vs_separated_1M_kJ"]
    return {"level": key.replace("_", "/"), "ts_irc": irc(ts), "max_irc": irc(k), "E_ts": E[ts],
            "E_max": E[k], "dE": d, "s_max": man["frames"][k]["s"], "dG_ts": g_ts, "dG_max": g_ts + d,
            "exp": exp, "exp_label": mb.get("label", xr.get("label", sm.get("label", "Experiment")))}


def _curve_max_section(cm: dict, note: str | None) -> list[str]:
    sg = lambda x: f"{x:+.1f} kJ/mol ({x / 4.184:+.1f} kcal/mol)"
    dif = lambda g: sg(g - cm["exp"]) if cm["exp"] is not None else "—"
    out = [f"### Barrier from the highest point of the {cm['level']} curve", "",
           f"The barrier quoted first is the one at the working-level transition state (IRC frame "
           f"{cm['ts_irc']}), with its thermal correction. On the {cm['level']} curve the highest point is "
           f"IRC frame {cm['max_irc']} (s = {cm['s_max']:+.2f} bohr·amu½), {fmt_kj(cm['dE'])} above the TS "
           f"frame ({fmt_kj(cm['E_max'])} against {fmt_kj(cm['E_ts'])} relative to {cm.get('zero', 'separated reactants')}, "
           "electronic). The barrier from that maximum adds this difference and keeps the thermal "
           "correction of the TS frame (there are no frequencies at the maximum). Taking the barrier "
           "from the highest point of high-level single points along a lower-level path is the IRCMax "
           "practice of Malick, Petersson and Montgomery (J. Chem. Phys. 1998, 108, 5704).", "",
           f"| ΔG‡ (298 K, 1 M), {cm['level']} | Value | Computed − measured |", "|---|---|---|",
           f"| From the working-level TS frame (IRC frame {cm['ts_irc']}) | {fmt_kj(cm['dG_ts'])} | "
           f"{dif(cm['dG_ts'])} |",
           f"| From the curve maximum (IRC frame {cm['max_irc']}; thermal correction from the TS frame) | "
           f"{fmt_kj(cm['dG_max'])} | {dif(cm['dG_max'])} |"]
    if cm["exp"] is not None:
        out.append(f"| **{cm['exp_label']}** | {fmt_kj(cm['exp'])} | — |")
    out.append("")
    if note:
        out += [note, ""]
    return out


def _path_quality(rx: Reaction) -> str:
    """Rule decided by J. Picione, 2026-11-11 (pipeline/path_quality.py)."""
    from .path_quality import card_text, quality
    return card_text(quality(rx))


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
    p = read_xyz(Path(f"{rx.ts_stem()}.xyz")).coords[list(atoms)]
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
                   + (f"**{fmt_kj(row['Ea_T_kJ'])}** |" if row["phase"] == "gas phase"
                      else f"{fmt_kj(row['Ea_T_kJ'])} |"))
    if gp.get("experimental_Ea_kJ"):
        x = gp["experimental_Ea_kJ"]
        out += ["", f"Experiment (gas phase): Ea ≈ {fmt_kj(x)}."]
        g = [r for r in gp["rows"] if r["phase"] == "gas phase"]
        so = [r for r in gp["rows"] if r["phase"] != "gas phase" and r["level"] == (g[-1]["level"] if g else None)]
        if g:      # decided by J. Picione, 2026-11-11: gas-phase measurement against gas-phase Ea
            b = g[-1]
            out += [f"Comparison with the gas-phase measurement (gas-phase computed Ea, bold above): "
                    f"{b['level'].replace('_', '/')} gives {fmt_kj(b['Ea_T_kJ'])}, "
                    f"{b['Ea_T_kJ'] - x:+.1f} kJ/mol ({(b['Ea_T_kJ'] - x) / 4.184:+.1f} kcal/mol) from experiment."
                    + (f" The {so[0]['phase']} value at the same level, {fmt_kj(so[0]['Ea_T_kJ'])}, is given "
                       "for reference only, not as the comparison." if so else "")]
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


def _two_way(cm: dict) -> str:
    a, b = cm["dG_ts"] - cm["exp"], cm["dG_max"] - cm["exp"]
    return (f" Against experiment the {cm['level']} ΔG‡ is {a:+.1f} kJ/mol ({a / 4.184:+.1f} kcal/mol) from the "
            f"TS frame and {b:+.1f} kJ/mol ({b / 4.184:+.1f} kcal/mol) from the curve maximum (IRC frame "
            f"{cm['max_irc']}; table above)." + (" " + cm["note"] if cm.get("note") else ""))


def _auto_limitations(rx, tz, svp, ts, diag, tr, has_exp: bool = True, gas: bool = False,
                      exp_dg: float | None = None, curve_max: dict | None = None, one: bool = False) -> list[str]:
    """Limitations written from the numbers, so they stay true when results change."""
    out = []
    wb, sb = rx.basis(), rx.basis(sp=True)
    rows = {r["level"]: r for r in diag["rows"]} if diag else {}
    # reference = the coupled-cluster row (2026-11-26: B3LYP is now the last row, and the cards built on 2026-11-25
    # compared against it while calling it the coupled-cluster number)
    best = next((r for r in diag["rows"] if "CCSD(T)" in r["level"]), diag["rows"][-1]) if diag else None
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
        func = rx.level().split("/")[0].replace("-D3(0)", "").replace("-D3(BJ)", "")      # "M06-2X", not "M06"
        out.append(f"{func} {'underestimates' if gap > 0 else 'overestimates'} this barrier: ΔG‡ = {fmt_kj(proto['dG_TS_vs_separated_1M_kJ'])} "
                   f"at {rx.level(sp=True)} versus {fmt_kj(best['dG_TS_vs_separated_1M_kJ'])} at "
                   f"{best['level'].replace('_', '/')} ({abs(gap):.0f} kJ/mol {'lower' if gap > 0 else 'higher'})."
                   + (" The comparison with experiment is the gas-phase Ea table above." if gas else
                      _two_way(curve_max) if curve_max and curve_max.get("exp") is not None and has_exp else
                      " The measured value is in the table above, with the differences below it."
                      if exp_dg is not None and abs(best["dG_TS_vs_separated_1M_kJ"] - exp_dg) > 4 else
                      "" if one else      # one molecule: the experiment section makes its own comparison
                      FORTUITOUS if has_exp else ""))
    rc_svp, rc_tz = svp["reactant_complex"], tz["reactant_complex"]
    _pos = lambda x: f"{fmt_kj(abs(x))} {'below' if x < 0 else 'above'}"     # sign in the word, not the number
    charged = any(s.get("smiles", "").count("-]") or s.get("smiles", "").count("+]")
                  for _, s in rx.species())
    txt = (f"{'Ion–dipole wells: the reactant complex' if charged else 'Reactant complex (a weak, dispersion-bound pair of the neutral reactants)'} lies {_pos(rc_svp)} separated reactants "
           f"on the working-level curve {'but only' if abs(rc_tz) < abs(rc_svp) and rc_tz * rc_svp > 0 else 'and'} "
           f"{_pos(rc_tz)} at {rx.basis(sp=True)}")
    if best:
        rc_best = best["dE_TS_vs_separated_kJ"] - best["dE_TS_vs_complex_kJ"]
        txt += (f", and {fmt_kj(abs(rc_best))} {'above' if rc_best > 0 else 'below'} at "
                f"{best['level'].replace('_', '/')}")
    txt += (f". In {rx.cfg['solvent']} the deep double well of gas-phase ion–molecule reactions is "
            "largely washed out; what remains here is a shallow, partly artificial well (a teaching "
            "point, not a failure)." if charged else
            "; the barriers above are also given from separated reactants.")
    if not one:                 # one molecule: there is no reactant complex
        out.append(txt)
    if not wb.startswith(("ma-", "aug-")):
        out.append(f"The {wb} and {sb} basis sets have no diffuse functions; small anions "
                   f"({_anions(rx)}) are described less well than neutral molecules.")
    kb = rx.cfg.get("key_bonds")      # the SN2 pilot has none: its key_atoms define Nu–C–LG
    breaks = any(b.get("role") == "breaks" for b in kb) if kb else "leaving_group" in rx.cfg.get("key_atoms", {})
    partners = "the base/nucleophile or leaving group" if breaks else "the nucleophile or the product anion"
    out.append(f"Implicit solvent only: no explicit {rx.cfg['solvent']} molecules"
               + (f", so no specific hydrogen bonds to {partners}." if charged and not one else "."))
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
        # one molecule on that side (the cyanide alkoxide, the cations): no fragments to drift apart (2026-11-26)
        single = one or len([s for role, s in rx.species() if role == f"{t['end']}s"]) == 1
        out.append(f"The {t['end']} end of the path was truncated where the energy had flattened: the "
                   f"IRC was stopped {who} at s = {t['s_bohr_amu05']:.2f} bohr·amu½ (IRC step "
                   f"{t['last_step']}), after the bonds had finished changing; "
                   + (f"the {t['end']}-end minimum was optimized from that last point." if single else
                      f"beyond it only the separating fragments drift apart. The {t['end']} complex was optimized "
                      f"from that last point."))
    if "key_bonds" in rx.cfg and not one:     # one molecule: its conformers are stated in card_notes.toml
        out.append("One conformer only: the path starts from a single hand-built TS guess; other "
                   "rotamers of the reactants were not searched.")
    hom = tr["canonical"]["HOMO"]["overlap_with_previous"]
    if min(hom) < 0.5:
        k = int(min(range(len(hom)), key=lambda i: hom[i]))
        out.append(f"The canonical HOMO changes character abruptly along the path (near frame {k}); "
                   f"canonical orbitals are delocalized and re-ordered by energy, which is why the "
                   f"localized IBOs are the default view.")
    return out
