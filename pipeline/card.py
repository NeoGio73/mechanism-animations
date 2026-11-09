"""Per-reaction card (card.md) filled from the computed results.

Free-text fields that need the chemist (course mapping, extra limitations) are read from
00_inputs/card_notes.toml when present, so a rerun never overwrites his wording."""
from __future__ import annotations

import tomllib

from .steps import Reaction, _read_json
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

    tz, svp = bar["def2-TZVP"], bar["def2-SVP"]
    changing = [t for t in tr["tracked"] if t["changes"]]
    orb_txt = "; ".join(f"{t['label_start']} → {t['label_end']}" for t in changing)
    min_ovl = min(t["min_overlap"] for t in changing) if changing else float("nan")

    lines = [
        f"# Reaction card: {rx.cfg['title']}", "",
        "| | |", "|---|---|",
        f"| Reaction | {rx.cfg['title']} |",
        f"| Model system | {rx.cfg['equation']} |",
        f"| Level of theory | {_level(rx.cfg['method'])}; barrier at {_level(rx.cfg['sp_method'])} "
        f"(single points on the working-level geometries) |",
        f"| Solvent model | CPCM/SMD, {rx.cfg['solvent']} |",
        f"| Barrier (def2-TZVP, electronic) | **{fmt_kj(tz['barrier_from_reactant_complex'])}** relative to "
        f"the reactant complex; {fmt_kj(tz['barrier_from_separated_reactants'])} relative to separated reactants |",
        f"| Barrier (def2-SVP, electronic, the path shown) | {fmt_kj(svp['barrier_from_reactant_complex'])} "
        f"relative to the reactant complex; {fmt_kj(svp['barrier_from_separated_reactants'])} relative to "
        f"separated reactants |",
        f"| Reaction energy (def2-TZVP, electronic) | {fmt_kj(tz['reaction_energy'])} |",
        f"| Imaginary frequency at TS | {abs(ts['imaginary_cm'][0]):.0f}i cm⁻¹ (exactly one) |",
        f"| TS geometry | {_bond(rx, 'nucleophile')} {ts['d_nu_center_A']:.2f} Å, "
        f"{_bond(rx, 'leaving_group')} {ts['d_center_lg_A']:.2f} Å, angle {ts['angle_nu_center_lg_deg']:.1f}° |",
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
    tok = method.split()
    disp = {"D3BJ": "-D3(BJ)", "D3ZERO": "-D3(0)", "D4": "-D4"}
    func = tok[0] + "".join(disp[t] for t in tok if t in disp)
    basis = next(t for t in tok if "def2-" in t and "/" not in t)
    return f"{func}/{basis}" + (", RIJCOSX" if "RIJCOSX" in tok else "")


def _bond(rx: Reaction, end: str) -> str:
    sy = _read_json(rx.root / "04_meshes/manifest.json")["symbols"]
    ka = rx.cfg["key_atoms"]
    a, b = (ka["nucleophile"], ka["center"]) if end == "nucleophile" else (ka["center"], ka["leaving_group"])
    return f"{sy[b] if end == 'nucleophile' else sy[a]}–{sy[a] if end == 'nucleophile' else sy[b]}"


def _auto_limitations(rx, tz, svp, ts, diag, tr) -> list[str]:
    """Limitations written from the numbers, so they stay true when results change."""
    out = []
    rows = {r["level"]: r for r in diag["rows"]} if diag else {}
    best = diag["rows"][-1] if diag else None          # last level = highest (reference)
    proto = next((r for r in rows.values() if "protocol" in r["level"]), None)
    if svp["barrier_from_separated_reactants"] < 0:
        out.append(f"On the working-level curve the TS lies {fmt_kj(-svp['barrier_from_separated_reactants'])} "
                   f"*below* separated reactants, which would suggest no barrier from separated "
                   f"reactants at all. This is mostly the small basis set describing HO⁻ poorly "
                   f"(at def2-TZVP the TS is {fmt_kj(tz['barrier_from_separated_reactants'])} above). "
                   + (f"The viewer and the main video therefore show the same path with "
                      f"{rx.cfg['default_curve'].replace('_', '/')} energies; the working-level "
                      f"curve stays selectable in the viewer."
                      if rx.cfg.get("default_curve", "working") != "working" else
                      "The viewer can show the same path with energies at a higher level."))
    if best and proto and best["dG_TS_vs_separated_1M_kJ"] and proto["dG_TS_vs_separated_1M_kJ"]:
        gap = best["dG_TS_vs_separated_1M_kJ"] - proto["dG_TS_vs_separated_1M_kJ"]
        out.append(f"B3LYP underestimates this barrier: ΔG‡ = {fmt_kj(proto['dG_TS_vs_separated_1M_kJ'])} "
                   f"at the protocol level versus {fmt_kj(best['dG_TS_vs_separated_1M_kJ'])} at "
                   f"{best['level'].replace('_', '/')} ({gap:.0f} kJ/mol lower). The experimental value "
                   f"is given above; the close agreement of the coupled-cluster number is partly "
                   f"fortuitous (implicit solvent, working-level geometries and thermal corrections).")
    rc_svp, rc_tz = svp["reactant_complex"], tz["reactant_complex"]
    txt = (f"Ion–dipole wells: the reactant complex lies {fmt_kj(-rc_svp)} below separated reactants "
           f"on the working-level curve but only {fmt_kj(-rc_tz)} below at def2-TZVP")
    if best:
        rc_best = best["dE_TS_vs_separated_kJ"] - best["dE_TS_vs_complex_kJ"]
        txt += (f", and {fmt_kj(abs(rc_best))} {'above' if rc_best > 0 else 'below'} at "
                f"{best['level'].replace('_', '/')}")
    txt += (". In water the double well of the gas-phase SN2 profile is largely washed out; "
            "what remains here is a shallow, partly artificial well (a teaching point, not a failure).")
    out.append(txt)
    out.append("The def2-SVP and def2-TZVP basis sets have no diffuse functions; small anions "
               "(HO⁻) are described less well than neutral molecules.")
    out.append("Implicit solvent only: no explicit water molecules, so no specific hydrogen "
               "bonds to the nucleophile or leaving group.")
    hom = tr["canonical"]["HOMO"]["overlap_with_previous"]
    if min(hom) < 0.5:
        k = int(min(range(len(hom)), key=lambda i: hom[i]))
        out.append(f"The canonical HOMO changes character abruptly along the path (near frame {k}); "
                   f"canonical orbitals are delocalized and re-ordered by energy, which is why the "
                   f"localized IBOs are the default view.")
    return out
