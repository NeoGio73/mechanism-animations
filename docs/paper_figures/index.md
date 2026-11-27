# Figure stills for the manuscript

PNG only, 2244 px wide = 190 mm at 300 dpi. One composite per reaction, plus one viewer screenshot. Each figure has a
caption file of the same name (`.md`) with the frames, energies, levels and isovalue. Made from the data on disk; no
calculation (`pipeline/paper_figure.py`; Blender EEVEE for the panels, matplotlib for the layout).

Each composite has a row of four IRC frames rendered with the video's styling and the same meshes as the viewer and video:
- the first frame;
- the reactant-side frame at half the barrier height (before the TS);
- the TS;
- the last frame.

The tracked orbitals that change are drawn on the frames, with forming and breaking bonds dashed and the key atoms
labelled. The default-level energy curve (DLPNO-CCSD(T)/ma-def2-TZVP) is below the frames, with every panel marked on it.

| File | Reaction | What is in it |
|---|---|---|
| `sn2_hydroxide_bromomethane.png` | SN2, hydroxide + bromomethane | four frames, two orbitals |
| `sn2_hydroxide_bromomethane_viewer.png` | SN2, hydroxide + bromomethane | viewer screenshot at the TS: energy diagram, orbital switches, camera and playback controls (headless Chrome, production data `?r=sn2_hydroxide_bromomethane_m062x`) |
| `e2_ethoxide_2_bromopropane.png` | E2, ethoxide + 2-bromopropane | four frames, three orbitals; second row: the first three frames in the Newman view along Cβ–Cα (Cβ in front) |
| `nucadd_cyanide_acetone.png` | Cyanide + acetone | four frames, two orbitals |
| `hshift_2_3_dimethyl_2_butyl.png` | Tertiary hydride shift | four frames with the tracked orbital; second row: the canonical LUMO (not an IBO) at the same four frames |
| `hshift_3_methyl_2_butyl.png` | Secondary hydride shift | three panels: the bridged eclipsed form (path start), the C2–C3 rotation TS and the tertiary cation (path end), each with the three tracked orbitals |
| `da_butadiene_ethylene.png` | Diels–Alder, butadiene + ethylene | four frames, three orbitals; second row: IRC frames 188 and 189, where the computed localization switches, as displayed (frame 188 shows the display-only mirror image) and as computed |

Notes:
- The energy curve is below the frames, not beside them, so that the four panels keep the full 190 mm width.
- Labels differ from the viewer and video in two places:
  - **Tertiary shift:** the zero is labelled as the cation (reactant = product), and there is no "separated products"
    line. The viewer and video data still say "separated reactants" and "separated products" for this one-ion
    degenerate shift.
  - **Cyanide:** the product line is labelled as the alkoxide, a single ion; the viewer and video say "separated
    products".
- **Diels–Alder "as computed" panels:** these use a separate mesh build without the display-only mirror frames, in
  `reactions/da_butadiene_ethylene_m062x/05_render/paper_figure/meshes_computed/`. The viewer data in `04_meshes/`
  were not changed.
- **Panel renders and render settings:** in `reactions/<id>_m062x/05_render/paper_figure/`.
