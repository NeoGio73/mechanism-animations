# MECHANISM_ANIMATIONS_HANDOFF.md

Handoff for Claude Code — "Computed Mechanism Animations" pilot (SN2 end to end).
Written 2026-11-08 from the design chat. The companion design document is the Claude Doc
"Computed Mechanism Animations — Project Brief"; this file is the build instruction set.

---

## 0. Working agreement

- The user (Prof. John Picione) is an organic chemist who does not code. **Claude writes
  all code and runs all jobs.** He reviews results, judges the chemistry, and relays
  messages between this build and the design chat by copy and paste. Keep status messages
  short and in plain language; one paragraph per milestone is right.
- He is a morning reviewer. Leave results in a state he can open and judge (a rendered
  frame, a short video, a number with units) rather than a log to read.
- Every number in any message or file carries a unit (kJ/mol, Å, s, MB). No bare numbers.
- In any 2D structure drawing produced for the project (CDXML, SVG), draw formal charges
  with a circle around the sign. This is his convention; follow it everywhere.
- Nothing in this project is private: inputs, outputs, scripts and models are all meant to
  be published openly (CC-BY for content, MIT for code). Write accordingly from day one.

## 1. Goal of the pilot

Produce one reaction, SN2 (HO⁻ + CH₃Br → CH₃OH + Br⁻), all the way through:

1. transition state and IRC computed in ORCA with an implicit solvent,
2. intrinsic bond orbitals (IBOs) computed at every IRC frame and tracked smoothly,
3. a rendered video (Blender) and an interactive three.js page whose scrubber is the
   energy diagram: molecule, orbitals and the energy marker move together,
4. a per-reaction "card" stating level of theory, solvent model, barrier, and known
   limitations.

If the pilot works, the same pipeline runs for 6–7 more textbook reactions (list in the
Brief). Build the pilot as the reusable pipeline, not as a one-off.

## 2. Environment (Omen, Windows 11)

- CPU AMD Ryzen 7 9800X3D (8 cores), 32 GB RAM, RTX 5080 (not used by ORCA).
- **ORCA 6.1.1** installed with the "Complete" option at `F:\ORCA\`. Verify the path on
  first run. ORCA parallel runs on Windows require **Microsoft MPI** and calling
  `orca.exe` by its full path; if MS-MPI is missing, tell the user: it is a free
  Microsoft download, and ask him to install it before continuing.
- Python 3 is available (check version); create a project venv. Expected packages:
  `rdkit`, `numpy`, `scipy`, `scikit-image` (marching cubes), `trimesh`,
  `open3d` or `pymeshlab` (decimation), `cclib` (ORCA output parsing, optional).
- Blender is installed; check version and whether the **Molecular Nodes** add-on is
  present. Install it if not (it reads XYZ trajectories directly).
- ChemOffice 2016 is present for 2D panels (CDXML in, PNG/SVG out) — optional for the pilot.

### First-run checks (do these before anything else and report results)

1. `F:\ORCA\orca.exe` on a water single point: `! B3LYP def2-SVP` — should finish in a
   few seconds. Report wall time and final SCF energy in Eh.
2. The same job with `%pal nprocs 8 end` — confirms MS-MPI works. Report the speedup.
3. `F:\ORCA\orca_plot.exe` runs and shows its menu.
4. Blender headless: `blender -b --python-expr "import bpy; print(bpy.app.version)"`.

## 3. Repository and layout

Create a new public GitHub repository `NeoGio73/mechanism-animations` (ask the user to
create it if you lack permission) with:

```
mechanism-animations/
  README.md                 project description, license, how to cite
  LICENSE                   MIT (code); note CC-BY 4.0 for content in README
  pipeline/                 Python package: geometry → ORCA → IRC → IBO → meshes → web
  reactions/
    sn2_hydroxide_bromomethane/
      00_inputs/            ORCA .inp files, starting XYZ, SMILES
      01_ts/                TS opt output, Hessian, frequency check
      02_irc/               IRC trajectory (xyz), energies (csv: frame, s, E_kJ_per_mol)
      03_orbitals/          per-frame IBO cubes (compressed), tracking log
      04_meshes/            decimated isosurfaces per frame (glb or binary json)
      05_render/            Blender scene + rendered mp4
      card.md               the per-reaction card (see §7)
  web/
    index.html              landing page listing reactions
    viewer/                 three.js viewer (one page per reaction, data-driven)
  docs/                     protocol notes for the paper
```

Large binary outputs (cube files, raw ORCA outputs) are archived to Zenodo later, not
committed; keep the repo under a few hundred MB. Commit decimated meshes and energies.

## 4. Computational protocol

Working level (fast, defensible for teaching pictures):

- `! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF` with CPCM implicit solvent; use SMD
  radii/parameters for water (`%cpcm smd true SMDsolvent "water" end` or the ORCA 6.1
  equivalent — **check the 6.1 manual for exact syntax before running**).
- Reported barrier: single points on the IRC endpoints and TS at `def2-TZVP`, same
  functional and solvent. Report both numbers; the pictures use the def2-SVP path.
- Energies in **kJ/mol** relative to the separated reactants (kcal/mol in parentheses
  wherever a human reads them).

Steps for SN2:

1. **Geometries.** Build HO⁻, CH₃Br, CH₃OH, Br⁻ from SMILES with RDKit; also a
   reactant complex and product complex. Pre-optimize with GFN2-xTB if xtb is
   available (optional), else go straight to DFT.
2. **TS guess.** Collinear O···C···Br, C–O ≈ 2.0 Å, C–Br ≈ 2.3 Å, planar CH₃ at C.
   Alternative if the guess fails: `! NEB-TS` between reactant and product complexes.
3. **TS optimization.** `! OptTS Freq` with `%geom Calc_Hess true Recalc_Hess 5 end`.
   Accept only with exactly **one imaginary frequency** whose mode is the O–C–Br motion.
   Report the frequency in cm⁻¹.
4. **IRC.** `! IRC` both directions from the TS, reading the TS Hessian
   (`%irc InitHess Read Hess_Filename "ts.hess" ... end`); target 40–80 frames total.
   Save `irc_full_trj.xyz` and the energy profile as CSV.
5. **Sanity checks.** The path must connect the reactant complex to the product complex;
   the Walden inversion must be visible in the frames; barrier should be in the expected
   range for solvated hydroxide + bromomethane (tens of kJ/mol, not hundreds). If the
   profile shows the gas-phase double well despite solvent, say so — it is a teaching
   point, not a failure.

Keep a `run_log.md` per reaction: every job, its wall time, and what was changed when a
job failed. Never silently re-run with different settings.

## 5. Orbitals: IBOs tracked along the IRC

Goal: at each IRC frame, localized occupied orbitals that map onto the curved arrows
(nucleophile lone pair → forming C–O bond; breaking C–Br bond → bromide lone pair).

1. For each frame geometry run a single point at the working level with
   `%loc LocMet IAOIBO end` (occupied orbitals only). Consult the ORCA 6.1 manual for the
   current keyword names and the output file (`.loc`) and how `orca_plot` plots localized
   orbitals to Gaussian cube files. Grid: 60–80 points per axis is enough for teaching
   isosurfaces; use the same box for every frame.
2. **Tracking.** Localized orbital order is arbitrary from frame to frame. Match orbitals
   between consecutive frames by maximum overlap (compute overlaps numerically from the
   cubes, or from coefficients if available) and relabel so that "orbital k" is the same
   chemical object across the path. Fix the arbitrary sign so the phase does not flip
   between frames. Log the assignment; flag any frame where the best match overlap drops
   below ~0.7 (that is where a bond changes identity — expected near the TS).
3. Keep only the orbitals that change: the O lone pair that becomes σ(C–O), the σ(C–Br)
   that becomes a Br lone pair, and optionally the σ*(C–Br) (needs virtuals; canonical
   LUMO is acceptable for that one). Everything else stays hidden by default.
4. Isosurfaces at ±0.05 e/bohr³ (adjust once visually), via marching cubes, decimated
   to ≤ 30 kB per orbital per frame, exported as glb or compact binary JSON. Record the
   isovalue in the card.

Fallback if IBOs prove awkward in this ORCA version: Pipek–Mezey or Foster–Boys localized
orbitals, same tracking logic. Canonical HOMO/LUMO as a user-selectable toggle in the
viewer, not the default.

## 6. Deliverables of the pilot

1. `02_irc/energies.csv` and `irc_full_trj.xyz`, plus a PNG of the energy profile.
2. A 15–30 s rendered video (Blender, Molecular Nodes; ball-and-stick, orbitals as
   two-color translucent lobes, white background, 1920 × 1080, 30 fps) with the energy
   diagram inset and a moving marker.
3. The three.js viewer page for SN2: energy diagram as scrubber (drag or arrow keys),
   molecule + tracked orbitals follow; toggles for orbitals on/off, atom labels, canonical
   vs localized; camera presets (side view down the O–C–Br axis is the money shot);
   keyboard accessible; works on GitHub Pages with no build step beyond static files.
4. `card.md` for SN2 (§7).
5. A short `docs/pipeline.md` describing how to add a reaction: which files to provide
   (SMILES of reactants/products, TS guess or NEB endpoints, solvent) and which commands
   to run.

## 7. The per-reaction card (template)

```
Reaction: SN2, hydroxide + bromomethane
Model system: HO⁻ + CH₃Br → CH₃OH + Br⁻
Level of theory: B3LYP-D3(BJ)/def2-SVP, RIJCOSX; barrier at def2-TZVP (same functional)
Solvent model: CPCM/SMD, water
Barrier (def2-TZVP, electronic): __ kJ/mol (__ kcal/mol) relative to reactant complex
Imaginary frequency at TS: __ cm⁻¹
IRC frames: __; reaction coordinate range: __ to __ bohr·amu^1/2
Orbitals shown: IBOs (isovalue ±0.05); tracked by overlap; canonical HOMO/LUMO optional
What this is: the minimum-energy path, not a trajectory; no temperature, no dynamics
Known artifacts/limitations: (fill in — e.g., residual ion–dipole wells, solvent limits)
Course mapping: McMurry chapter __ (user fills in)
```

## 8. Order of work and check-ins

1. First-run checks (§2) → report.
2. SN2 TS + IRC → report barrier, imaginary frequency, and a PNG of the profile.
   **Stop for the user's judgment before spending time on orbitals.**
3. IBOs + tracking → report with 3–4 still frames (reactant, before TS, TS, product).
4. Video + viewer → report with the link to the local page and the mp4.
5. Card + pipeline docs → report.

Ask the user only for decisions a chemist should make (is this picture faithful to what
the course teaches? kJ/mol or kcal/mol first?). Do not ask him to run commands.

## 9. Known risks and what to do

- **TS search fails repeatedly:** switch to NEB-TS from the reactant and product
  complexes; if still failing, tighten the guess from a relaxed scan of the C–O distance.
- **Gas-phase-like double well with solvent on:** report it; it is expected behavior for
  this reaction class and will appear in the card and the paper as a teaching point.
- **Orbital flicker:** always a tracking or phase problem, never "the chemistry". Fix the
  matching; never hand-edit frames.
- **File size:** per-frame meshes must stay small enough that one reaction page loads in
  under 10 MB total. Decimate harder or drop frames before compromising on that.
- **Blender automation fights back:** render from a Python script in headless mode;
  do not depend on the user clicking anything in Blender.
