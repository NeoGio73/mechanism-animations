# How to add a reaction

The pipeline turns one settings file into: transition state, IRC, energies, tracked orbitals,
web viewer data and a rendered video. Every step reads the previous step's files from disk, so
steps can be rerun one at a time, and a finished ORCA job with an unchanged input is reused
rather than recomputed.

## 1. What to provide

Make a folder `reactions/<reaction_id>/00_inputs/` with two files.

**Standard level of theory (decided by J. Picione, 2026-11-09; details in `docs/protocol.md`):**
geometries, TS, IRC and orbitals at B3LYP-D3(BJ)/**ma-def2-SVP** in SMD solvent (diffuse
functions on the heavy atoms, because most of these systems are anions); energies along the path
recomputed at DLPNO-CCSD(T)/ma-def2-TZVP, and that curve shown by default. Every new reaction
uses `method = "B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF"`,
`sp_method = "B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF"`, `path_levels` with the
DLPNO-CCSD(T) level, and `default_curve = "DLPNO-CCSD(T)_ma-def2-TZVP"`. (The SN2 pilot was
computed at def2-SVP and is due to be rerun at ma-def2-SVP; see `docs/tasks.md`.)

**`reaction.toml`** (copy `reactions/sn2_hydroxide_bromomethane/00_inputs/reaction.toml`).
TOML rule: put all plain settings before the first `[section]` line; anything after a
section header belongs to that section. Optionally add `card_notes.toml` (course mapping,
experimental comparison, extra limitations) for the hand-written parts of the card. Cite the
primary paper with its DOI and add every source to `docs/references.md`. Leave
`experimental_confirmed = false` (the card then says "Value to be confirmed against the primary
paper by J. Picione") until he has checked the value.

| Field | Meaning |
|---|---|
| `id`, `title`, `equation` | names used on the web page and in the video |
| `charge`, `mult` | total charge and spin multiplicity of the reacting system |
| `solvent` | SMD solvent name as ORCA spells it (e.g. `water`, `methanol`, `dmso`) |
| `method` | working level for geometries, TS, IRC and orbitals |
| `sp_method` | larger-basis level for the reported barrier |
| `nprocs`, `maxcore_mb` | cores per ORCA job and memory per core (MB) |
| `ts_guess` | file name of the transition-state guess (XYZ, Å) |
| `irc_scale_displ_sd` | optional IRC step; smaller gives more frames (0.05 gave 60 for SN2) |
| `irc_stop_dE_kJ`, `irc_stop_steps` | optional: IRC stopping rule, default 0.05 kJ/mol per step for 20 steps after the last bond change |
| `isovalue` | orbital isosurface value, e^1/2 bohr^-3/2 (0.08 for SN2) |
| `path_levels` | optional: levels (labels from `[[diagnostic_levels]]`) to recompute at every IRC frame |
| `video_timing` | optional: `"displacement"` (default) moves the video at constant summed atomic displacement between IRC points, so bond changes get the time and a drifting-fragment tail passes quickly; `"path_length"` moves at constant s (SN2 pilot) |
| `default_curve` | energy curve shown first in the viewer and in `<id>.mp4`: `"working"` or one of `path_levels` (SN2: DLPNO-CCSD(T)) |
| `ts_preopt` | optional: relax the guess with every changing bond frozen before the TS search |
| `view_axis`, `view_axis_label` | optional: the "down the axis" camera looks along this bond instead (E2: Newman view along Cβ–Cα) |
| `[key_atoms]` | 0-based indices in the TS guess: `nucleophile`, `center`, `leaving_group` (also set the camera: nucleophile/base left, leaving group right) |
| `[[key_bonds]]` | optional: every bond that changes, `atoms`, `role` (`forms`, `breaks`, `shortens` for C–C → C=C, `lengthens` for C=C → C–C), `label`; without it, nucleophile–center forms and center–leaving group breaks |
| `[ts_dihedrals]` | optional: dihedrals to report at the TS (label = four atom indices) |
| `[ts_angles]` | optional: angles to report at the TS (label = three atom indices, middle atom at the vertex; cyanide + acetone: the Bürgi–Dunitz angle C(N)···C=O); reported, never constrained |
| `[ts_mode]` | optional: `dominant_atom` must move most in the imaginary mode; `reject_atom` moving along `reject_axis` marks the TS as the SN2 competitor and rejects it; `synchronous` (labels of bonds that must form together) rejects the TS unless their lengths agree within `synchronous_max_spread_A` (0.05 Å) and their stretch rates in the mode within `synchronous_min_rate_ratio` (0.8) |
| `[ts_expected]` | optional: expected ranges per bond (Å), dihedral or angle (°); outside gives a warning |
| `[ts_required]` | optional: same form, but outside the range rejects the TS (Diels–Alder: both forming bonds 2.15–2.35 Å) |
| `[[reactants]]`, `[[products]]` | `name` + `smiles` of each separated species; optional `start_xyz` (file in `00_inputs`) fixes the starting conformer. A species that optimizes onto a saddle point (e.g. an eclipsed methyl) is displaced 0.1 Å along the imaginary mode and re-optimized once with the same keywords (`<name>_min`) |
| `[[references]]` | optional: species optimized and reported but not part of the barrier (Diels–Alder: s-trans butadiene, the ground-state conformer) |
| `[[diagnostic_levels]]` | optional extra single points reported next to the protocol numbers |

`key_atoms` drive the sanity checks (is the imaginary mode the right motion, which IRC end
is the reactant side), the dashed partial bonds, and the distance readouts. For reactions that
are not substitutions, use the two bonds that change most: `nucleophile`–`center` is the bond
that forms, `center`–`leaving_group` the bond that breaks.

**TS guess (`ts_guess.xyz`)**: a reasonable starting structure with the forming and breaking
bonds stretched (for SN2: collinear Nu···C···LG, C–Nu ≈ 2.0 Å, C–LG ≈ 2.3 Å, planar CH₃).

## 2. Commands

All commands run from the repository folder with the project Python:

```
.venv\Scripts\python.exe -m pipeline reactions\<reaction_id> <step> [<step> ...]
```

| Step | What it does | Output |
|---|---|---|
| `fragments` | optimizes each separated reactant and product (Opt Freq) | `01_ts/fragments/`, `01_ts/fragments.json` |
| `ts` | OptTS + Freq from the guess; checks for exactly one imaginary frequency and that it is the bond-forming/breaking motion | `01_ts/ts/`, `01_ts/ts.json` |
| `irc` | IRC from the TS, reading the TS Hessian: one job per direction (`irc_f`, `irc_b`), each ended by the IRC stopping rule (docs/protocol.md), then joined (backward reversed + TS + forward). Older one-job runs (`irc.out`) are reused; one stopped by decision is joined from `02_irc/irc_stop.json` | `02_irc/irc_IRC_Full_trj.xyz`, `02_irc/irc_IRC_F.xyz`, `02_irc/irc_IRC_B.xyz`, `02_irc/irc_path.json` |
| `profile` | orders frames reactant → product, mass-weighted path length s, energies relative to separated reactants; flags any point within two steps of the TS that lies above it (kept in the data, left out of viewer and video, stated on the card) | `02_irc/energies.csv`, `02_irc/irc_full_trj.xyz`, `02_irc/energy_profile.png` |
| `endpoints` | optimizes both IRC ends into the reactant and product complexes (Opt Freq) | `02_irc/endpoints/`, `02_irc/endpoints.json` |
| `singlepoints` | `sp_method` single points on all stationary points | `01_ts/sp_tzvp.json` |
| `barriers` | barrier and reaction energies at both levels | `02_irc/barriers.json` |
| `diagnostics` | extra levels from `[[diagnostic_levels]]`, plus a free-energy barrier | `01_ts/diagnostics.md` |
| `path_singlepoints` | energy at every IRC frame at the levels in `path_levels` (labels from `[[diagnostic_levels]]`); two frames at a time | `02_irc/energies_<level>.csv` |
| `profile_levels` | one plot of the energy along the path at every level | `02_irc/energy_profile_levels.png` |
| `orbital_sps` | single point with IBO localization at every IRC frame | `03_orbitals/frames/` |
| `orbital_track` | matches orbitals frame to frame by overlap, fixes signs, finds the orbitals that change | `03_orbitals/tracking.json`, `03_orbitals/tracking_log.md` |
| `orbital_cubes` | cube files (gzipped) for the changing IBOs and canonical HOMO/LUMO | `03_orbitals/cubes/` |
| `meshes` | isosurfaces (± `isovalue`), decimated, packed for the web | `04_meshes/manifest.json`, `04_meshes/meshes.bin` |
| `video_preview` | four half-resolution test frames | `05_render/preview.mp4` |
| `video` | full 1920 × 1080, 30 fps Blender render, then `video_composite` | `05_render/png/` |
| `video_composite` | energy inset and labels added to the rendered frames; one video per energy curve (no re-render, about 15 s each) | `05_render/<id>.mp4` (default curve), `<id>_<level>.mp4`, `stills*/` |
| `card` | writes the reaction card from the computed numbers | `card.md` |

A full run for a new reaction:

```
.venv\Scripts\python.exe -m pipeline reactions\<id> fragments ts irc profile endpoints singlepoints barriers diagnostics path_singlepoints profile_levels orbital_sps orbital_track orbital_cubes meshes video card
```

Then add one line for the reaction to `REACTIONS` in `web/index.html`. The viewer page is
`web/viewer/index.html?r=<id>`; it needs no other changes.

## 3. What to check by eye

- `01_ts/ts.json`: `n_imaginary` must be 1, and `mode_check.central_atom_share` should be large
  (the central atom carries most of the motion).
- `02_irc/energy_profile.png`: the path should run from a reactant-side minimum through the TS to
  a product-side minimum.
- `03_orbitals/tracking_log.md`: the "changes" column should name the orbitals that the curved
  arrows describe. Frames flagged below 0.7 overlap are expected near the TS; flicker anywhere
  else means a tracking problem (never edit frames by hand).
- `run_log.md`: every ORCA job with its wall time, and a note for every rerun.

## 4. Software

ORCA 6.1.1 (`F:\ORCA`, change with the `ORCA_DIR` environment variable), Microsoft MPI 10.1,
Python 3.15 venv (`rdkit numpy scipy matplotlib scikit-image trimesh fast-simplification
cclib imageio-ffmpeg pillow`), Blender 5.2 (EEVEE, headless; no add-ons needed).
On this machine, McAfee real-time scanning made ORCA up to 10 × slower; see
`docs/environment_checks.md`.
