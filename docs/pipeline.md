# How to add a reaction

The pipeline turns one settings file into: transition state, IRC, energies, tracked orbitals,
web viewer data and a rendered video. Every step reads the previous step's files from disk, so
steps can be rerun one at a time, and a finished ORCA job with an unchanged input is reused
rather than recomputed.

## 1. What to provide

Make a folder `reactions/<reaction_id>/00_inputs/` with two files.

**`reaction.toml`** (copy `reactions/sn2_hydroxide_bromomethane/00_inputs/reaction.toml`).
TOML rule: put all plain settings before the first `[section]` line; anything after a
section header belongs to that section. Optionally add `card_notes.toml` (course mapping,
experimental comparison, extra limitations) for the hand-written parts of the card.

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
| `isovalue` | orbital isosurface value, e^1/2 bohr^-3/2 (0.08 for SN2) |
| `path_levels` | optional: levels (labels from `[[diagnostic_levels]]`) to recompute at every IRC frame |
| `default_curve` | energy curve shown first in the viewer and in `<id>.mp4`: `"working"` or one of `path_levels` (SN2: DLPNO-CCSD(T)) |
| `[key_atoms]` | 0-based indices in the TS guess: `nucleophile`, `center`, `leaving_group` |
| `[[reactants]]`, `[[products]]` | `name` + `smiles` of each separated species |
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
| `irc` | IRC in both directions from the TS, reading the TS Hessian | `02_irc/irc*.xyz`, `02_irc/irc.out` |
| `profile` | orders frames reactant → product, mass-weighted path length s, energies relative to separated reactants | `02_irc/energies.csv`, `02_irc/irc_full_trj.xyz`, `02_irc/energy_profile.png` |
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
