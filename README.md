# Computed Mechanism Animations

Textbook organic reaction mechanisms, computed rather than drawn. For each reaction:

1. the transition state and intrinsic reaction coordinate (IRC) are computed with DFT in
   implicit solvent (ORCA 6.1.1, SMD);
2. intrinsic bond orbitals (IBOs) are computed at every frame of the path and tracked from
   frame to frame, so the orbitals that the curved arrows describe can be watched as they change;
3. the result is shown as an interactive web page, where the energy diagram is the scrubber,
   and as a rendered video;
4. each reaction has a card stating the level of theory, solvent model, barrier and known
   limitations.

Pilot reaction: SN2, HO⁻ + CH₃Br → CH₃OH + Br⁻ (`reactions/sn2_hydroxide_bromomethane/`).

## Viewing

The site is static files only. Locally:

```
.venv\Scripts\python.exe -m http.server 8765
```

then open <http://localhost:8765/web/index.html>. On GitHub Pages, serve the repository root.

## Layout

```
pipeline/     Python package: geometry → ORCA → IRC → IBO → meshes → web / video
reactions/    one folder per reaction (inputs, TS, IRC, orbitals, meshes, render, card.md, run_log.md)
web/          landing page and the data-driven three.js viewer
docs/         pipeline instructions (docs/pipeline.md) and protocol notes
```

Large raw outputs (wavefunction files, cube files, full ORCA outputs) are archived separately
(Zenodo) rather than committed; decimated meshes, energies and inputs are in the repository.

## How to add a reaction

See [docs/pipeline.md](docs/pipeline.md).

## License and citation

Code: MIT (see `LICENSE`). Content (structures, energies, meshes, images, videos, text):
CC-BY 4.0. Please cite this repository and ORCA:

- F. Neese, "Software update: The ORCA program system—Version 6.0", *WIREs Comput. Mol. Sci.*
  2025, 15, e70019.
- G. Knizia, "Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical
  Concepts", *J. Chem. Theory Comput.* 2013, 9, 4834–4843 (IBOs).

Project lead: John Picione.
