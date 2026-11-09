# Environment checks (Omen, Windows 11) — 2026-11-08

Machine: AMD Ryzen 7 9800X3D (8 cores / 16 threads), 32 GB RAM.

| Component | Version / location | Status |
|---|---|---|
| ORCA | 6.1.1, `F:\ORCA\` | works |
| Microsoft MPI | 10.1.12498.52 (runtime + SDK) | installed, parallel runs complete |
| orca_plot | `F:\ORCA\orca_plot.exe` | works (interactive menu, cube output) |
| Blender | 5.2.0, `C:\Program Files\Blender Foundation\Blender 5.2\` | headless works (2.0 s start); Molecular Nodes not installed |
| Python | 3.15.0, project venv `.venv\` | rdkit, numpy, scipy, scikit-image, trimesh, cclib, fast-simplification available; pymeshlab and open3d have no 3.15 wheels |
| xtb | — | not installed (optional) |

## Timings

Water, B3LYP/def2-SVP single point (final energy −76.321119597 Eh in every run):

| Cores | Wall time |
|---|---|
| 1 | 6.6 s (first ever launch: 34 s) |
| 8 | 40 s (first launch: 66 s) |

SN2 TS guess (HO⁻···CH₃···Br, charge −1), working level
B3LYP-D3(BJ)/def2-SVP, RIJCOSX, TightSCF, CPCM/SMD water, energy + gradient
(final energy −2689.356855334 Eh in every run):

| Cores | Wall time |
|---|---|
| 1 | 21.3 s |
| 4 | 34.3 s |
| 8 | 53.1 s |

Parallel runs give identical energies but are slower than serial. During an
8-process run the processes were mostly idle (total CPU use 18–20 %), and the
gradient step took the same time on 1 and 8 cores (9.1 s vs 9.3 s). Windows
Defender is not running; McAfee real-time protection was active.

### After turning off McAfee real-time scanning

Same SN2 job, same input, identical energies:

| Cores | Wall time | SCF | Gradient |
|---|---|---|---|
| 1 | 8.7 s | 5.5 s | 2.7 s |
| 4 | 4.6 s | 2.8 s | 1.0 s |
| 8 | 5.3 s | 3.2 s | 0.9 s |

McAfee real-time scanning was the cause: it made serial jobs about 2.4× slower
and parallel jobs about 10× slower. With it off, 4 cores give a 1.9× speedup on
this 7-atom system; 8 cores give no further gain because the molecule is too
small to split eight ways. Use `nprocs 4` for SN2; re-test 8 for the larger
reactions.

Note: a permanent fix that keeps protection on is a McAfee scan exclusion for
`F:\ORCA\`, `C:\Program Files\Microsoft MPI\Bin\` and the project folder.
