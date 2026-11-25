# Zenodo preparation: what is not in the repository (2026-11-25 09:30)

Measured on the working copy with the IBO lanes still running (E2 and tertiary IBO frames partly done).
GB = 10⁹ bytes. Python environments (`.venv`, `.venv-pysis`) are not counted.

## In the repository (tracked)

450 MB in 10 297 files. A GitHub release archives only this, as one zip file.

## Not in the repository: 106.8 GB in 266 955 files

| Part | Size now | Worth archiving? |
|---|---|---|
| ORCA scratch left by stopped jobs (`*.tmp`, `*.tmp.0`–`.3` and similar) | 37.1 GB | No: intermediate integral files, no results |
| Six final reactions (`*_m062x`), kept results | 43.7 GB | Yes |
| Six final reactions: attempts, stopped, crashed and level-test folders | 5.3 GB | Optional (the run logs describe them) |
| Earlier folders (pilot level, superseded: `e2_ethoxide_2_bromopropane`, `sn2_hydroxide_bromomethane_masvp`, …) | 20.8 GB | Optional |

### Six final reactions, kept results, by kind (largest first)

| Kind | Size |
|---|---|
| IRC runs (`02_irc/true_irc`): `.gbw` 7.6, `.densities` 2.4, `.cpcm` 1.3, `.out` 0.3 | 11.6 GB |
| DLPNO-CCSD(T) path single points: `.gbw` 5.3, `.loc` 5.3, `.densities` 2.2, `.cpcm` 0.4, `.out` 0.2 | 13.4 GB |
| M06-2X/ma-def2-TZVP path single points: `.gbw` 4.1, `.densities` 2.8, `.cpcm` 0.5, `.out` 0.2 | 7.6 GB |
| B3LYP/ma-def2-TZVP path single points: `.gbw` 3.6, `.densities` 2.4, `.cpcm` 0.5, `.out` 0.2 | 6.7 GB |
| IBO frames so far (`03_orbitals/frames`; 539 of 2379 frames) | 2.8 GB |
| TS, stationary points, other | about 1.6 GB |

**When the IBO step ends:** IBO frames take about 5 MB per frame (E2 5.7 MB, tertiary 4.7 MB), so
2379 frames are about 12 GB, plus the orbital cubes (about 3–5 GB, already gzipped). The kept results of
the six final reactions then total about **55–57 GB uncompressed**.

**Compression (measured on E2 frame 134):** `.gbw` and `.loc` gzip to 27 % of their size; `.densities`
to 95 %. `.densities` and `.cpcm` can be regenerated from `.gbw` and the input.

| Archive set | Estimated size |
|---|---|
| Six final reactions, all kept files, gzipped | about 25–30 GB |
| Same without `.densities` and `.cpcm` | about 15–20 GB |
| Plus attempts and earlier folders | about 35–40 GB |

Zenodo's limit is 50 GB and 100 files per record, so the raw output must be bundled into archives
(for example one per reaction and kind, about 30 files).
