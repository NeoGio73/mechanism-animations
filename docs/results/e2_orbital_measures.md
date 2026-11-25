# E2: when do the three electron pairs move? Two orbital measures and the distances (2026-11-25)

s in bohr·amu½; the TS is frame 132 (s = 0); 0.05 bohr·amu½ per frame near the TS. All values are from files
already on disk; no new calculation.

**Measure 1: Mulliken populations** (`pipeline/orbital_order.py`, `03_orbitals/orbital_order.json`). These are the
Mulliken atomic populations of each IBO in the full ma-def2-SVP basis (C·SC summed over each atom's basis
functions), not IAO charges. They can be slightly negative: C1 in the bromide lone pair at the product end is
−0.008. Rate = sum over atoms of |dp/ds|.

**Measure 2: orbital centroid** (`pipeline/tracking_analysis.py centroids`, `03_orbitals/centroids.json`). This
is the centroid ⟨r⟩ of each IBO computed from its cube file (weighted by φ², normalized on the grid). It does not
depend on how charge is shared between basis functions. One coordinate per orbital: f, the projection of the
centroid on the line from the start site to the end site, in that frame's geometry (f = 0 at the start site,
1 at the end site). The sites are:

| Orbital | Start site | End site | f at the reactant end | f at the product end |
|---|---|---|---|---|
| O lone pair → σ(O–H) | O | midpoint of O–H | 0.08 | 1.09 |
| σ(Cβ–H) → π(C=C) | midpoint of Cβ–H | midpoint of Cα–Cβ | −0.10 | 0.96 |
| σ(Cα–Br) → Br lone pair | midpoint of Cα–Br | Br | 0.11 | 0.88 |

"Half-change" is the frame where f, or the populations, have covered half of their total change.

| Orbital | Mulliken: fastest | Mulliken: half-change | Mulliken: fast stretch | Centroid: fastest | Centroid: half-change | Centroid: fast stretch |
|---|---|---|---|---|---|---|
| O lone pair → σ(O–H) | +0.25 | +0.30 | −0.30 to +0.55 | **−0.05** | **−0.20** | −0.40 to +0.45 |
| σ(Cβ–H) → π(C=C) | **−0.10** | **−0.10** | −0.40 to +0.75 | +0.20 | +0.40 | −0.35 to +0.65 |
| σ(Cα–Br) → Br lone pair | +0.40 | +0.25 | −0.45 to +1.00 | +0.50 | +0.80 | −0.40 to +2.00 |

The fast stretch is where the rate is above half its peak.

| Distance | Fastest change |
|---|---|
| O–H | 0.00 (TS) |
| Cβ–H | +0.20 |
| Cα–Cβ | +0.55 |
| Cα–Br | +1.60 (slow over a long stretch) |

**Reading:**
- **Mulliken:** the C–H pair goes first, and the O lone pair and the C–Br pair follow almost together; the last
  two swap order between the fastest and half-change points.
- **Centroid:** the order is O lone pair, then C–H → π, then C–Br, by both the fastest point and the half-change
  point. That is the same order as the distances. The three fastest points span 0.55 bohr·amu½ (−0.05 to +0.50)
  and the three half-change points span 1.0 bohr·amu½ (−0.20 to +0.80). The C–Br pair's centroid has the
  broadest fast stretch, out to +2.0.
- **So the two orbital measures disagree on the order.** Only the distances and the centroid agree.

## The C–Br orbital along the path

| s (bohr·amu½) | Frame | Mulliken: share on Br | Share on Cα | Mulliken: fraction of the total change done | Centroid: fraction of the total change done |
|---|---|---|---|---|---|
| reactant end | 0 | 65.1 % | 35.3 % | 0 % | 0 % |
| 0.00 (TS) | 132 | 78.9 % | 20.3 % | 40 % | 30 % |
| +0.40 | 140 | 84.1 % | 15.0 % | 55 % | 40 % |
| +1.00 | 152 | 90.2 % | 8.8 % | 72 % | 55 % |
| +1.60 | 164 | 93.8 % | 5.5 % | 83 % | 67 % |
| product end | 531 | 99.8 % | −0.8 % | 100 % | 100 % |

## A third measure (IAO charges of single orbitals): not available (2026-11-25)

ORCA's localization program (orca_loc, ORCA 6.1.1) was tested on frame 132 at print levels 3, 5, 10 and 100. It
prints IAO partial charges only per atom, summed over all localized orbitals (with ORCA's warning that these mean
something only when every occupied orbital is included; the core orbitals were not localized), and per-orbital
compositions only as Mulliken populations. It cannot print the IAO charges of single orbitals. As J. Picione
decided, nothing further was run (no option B). **The order in which the three electron pairs move is not resolved
by these calculations:** the two available orbital measures disagree, and the card describes the atoms only.
