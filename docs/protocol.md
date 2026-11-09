# Computational protocol (draft methods text)

## Standard for the reaction set (decided by J. Picione, 2026-11-09)

From E2 onward, every reaction uses one working level for geometries, the TS, the IRC and the
orbitals: **B3LYP-D3(BJ)/ma-def2-SVP**, RIJCOSX (def2/J), TightSCF, SMD solvent. The minimally
augmented basis adds diffuse functions to the heavy atoms, which the anions in most of these
reactions need. Barrier single points use B3LYP-D3(BJ)/ma-def2-TZVP. Energies at every IRC frame are
recomputed at DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO), and that curve is the default in the viewer and
the main video; M06-2X-D3(0)/ma-def2-TZVP is also computed along the path for comparison.

Why: in the SN2 pilot the def2-SVP working level (no diffuse functions) put the TS 10.9 kJ/mol
(2.6 kcal/mol) below separated reactants. The barrier needed rescuing with coupled-cluster single
points, which gave ΔG‡ = 95.6 kJ/mol (22.9 kcal/mol) against 95.0 kJ/mol (22.7 kcal/mol) measured.
The SN2 pilot below was computed at def2-SVP. It will be rerun at ma-def2-SVP so the whole set shares
one working level (`docs/tasks.md`); the def2-SVP results stay until then.

### Points next to the TS that lie above it (decided by J. Picione, 2026-11-09)

ORCA's first IRC step leaves the TS along a straight line. On a strongly curved path (a hydrogen
transfer) that point can lie slightly above the TS. Any IRC point within two steps of s = 0 that
lies above the TS energy is flagged by `profile`, kept in the data files, left out of the viewer
path and the video so the displayed maximum is the TS, and stated on the card with its energy.
E2: the first reactant-side point is 1.3 kJ/mol (0.3 kcal/mol) above the TS.

### IRC stopping rule (decided by J. Picione, 2026-11-09)

Each IRC direction is ended when the energy has changed by less than **0.05 kJ/mol
(0.012 kcal/mol) per step for 20 consecutive steps** after the last bond-forming or
bond-breaking event. A forming or breaking key bond counts as bonded at ≤ 1.25 × the sum of the
covalent radii and as not bonded at ≥ 1.5 ×. The last event is the first point from which every
such bond keeps its final state to the end. Past that point the fragments only drift apart
through the solvent, and ORCA's gradient criteria can take hundreds of steps to be met. The pipeline
runs the two directions as separate jobs, checks the rule every 30 s, cuts the path at the step
where the rule was met, and logs the stop in `run_log.md`. The end complex is then optimized from
that last point (Opt Freq, TightOpt if an imaginary mode remains), and the card states the s value
where the path was cut.

First use: E2 (2026-11-09). The product-side IRC was stopped by hand at step 115, before the rule
existed. Products had formed by step 56 (O–H 0.98 Å, C=C 1.34 Å). Each step changed the energy by
less than 0.05 kJ/mol only from step 79 on, and the rule would have stopped that direction at
step 98.

## SN2 pilot

Settings as actually used for the SN2 pilot (`reactions/sn2_hydroxide_bromomethane`). Every
job, its wall time and every rerun are in that folder's `run_log.md`.

## Electronic structure

All calculations used ORCA 6.1.1. The working level, used for all geometries, the transition
state (TS), the intrinsic reaction coordinate (IRC) and the orbitals, was
B3LYP-D3(BJ)/def2-SVP with the RIJCOSX approximation (def2/J auxiliary basis) and TightSCF
convergence, in implicit water: C-PCM with SMD parameters (`! CPCM(water)` with
`%cpcm smd true SMDsolvent "WATER" end`). ORCA's final energies include the SMD CDS term
and the D3(BJ) correction.

Separated reactants and products were built from SMILES (RDKit ETKDG + MMFF) and optimized
with frequency calculations at the working level (bromide: single point).

The TS was optimized (`OptTS`) from a collinear guess (C–O 2.0 Å, C–Br 2.3 Å, planar CH₃,
C–O–H 105°) with the exact Hessian computed at the start and every 5 steps
(`Calc_Hess true`, `Recalc_Hess 5`), followed by a frequency calculation. The TS was
accepted only with exactly one imaginary frequency whose normal mode is dominated by the
central carbon moving along the O···Br axis (share of the mode on C and alignment with the
axis are recorded in `01_ts/ts.json`).

The IRC was followed in both directions from the TS using the TS Hessian (`InitHess read`),
with `Scale_Displ_SD 0.05` to obtain 60 frames (the default step gave 37, with only 7 across
the steep central region). Frames are ordered reactant → product, and the mass-weighted path
length s (bohr·amu^1/2) is accumulated from consecutive frames, with s = 0 at the TS. Both IRC
end geometries were reoptimized (Opt Freq) into the reactant and product complexes; a
minimum showing a small imaginary mode was displaced 0.1 Å along that mode and reoptimized
with TightOpt.

Barrier heights were recomputed by single points on the working-level geometries at
B3LYP-D3(BJ)/def2-TZVP (protocol), B3LYP-D3(BJ)/ma-def2-TZVP, M06-2X-D3(0)/ma-def2-TZVP and
DLPNO-CCSD(T)/ma-def2-TZVP (TightPNO, AutoAux), all in SMD water. Free-energy barriers add
the working-level thermal corrections (298.15 K, rigid rotor/harmonic oscillator as
implemented in ORCA) and the 1 atm → 1 M standard-state correction (RT ln 24.47 =
7.92 kJ/mol per species). The M06-2X and DLPNO-CCSD(T) energies were also computed at every IRC
frame, giving energy curves at those levels along the working-level path. The viewer and the
main video show the DLPNO-CCSD(T) curve by default; the working-level and M06-2X curves are
selectable in the viewer and rendered as separate videos.

## Orbitals

At every IRC frame a working-level single point was followed by intrinsic bond orbital
localization of the occupied valence orbitals (`%loc LocMet IAOIBO end`; ORCA localizes
the orbitals above the frozen core). Coefficients and the AO overlap matrix were exported
with `orca_2json`.

Tracking: orbitals of neighbouring frames k and k+1 are matched by maximizing |O| with
O = C_k^T S C_{k+1}, S = (S_k + S_{k+1})/2 (Hungarian assignment), propagating from the
reactant end. Each match's sign is fixed so that the overlap is positive, so phases never
flip between frames. Matches with |O| < 0.7 are flagged (none occurred for SN2; the two
orbitals shown never fell below 0.99). Each tracked orbital is labelled from its Mulliken
atomic populations (atoms above 0.15 e) at both ends of the path. The orbitals whose set of
atoms changes between the ends are the ones shown by default: for SN2, an O lone pair →
σ(C–O), and σ(C–Br) → a Br lone pair. The overall sign of each orbital is chosen once, at
the first frame, so that its larger lobe is the + phase. Canonical HOMO and LUMO are tracked
for sign only (their identity is fixed by energy order).

Cube files (70 points per axis, orca_plot) were contoured at ±0.08 e^1/2 bohr^-3/2 by marching
cubes (scikit-image) and decimated to at most 900 triangles per lobe (fast-simplification).
Vertex positions are stored as 16-bit integers in a common box, with 16-bit triangle indices:
about 15 kB per orbital per frame, 3.7 MB for the whole SN2 page.

## Rendering

Web: three.js r160, static files; ball-and-stick model linearly interpolated between IRC
frames, orbitals shown per nearest frame. The forming and breaking bonds are drawn solid
below 1.15 × their equilibrium length (taken from the path ends) and dashed up to 1.6 ×.

Video: Blender 5.2 (EEVEE, headless Python), 1920 × 1080, 30 fps, 20 s. Playback moves at
constant speed in s for the SN2 pilot; from E2 on it moves at constant summed atomic
displacement between IRC points (decided by J. Picione, 2026-11-09), so the bond-changing region
gets the time and a long drifting-fragment tail passes quickly; the interactive viewer keeps every
frame. Start and stop are eased; orbital surfaces cross-fade between
neighbouring IRC frames. The energy inset and labels are composited in Python (Pillow,
matplotlib) and encoded with H.264 (ffmpeg, CRF 18).
