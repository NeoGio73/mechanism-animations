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
E2: the first reactant-side point is 1.3 kJ/mol (0.3 kcal/mol) above the TS. On E2 the M06-2X/ma-def2-TZVP
curve peaks one frame from the working-level TS (that same point, 0.6 kJ/mol (0.1 kcal/mol) above the TS
frame) because the single points sit on a path optimized at another level, while the
DLPNO-CCSD(T)/ma-def2-TZVP maximum falls on the TS frame.

### Implicit solvent and anionic products (noted by J. Picione, 2026-11-10)

SMD represents water as a continuum, with no explicit water molecules. It therefore misses the specific hydrogen bonds that water donates to a localized anion, and it under-stabilizes an alkoxide product. A computed step that forms an alkoxide (cyanide + acetone) is therefore likely too endergonic. The error comes from the solvent model, not the electronic-structure level: for cyanide + acetone the alkoxide lies +34.4 kJ/mol (+8.2 kcal/mol) above separated reactants in free energy at the working level and +36.3 kJ/mol (+8.7 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (298 K, 1 M). No explicit water molecules are added; the card states the limitation.

### Result: the 3-methyl-2-butyl cation has no classical secondary minimum to shift from (accepted by J. Picione, 2026-11-10)

The first choice for the carbocation 1,2-hydride shift was 3-methyl-2-butyl → 2-methyl-2-butyl. The rule set beforehand was that the secondary cation must be a true minimum (all real frequencies, not hydrogen-bridged) with a barrier to the shift; otherwise switch to the degenerate shift in the 2,3-dimethyl-2-butyl cation.

The secondary cation was optimized at the working level, in SMD dichloromethane and in the gas phase, from three rotamers about C2–C3 (`reactions/hshift_3_methyl_2_butyl`, run log there). Two findings came out:

- **No barrier for the shift itself.** With the C3–H in line with the empty p orbital on C2, the optimization runs downhill without a barrier to the tertiary 2-methyl-2-butyl cation, in both phases.
- **The minimum that does exist is partly methyl-bridged, not classical.** With the C3–H turned 64–78° away from the p axis, the optimization finds a true minimum: all frequencies real, the hydrogen still on C3. But one C3–CH₃ bond is stretched to 1.60–1.61 Å in DCM and 1.66–1.67 Å in the gas phase. That methyl leans toward C2, with a C2–C3–CH₃ angle of 92–95° in DCM and 82–85° in the gas phase. This is strong C–C hyperconjugation. Only rotation about C2–C3 keeps the hydrogen from shifting.

The tertiary cation lies 56.9 kJ/mol (13.6 kcal/mol) below this minimum in DCM and 54.2 kJ/mol (13.0 kcal/mol) below in the gas phase (ΔE).

This is a chemical result, not a failure of the calculation. A secondary cation next to a tertiary centre does not survive as a classical ion at this level; it is held only by its conformation. The degenerate tertiary → tertiary shift in the 2,3-dimethyl-2-butyl cation is used instead.

The literature on the simplest secondary case, the 2-butyl cation, points the same way. Experimentally, its hydrogen scrambling is not frozen out even at −190 °C in the solid (barrier below 2.4 kcal/mol, 10 kJ/mol), and isotope experiments indicate two forms in rapid equilibrium, a hydrogen-bridged form about 0.4 kcal/mol (1.7 kJ/mol) more stable (Saunders and Jiménez-Vázquez, *Chem. Rev.* 1991, 91, 375–397). Ab initio, Carneiro, Schleyer, Koch and Raghavachari (*J. Am. Chem. Soc.* 1990, 112, 4064–4066) found for the 2-butyl cation:

- the open classical form is a transition structure (one imaginary frequency), 2.8 kcal/mol (11.7 kJ/mol) above the methyl-bridged form at MP4(sdtq)/6-311G**//MP2(fu)/6-31G** with zero-point energy;
- the methyl-bridged and hydrogen-bridged forms are both minima, within 0.8 kcal/mol (3.3 kJ/mol); with zero-point energy the hydrogen-bridged form is lower by 0.3 kcal/mol (1.3 kJ/mol).

The partially methyl-bridged minimum found here for the 3-methyl-2-butyl cation is the same kind of structure reported for the 2-butyl cation. The bridging angle depends on the method in the same way: Hartree–Fock and B3LYP give a wider angle, MP2 and M06-2X a narrower one.

| Cation | Level | C⁺–C–CH₃ angle of the methyl-bridged form |
|---|---|---|
| 2-butyl (Carneiro et al. 1990) | HF/6-31G* | 102.4° |
| 2-butyl (Carneiro et al. 1990) | MP2(fu)/6-31G* | 77.4° |
| 2-butyl (Carneiro et al. 1990) | MP2/6-31G** | 76.9° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | B3LYP-D3(BJ)/ma-def2-SVP | 81.6° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | M06-2X/ma-def2-SVP | 71.3° |
| 3-methyl-2-butyl, anti minimum, gas phase (this work) | RI-MP2/ma-def2-SVP | 71.6° |

One difference: in the 2-butyl cation the hydrogen-bridged form is a minimum; in the 3-methyl-2-butyl cation it is not, because the hydride shift from it leads downhill to the tertiary 2-methyl-2-butyl cation. Full references in docs/references.md.

Two follow-ups were asked for:
- whether the partial methyl bridging survives at M06-2X and at a wavefunction level;
- whether a C2–C3 rotation transition state connects this minimum to the tertiary cation.

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
the working-level thermal corrections (298.15 K, 1 atm, ideal gas, rigid rotor/harmonic
oscillator with Grimme's quasi-RRHO vibrational entropy, reference frequency 100 cm⁻¹, ORCA's
default) and the 1 atm → 1 M standard-state correction (RT ln 24.47 =
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
atomic populations (atoms above 0.15 e) at both ends of the path; a two-atom orbital with
less than 0.10 of its population in s-type basis functions is labelled π, otherwise σ, and a
multiple bond whose orbitals are all above 0.10 is labelled as bent bonds (IBO localization
sometimes mixes σ and π into two bent bonds, not necessarily equivalent). The orbitals whose set of
atoms changes between the ends are the ones shown by default: for SN2, an O lone pair →
σ(C–O), and σ(C–Br) → a Br lone pair. The overall sign of each orbital is chosen once, at
the first frame, so that its larger lobe is the + phase. Canonical HOMO and LUMO are tracked
for sign only (their identity is fixed by energy order).

Abrupt switch of the localization (Diels–Alder, tested 2026-11-10): the three shown orbitals
have their lowest frame-to-frame overlap at IRC step 59→60, one to two steps after the TS
(0.827–0.838). IBO single points at ¼ and ½ of that step (Cartesian interpolation between the two
frames, same input as the frames; `03_orbitals/step_test_59_60/`) show that the localized
orbitals switch abruptly within the first quarter of the step: overlap 0.79 from frame 59 to the ¼
point, then 0.996–0.997 on to frame 60. More IRC frames would therefore not spread the change, and
none were added. SN2 and E2 show no such switch: every frame-to-frame overlap of their shown
orbitals is at least 0.99 (lowest 0.993 and 0.990); nor does the cyanide + acetone addition
(lowest 0.992).

Display rotation of bent bonds (decided by J. Picione, 2026-11-10; first use: Diels–Alder): where
IBO localization gives an isolated double bond as two bent bonds (for ethylene at the Diels–Alder
reactant end, two non-equivalent ones, s-character 0.18 and 0.24), the pair is rotated within itself,
for display only, into the σ + π pair that spans the same two orbitals, and the π orbital is shown.
This applies from the reactant end at every frame where the computed orbital is more than 0.5° from
that split (π = the pair's combination with the least s-character), so the rotated frames join the
unrotated ones smoothly (Diels–Alder: IRC frames 0–47, join overlap 0.9998); a rotation inside an
occupied pair changes no energy or density, and the unrotated orbitals stay in the data files
(`03_orbitals/display_rotation.json` records the frames, angles and overlaps).

Arrow direction in pericyclic reactions (decided by J. Picione, 2026-11-10): in a pericyclic
reaction the clockwise and counter-clockwise curved-arrow drawings are mirror images of the
same mechanism, and a mirror-symmetric transition state favours neither. The direction shown
is the one the orbital localization settled on, not a physical preference; `reaction.toml`
lists that assignment, and where `[mirror]` is set the viewer can show the mirror-image set,
made by reflecting the computed orbital surfaces through the TS mirror plane at every frame
(no new calculation).

Cube files (70 points per axis, orca_plot) were contoured at ±0.08 e^1/2 bohr^-3/2 by marching
cubes (scikit-image) and decimated to at most 900 triangles per lobe (fast-simplification).
Vertex positions are stored as 16-bit integers in a common box, with 16-bit triangle indices:
about 15 kB per orbital per frame, 3.7 MB for the whole SN2 page.

## Rendering

Web: three.js r160, static files; ball-and-stick model linearly interpolated between IRC
frames, orbitals shown per nearest frame. The forming and breaking bonds are drawn solid
below 1.15 × their equilibrium length (taken from the path ends) and dashed up to 1.6 ×.

Video: Blender 5.2 (EEVEE, headless Python), 1920 × 1080, 30 fps, 20 s. EEVEE renders on the GPU
(NVIDIA GeForce RTX 5080, OpenGL; E2: 0.33 s per frame), and every render log records the engine,
graphics backend and GPU (`RENDER_DEVICE` line in `05_render/blender.log`). Playback moves at
constant speed in s for the SN2 pilot. From E2 on, playback time is allotted in proportion to the
summed change in the key bond lengths (for E2: Cβ–H, O–H, Cα–Br, Cα–Cβ), with a floor of about
0.3 s per bohr·amu^1/2 of path and the tail after the last bond change compressed to 1 s (decided
by J. Picione, 2026-11-09), so the bond-changing region gets the time and a long drifting-fragment
tail passes quickly. Timing by summed atomic displacement is kept as a selectable alternative
(`video_timing` in `reaction.toml`); the interactive viewer keeps every frame. Start and stop are eased; orbital surfaces cross-fade between
neighbouring IRC frames. The energy inset and labels are composited in Python (Pillow,
matplotlib) and encoded with H.264 (ffmpeg, CRF 18). Rendered frames carry a SHA-256 fingerprint
of everything that decides their look (orbital meshes, geometries, timeline, render settings and
the render script; not the energy curves, which are drawn by the compositor). The renderer
deletes all earlier frames when the fingerprint changes, and the compositor refuses frames whose
fingerprint does not match (added 2026-11-10 after a Diels–Alder render reused frames made before
the orbitals were localized; SN2 and E2 were checked and were not affected).
