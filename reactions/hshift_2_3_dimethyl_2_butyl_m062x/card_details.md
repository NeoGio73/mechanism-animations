# Card details: 1,2-Hydride shift, 2,3-dimethyl-2-butyl cation (degenerate)

## Energies at every level (same M06-2X-D3(0)/ma-def2-SVP geometries, SMD dichloromethane)

| Level | ΔE‡ from the optimized cation (`dimethylbutyl_cation`, lowest in free energy) | ΔE‡ from the reactant-end minimum (optimized from the start of the path) | Reaction energy ΔrE | ΔG‡ from the optimized cation (`dimethylbutyl_cation`, lowest in free energy) (298 K) | ΔG‡ from the reactant-end minimum (optimized from the start of the path) (298 K) | Reaction free energy ΔrG (298 K) |
|---|---|---|---|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (geometries) | 15.3 kJ/mol (3.7 kcal/mol) | 15.5 kJ/mol (3.7 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) | 14.1 kJ/mol (3.4 kcal/mol) | 13.0 kJ/mol (3.1 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP | 19.2 kJ/mol (4.6 kcal/mol) | 19.4 kJ/mol (4.6 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) | 18.0 kJ/mol (4.3 kcal/mol) | 16.8 kJ/mol (4.0 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP (default curve) | 18.2 kJ/mol (4.3 kcal/mol) | 18.4 kJ/mol (4.4 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) | 17.0 kJ/mol (4.1 kcal/mol) | 15.8 kJ/mol (3.8 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP | 28.5 kJ/mol (6.8 kcal/mol) | 28.7 kJ/mol (6.9 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) | 27.3 kJ/mol (6.5 kcal/mol) | 26.2 kJ/mol (6.3 kcal/mol) | 0.0 kJ/mol (0.0 kcal/mol) |

ΔG adds the M06-2X-D3(0)/ma-def2-SVP thermal corrections (298.15 K, quasi-RRHO) to each level's electronic energies; there is one molecule throughout, so no standard-state correction applies.

## Comparison with the measurement at 135 K

ΔG‡ at 135 K, from the lowest conformer of the cation (`01_ts/rotamer_zero.json`; same thermal model as at 298.15 K). Measured: 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol), from ¹³C NMR line broadening in a superacid solution (SO₂ClF/SO₂F₂), 135 K.

| Level | ΔG‡ (135 K) | Computed − measured |
|---|---|---|
| M06-2X-D3(0)/ma-def2-SVP (geometries) | 13.2 kJ/mol (3.2 kcal/mol) | +0.2 kJ/mol (+0.04 kcal/mol) |
| M06-2X-D3(0)/ma-def2-TZVP | 17.1 kJ/mol (4.1 kcal/mol) | +4.1 kJ/mol (+1.0 kcal/mol) |
| DLPNO-CCSD(T)/ma-def2-TZVP (default curve) | 16.1 kJ/mol (3.8 kcal/mol) | +3.1 kJ/mol (+0.7 kcal/mol) |
| B3LYP-D3(BJ)/ma-def2-TZVP | 26.4 kJ/mol (6.3 kcal/mol) | +13.4 kJ/mol (+3.2 kcal/mol) |

## The measured value: source and derivation

Experiment: Saunders and Kates 1978 (M. Saunders and M. R. Kates, J. Am. Chem. Soc. 1978, 100, 7082–7083, doi:10.1021/ja00490a060), compound I, the degenerate 1,2-hydride shift in the 2,3-dimethyl-2-butyl cation: ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol) at −138 °C (135 K); k = 3.1 × 10⁷ s⁻¹ at 135 K, 4.6 × 10⁷ s⁻¹ at −128 °C (145 K) and 1.4 × 10⁸ s⁻¹ at −119 °C (154 K). Solvent 50:50 SO₂ClF/SO₂F₂ (superacid medium); ¹³C NMR line broadening at 67.9 MHz in the fast-exchange limit; temperatures ±2 °C. The authors did not fit an Arrhenius equation, so there is no measured activation energy.

The comparison is made at 135 K, not 298 K, and in a different medium: the experiment is in a superacid medium, the calculation in SMD dichloromethane (and in the gas phase). Computed ΔG‡ is measured from the lowest cation conformer. At 135 K, RT is only 1.12 kJ/mol, so a small error in the barrier is a large error in rate. The working-level ΔG‡ at 135 K in dichloromethane, 13.2 kJ/mol (3.2 kcal/mol), is 0.2 kJ/mol above the reported 13.0 kJ/mol, a factor of about 1.2 in rate; the DLPNO-CCSD(T)/ma-def2-TZVP value, 16.1 kJ/mol (3.8 kcal/mol), is 3.1 kJ/mol (0.7 kcal/mol) above it, a factor of about 16. The rate ratio therefore looks much worse than the agreement in energy.

The authors infer from their data that hyperconjugation is weaker at the transition state than in the cation.

At M06-2X/ma-def2-SVP (DefGrid3; the working level for the two cation reactions) the top of the path is flat within 0.06 kJ/mol (0.014 kcal/mol): the transition state between the open cation and the hydrogen-bridged structure lies 0.010 kJ/mol above the bridged structure in dichloromethane and 0.056 kJ/mol in the gas phase. The bridged well is shallower than the zero-point energy of the transfer mode (0.38 kJ/mol in dichloromethane, 0.54 kJ/mol in the gas phase). The production transition state (115.2i cm⁻¹) gives ΔG‡ = 13.2 kJ/mol (3.2 kcal/mol) at 135 K in dichloromethane at this level, against the measured 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol).

## Path and orbitals

- IRC frames: 543; reaction coordinate s from −13.29 to 13.97 bohr·amu½ (TS at s = 0).
- Path quality: bond-changing region IRC frames 55–487, median step there 0.050 bohr·amu½; steps longer than 2 × that median: frames 271→272 (s +0.28 → +0.44), 0.162 bohr·amu½ = 3.2 × median.
- Orbitals: IBOs (IAO/IBO localization, occupied valence), isovalue ±0.08 e^1/2 bohr^−3/2 (orbital amplitude); tracked frame to frame by maximum overlap with fixed phase (lowest frame-to-frame overlap of each shown orbital: σ(C3–H) → σ(C2–H) (three-centre at the TS) 1.000 (IRC frames 271→272); for comparison, lowest overlap of the shown orbitals in Diels–Alder 0.714, E2 1.000, Secondary hydride shift 1.000, Cyanide + acetone 1.000, SN2 1.000); canonical HOMO/LUMO optional in the viewer.

## Technical notes

- M06-2X overestimates this barrier: ΔG‡ = 18.0 kJ/mol (4.3 kcal/mol) at M06-2X-D3(0)/ma-def2-TZVP versus 17.0 kJ/mol (4.1 kcal/mol) at DLPNO-CCSD(T)/ma-def2-TZVP (1 kJ/mol higher).
- Implicit solvent only: no explicit dichloromethane molecules.
- The product end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = 13.97 bohr·amu½ (IRC step 129), after the bonds had finished changing; the product-end minimum was optimized from that last point.
- The reactant end of the path was truncated where the energy had flattened: the IRC was stopped by the IRC stopping rule (docs/protocol.md) at s = −13.29 bohr·amu½ (IRC step 264), after the bonds had finished changing; the reactant-end minimum was optimized from that last point.
- Solvent: the experiments on these cations are in superacid media, which the SMD solvent model has no parameters for. Dichloromethane (SMD) is used as a polar, non-nucleophilic stand-in; the gas phase is computed alongside for comparison.
- Tunnelling by the migrating hydrogen is not included. It may matter at the low temperatures of the NMR experiments (135–154 K). The Wigner factor at 135 K is about 1.06 (imaginary frequency 115.2i cm⁻¹ at the working level in DCM, hν‡/k_BT = 1.23). It is small because the top of this barrier is very flat; a flat-topped, wide barrier gives little tunnelling, so the true correction is probably also small, but the Wigner formula is only a rough indicator.
- Only half of the path is computed: the IRC runs from the open cation (H on C3) through the TS to the C₂-symmetric bridged structure. The second half, bridged → TS′ → open cation (H on C2), is its mirror image under the C2 ⇄ C3 swap; those frames are labelled as mirror images in the viewer and are not IRC points. The bridged side is limited by SCF noise: the bridged structure lies only about 0.01 kJ/mol below the TS in DCM, below the SCF noise floor of 0.034 kJ/mol.
- At this level (M06-2X/ma-def2-SVP, SMD dichloromethane) the bridged structure and the transition state are not distinguishable within numerical precision: the bridged structure lies 0.012 kJ/mol below the TS, below the SCF noise floor of 0.034 kJ/mol and the optimizer's energy tolerance. The top of the barrier is therefore a flat plateau, not a sharp maximum.

## Files

- `02_irc/energies.csv`: frame, s (bohr·amu½), E (Eh, kJ/mol, kcal/mol), key distances (Å)
- `02_irc/irc_full_trj.xyz`: IRC frames, reactant → product
- `02_irc/energy_profile.png`: energy profile
- `03_orbitals/tracking_log.md`: orbital matching along the path
- `05_render/hshift_2_3_dimethyl_2_butyl.mp4`: rendered video, energy inset at DLPNO-CCSD(T)/ma-def2-TZVP (`hshift_2_3_dimethyl_2_butyl_<level>.mp4`: the same video with the other energy curves)
- `web/viewer/index.html?r=hshift_2_3_dimethyl_2_butyl`: interactive viewer
- `run_log.md`: every ORCA job and rerun

Content CC-BY 4.0; code MIT. Computed with ORCA 6.1.1.
