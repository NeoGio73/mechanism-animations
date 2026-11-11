# Tertiary → tertiary 1,2-hydride shift (2,3-dimethyl-2-butyl cation): M06-2X geometry check

Computed 2026-11-11 (queue task `hshift_fallback`, pipeline/level_test.py). Numbers from
`summary.json` (this folder's parent), `../../low_T_135K.json`.

## Main finding: at M06-2X the bridged ion is a minimum, not a transition state

The M06-2X-D3(0)/ma-def2-SVP (SMD dichloromethane, DefGrid3) TS search started from the B3LYP TS. The first
M06-2X Hessian there had one imaginary mode, 214i cm⁻¹. The optimization converged to the symmetric
hydrogen-bridged structure with **all frequencies real** (lowest 64 cm⁻¹; the bridging H moves mainly in
a real mode near 1165 cm⁻¹). At B3LYP-D3(BJ)/ma-def2-SVP the same symmetric structure is the TS (563i cm⁻¹).

So M06-2X describes the shift as two steps through a shallow bridged intermediate: open → TS → bridged →
TS → open. The two transition states on either side of the bridged ion have not been located. The energy of the
bridged ion is therefore a lower bound to the M06-2X barrier, not the barrier itself.

| Geometry | C–H (both) | C–C (shift bond) | Character |
|---|---|---|---|
| B3LYP-D3(BJ)/ma-def2-SVP | 1.333 Å | 1.414 Å | TS, 563i cm⁻¹ |
| M06-2X-D3(0)/ma-def2-SVP | 1.320 / 1.321 Å | 1.408 Å | minimum, all real |

## Energies, SMD dichloromethane, relative to the open tertiary cation (lowest conformer)

ΔG at 135 K (the temperature of the measurement) and 298 K. Thermal corrections are from the
geometry's own level: B3LYP for the B3LYP geometry (a TS, with the transfer mode removed) and M06-2X for
the M06-2X geometry (a minimum, with the transfer mode kept as a vibration). This is why the thermal correction is −7.6 kJ/mol
on the B3LYP geometry and −2.9 kJ/mol on the M06-2X geometry at 135 K.

| Geometry | Energy level | ΔE | ΔG (135 K) | ΔG (298 K) | ΔG(135 K) − exp |
|---|---|---|---|---|---|
| B3LYP (TS) | B3LYP/ma-def2-SVP (own) | 26.0 (6.2) | 18.4 (4.4) | 17.9 (4.3) | +5.4 |
| B3LYP (TS) | B3LYP/ma-def2-TZVP | 29.5 (7.0) | 21.9 (5.2) | 21.4 (5.1) | +8.9 |
| B3LYP (TS) | M06-2X/ma-def2-TZVP | 19.6 (4.7) | 12.0 (2.9) | 11.5 (2.8) | −1.0 |
| B3LYP (TS) | DLPNO-CCSD(T)/ma-def2-TZVP | 18.4 (4.4) | 10.8 (2.6) | 10.3 (2.5) | −2.2 |
| M06-2X (bridged minimum) | M06-2X/ma-def2-SVP (own) | 15.5 (3.7) | 12.7 (3.0) | 12.2 (2.9) | (lower bound) |
| M06-2X (bridged minimum) | M06-2X/ma-def2-TZVP | 19.2 (4.6) | 16.3 (3.9) | 15.9 (3.8) | (lower bound) |
| M06-2X (bridged minimum) | DLPNO-CCSD(T)/ma-def2-TZVP | 18.4 (4.4) | 15.5 (3.7) | 15.1 (3.6) | (lower bound) |

kJ/mol (kcal/mol). Experiment (Saunders and Kates 1978, superacid SO₂ClF/SO₂F₂, 135 K): ΔG‡ = 13.0 ± 0.4 kJ/mol
(3.1 ± 0.1 kcal/mol). On the M06-2X geometry the "ΔG" is the free energy of the bridged intermediate, not a barrier.
The minimum on the M06-2X side is the M06-2X re-optimization of the B3LYP lowest conformer.

## Gas phase, B3LYP geometry (TS), 135 K

| Energy level | ΔE | ΔG (135 K) | − exp |
|---|---|---|---|
| B3LYP/ma-def2-SVP (own) | 26.6 (6.3) | 19.3 (4.6) | +6.3 |
| B3LYP/ma-def2-TZVP | 30.3 (7.2) | 23.0 (5.5) | +10.0 |
| M06-2X/ma-def2-TZVP | 21.1 (5.0) | 13.8 (3.3) | +0.8 |
| DLPNO-CCSD(T)/ma-def2-TZVP | 20.1 (4.8) | 12.8 (3.1) | −0.2 |

B3LYP gas-phase thermal correction −7.3 kJ/mol at 135 K. No M06-2X geometries were computed in the gas phase.

## Reading

- The coupled-cluster energy of the bridged structure is the same on both geometries: 18.37 and 18.44 kJ/mol
  above the open cation (DCM). The two geometries differ by only 0.013 Å in C–H, and the CC surface is flat there.
- The ΔG difference between the two geometries (10.8 against 15.5 kJ/mol at CC, 135 K) comes almost entirely
  from the thermal correction, i.e. from whether the bridged structure is treated as a TS (B3LYP) or a minimum
  (M06-2X). It does not come from the electronic energy.
- Whether the bridged form is a minimum or a TS at the coupled-cluster level is not known from single points.
  The flanking M06-2X transition states have not been located, so no M06-2X barrier is available yet.
- `summary.json` lists `imag_cm: []` for the M06-2X "ts"; the task did not flag it. Found on reading the output.
