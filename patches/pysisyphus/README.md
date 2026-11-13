# Modifications to pysisyphus 1.0.0

J. Picione, 2026-11-13: every modification to pysisyphus is saved as a patch file, with one line on
what it changes and one on the test that shows it works, so the paper can state exactly what was
modified.

**Base:** pysisyphus 1.0.0 (GPL-3), installed unmodified in `.venv-pysis`. All 1014 package files match
their install record hashes (`tests/test_pysisyphus_patches.py`). The original files that are
patched, by SHA-256:

| File | SHA-256 |
|---|---|
| `pysisyphus/irc/GonzalezSchlegel.py` | `cac9e897a56b427e7a7c982a67a76c5cb06a809da78603650caeb0fd84d09cb2` |
| `pysisyphus/irc/IRC.py` | `5eee231ccafa2138d00cd9d36f068cfcb4e9d71dba5d544b9b61ff917b9f2855` |
| `pysisyphus/irc/EulerPC.py` | `1a1f6f98509c197a1f858eb411248bdd933222784ca58c8c0fe9e4ab2459a7a6` |

**How the runs use them:** the runs do not edit the package. They use the same changes as subclass
overrides in `pipeline/true_irc2.py`:
- `GSmw` carries patches 0001–0004.
- `OwnStop` carries 0005.
- `EulerPCmw` carries 0006–0008.

The patch files state the same changes at the source level. The test shows that the patched package
and the overrides give identical paths.

## The patches (apply in order, `git apply -p1` inside `site-packages`)

| Patch | What it changes | Test that shows it works |
|---|---|---|
| 0001-gs-bofill-update | Gonzalez–Schlegel micro-steps update the Hessian with Bofill instead of BFGS; no update for a zero step. | Replay of the recorded first step of tail_gs_reactant: BFGS turns the lowest curvature from −0.0152 to −0.621 au, as in the failed run; Bofill (with 0002) gives −0.0027 au. |
| 0002-gs-project-translation-rotation | Gonzalez–Schlegel removes the three translations and three rotations of the current frame from the gradient, the update differences and the Hessian used for the step. | Same replay. Patched package and run code agree on a model potential to 1e-10 bohr. |
| 0003-gs-bracketed-lambda | λ found by a bracketed search (brentq) below the lowest eigenvalue, where the root is unique, instead of Newton's method from a guess. | Same replay: Newton fails on the recorded data; the bracketed search lands on the hypersphere (radius error 3e-16). The SN2 validation path had no failed search in 939 gradients. |
| 0004-gs-failed-search-is-failure | A constrained search that does not converge raises `IRCIntegratorFailure`; 1.0.0 marked the IRC converged. | tail_gs_reactant before the fixes "converged" after one step, far from the minimum. Now such a side is recorded as incomplete, never as a result. |
| 0005-irc-external-stop-rule | `IRC(stop_func=...)`: a rule checked after every step ends the direction. The runs pass the project stopping rule and switch pysisyphus's own stops off. | SN2 validation path: both sides ended by this rule, ends −0.005 and +0.322 kJ/mol from the optimized minima. |
| 0006-eulerpc-mass-weighted-step-length | EulerPC measures the step in mass-weighted coordinates (bohr·amu½), not plain bohr. | Measured steps on every EulerPC run: median 0.1000 and 0.0500 bohr·amu½ for steps of 0.10 and 0.05. With 1.0.0, a "step" of 0.10 gave 0.35. |
| 0007-eulerpc-fixed-euler-substep | Fixed conversion factor 2 for the Euler predictor substep (step × 2 / max_pred_steps). | With max_pred_steps = 5000 the predictor reaches the full step on every EulerPC run. |
| 0008-eulerpc-corrector-fallback | If the mBS corrector returns less than half a step (oscillation, down to zero), the predictor point is kept and counted. | Before the fix both EulerPC tail runs sent all-NaN coordinates to ORCA. After it the four Hessian-every-5 tail runs had 24, 40, 35 and 19 fallbacks and no NaN. |

The first lines of each `.patch` file repeat its "What" and "Test" in full.

**Which runs used which patches:**
- Production paths (all reactions, `pipeline/true_irc_prod.py`) and the SN2 validation path: Gonzalez–Schlegel with 0001–0005.
- EulerPC (0005–0008) was used for these tests: the tail tests, the isotope control (h_eulerpc, d4_eulerpc), and the coupled-cluster maximum on h_eulerpc.

## Settings (parameters, not code changes)

**Production:**
- Step and first displacement: `step_length = 0.05`, `displ = "length"`, `displ_length = 0.05`.
- Hessian: `hessian_init = "calc"`, which is ORCA's exact TS Hessian, then updates only (`hessian_recalc = None`).
- pysisyphus's own stops switched off: `rms_grad_thresh = 1e-12`, `energy_thresh = −1`.
- `max_cycles = 3000`.

**EulerPC tests:**
- `max_pred_steps = 5000` (default 500).
- Corrector `"mbs"`.
- `hessian_recalc = 5` for the Hessian-every-5 option.

## Bookkeeping without effect on the path

The overrides also do the following. None of it changes a path:
- record each micro-step's λ and lowest eigenvalue (`micro_log`);
- record cycles that hit the micro-cycle limit;
- record the corrector fallbacks;
- stop with an error before ORCA is sent non-finite coordinates.

## Reproduce

```
.venv-pysis/Scripts/python.exe -X utf8 patches/pysisyphus/make_patches.py         # writes the 8 patch files
.venv-pysis/Scripts/python.exe -X utf8 tests/test_pysisyphus_patches.py           # 6 checks
```

The tests do four things:
1. They check the installed package against its record hashes.
2. They apply the 8 patches to a copy of it.
3. They compare the patched Gonzalez–Schlegel and EulerPC with the run code on a model potential (Morse pairs, four atoms).
4. They replay the recorded failed step with the original and the patched code.

Result on 2026-11-13: 6 of 6 passed.
