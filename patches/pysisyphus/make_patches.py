"""Write the patch files against pysisyphus 1.0.0 (J. Picione, 2026-11-13: every modification to
pysisyphus as a patch file, so the paper can state exactly what was modified).

The runs themselves use the same changes as subclass overrides in pipeline/true_irc2.py (GSmw,
EulerPCmw, OwnStop); the installed package is unmodified (checked against its RECORD hashes).
tests/test_pysisyphus_patches.py applies these patches to a copy of the installed package and checks
that the patched classes give the same numbers as the overrides.

Each patch is one change, applied in order (0001 … 0008); each step below edits the source as left by
the previous one. Usage: .venv-pysis/Scripts/python.exe patches/pysisyphus/make_patches.py
"""
from __future__ import annotations

import difflib
from pathlib import Path

import pysisyphus

HERE = Path(__file__).resolve().parent
PKG = Path(pysisyphus.__file__).resolve().parent
assert pysisyphus.__version__ == "1.0.0", pysisyphus.__version__


def sub(text: str, old: str, new: str) -> str:
    assert text.count(old) == 1, old[:80]
    return text.replace(old, new)


GS, IRCF, EPC = "irc/GonzalezSchlegel.py", "irc/IRC.py", "irc/EulerPC.py"
src = {f: (PKG / f).read_text(encoding="utf-8") for f in (GS, IRCF, EPC)}
steps = []


def step(name: str, f: str, what: str, test: str, edits: list[tuple[str, str]]):
    old = src[f]
    new = old
    for a, b in edits:
        new = sub(new, a, b)
    src[f] = new
    steps.append((name, f, old, new, what, test))


# 0001 ------------------------------------------------------------------------------------------
step("0001-gs-bofill-update", GS,
     "Gonzalez–Schlegel micro-steps update the Hessian with Bofill instead of BFGS (BFGS is unstable "
     "where a curvature is negative, as in the flat tails); the update is skipped for a zero step.",
     "Replay of the recorded first micro-step of tail_gs_reactant (before the fixes): BFGS turned the lowest "
     "curvature from −0.015 to −0.62 au; with Bofill it stays near −0.015 au "
     "(tests/test_pysisyphus_patches.py::test_replay_recorded_gs_step).",
     [("from pysisyphus.optimizers.hessian_updates import bfgs_update\n",
       "from pysisyphus.optimizers.hessian_updates import bofill_update\n"),
      ("""        # Or update Hessian
        else:
            dH, _ = bfgs_update(self.mw_hessian, coords_diff, gradient_diff)
            self.mw_hessian += dH
""",
       """        # Or update Hessian (Bofill: stable where a curvature is negative; no update for a zero step)
        elif np.linalg.norm(coords_diff) > 1e-10:
            dH, _ = bofill_update(self.mw_hessian, coords_diff, gradient_diff)
            self.mw_hessian += dH
""")])

# 0002 ------------------------------------------------------------------------------------------
step("0002-gs-project-translation-rotation", GS,
     "Gonzalez–Schlegel removes the three translations and three rotations of the current frame (about its "
     "centre of mass, mass-weighted) from the gradient, the step and gradient differences used in the update, "
     "and the Hessian used for the step.",
     "With the projection, the six zero modes stay at zero through the updates, so the lowest eigenvalue is a "
     "true curvature (replay test above); the validation path keeps the TS orientation (RMSD to the minima "
     "0.027 and 0.080 Å after alignment).",
     [("""class GonzalezSchlegel(IRC):""",
       """def trans_rot_projector(cart_bohr, masses):
    \"\"\"I − QQᵀ in mass-weighted Cartesians of the current frame; Q = orthonormal basis of the three
    translations and three rotations about the centre of mass.\"\"\"
    x = np.asarray(cart_bohr).reshape(-1, 3)
    m = np.asarray(masses)
    r = x - (m[:, None] * x).sum(0) / m.sum()
    sq = np.sqrt(m)[:, None]
    vecs = []
    for k in range(3):
        e = np.zeros(3)
        e[k] = 1.0
        vecs.append((sq * e).ravel())
        vecs.append((sq * np.cross(e, r)).ravel())
    U, s, _ = np.linalg.svd(np.array(vecs).T, full_matrices=False)
    Q = U[:, s > 1e-6 * s.max()]
    return np.eye(x.size) - Q @ Q.T


class GonzalezSchlegel(IRC):"""),
      ("""        # Calculate gradient at current coordinates
        gradient = self.mw_gradient
""",
       """        # Calculate gradient at current coordinates; translations and rotations projected out
        P = trans_rot_projector(self.geometry.cart_coords, self.geometry.masses)
        gradient = P @ self.mw_gradient
"""),
      ("""        gradient_diff = gradient - self.prev_grad
        coords_diff = self.mw_coords - self.prev_coords
""",
       """        gradient_diff = gradient - P @ self.prev_grad
        coords_diff = P @ (self.mw_coords - self.prev_coords)
"""),
      ("""        eigvals, eigvecs = np.linalg.eigh(self.mw_hessian)
""",
       """        H = P @ self.mw_hessian @ P
        H = 0.5 * (H + H.T)
        eigvals, eigvecs = np.linalg.eigh(H)
""")])

# 0003 ------------------------------------------------------------------------------------------
step("0003-gs-bracketed-lambda", GS,
     "Gonzalez–Schlegel finds λ by a bracketed search (scipy brentq) below the lowest eigenvalue, where the "
     "constraint function rises monotonically from −r² to +∞ so the root is unique, instead of Newton's "
     "method from a guess; no root raises IRCIntegratorFailure.",
     "Replay of the recorded micro-step: Newton diverges (pysisyphus 1.0.0), the bracketed search returns a "
     "point on the hypersphere (|p| = r to 1e-10) (test_replay_recorded_gs_step); the tail tests and the "
     "validation path ran 939 gradients with no failed search.",
     [("from scipy.optimize import newton\n", "from scipy.optimize import brentq\n"),
      ("""class GonzalezSchlegel(IRC):""",
       """class IRCIntegratorFailure(Exception):
    \"\"\"A failed integration step; never reported as convergence.\"\"\"


class GonzalezSchlegel(IRC):"""),
      ("""        # Initial guess for λ.
        # λ must be smaller then the smallest eigenvector
        lambda_0 = big_eigvals[0]
        lambda_0 *= 1.5 if (lambda_0 < 0) else 0.5
        self.log(f"\\tSmallest eigenvalue is {big_eigvals[0]:.4f}, λ_0={lambda_0:.4f}.")
        # Find the root with scipy
        lambda_ = newton(on_sphere, lambda_0, maxiter=500)
        self.log(f"\\tDetermined λ={lambda_0:.4f} from Newtons method.")
""",
       """        # λ below the smallest eigenvalue: there on_sphere rises monotonically from −r² to +∞, so a
        # bracketed search finds the one root.
        b0 = big_eigvals[0]
        hi = b0 - 1e-10 * max(1.0, abs(b0))
        if not on_sphere(hi) > 0:
            raise IRCIntegratorFailure(f"GS micro-step {counter}: no root below the lowest eigenvalue ({b0:.3e})")
        lo, k = b0 - 1.0, 0
        while on_sphere(lo) >= 0:
            lo, k = b0 - 2 * (b0 - lo), k + 1
            if k > 200:
                raise IRCIntegratorFailure(f"GS micro-step {counter}: λ bracket not found")
        lambda_ = brentq(on_sphere, lo, hi, xtol=1e-14, maxiter=500)
        self.log(f"\\tSmallest eigenvalue is {b0:.4f}, λ={lambda_:.6f} from a bracketed search.")
""")])

# 0004 ------------------------------------------------------------------------------------------
step("0004-gs-failed-search-is-failure", GS,
     "A Gonzalez–Schlegel constrained search that does not converge raises IRCIntegratorFailure; "
     "pysisyphus 1.0.0 printed a message and marked the IRC converged.",
     "tail_gs_reactant before the fixes 'converged' after one step at RMS gradient 1e-4 Eh/bohr, far from the "
     "minimum (failed_before_fixes_1113); with the change such a side ends as a failure, never as a result "
     "(the driver records it as INCOMPLETE).",
     [("""            except RuntimeError:
                print("Constrained search did not converge!")
                self.converged = True
                return
""",
       """            except RuntimeError as err:
                raise IRCIntegratorFailure(
                    f"GS constrained search did not converge in micro-step {i}: {err}"
                ) from err
""")])

# 0005 ------------------------------------------------------------------------------------------
step("0005-irc-external-stop-rule", IRCF,
     "IRC accepts stop_func: called after every step; when it returns a true value the direction ends as "
     "converged. The runs pass the project stopping rule (|ΔE| < 0.05 kJ/mol per step for 20 steps after the "
     "last bond event, pipeline/steps.irc_stop_index) and switch pysisyphus's own stops off "
     "(rms_grad_thresh=1e-12, energy_thresh=−1).",
     "SN2 validation path: both sides ended by this rule, ends −0.005 and +0.322 kJ/mol from the optimized "
     "minima (true_irc2/final_final_gs_s005.md).",
     [("""        dump_every=5,
    ):
        \"\"\"Base class for IRC calculations.""",
       """        dump_every=5,
        stop_func=None,
    ):
        \"\"\"Base class for IRC calculations."""),
      ("""        self.max_cycles = max_cycles
        self.downhill = downhill
""",
       """        self.max_cycles = max_cycles
        self.stop_func = stop_func
        self.downhill = downhill
"""),
      ("""        \"\"\"Base class for IRC calculations.
""",
       """        \"\"\"Base class for IRC calculations.

        stop_func: optional callable(irc) -> bool, checked after every step; True ends the direction
        as converged ("Integrator indicated convergence!").
"""),
      ("""            last_energy = self.irc_energies[-2]
            this_energy = self.irc_energies[-1]
""",
       """            if self.stop_func is not None and self.stop_func(self):
                self.converged = True

            last_energy = self.irc_energies[-2]
            this_energy = self.irc_energies[-1]
""")])

# 0006 ------------------------------------------------------------------------------------------
step("0006-eulerpc-mass-weighted-step-length", EPC,
     "EulerPC measures the integration (step) length in mass-weighted coordinates, bohr·amu½; pysisyphus "
     "1.0.0 converts it back to plain Cartesian bohr.",
     "Measured mass-weighted steps on the EulerPC runs equal the set step (0.10 or 0.05 bohr·amu½) instead of "
     "0.18–0.35 bohr·amu½ for a 'step' given in plain bohr (protocol, attempt1_step_in_plain_bohr).",
     [("""            \"\"\"Returns length of integration path done in mass-weighted coordinates
            in un-mass-weighted coordinates.\"\"\"
            return np.linalg.norm((cur_mw_coords - init_mw_coords) / self.m_sqrt)
""",
       """            \"\"\"Returns length of integration path in mass-weighted coordinates
            (bohr·amu½).\"\"\"
            return np.linalg.norm(cur_mw_coords - init_mw_coords)
""")])

# 0007 ------------------------------------------------------------------------------------------
step("0007-eulerpc-fixed-euler-substep", EPC,
     "EulerPC uses a fixed conversion factor of 2 for the predictor Euler substep (substep = step × 2 / "
     "max_pred_steps; with max_pred_steps = 5000, step/2500), since the step is now mass-weighted (0006).",
     "With 5000 substeps the predictor reaches the full step on every EulerPC run (no 'predictor did not "
     "reach' stall, which the default 500 gave in the flat tail on 2026-11-12).",
     [("""    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)
""",
       """    def get_conv_fact(self, mw_grad, min_fact=2.0):
        # the step length is mass-weighted (see get_integration_length_func): no unit conversion
        return 2.0

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)
""")])

# 0008 ------------------------------------------------------------------------------------------
step("0008-eulerpc-corrector-fallback", EPC,
     "If the mBS corrector returns a point less than half the step length from the start (it does so on an "
     "oscillation of the interpolated surface, down to a zero step), EulerPC keeps the predictor point and "
     "counts the event in self.fallbacks; a zero step never reaches the Hessian update (which then filled "
     "the Hessian with NaN).",
     "Both EulerPC tail runs before the fix sent all-NaN coordinates to ORCA; after it the Hessian-every-5 "
     "runs completed with 35 and 19 fallbacks and no NaN (true_irc2/tails.md).",
     [("""        self.corr_func = corr_funcs[corr_func]
""",
       """        self.corr_func = corr_funcs[corr_func]
        self.fallbacks = list()
"""),
      ("""    def corrector_step(self, init_mw_coords, step_length, dwi):
        self.log("Corrector step using mBS integration")
""",
       """    def corrector_step(self, init_mw_coords, step_length, dwi):
        # self.mw_coords holds the predictor point when the corrector is called (see step)
        predictor = self.mw_coords.copy()
        x = self._mbs_corrector_step(init_mw_coords, step_length, dwi)
        got = float(np.linalg.norm(x - init_mw_coords))
        if got < 0.5 * step_length:
            self.fallbacks.append({"cycle": int(self.cur_cycle), "corrected_length": got})
            self.log(f"Corrected step {got:.4f} < half the step: predictor point kept")
            return predictor
        return x

    def _mbs_corrector_step(self, init_mw_coords, step_length, dwi):
        self.log("Corrector step using mBS integration")
""")])


# The "Test:" line of each patch, with the numbers measured by tests/test_pysisyphus_patches.py and
# on the recorded runs (2026-11-13). These replace the provisional lines given with each step above.
TESTS = {
    "0001-gs-bofill-update":
        "Replay of the recorded first step of tail_gs_reactant (before the fixes; tests/test_pysisyphus_patches.py): "
        "BFGS turns the lowest curvature from −0.0152 to −0.621 au, as in the run; Bofill (with 0002) gives −0.0027 au.",
    "0002-gs-project-translation-rotation":
        "Same replay: with the projection the six translation and rotation modes stay out of the step; patched "
        "package = run code on a model potential to 1e-10 bohr (test_patched_equals_override).",
    "0003-gs-bracketed-lambda":
        "Same replay: Newton's search fails on the BFGS Hessian (pysisyphus 1.0.0 then calls the IRC converged); the "
        "bracketed search returns a point on the hypersphere (radius error 3e-16); the SN2 validation path "
        "(939 gradients) had no failed search.",
    "0004-gs-failed-search-is-failure":
        "tail_gs_reactant before the fixes 'converged' after one step, far from the minimum (failed_before_fixes_1113); "
        "with the change such a side ends as a failure and the driver records it as incomplete, never as a result.",
    "0005-irc-external-stop-rule":
        "SN2 validation path: both sides ended by this rule, ends −0.005 and +0.322 kJ/mol from the optimized minima "
        "(true_irc2/final_final_gs_s005.md); patched package = run code on a model potential (stand-in rule).",
    "0006-eulerpc-mass-weighted-step-length":
        "Measured mass-weighted steps on every EulerPC run: median 0.1000 (0.10 runs) and 0.0500 (0.05 runs) "
        "bohr·amu½; with the 1.0.0 length a 'step' of 0.10 gave 0.35 bohr·amu½ (attempt1_step_in_plain_bohr).",
    "0007-eulerpc-fixed-euler-substep":
        "With max_pred_steps = 5000 the predictor reaches the full step on every EulerPC run (medians above); "
        "patched package = run code on a model potential to 1e-10 bohr.",
    "0008-eulerpc-corrector-fallback":
        "Before the fix both EulerPC tail runs sent all-NaN coordinates to ORCA (g0066, g0072); after it the four "
        "Hessian-every-5 tail runs had 24, 40, 35 and 19 fallbacks and no NaN (true_irc2/tails.md).",
}


if __name__ == "__main__":
    for old in HERE.glob("0*.patch"):
        old.unlink()
    for name, f, old, new, what, test in steps:
        test = TESTS[name]
        d = difflib.unified_diff(old.splitlines(keepends=True), new.splitlines(keepends=True),
                                 fromfile=f"a/pysisyphus/{f}", tofile=f"b/pysisyphus/{f}", n=3)
        body = "".join(d)
        head = f"pysisyphus 1.0.0: {name}\n\nWhat: {what}\nTest: {test}\n\n"
        (HERE / f"{name}.patch").write_text(head + body, encoding="utf-8", newline="\n")
        print(f"{name}.patch  ({f}, {body.count(chr(10))} lines)")
