"""The pysisyphus patches (patches/pysisyphus/*.patch) state exactly what the runs changed.

1. They apply cleanly, in order, to a copy of the installed pysisyphus 1.0.0, whose files match its
   install record (no other modification).
2. The patched classes give the same path as the overrides the runs use (pipeline/true_irc2.GSmw and
   EulerPCmw) on a model potential: same points to 1e-10 bohr, same energies.
3. Replay of the recorded first step of tail_gs_reactant (2026-11-13, before the fixes): pysisyphus
   1.0.0 reproduces the recorded run and fails in micro-step 1 (lowest curvature driven to about
   −0.62 au by BFGS, Newton's λ search diverges); the patched micro-step on the same data keeps the
   curvature small and returns a point on the hypersphere.

Run: .venv-pysis/Scripts/python.exe -X utf8 tests/test_pysisyphus_patches.py   (plain runner: the
pysisyphus environment has no pytest, and its package list stays as frozen)
"""
from __future__ import annotations

import base64
import csv
import hashlib
import importlib.metadata as md
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATCHES = sorted((ROOT / "patches/pysisyphus").glob("0*.patch"))
CASES = ROOT / "tests/pysis_patch_cases.py"


def _installed():
    d = md.distribution("pysisyphus")
    return d, Path(d.locate_file(""))


def test_installed_package_unmodified():
    d, root = _installed()
    assert d.version == "1.0.0"
    bad = []
    for path, h, _ in csv.reader(open(root / f"pysisyphus-{d.version}.dist-info/RECORD", encoding="utf-8")):
        if h and path.startswith("pysisyphus/"):
            algo, want = h.split("=", 1)
            got = base64.urlsafe_b64encode(hashlib.new(algo, (root / path).read_bytes()).digest()).rstrip(b"=").decode()
            if got != want:
                bad.append(path)
    assert not bad, bad


def make_patched(base: Path) -> Path:
    """A copy of the installed package with all patches applied in order (git apply)."""
    _, root = _installed()
    dst = base / "pysis_patched"
    dst.mkdir()
    shutil.copytree(root / "pysisyphus", dst / "pysisyphus", ignore=shutil.ignore_patterns("__pycache__"))
    assert len(PATCHES) == 8
    for p in PATCHES:
        r = subprocess.run(["git", "apply", "--check", "-p1", str(p)], cwd=dst, capture_output=True, text=True)
        assert r.returncode == 0, f"{p.name}: {r.stderr}"
        r = subprocess.run(["git", "apply", "-p1", str(p)], cwd=dst, capture_output=True, text=True)
        assert r.returncode == 0, f"{p.name}: {r.stderr}"
    return dst


def run_case(case, pythonpath=None):
    env = dict(os.environ)
    if pythonpath:
        env["PYTHONPATH"] = str(pythonpath)
    r = subprocess.run([sys.executable, "-X", "utf8", "-W", "ignore", str(CASES), case], capture_output=True,
                       text=True, env=env, cwd=ROOT)
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT ")), None)
    assert line, f"{case}: no result\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
    return json.loads(line[7:])


def test_patched_import_is_the_copy(patched):
    r = subprocess.run([sys.executable, "-c", "import importlib; print(importlib.import_module('pysisyphus.irc.GonzalezSchlegel').__file__)"],
                       capture_output=True, text=True, env={**os.environ, "PYTHONPATH": str(patched)})
    assert str(patched) in r.stdout, r.stdout + r.stderr[-500:]


def test_patched_equals_override(patched, kind):
    a = run_case(f"toy_{kind}_override")
    b = run_case(f"toy_{kind}_patched", patched)
    assert a["n"] == b["n"] and a["n"] >= 5
    assert a["converged"] and b["converged"]
    assert np.max(np.abs(np.array(a["coords"]) - np.array(b["coords"]))) < 1e-10
    assert np.max(np.abs(np.array(a["energies"]) - np.array(b["energies"]))) < 1e-12
    if kind == "epc":
        assert a["fallbacks"] == b["fallbacks"]


def test_replay_recorded_gs_step(patched):
    o = run_case("replay_original")
    assert o["micro0_lands_on"] == "g0003", o           # the replay reproduces the recorded run
    assert o["lowest_after_update"] < -0.5               # BFGS: about −0.62 au
    assert o["micro1"].startswith("RuntimeError")        # Newton's λ search fails (1.0.0 then calls it converged)
    p = run_case("replay_patched", patched)
    assert p["on_sphere_err"] < 1e-9, p                  # bracketed λ: on the hypersphere (internal space)
    assert p["lowest_after_update"] > -0.1, p            # Bofill + projection: the curvature stays small


if __name__ == "__main__":
    results = []

    def check(name, fn, *a):
        try:
            fn(*a)
            results.append((name, "ok"))
        except AssertionError as e:
            results.append((name, f"FAILED: {e}"))
        print(f"{name}: {results[-1][1]}", flush=True)

    check("installed package unmodified", test_installed_package_unmodified)
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        try:
            patched = make_patched(Path(tmp))
            results.append(("8 patches apply in order", "ok"))
            print("8 patches apply in order: ok", flush=True)
        except AssertionError as e:
            print(f"patches apply: FAILED: {e}")
            sys.exit(1)
        check("patched import is the copy", test_patched_import_is_the_copy, patched)
        check("Gonzalez-Schlegel: patched = override", test_patched_equals_override, patched, "gs")
        check("EulerPC: patched = override", test_patched_equals_override, patched, "epc")
        check("replay of the recorded GS step", test_replay_recorded_gs_step, patched)
    print(f"{sum(r == 'ok' for _, r in results)} of {len(results)} passed")
    sys.exit(0 if all(r == "ok" for _, r in results) else 1)
