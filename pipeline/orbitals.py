"""Intrinsic bond orbitals (IBOs) along the IRC: compute, track, and plot.

Per frame: one ORCA single point with IAO/IBO localization (writes .loc), export of
coefficients + AO overlap matrix through orca_2json, and cube files through orca_plot.

Tracking: orbitals of neighbouring frames are matched by maximum |overlap|,
    O(k, k+1) = C_k^T S C_{k+1},  S = (S_k + S_{k+1}) / 2,
which is accurate for the small steps of an IRC. Signs are fixed so that every
matched pair has positive overlap, so phases never flip between frames.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from . import orca
from .geometry import Geometry


@dataclass
class FrameOrbitals:
    C: np.ndarray            # (nbf, nmo) MO coefficients
    occ: np.ndarray          # (nmo,)
    S: np.ndarray            # (nbf, nbf) AO overlap
    ao_atom: np.ndarray      # (nbf,) atom index of each AO
    symbols: list[str]

    def populations(self, idx: list[int]) -> np.ndarray:
        """Mulliken atomic populations (n_orb, n_atoms) of selected orbitals."""
        C = self.C[:, idx]
        pa = C * (self.S @ C)
        nat = len(self.symbols)
        return np.stack([pa[self.ao_atom == a].sum(axis=0) for a in range(nat)], axis=1)


def _env() -> dict:
    env = os.environ.copy()
    env["PATH"] = f"{orca.MSMPI_BIN};{orca.ORCA_DIR};" + env.get("PATH", "")
    return env


def export_json(wavefile: Path) -> FrameOrbitals:
    """orca_2json on a .gbw or .loc file; returns coefficients and overlap."""
    wd = wavefile.parent
    stem = wavefile.stem + ("_loc" if wavefile.suffix == ".loc" else "_can")
    gbw = wd / f"{stem}.gbw"
    shutil.copy(wavefile, gbw)
    (wd / f"{stem}.json.conf").write_text(json.dumps(
        {"MOCoefficients": True, "Basisset": True, "1elIntegrals": ["S"], "JSONFormats": ["json"]}))
    subprocess.run([str(orca.ORCA_DIR / "orca_2json.exe"), gbw.name, "-json"], cwd=wd,
                   env=_env(), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
    mol = json.loads((wd / f"{stem}.json").read_text())["Molecule"]
    mos = mol["MolecularOrbitals"]["MOs"]
    C = np.array([m["MOCoefficients"] for m in mos]).T
    occ = np.array([m["Occupancy"] for m in mos])
    S = np.array(mol["S-Matrix"])
    labels = mol["MolecularOrbitals"]["OrbitalLabels"]
    ao_atom = np.array([int("".join(ch for ch in lab.split()[0] if ch.isdigit())) for lab in labels])
    symbols = [a["ElementLabel"] for a in mol["Atoms"]]
    # keep the json small on disk: drop the big intermediate gbw copy
    gbw.unlink()
    np.savez_compressed(wd / f"{stem}.npz", C=C, occ=occ, S=S, ao_atom=ao_atom)
    (wd / f"{stem}.json").unlink()
    return FrameOrbitals(C, occ, S, ao_atom, symbols)


def load_npz(path: Path, symbols: list[str]) -> FrameOrbitals:
    z = np.load(path)
    return FrameOrbitals(z["C"], z["occ"], z["S"], z["ao_atom"], symbols)


def localized_range(out_text: str) -> tuple[int, int]:
    """0-based inclusive orbital range ORCA localized ('Orbital range for localization ... a to b')."""
    for ln in out_text.splitlines():
        if "Orbital range for localization" in ln:
            p = ln.split("...")[1].split()
            return int(p[0]), int(p[2])
    raise ValueError("no localization range in output")


def overlap(a: FrameOrbitals, b: FrameOrbitals, ia: list[int], ib: list[int]) -> np.ndarray:
    S = 0.5 * (a.S + b.S)
    return a.C[:, ia].T @ S @ b.C[:, ib]


def track(frames: list[FrameOrbitals], idx: list[list[int]], ref: int = 0):
    """Follow orbitals from frame `ref` outward in both directions.

    idx[k] lists the orbital indices (in frame k) that take part in the matching.
    Returns perm[k][j] = orbital index in frame k for tracked orbital j,
            sign[k][j] = +1/-1 to apply, best[k][j] = |overlap| with the previous frame.
    """
    n = len(frames)
    m = len(idx[ref])
    perm = [None] * n
    sign = [None] * n
    best = [None] * n
    perm[ref] = list(idx[ref])
    sign[ref] = [1.0] * m
    best[ref] = [1.0] * m
    for order in (range(ref + 1, n), range(ref - 1, -1, -1)):
        for k in order:
            p = k - 1 if k > ref else k + 1
            O = overlap(frames[p], frames[k], perm[p], idx[k])
            O = O * np.array(sign[p])[:, None]
            r, c = linear_sum_assignment(-np.abs(O))
            perm[k] = [idx[k][c[j]] for j in np.argsort(r)]
            sel = O[r, c][np.argsort(r)]
            sign[k] = [1.0 if v >= 0 else -1.0 for v in sel]
            best[k] = [float(abs(v)) for v in sel]
    return perm, sign, best


def plot_cubes(wavefile: Path, orbitals: list[int], grid: int = 70) -> list[Path]:
    """Cube files for the listed orbitals (0-based) via orca_plot's interactive menu."""
    wd = wavefile.parent
    before = set(wd.glob("*.cube"))
    cmds = ["1", "1", "4", str(grid), "5", "7"]
    for o in orbitals:
        cmds += ["2", str(o), "11"]
    cmds += ["12"]
    # bytes, not text: orca_plot rejects some console encodings of the menu answers
    subprocess.run([str(orca.ORCA_DIR / "orca_plot.exe"), wavefile.name, "-i"], cwd=wd,
                   input=("\n".join(cmds) + "\n").encode("ascii"), env=_env(),
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, timeout=600)
    return sorted(set(wd.glob("*.cube")) - before)


def read_cube(path: Path):
    """Returns (values (nx,ny,nz) in e^1/2/bohr^3/2, origin (bohr), axes (3x3 bohr), symbols, coords bohr)."""
    lines = Path(path).read_text().splitlines()
    nat, *org = lines[2].split()[:4]
    nat = int(nat)
    origin = np.array([float(v) for v in org])
    dims, axes = [], []
    for i in range(3):
        p = lines[3 + i].split()
        dims.append(int(p[0]))
        axes.append([float(v) for v in p[1:4]])
    natoms = abs(nat)
    start = 6 + natoms + (1 if nat < 0 else 0)
    vals = np.array(" ".join(lines[start:]).split(), dtype=np.float32)
    return vals.reshape(dims), origin, np.array(axes)
