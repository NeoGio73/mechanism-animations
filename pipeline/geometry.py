"""XYZ input/output and 3D structures from SMILES."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .units import ANGSTROM_TO_BOHR, MASSES


@dataclass
class Geometry:
    symbols: list[str]
    coords: np.ndarray          # (N, 3) in Å
    comment: str = ""

    def distance(self, i: int, j: int) -> float:
        return float(np.linalg.norm(self.coords[i] - self.coords[j]))

    def angle(self, i: int, j: int, k: int) -> float:
        a = self.coords[i] - self.coords[j]
        b = self.coords[k] - self.coords[j]
        c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
        return float(np.degrees(np.arccos(np.clip(c, -1, 1))))

    def masses(self) -> np.ndarray:
        return np.array([MASSES[s] for s in self.symbols])

    def xyz_block(self) -> str:
        return "\n".join(f"{s:<3}{x:14.8f}{y:14.8f}{z:14.8f}"
                         for s, (x, y, z) in zip(self.symbols, self.coords))

    def to_xyz(self, comment: str | None = None) -> str:
        c = self.comment if comment is None else comment
        return f"{len(self.symbols)}\n{c}\n{self.xyz_block()}\n"

    def write(self, path: Path, comment: str | None = None) -> None:
        Path(path).write_text(self.to_xyz(comment))


def read_xyz_frames(path: Path) -> list[Geometry]:
    lines = Path(path).read_text().splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1]
        syms, xyz = [], []
        for ln in lines[i + 2:i + 2 + n]:
            p = ln.split()
            syms.append(p[0])
            xyz.append([float(v) for v in p[1:4]])
        frames.append(Geometry(syms, np.array(xyz), comment))
        i += 2 + n
    return frames


def read_xyz(path: Path) -> Geometry:
    return read_xyz_frames(path)[-1]


def mass_weighted_step(a: Geometry, b: Geometry) -> float:
    """Mass-weighted distance between two frames, in bohr·amu^1/2."""
    d = (b.coords - a.coords) * ANGSTROM_TO_BOHR
    return float(np.sqrt((a.masses()[:, None] * d ** 2).sum()))


def from_smiles(smiles: str, seed: int = 7) -> Geometry:
    """3D structure from SMILES (RDKit ETKDG + MMFF where possible)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol.GetNumAtoms() == 1:
        return Geometry([mol.GetAtomWithIdx(0).GetSymbol()], np.zeros((1, 3)), smiles)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass  # MMFF lacks some ions; DFT optimizes anyway
    conf = mol.GetConformer()
    xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return Geometry([a.GetSymbol() for a in mol.GetAtoms()], xyz, smiles)


def formal_charge(smiles: str) -> int:
    from rdkit import Chem
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
