"""Build the E2 transition-state guess: EtO⁻ + 2-bromopropane, anti-periplanar.

Run from the repository folder:  .venv\\Scripts\\python.exe reactions\\e2_ethoxide_2_bromopropane\\00_inputs\\make_ts_guess.py
Writes ts_guess.xyz next to this file and prints the atom indices used in reaction.toml.

Recipe (J. Picione, 2026-11-09): MMFF geometry of 2-bromopropane; methyl H anti to Br (H–Cβ–Cα–Br near 180°);
Cα–Br 2.30 Å, Cβ–H 1.30 Å, Cα–Cβ 1.45 Å; ethoxide O on the Cβ–H line, O···H 1.30 Å,
ethyl group pointing away from the substrate. The pipeline then relaxes
everything else with these four bonds frozen (ts_preopt) before the TS search.
"""
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def embed(smiles, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    return m, m.GetConformer().GetPositions()


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def unit(v):
    return v / np.linalg.norm(v)


sub, X = embed("CC(C)Br")
sym = [a.GetSymbol() for a in sub.GetAtoms()]
br = sym.index("Br")
ca = sub.GetAtomWithIdx(br).GetNeighbors()[0].GetIdx()
cb = next(n.GetIdx() for n in sub.GetAtomWithIdx(ca).GetNeighbors() if n.GetSymbol() == "C")
hs = [n.GetIdx() for n in sub.GetAtomWithIdx(cb).GetNeighbors() if n.GetSymbol() == "H"]
h = max(hs, key=lambda i: abs(dihedral(X[i], X[cb], X[ca], X[br])))
print(f"anti H: dihedral H–Cβ–Cα–Br {dihedral(X[h], X[cb], X[ca], X[br]):.1f}°")

X = X.copy()
X[br] = X[ca] + unit(X[br] - X[ca]) * 2.30          # Cα–Br stretched
dHB = unit(X[h] - X[cb])
# Cα–Cβ shortened to 1.45 Å: move Cβ and its other two H together toward Cα
cb_group = [cb] + [i for i in hs if i != h]
shift = unit(X[ca] - X[cb]) * (np.linalg.norm(X[ca] - X[cb]) - 1.45)
X[cb_group] += shift
X[h] = X[cb] + dHB * 1.30                           # Cβ–H stretched
o_pos = X[h] + dHB * 1.30                           # O on the Cβ–H line, O···H 1.30 Å

eto, Y = embed("CC[O-]")
esym = [a.GetSymbol() for a in eto.GetAtoms()]
o = esym.index("O")
c1 = eto.GetAtomWithIdx(o).GetNeighbors()[0].GetIdx()
# orient ethoxide: O at o_pos, O→C1 at ~110° to O→H, ethyl pointing away from Br and the substrate
Y = Y - Y[o]
toward_h = -dHB
away = unit(X[ca] - X[br] + (o_pos - X[[ca, cb]].mean(0)))
perp = unit(away - np.dot(away, toward_h) * toward_h)
target = unit(np.cos(np.radians(110)) * toward_h + np.sin(np.radians(110)) * perp)


def rot_to(a, b):
    a, b = unit(a), unit(b)
    v, c = np.cross(a, b), np.dot(a, b)
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


Y = Y @ rot_to(Y[c1], target).T
Y = Y + o_pos

symbols = sym + esym
coords = np.vstack([X, Y])
n = len(sym)
dmin = min(np.linalg.norm(coords[i] - coords[j]) for i in range(n) for j in range(n, len(symbols))
           if not (i == h and j == n + o))
print(f"closest contact between the fragments (excluding H···O): {dmin:.2f} Å")
out = Path(__file__).with_name("ts_guess.xyz")
lines = [str(len(symbols)), "E2 TS guess: EtO- + 2-bromopropane, anti-periplanar (make_ts_guess.py)"]
lines += [f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(symbols, coords)]
out.write_text("\n".join(lines) + "\n")
print(f"wrote {out}")
print(f"indices: base O = {n + o}, H = {h}, Cβ = {cb}, Cα = {ca}, Br = {br}")
