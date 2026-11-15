# End check of the open side (backward) against open_cation, with permutation-aware heavy-atom RMSD.
import itertools, json, numpy as np
from pathlib import Path
from pipeline import orca
from pipeline.steps import Reaction
from pipeline.geometry import Geometry as PGeom, read_xyz
from pipeline.path_angles import _minima, _xyz_frames, rmsd, HARTREE_KJ
from pipeline.end_check import _ref_opt_keyword
R = Path('reactions/hshift_2_3_dimethyl_2_butyl_m062x').resolve()
rx = Reaction.load(R)
mins = {m['name']: m for m in _minima(R)}
oc = mins['open_cation']
sym, x = _xyz_frames(R / '02_irc/true_irc/backward/path.xyz')[-1]
pj = json.loads((R / '02_irc/true_irc/backward/path.json').read_text(encoding='utf-8'))
Eend = pj['energies'][-1]
C = [i for i, s in enumerate(sym) if s == 'C']
def best_heavy(a, b):
    # carbons bonded to exactly one other carbon are methyls; try every relabelling of the carbons that keeps the
    # C–C bond graph (covers the C2/C3 swap and methyl exchanges), heavy atoms only
    def cc(xx): return {(i, j) for i in C for j in C if i < j and np.linalg.norm(xx[i] - xx[j]) < 1.75}
    ga = cc(a); best = 9e9
    for p in itertools.permutations(C):
        m = dict(zip(C, p))
        if {tuple(sorted((m[i], m[j]))) for i, j in ga} != cc(b): continue
        best = min(best, rmsd(a[list(C)], b[[m[i] for i in C]]))
    return best
print('end - open_cation: %.3f kJ/mol; heavy RMSD (perm-aware) %.3f A' % ((Eend - oc['E_eh']) * HARTREE_KJ, best_heavy(x, oc['xyz'])))
kw = _ref_opt_keyword(R, 'open_cation')
jd = R / '02_irc/true_irc/end_check/backward'; jd.mkdir(parents=True, exist_ok=False)
r = orca.run(jd, 'opt', rx.inp(f'{kw} Freq', PGeom(sym, x, 'last frame of the backward (open) side'), blocks='%geom MaxIter 200 end'),
             rx.log, note=f'end check, open side: {kw} Freq from the last frame of the backward side (reference open_cation; permutation-aware comparison)')
g = read_xyz(jd / 'opt.xyz'); xo = np.asarray(g.coords)
out = {'keyword': kw, 'ok': r.ok, 'imaginary_cm': r.imaginary, 'end_minus_open_kJ': (Eend - oc['E_eh']) * HARTREE_KJ,
       'end_heavy_rmsd_perm_A': best_heavy(x, oc['xyz']), 'opt_minus_open_kJ': (r.energy_eh - oc['E_eh']) * HARTREE_KJ,
       'opt_heavy_rmsd_perm_A': best_heavy(xo, oc['xyz'])}
print(json.dumps(out, indent=1))
(R / '02_irc/true_irc/end_check_open_side.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
