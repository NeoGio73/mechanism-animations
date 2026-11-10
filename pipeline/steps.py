"""Pipeline steps for one reaction. Each step reads earlier results from disk,
so steps can be rerun one at a time; finished ORCA jobs with unchanged input are reused."""
from __future__ import annotations

import json
import shutil
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import orca
from .geometry import Geometry, formal_charge, from_smiles, read_xyz
from .units import fmt_kj, kj


@dataclass
class Reaction:
    root: Path
    cfg: dict

    @classmethod
    def load(cls, root: Path) -> "Reaction":
        root = Path(root).resolve()
        cfg = tomllib.loads((root / "00_inputs" / "reaction.toml").read_text(encoding="utf-8"))
        return cls(root, cfg)

    @property
    def log(self) -> orca.RunLog:
        return orca.RunLog(self.root / "run_log.md")

    def d(self, sub: str) -> Path:
        p = self.root / sub
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def solvent_block(self) -> str:
        return f'%cpcm smd true SMDsolvent "{self.cfg["solvent"].upper()}" end'

    def keywords(self, job: str, sp: bool = False) -> str:
        m = self.cfg["sp_method"] if sp else self.cfg["method"]
        return f'{m} CPCM({self.cfg["solvent"]}) {job}'.strip()

    def inp(self, job: str, geom: Geometry, charge=None, mult=None, blocks="", sp=False,
            nprocs=None) -> str:
        return orca.make_input(
            self.keywords(job, sp), geom,
            self.cfg["charge"] if charge is None else charge,
            self.cfg["mult"] if mult is None else mult,
            blocks=self.solvent_block + "\n" + blocks,
            nprocs=self.cfg["nprocs"] if nprocs is None else nprocs,
            maxcore_mb=self.cfg["maxcore_mb"])

    def level(self, sp: bool = False) -> str:
        """'B3LYP-D3(BJ)/ma-def2-SVP' for the working level (sp=True: the single-point level)."""
        return level_name(self.cfg["sp_method"] if sp else self.cfg["method"])

    def basis(self, sp: bool = False) -> str:
        return self.level(sp).split("/")[1]

    def species(self):
        """Reactants, products and optional `[[references]]` (optimized and reported, never part
        of the barrier, e.g. the s-trans ground-state conformer of a diene that reacts s-cis)."""
        for role in ("reactants", "products", "references"):
            for s in self.cfg.get(role, []):
                yield role, s

    def key(self, name: str) -> int:
        return self.cfg["key_atoms"][name]

    def bonds(self, g: Geometry) -> tuple[float, float]:
        """(nucleophile–center, center–leaving group) distances in Å."""
        return (g.distance(self.key("nucleophile"), self.key("center")),
                g.distance(self.key("center"), self.key("leaving_group")))

    def key_bonds(self, symbols: list[str] | None = None) -> list[dict]:
        """Bonds that change along the reaction: [{atoms: [a, b], role, label}].

        role: "forms" (long → bonded), "breaks" (bonded → long), "shortens" (stays bonded,
        bond order rises, e.g. C–C → C=C in an elimination; never drawn dashed) or "lengthens"
        (stays bonded, bond order falls, e.g. C=C → C–C in a cycloaddition; never drawn dashed).
        From `[[key_bonds]]` in reaction.toml; without it, nucleophile–center forms and
        center–leaving_group breaks (SN2)."""
        if "key_bonds" in self.cfg:
            kb = [{"atoms": list(b["atoms"]), "role": b["role"], "label": b.get("label")}
                  for b in self.cfg["key_bonds"]]
        else:
            kb = [{"atoms": [self.key("center"), self.key("nucleophile")], "role": "forms", "label": None},
                  {"atoms": [self.key("center"), self.key("leaving_group")], "role": "breaks", "label": None}]
        if symbols is None:
            symbols = read_xyz(self.root / "00_inputs" / self.cfg["ts_guess"]).symbols
        for b in kb:
            if not b["label"]:
                b["label"] = bond_label(symbols, *b["atoms"])
        return kb

    def bond_lengths(self, g: Geometry) -> list[float]:
        return [g.distance(*b["atoms"]) for b in self.key_bonds(g.symbols)]

    def progress(self, g: Geometry) -> float:
        """Larger = further toward products: breaking bonds long, forming bonds short (Å)."""
        return sum((1 if b["role"] == "breaks" else -1) * g.distance(*b["atoms"])
                   for b in self.key_bonds(g.symbols) if b["role"] in ("forms", "breaks"))


def level_name(method: str) -> str:
    """'B3LYP D3BJ ma-def2-SVP def2/J RIJCOSX TightSCF' → 'B3LYP-D3(BJ)/ma-def2-SVP'"""
    tok = method.split()
    disp = {"D3BJ": "-D3(BJ)", "D3ZERO": "-D3(0)", "D4": "-D4"}
    func = tok[0] + "".join(disp[t] for t in tok if t in disp)
    basis = next(t for t in tok if "def2-" in t and "/" not in t)
    return f"{func}/{basis}"


def bond_label(symbols, a: int, b: int) -> str:
    """'C–O', 'C–Br', 'O–H': carbon first, hydrogen last, otherwise as given."""
    pa, pb = symbols[a], symbols[b]
    rank = lambda s: 0 if s == "C" else (2 if s == "H" else 1)
    if rank(pb) < rank(pa):
        pa, pb = pb, pa
    return f"{pa}–{pb}"


def _write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------- separated species
def fragments(rx: Reaction) -> None:
    """Optimize each separated reactant and product (Opt Freq), working level."""
    res = {}
    for role, s in rx.species():
        wd = rx.d(f"01_ts/fragments/{s['name']}")
        # optional fixed starting geometry (e.g. a chosen conformer), else built from SMILES
        g = read_xyz(rx.root / "00_inputs" / s["start_xyz"]) if "start_xyz" in s else from_smiles(s["smiles"])
        g.write(wd / "start.xyz")
        q, m = formal_charge(s["smiles"]), s.get("mult", 1)
        job = "Opt Freq" if len(g.symbols) > 1 else ""
        r = orca.run(wd, s["name"], rx.inp(job, g, q, m, nprocs=1), rx.log,
                     note=f"{role[:-1]} {s['smiles']}")
        if not r.ok or r.opt_converged is False:
            raise orca.OrcaError(f"{s['name']} failed; see {r.out}")
        res[s["name"]] = {"role": role, "smiles": s["smiles"], "charge": q,
                          "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
                          "imaginary_cm": r.imaginary}
        if r.imaginary:
            print(f"WARNING {s['name']} has imaginary frequencies {r.imaginary}")
    _write_json(rx.d("01_ts") / "fragments.json", res)


# ---------------------------------------------------------------- transition state
def check_ts_mode(rx: Reaction, hess_path: Path, geom: Geometry) -> dict:
    """Is the imaginary mode the Nu–C–LG motion? Measures how much of the mode is the
    central atom moving along the nucleophile→leaving-group axis."""
    h = orca.read_hess(hess_path)
    k = int(np.argmin(h["freqs"]))
    v = h["modes"][:, k].reshape(-1, 3)
    nu, c, lg = rx.key("nucleophile"), rx.key("center"), rx.key("leaving_group")
    axis = geom.coords[lg] - geom.coords[nu]
    axis /= np.linalg.norm(axis)
    share = np.linalg.norm(v, axis=1) ** 2 / (v ** 2).sum()
    # relative motion along the axis: C moves toward one end while the ends move apart less
    along_c = abs(np.dot(v[c], axis)) / (np.linalg.norm(v[c]) + 1e-12)
    # General test, any number of changing bonds: how fast each key bond stretches along the
    # mode. A correct mode shortens every forming bond while it lengthens every breaking bond.
    rates = []
    for b in rx.key_bonds(geom.symbols):
        a1, a2 = b["atoms"]
        u = geom.coords[a1] - geom.coords[a2]
        u /= np.linalg.norm(u)
        rates.append({"bond": b["label"], "role": b["role"], "atoms": b["atoms"],
                      "stretch_rate": float(np.dot(v[a1] - v[a2], u))})
    big = max(abs(x["stretch_rate"]) for x in rates)
    signs = {x["role"]: set() for x in rates}
    for x in rates:
        if abs(x["stretch_rate"]) > 0.15 * big:
            signs[x["role"]].add(np.sign(x["stretch_rate"]))
    f, br = signs.get("forms", set()), signs.get("breaks", set())
    consistent = len(f) <= 1 and len(br) <= 1 and (not f or not br or f != br)
    return {"freq_cm": float(h["freqs"][k]),
            "central_atom_share": float(share[c]),
            "central_atom_alignment_with_axis": float(along_c),
            "largest_atom": geom.symbols[int(np.argmax(share))] + str(int(np.argmax(share))),
            "key_bond_stretch_rates": rates,
            "forming_and_breaking_bonds_move_together": bool(consistent)}


def dihedral(g: Geometry, atoms) -> float:
    """Dihedral angle a–b–c–d in degrees (−180 … 180)."""
    p = g.coords[list(atoms)]
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _accept_ts(rx: Reaction, g: Geometry, out: dict) -> dict:
    """Accept or reject the TS from its imaginary mode and geometry.

    Always: exactly one imaginary frequency; forming bonds shorten while breaking bonds
    lengthen. Optional `[ts_mode]` in reaction.toml:
      dominant_atom = i        the atom that must move most (e.g. the transferring H in E2)
      reject_axis = [a, b]     with reject_atom = c: reject when c moves mainly along a···b
                               (the SN2 motion of the central carbon), classified as SN2
      synchronous = [labels]   key bonds that must form together: lengths within
                               synchronous_max_spread_A (0.05 Å) and stretch rates in the mode
                               within synchronous_min_rate_ratio (0.8); else "asynchronous"
    Optional `[ts_expected]`: label = [low, high] for bonds (Å) and dihedrals (°, compared
    as |angle|); outside the range gives a warning, not a rejection. `[ts_required]`: same
    form, but outside the range rejects the TS."""
    m = out["mode_check"]
    cfg = rx.cfg.get("ts_mode", {})
    reasons, cls = [], "expected TS"
    if out["n_imaginary"] != 1:
        reasons.append(f"{out['n_imaginary']} imaginary frequencies, need exactly 1")
    rates = {x["bond"]: x for x in m["key_bond_stretch_rates"]}
    forms = [x["stretch_rate"] for x in rates.values() if x["role"] == "forms"]
    ref = np.sign(forms[0]) if forms else 0
    for x in rates.values():
        want = {"forms": ref, "shortens": ref, "breaks": -ref, "lengthens": -ref}[x["role"]]
        if ref and np.sign(x["stretch_rate"]) != want:
            reasons.append(f"{x['bond']} moves the wrong way in the mode ({x['role']})")
    if "dominant_atom" in cfg:
        h = orca.read_hess(rx.root / "01_ts/ts/ts.hess")
        v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
        share = np.linalg.norm(v, axis=1) ** 2 / (v ** 2).sum()
        m["atom_shares"] = {f"{g.symbols[i]}{i}": round(float(share[i]), 3)
                            for i in np.argsort(-share)[:5]}
        if int(np.argmax(share)) != cfg["dominant_atom"]:
            reasons.append(f"largest motion is {m['largest_atom']}, not "
                           f"{g.symbols[cfg['dominant_atom']]}{cfg['dominant_atom']}")
        if "reject_axis" in cfg:
            a, b = cfg["reject_axis"]
            c = cfg["reject_atom"]
            ax = g.coords[b] - g.coords[a]
            ax /= np.linalg.norm(ax)
            along = abs(np.dot(v[c], ax)) / (np.linalg.norm(v[c]) + 1e-12)
            m["reject_atom_alignment"] = float(along)
            if along > 0.7 and share[c] >= share[cfg["dominant_atom"]]:
                cls = (f"SN2-like: {g.symbols[c]}{c} moves along the "
                       f"{g.symbols[a]}{a}···{g.symbols[b]}{b} axis")
    if "synchronous" in cfg:
        # bonds that must form together: equal lengths and equal stretch rates in the mode
        sb = cfg["synchronous"]
        d = [out["bonds_A"][lab] for lab in sb]
        r = [abs(rates[lab]["stretch_rate"]) for lab in sb]
        m["synchronous_length_spread_A"] = float(max(d) - min(d))
        m["synchronous_rate_ratio"] = float(min(r) / max(r))
        if max(d) - min(d) > cfg.get("synchronous_max_spread_A", 0.05):
            cls = f"asynchronous: {', '.join(f'{lab} {x:.2f} Å' for lab, x in zip(sb, d))}"
        elif min(r) / max(r) < cfg.get("synchronous_min_rate_ratio", 0.8):
            cls = f"asynchronous mode: stretch-rate ratio {min(r) / max(r):.2f}"
    for label, (lo, hi) in rx.cfg.get("ts_required", {}).items():
        val = out["bonds_A"].get(label, out["dihedrals_deg"].get(label, out.get("angles_deg", {}).get(label)))
        x = abs(val) if label in out["dihedrals_deg"] else val
        if not lo <= x <= hi:
            reasons.append(f"{label} {x:.2f} outside the required {lo}–{hi}")
    if reasons and cls == "expected TS":
        cls = "unexpected mode"
    warnings = []
    for label, (lo, hi) in rx.cfg.get("ts_expected", {}).items():
        val = out["bonds_A"].get(label, out["dihedrals_deg"].get(label, out.get("angles_deg", {}).get(label)))
        if val is None:
            continue
        x = abs(val) if label in out["dihedrals_deg"] else val
        if not lo <= x <= hi:
            warnings.append(f"{label} {x:.2f} outside the expected {lo}–{hi}")
    for w in warnings:
        print("WARNING", w)
    return {"accepted": not reasons and cls == "expected TS", "classification": cls,
            "reasons": reasons, "range_warnings": warnings}


def _constraint_block(rx: Reaction) -> str:
    """Freeze every changing bond at its guess value (pre-optimization)."""
    lines = [f"    {{B {b['atoms'][0]} {b['atoms'][1]} C}}" for b in rx.key_bonds()]
    return "%geom\n  Constraints\n" + "\n".join(lines) + "\n  end\nend"


def ts(rx: Reaction) -> None:
    """OptTS + Freq from the TS guess, exact Hessian at start and every 5 steps.

    With `ts_preopt = true` in reaction.toml, the guess is first optimized with the forming
    and breaking bonds frozen, so everything else relaxes before the TS search."""
    wd = rx.d("01_ts/ts")
    g = read_xyz(rx.root / "00_inputs" / rx.cfg["ts_guess"])
    if rx.cfg.get("ts_preopt"):
        pw = rx.d("01_ts/preopt")
        r = orca.run(pw, "preopt", rx.inp("Opt", g, blocks=_constraint_block(rx)), rx.log,
                     note="TS guess relaxed with forming/breaking bonds frozen")
        if not r.ok:
            raise orca.OrcaError(f"constrained pre-optimization failed; see {r.out}")
        if not r.opt_converged:
            # Only a starting guess: accept it when the energy has settled (flat modes such as a
            # rotating ethyl group can keep the step criteria from being met).
            E = [float(ln.split()[-1]) for ln in r.out.read_text(errors="replace").splitlines()
                 if ln.startswith("FINAL SINGLE POINT ENERGY")][-10:]
            spread = kj(max(E) - min(E))
            if len(E) < 10 or spread > 0.5:
                raise orca.OrcaError(f"constrained pre-optimization did not settle "
                                     f"({spread:.2f} kJ/mol over the last cycles); see {r.out}")
            note = (f"Constrained pre-optimization reached its cycle limit without meeting the step "
                    f"criteria, but its energy varied by only {spread:.2f} kJ/mol over the last 10 "
                    f"cycles (flat surface); accepted as the TS starting guess.")
            if note not in rx.log.path.read_text(encoding="utf-8"):
                rx.log.note(note)
            print("NOTE", note)
        g = read_xyz(pw / "preopt.xyz")
    blocks = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end"
    r = orca.run(wd, "ts", rx.inp("OptTS Freq", g, blocks=blocks), rx.log, note="TS search from guess")
    if not r.ok or not r.opt_converged:
        raise orca.OrcaError(f"TS optimization failed; see {r.out}")
    gts = read_xyz(wd / "ts.xyz")
    mode = check_ts_mode(rx, wd / "ts.hess", gts)
    d_nu, d_lg = rx.bonds(gts)
    bonds = {b["label"]: gts.distance(*b["atoms"]) for b in rx.key_bonds(gts.symbols)}
    dihedrals = {name: dihedral(gts, atoms) for name, atoms in rx.cfg.get("ts_dihedrals", {}).items()}
    # optional `[ts_angles]`: label = [a, b, c], angle a–b–c in degrees (e.g. Bürgi–Dunitz Nu···C=O)
    angles = {name: gts.angle(*atoms) for name, atoms in rx.cfg.get("ts_angles", {}).items()}
    out = {"E_eh": r.energy_eh, "G_eh": r.gibbs_eh, "imaginary_cm": r.imaginary,
           "n_imaginary": len(r.imaginary), "mode_check": mode,
           "bonds_A": bonds, "dihedrals_deg": dihedrals, "angles_deg": angles,
           "d_nu_center_A": d_nu, "d_center_lg_A": d_lg,
           "angle_nu_center_lg_deg": gts.angle(rx.key("nucleophile"), rx.key("center"),
                                               rx.key("leaving_group"))}
    verdict = _accept_ts(rx, gts, out)
    out.update(verdict)
    _write_json(wd.parent / "ts.json", out)
    print(json.dumps(out, indent=2, ensure_ascii=False))
    if not verdict["accepted"]:
        # keep the rejected TS for the record and stop; the guess must be corrected
        k = len(list(rx.d("01_ts/attempts").glob("*"))) + 1
        dest = rx.root / f"01_ts/attempts/{k}_{verdict['classification'].split(':')[0].replace(' ', '_')}"
        shutil.copytree(wd, dest)
        shutil.copy(wd.parent / "ts.json", dest / "ts.json")
        rx.log.note(f"TS rejected ({verdict['classification']}); kept in `{dest.relative_to(rx.root).as_posix()}`. "
                    f"Reasons: {'; '.join(verdict['reasons'])}")
        raise orca.OrcaError(f"TS rejected: {verdict['classification']}: {verdict['reasons']}")


# ---------------------------------------------------------------- IRC
COV_RADII_A = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
               "Cl": 1.02, "Br": 1.20, "I": 1.39}


def _bond_states(rx: Reaction, g: Geometry) -> tuple:
    """For each forming/breaking key bond: 1 = bonded (≤ 1.25 × sum of covalent radii),
    0 = not bonded (≥ 1.5 ×), None = in between."""
    out = []
    for b in rx.key_bonds(g.symbols):
        if b["role"] not in ("forms", "breaks"):
            continue
        a, c = b["atoms"]
        r = COV_RADII_A[g.symbols[a]] + COV_RADII_A[g.symbols[c]]
        d = g.distance(a, c)
        out.append(1 if d <= 1.25 * r else 0 if d >= 1.5 * r else None)
    return tuple(out)


def _trj_energies_eh(frames) -> list[float]:
    """Energies from ORCA IRC trajectory comment lines ('Coordinates from ORCA-job ... E <Eh>')."""
    return [float(f.comment.split(" E ")[-1].split()[0]) for f in frames]


def irc_stop_index(rx: Reaction, frames) -> tuple[int, int] | None:
    """IRC stopping rule (docs/protocol.md): stop once the energy has changed by less than
    irc_stop_dE_kJ (default 0.05 kJ/mol) per step for irc_stop_steps (default 20) consecutive
    steps after the last bond-forming or bond-breaking event. The last event is the first frame
    from which every forming/breaking key bond keeps its final state (bonded / not bonded).
    Returns (last_event_frame, stop_frame) or None."""
    dE_max = float(rx.cfg.get("irc_stop_dE_kJ", 0.05))
    n_need = int(rx.cfg.get("irc_stop_steps", 20))
    if len(frames) < n_need + 2:
        return None
    states = [_bond_states(rx, g) for g in frames]
    end = states[-1]
    if None in end:
        return None
    k0 = max((k for k, s in enumerate(states) if s != end), default=-1) + 1
    E = [kj(e) for e in _trj_energies_eh(frames)]
    run = 0
    for k in range(k0 + 1, len(E)):
        run = run + 1 if abs(E[k] - E[k - 1]) < dE_max else 0
        if run >= n_need:
            return k0, k
    return None


def _irc_blocks(rx: Reaction, direction: str) -> str:
    step = rx.cfg.get("irc_scale_displ_sd")       # smaller → more frames along the path
    extra = f"  Scale_Displ_SD {step}\n" if step else ""
    return (f'%irc\n  MaxIter 300\n  Direction {direction}\n  InitHess read\n'
            f'  Hess_Filename "ts.hess"\n{extra}  PrintLevel 1\nend')


def irc(rx: Reaction) -> None:
    """IRC from the optimized TS, reading the TS Hessian.

    New runs: one job per direction (irc_f, irc_b), each ended by the IRC stopping rule when it
    applies, then joined into irc_IRC_Full_trj.xyz + irc_path.json. Older runs made as one
    'Direction both' job (irc.out) are reused; if such a run was stopped by decision
    (02_irc/irc_stop.json), its partial trajectories are joined the same way."""
    wd = rx.d("02_irc")
    g = read_xyz(rx.root / "01_ts/ts/ts.xyz")
    step = rx.cfg.get("irc_scale_displ_sd")
    how = f", Scale_Displ_SD {step}" if step else ", default step"
    legacy_out = wd / "irc.out"
    if (wd / "irc_stop.json").exists() and legacy_out.exists():
        _join_irc(rx, "irc", "irc", _read_json(wd / "irc_stop.json"))
        return
    if legacy_out.exists() and "ORCA TERMINATED NORMALLY" in legacy_out.read_text(errors="replace"):
        shutil.copy(rx.root / "01_ts/ts/ts.hess", wd / "ts.hess")
        r = orca.run(wd, "irc", rx.inp("IRC", g, blocks=_irc_blocks(rx, "both")), rx.log,
                     note="IRC both directions" + how)
        if not r.ok:
            raise orca.OrcaError(f"IRC failed; see {r.out}")
        return
    shutil.copy(rx.root / "01_ts/ts/ts.hess", wd / "ts.hess")
    stops = {}
    for direction, tag in (("forward", "F"), ("backward", "B")):
        name = f"irc_{tag.lower()}"
        trj = wd / f"{name}_IRC_{tag}_trj.xyz"

        def watch(trj=trj):
            from .geometry import read_xyz_frames
            if not trj.exists():
                return None
            try:
                frames = read_xyz_frames(trj)
            except (ValueError, IndexError):        # file caught mid-write
                return None
            hit = irc_stop_index(rx, frames)
            if hit:
                return (f"IRC stopping rule met at step {hit[1]} (energy change < "
                        f"{rx.cfg.get('irc_stop_dE_kJ', 0.05)} kJ/mol per step for "
                        f"{rx.cfg.get('irc_stop_steps', 20)} steps after the last bond change "
                        f"at step {hit[0]})")
            return None

        r = orca.run(wd, name, rx.inp("IRC", g, blocks=_irc_blocks(rx, direction)), rx.log,
                     note=f"IRC {direction}{how}", watch=watch)
        if not r.ok and not r.stopped:
            raise orca.OrcaError(f"IRC {direction} failed; see {r.out}")
        if r.stopped:
            stops[tag] = {"by": "rule", "reason": r.stopped}
    _join_irc(rx, "irc_f", "irc_b", {"directions": stops} if stops else {})


def _join_irc(rx: Reaction, name_f: str, name_b: str, stop_info: dict) -> None:
    """Join backward (reversed) + TS + forward trajectories into irc_IRC_Full_trj.xyz (the same
    order ORCA uses), write the two end geometries and irc_path.json (energies, TS index, stops)."""
    from .geometry import read_xyz_frames
    wd = rx.root / "02_irc"
    ts_g = read_xyz(rx.root / "01_ts/ts/ts.xyz")
    out_f = (wd / f"{name_f}.out").read_text(errors="replace")
    e_ts = float(out_f.split("FINAL SINGLE POINT ENERGY", 1)[1].split()[0])   # TS point of the IRC job
    sides = {}
    for name, tag in ((name_f, "F"), (name_b, "B")):
        frames = read_xyz_frames(wd / f"{name}_IRC_{tag}_trj.xyz")
        E = _trj_energies_eh(frames)
        info = stop_info.get("directions", {}).get(tag)
        if info and info.get("by") == "rule":      # cut at the rule's frame (poll timing adds a few)
            hit = irc_stop_index(rx, frames)
            if hit:
                frames, E = frames[:hit[1] + 1], E[:hit[1] + 1]
                info["last_step"] = hit[1]
                info["last_event_step"] = hit[0]
        sides[tag] = (frames, E)
    fb, eb = sides["B"]
    ff, ef = sides["F"]
    frames = fb[::-1] + [ts_g] + ff
    E = eb[::-1] + [e_ts] + ef
    (wd / "irc_IRC_Full_trj.xyz").write_text(
        "".join(g.to_xyz(f"IRC frame {k} E {e:.12f}") for k, (g, e) in enumerate(zip(frames, E))))
    fb[-1].write(wd / "irc_IRC_B.xyz", f"IRC B end E {eb[-1]:.12f}")
    ff[-1].write(wd / "irc_IRC_F.xyz", f"IRC F end E {ef[-1]:.12f}")
    _write_json(wd / "irc_path.json", {
        "energies_eh": E, "ts_index": len(fb), "n_backward": len(fb), "n_forward": len(ff),
        "stops": stop_info.get("directions", {}),
        "note": "frames = backward (reversed) + TS + forward; energies from the IRC trajectories"})


def endpoints(rx: Reaction) -> None:
    """Optimize the two IRC end geometries into the reactant and product complexes."""
    wd = rx.d("02_irc")
    out = {}
    ends = {e: read_xyz(wd / f"irc_IRC_{e}.xyz") for e in ("F", "B")}
    reactant_end = min(ends, key=lambda e: rx.progress(ends[e]))
    for end, g in ends.items():
        label = "reactant_complex" if end == reactant_end else "product_complex"
        ew = rx.d(f"02_irc/endpoints/{label}")
        g.write(ew / "start.xyz", f"IRC {end} end")
        r = orca.run(ew, label, rx.inp("Opt Freq", g), rx.log, note=f"from IRC {end} end")
        name = label
        if r.ok and r.opt_converged is False:
            # Loose complexes (fragments drifting on a flat surface) can use up ORCA's default
            # 3N cycles: continue once from the last geometry with a larger cycle limit.
            prev = name
            name = f"{label}_cont"
            r = orca.run(ew, name, rx.inp("Opt Freq", read_xyz(ew / f"{prev}.xyz"),
                                          blocks="%geom\n  MaxIter 200\nend"), rx.log,
                         note=f"rerun: {prev} reached ORCA's default cycle limit without converging "
                              "(flat surface); continued from its last geometry with MaxIter 200")
        if not r.ok or not r.opt_converged:
            raise orca.OrcaError(f"{name} optimization failed; see {r.out}")
        if r.imaginary:
            # A minimum with a small imaginary mode (flat surface): push the structure along that
            # mode (largest atom moves 0.1 Å) and re-optimize with tighter criteria, once.
            h = orca.read_hess(ew / f"{name}.hess")
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            gp = read_xyz(ew / f"{name}.xyz")
            gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
            prev, name = name, f"{label}_tight"
            r = orca.run(ew, name, rx.inp("TightOpt Freq", gp), rx.log,
                         note=f"rerun: {prev} had imaginary mode(s) {[round(f, 1) for f in r.imaginary]} "
                              "cm-1; displaced 0.1 Å along it, TightOpt")
            if not r.ok or not r.opt_converged:
                raise orca.OrcaError(f"{name} optimization failed; see {r.out}")
        go = read_xyz(ew / f"{name}.xyz")
        d_nu2, d_lg2 = rx.bonds(go)
        out[label] = {"irc_end": end, "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
                      "imaginary_cm": r.imaginary, "xyz": f"02_irc/endpoints/{label}/{name}.xyz",
                      "bonds_A": {b["label"]: go.distance(*b["atoms"]) for b in rx.key_bonds(go.symbols)},
                      "d_nu_center_A": d_nu2, "d_center_lg_A": d_lg2}
        if r.imaginary:
            print(f"WARNING {label} still has imaginary frequencies {r.imaginary}")
    if set(out) != {"reactant_complex", "product_complex"}:
        raise orca.OrcaError(f"IRC ends do not give one reactant and one product side: {out}")
    _write_json(wd / "endpoints.json", out)
    print(json.dumps(out, indent=2))


# ---------------------------------------------------------------- big-basis single points
def singlepoints(rx: Reaction) -> None:
    """Single points at `sp_method` (larger basis, same solvent) on all stationary points."""
    out = {}
    for name, (xyz, q, m) in _stationary_points(rx).items():
        r = orca.run(rx.d(f"01_ts/sp_tzvp/{name}"), name,
                     rx.inp("", read_xyz(xyz), q, m, sp=True), rx.log, note=f"{rx.basis(sp=True)} single point")
        if not r.ok:
            raise orca.OrcaError(f"single point {name} failed; see {r.out}")
        out[name] = r.energy_eh
    _write_json(rx.d("01_ts") / "sp_tzvp.json", out)


# ---------------------------------------------------------------- summary of energies
def barriers(rx: Reaction) -> dict:
    frag = _read_json(rx.root / "01_ts/fragments.json")
    tsj = _read_json(rx.root / "01_ts/ts.json")
    ep = _read_json(rx.root / "02_irc/endpoints.json")
    sp = _read_json(rx.root / "01_ts/sp_tzvp.json")
    reac = [n for n, v in frag.items() if v["role"] == "reactants"]
    prod = [n for n, v in frag.items() if v["role"] == "products"]

    def table(E: dict) -> dict:
        r0 = sum(E[n] for n in reac)
        p0 = sum(E[n] for n in prod)
        rel = lambda e: kj(e - r0)
        return {"separated_reactants": 0.0,
                "reactant_complex": rel(E["reactant_complex"]),
                "ts": rel(E["ts"]),
                "product_complex": rel(E["product_complex"]),
                "separated_products": rel(p0),
                "barrier_from_reactant_complex": kj(E["ts"] - E["reactant_complex"]),
                "barrier_from_separated_reactants": rel(E["ts"]),
                "reaction_energy": rel(p0)}

    E_svp = {n: v["E_eh"] for n, v in frag.items()}
    E_svp.update(ts=tsj["E_eh"], reactant_complex=ep["reactant_complex"]["E_eh"],
                 product_complex=ep["product_complex"]["E_eh"])
    out = {"units": "kJ/mol, electronic energies relative to separated reactants",
           "labels": {"working": rx.level(), "sp": rx.level(sp=True)},
           "working": table(E_svp), "sp": table(sp)}
    G = {n: v["G_eh"] for n, v in frag.items()}
    G.update(ts=tsj["G_eh"], reactant_complex=ep["reactant_complex"]["G_eh"],
             product_complex=ep["product_complex"]["G_eh"])
    if all(v is not None for v in G.values()):
        out["working_gibbs_298K_1atm"] = table(G)
    _write_json(rx.root / "02_irc" / "barriers.json", out)
    for lvl in ("working", "sp"):
        t = out[lvl]
        print(f"{out['labels'][lvl]}: barrier from complex {fmt_kj(t['barrier_from_reactant_complex'])}, "
              f"from separated reactants {fmt_kj(t['barrier_from_separated_reactants'])}, "
              f"reaction energy {fmt_kj(t['reaction_energy'])}")
    return out


RT_KJ = 8.314462618e-3 * 298.15             # kJ/mol at 298.15 K
STD_STATE_KJ = RT_KJ * np.log(24.4654)      # 1 atm → 1 mol/L, per species (7.92 kJ/mol)


def _stationary_points(rx: Reaction) -> dict:
    """name -> (xyz, charge, mult) for every stationary point used in the energetics."""
    pts = {"ts": (rx.root / "01_ts/ts/ts.xyz", rx.cfg["charge"], rx.cfg["mult"])}
    ep = _read_json(rx.root / "02_irc/endpoints.json")
    for label in ("reactant_complex", "product_complex"):
        pts[label] = (rx.root / ep[label]["xyz"], rx.cfg["charge"], rx.cfg["mult"])
    for role, s in rx.species():
        p = rx.root / f"01_ts/fragments/{s['name']}"
        xyz = p / f"{s['name']}.xyz"
        pts[s["name"]] = (xyz if xyz.exists() else p / "start.xyz",
                          formal_charge(s["smiles"]), s.get("mult", 1))
    return pts


def diagnostics(rx: Reaction) -> None:
    """Single points at the extra levels in reaction.toml; table of barriers per level.

    Free-energy barriers use the working-level thermal corrections (G - E at the working level)
    and the 1 atm → 1 M standard-state correction for the bimolecular step."""
    frag = _read_json(rx.root / "01_ts/fragments.json")
    tsj = _read_json(rx.root / "01_ts/ts.json")
    ep = _read_json(rx.root / "02_irc/endpoints.json")
    reac = [n for n, v in frag.items() if v["role"] == "reactants"]
    prod = [n for n, v in frag.items() if v["role"] == "products"]
    gcorr = {n: v["G_eh"] - v["E_eh"] for n, v in frag.items() if v["G_eh"] is not None}
    gcorr["ts"] = tsj["G_eh"] - tsj["E_eh"]
    for k in ("reactant_complex", "product_complex"):
        gcorr[k] = ep[k]["G_eh"] - ep[k]["E_eh"]

    levels = [{"label": rx.level().replace("/", "_") + " (working)", "E": {
        **{n: v["E_eh"] for n, v in frag.items()}, "ts": tsj["E_eh"],
        "reactant_complex": ep["reactant_complex"]["E_eh"],
        "product_complex": ep["product_complex"]["E_eh"]}}]
    levels.append({"label": rx.level(sp=True).replace("/", "_") + " (protocol)",
                   "E": _read_json(rx.root / "01_ts/sp_tzvp.json")})
    for lv in rx.cfg.get("diagnostic_levels", []):
        E = {}
        for name, (xyz, q, m) in _stationary_points(rx).items():
            wd = rx.d(f"01_ts/diagnostics/{lv['label']}/{name}")
            inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', read_xyz(xyz), q, m,
                                  blocks=rx.solvent_block, nprocs=rx.cfg["nprocs"],
                                  maxcore_mb=rx.cfg["maxcore_mb"])
            r = orca.run(wd, name, inp, rx.log, note=f"diagnostic {lv['label']}")
            if not r.ok:
                print(f"WARNING diagnostic {lv['label']} {name} failed; level skipped")
                E = None
                break
            E[name] = r.energy_eh
        if E:
            levels.append({"label": lv["label"], "E": E})

    rows = []
    for lv in levels:
        E = lv["E"]
        r0 = sum(E[n] for n in reac)
        dE_sep = kj(E["ts"] - r0)
        dE_rc = kj(E["ts"] - E["reactant_complex"])
        dE_rxn = kj(sum(E[n] for n in prod) - r0)
        dG_sep = None
        if all(n in gcorr for n in reac):
            dG_sep = dE_sep + kj(gcorr["ts"] - sum(gcorr[n] for n in reac)) \
                - STD_STATE_KJ * (len(reac) - 1)
        rows.append({"level": lv["label"], "dE_TS_vs_separated_kJ": dE_sep,
                     "dE_TS_vs_complex_kJ": dE_rc, "dE_reaction_kJ": dE_rxn,
                     "dG_TS_vs_separated_1M_kJ": dG_sep})
    _write_json(rx.root / "01_ts" / "diagnostics.json",
                {"units": "kJ/mol", "rows": rows,
                 "energies_eh": {lv["label"]: lv["E"] for lv in levels[2:]}})
    md = [f"| Level (SMD {rx.cfg['solvent']}, {rx.basis()} geometries) | TS − separated reactants, ΔE | "
          "TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |",
          "|---|---|---|---|---|"]
    for r in rows:
        g = fmt_kj(r["dG_TS_vs_separated_1M_kJ"]) if r["dG_TS_vs_separated_1M_kJ"] is not None else "—"
        md.append(f"| {r['level']} | {fmt_kj(r['dE_TS_vs_separated_kJ'])} | "
                  f"{fmt_kj(r['dE_TS_vs_complex_kJ'])} | {fmt_kj(r['dE_reaction_kJ'])} | {g} |")
    (rx.root / "01_ts" / "diagnostics.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


def path_sp_frames(rx: Reaction, n: int, its: int) -> list[int]:
    """IRC frames for the path single points. Every frame, unless reaction.toml sets
    `path_sp_sampling = [a, b]`: every frame within a steps of the TS, every second frame out to
    b steps, every fourth beyond, and always both ends (decided by J. Picione, 2026-11-09, for
    paths where coupled cluster at every frame would take more than about 8 h)."""
    samp = rx.cfg.get("path_sp_sampling")
    if not samp:
        return list(range(n))
    a, b = samp
    keep = []
    for k in range(n):
        d = abs(k - its)
        if d <= a or (d <= b and d % 2 == 0) or (d > b and d % 4 == 0) or k in (0, n - 1):
            keep.append(k)
    return keep


def path_singlepoints(rx: Reaction) -> None:
    """Energy along the IRC at the levels named in `path_levels` (labels from
    [[diagnostic_levels]]), on the working-level geometries. Relative to separated reactants
    at the same level. Writes 02_irc/energies_<label>.csv. Two frames run at a time.
    With `path_sp_sampling` only some frames are computed; the others are filled by a cubic
    spline in s through the computed values, and the `computed` column marks which is which."""
    from concurrent.futures import ThreadPoolExecutor
    from .geometry import read_xyz_frames

    levels = {lv["label"]: lv for lv in rx.cfg.get("diagnostic_levels", [])}
    frames = read_xyz_frames(rx.root / "02_irc" / "irc_full_trj.xyz")
    frag = _read_json(rx.root / "01_ts/fragments.json")
    reac = [n for n, v in frag.items() if v["role"] == "reactants"]
    prod = [n for n, v in frag.items() if v["role"] == "products"]
    flog = orca.RunLog(rx.d("02_irc/path_sp") / "run_log_path_sp.md")
    if not rx.cfg.get("path_levels"):
        print("WARNING no top-level path_levels in reaction.toml; nothing to do")
    for label in rx.cfg.get("path_levels", []):
        lv = levels[label]
        ref = {}
        for n in reac + prod:
            r = orca.parse(rx.root / f"01_ts/diagnostics/{label}/{n}", n, float("nan"))
            if not r.ok:
                raise orca.OrcaError(f"run the diagnostics step first ({label}/{n} missing)")
            ref[n] = r.energy_eh
        e_reac = sum(ref[n] for n in reac)
        e_prod = sum(ref[n] for n in prod)

        def job(k):
            wd = rx.d(f"02_irc/path_sp/{label}/f{k:03d}")
            inp = orca.make_input(f'{lv["method"]} CPCM({rx.cfg["solvent"]})', frames[k],
                                  rx.cfg["charge"], rx.cfg["mult"], blocks=rx.solvent_block,
                                  nprocs=rx.cfg["nprocs"], maxcore_mb=rx.cfg["maxcore_mb"])
            r = orca.run(wd, "f", inp, flog, note=f"{label} IRC frame {k}")
            for tmp in wd.glob("*.tmp*"):
                tmp.unlink()
            if not r.ok:
                raise orca.OrcaError(f"{label} frame {k} failed; see {r.out}")
            return r

        meta = _read_json(rx.root / "02_irc/irc_meta.json")
        todo = path_sp_frames(rx, len(frames), meta["ts_frame"])
        with ThreadPoolExecutor(max_workers=2) as ex:
            res = list(ex.map(job, todo))
        Ec = [r.energy_eh for r in res]
        total = sum(r.wall_s for r in res if r.wall_s == r.wall_s)
        if len(todo) < len(frames):
            from scipy.interpolate import CubicSpline
            s = [float(ln.split(",")[1]) for ln in
                 (rx.root / "02_irc/energies.csv").read_text().splitlines()[1:]]
            spline = CubicSpline([s[k] for k in todo], Ec)
            E = [Ec[todo.index(k)] if k in todo else float(spline(s[k])) for k in range(len(frames))]
        else:
            E = Ec
        rows = ["frame,E_eh,E_kJ_per_mol,E_kcal_per_mol,computed"]
        for k, e in enumerate(E):
            rel = kj(e - e_reac)
            rows.append(f"{k},{e:.10f},{rel:.3f},{rel / 4.184:.3f},{1 if k in todo else 0}")
        (rx.root / "02_irc" / f"energies_{label}.csv").write_text("\n".join(rows) + "\n")
        _write_json(rx.root / "02_irc" / f"energies_{label}.json",
                    {"label": label, "method": lv["method"], "E_separated_products_kJ": kj(e_prod - e_reac),
                     "max_frame": int(np.argmax(E)), "E_max_kJ": kj(max(E) - e_reac),
                     "computed_frames": todo, "n_computed": len(todo), "n_frames": len(frames),
                     "interpolation": None if len(todo) == len(frames) else
                     "cubic spline in s through the computed values"})
        what = (f"all {len(frames)}" if len(todo) == len(frames) else
                f"{len(todo)} of {len(frames)} (path_sp_sampling {rx.cfg['path_sp_sampling']}; "
                f"others by cubic spline in s)")
        rx.log.note(f"Single points at {what} IRC frames at {label} "
                    f"(`02_irc/path_sp/{label}/`, log `02_irc/path_sp/run_log_path_sp.md`): "
                    f"{total:.0f} s of job time; highest frame {int(np.argmax(E))}, "
                    f"{kj(max(E) - e_reac):.1f} kJ/mol above separated reactants.")
        print(f"{label}: max at frame {int(np.argmax(E))}, {fmt_kj(kj(max(E) - e_reac))} "
              f"above separated reactants")


def _irc_summary_energies(out_text: str) -> list[float]:
    """Energies (Eh) from the final 'IRC PATH SUMMARY' table, in trajectory order."""
    i = out_text.rfind("IRC PATH SUMMARY")
    if i < 0:
        raise orca.OrcaError("no IRC PATH SUMMARY in output")
    E = []
    for ln in out_text[i:].splitlines()[1:]:
        p = ln.split()
        if len(p) >= 4 and p[0].isdigit():
            try:
                E.append(float(p[1]))
            except ValueError:
                pass
        elif E and not ln.strip():
            break
    return E


def profile(rx: Reaction) -> None:
    """IRC trajectory ordered reactant → product, mass-weighted path length s, energies.csv, PNG."""
    from .geometry import mass_weighted_step, read_xyz_frames

    wd = rx.root / "02_irc"
    frames = read_xyz_frames(wd / "irc_IRC_Full_trj.xyz")
    joined = _read_json(wd / "irc_path.json") if (wd / "irc_path.json").exists() else None
    if joined:                                  # path joined from per-direction / stopped runs
        E, its = joined["energies_eh"], joined["ts_index"]
    else:
        E = _irc_summary_energies((wd / "irc.out").read_text(errors="replace"))
        its = None
    if len(E) != len(frames):
        raise orca.OrcaError(f"{len(frames)} IRC frames but {len(E)} energies")
    reverse = rx.progress(frames[0]) > rx.progress(frames[-1])   # first frame is product side
    if reverse:
        frames, E = frames[::-1], E[::-1]
    if its is None:
        its = int(np.argmax(E))
    elif reverse:
        its = len(frames) - 1 - its
    # Points within two steps of the TS that lie above it (ORCA's straight-line first step off a
    # curved path) stay in the data files but are left out of the viewer path and the video.
    above = [k for k in range(max(0, its - 2), min(len(E), its + 3)) if k != its and E[k] > E[its]]
    for k in above:
        print(f"NOTE IRC frame {k} lies {kj(E[k] - E[its]):.2f} kJ/mol above the TS point; kept in "
              f"the data files, left out of the viewer path and the video")
    far = [k for k in range(len(E)) if E[k] > E[its] and k not in above]
    if far:
        print(f"WARNING IRC frames {far} lie above the TS point further than two steps from it")
    steps = [0.0] + [mass_weighted_step(a, b) for a, b in zip(frames[:-1], frames[1:])]
    s = np.cumsum(steps)
    s -= s[its]

    frag = _read_json(rx.root / "01_ts/fragments.json")
    e_reac = sum(v["E_eh"] for v in frag.values() if v["role"] == "reactants")
    e_prod = sum(v["E_eh"] for v in frag.values() if v["role"] == "products")
    rel = np.array([kj(e - e_reac) for e in E])

    kb = rx.key_bonds(frames[0].symbols)
    cols = [f"d_{b['label'].replace('–', '-')}_{b['atoms'][0]}_{b['atoms'][1]}_A" for b in kb]
    lines = ["frame,s_bohr_amu05,E_eh,E_kJ_per_mol,E_kcal_per_mol," + ",".join(cols)]
    traj = []
    for k, (g, e) in enumerate(zip(frames, E)):
        ds = ",".join(f"{d:.4f}" for d in rx.bond_lengths(g))
        lines.append(f"{k},{s[k]:.4f},{e:.10f},{rel[k]:.3f},{rel[k] / 4.184:.3f},{ds}")
        traj.append(g.to_xyz(f"frame {k} s={s[k]:.4f} bohr*amu^0.5 E={e:.10f} Eh "
                             f"E_rel={rel[k]:.3f} kJ/mol"))
    (wd / "energies.csv").write_text("\n".join(lines) + "\n")
    (wd / "irc_full_trj.xyz").write_text("".join(traj))
    meta = {"n_frames": len(frames), "ts_frame": its, "s_min": float(s[0]), "s_max": float(s[-1]),
            "E_separated_reactants_kJ": 0.0, "E_separated_products_kJ": kj(e_prod - e_reac),
            "units": {"s": "bohr*amu^0.5", "E": "kJ/mol relative to separated reactants"}}
    if above:
        meta["display_excluded_frames"] = [
            {"frame": k, "s_bohr_amu05": float(s[k]), "dE_above_TS_kJ": kj(E[k] - E[its]),
             "steps_from_TS": abs(k - its)} for k in above]
    if joined and joined.get("stops"):
        # which end of the ordered path each stopped direction is: Full = B reversed + TS + F
        b_first = not reverse
        meta["truncated"] = []
        for tag, st in joined["stops"].items():
            at_start = (tag == "B") == b_first
            meta["truncated"].append({
                "end": "reactant" if at_start else "product", "irc_direction": tag,
                "last_step": st.get("last_step"), "s_bohr_amu05": float(s[0] if at_start else s[-1]),
                "stopped_by": st.get("by"), "reason": st.get("reason", "")})
    _write_json(wd / "irc_meta.json", meta)
    plot_profile(rx, s, rel, its, kj(e_prod - e_reac), wd / "energy_profile.png", flagged=above)
    print(json.dumps(meta, indent=2))


def plot_profile(rx: Reaction, s, rel, its, e_prod_kj, path: Path, level=None, flagged=()) -> None:
    level = level or rx.level()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    ax.plot(s, rel, "-", color="#1f4e79", lw=2.2, zorder=3)
    ax.plot(s, rel, "o", color="#1f4e79", ms=3.2, zorder=4)
    for k in flagged:          # above the TS next to it: kept in the data, left out of the display
        ax.plot(s[k], rel[k], "o", ms=6, mfc="white", mec="#c0392b", mew=1.4, zorder=5)
        ax.annotate(f"first step off the TS, {rel[k] - rel[its]:.1f} kJ/mol above it\n"
                    "(not shown in viewer or video)", (s[k], rel[k]), textcoords="offset points",
                    xytext=(-10, 18), ha="right", fontsize=7, color="#c0392b")
    ax.axhline(0, color="#888", lw=1, ls="--")
    ax.axhline(e_prod_kj, color="#888", lw=1, ls=":")
    x1 = s[-1]
    ax.text(x1, 2, "separated reactants  0 kJ/mol", fontsize=8, color="#555", va="bottom", ha="right")
    # label at the right unless the curve runs along the line there (then at the left end)
    right = abs(rel[-1] - e_prod_kj) > 15
    ax.text(x1 if right else s[0], e_prod_kj + 2, f"separated products  {e_prod_kj:.0f} kJ/mol",
            fontsize=8, color="#555", va="bottom", ha="right" if right else "left")
    for k, lab in ((0, "reactant\ncomplex"), (its, "TS"), (len(s) - 1, "product\ncomplex")):
        ax.annotate(f"{lab}\n{rel[k]:.1f} kJ/mol\n({rel[k] / 4.184:.1f} kcal/mol)",
                    (s[k], rel[k]), textcoords="offset points",
                    xytext=(16, -6) if k == its else (0, -38),
                    ha="left" if k == its else "center", fontsize=8)
    ax.set_xlabel("reaction coordinate s (bohr·amu$^{1/2}$)")
    ax.set_ylabel("electronic energy (kJ/mol)")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / 4.184, lambda x: x * 4.184))
    sec.set_ylabel("(kcal/mol)")
    ax.set_title(f"{rx.cfg['equation']}   IRC, {level}, SMD({rx.cfg['solvent']})", fontsize=10)
    lo = min(rel.min(), e_prod_kj)
    ax.set_ylim(lo - 45, max(rel.max(), 0) + 35)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def profile_levels(rx: Reaction) -> None:
    """One plot: the energy along the same IRC path at the working level and every path level."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = (rx.root / "02_irc/energies.csv").read_text().splitlines()[1:]
    s = np.array([float(r.split(",")[1]) for r in rows])
    curves = [(f"{rx.level()} (working level, IRC)", [float(r.split(",")[3]) for r in rows], "#9aa4b1")]
    colors = ["#1f4e79", "#c0392b", "#2e8b57"]
    for i, label in enumerate(rx.cfg.get("path_levels", [])):
        csv = rx.root / "02_irc" / f"energies_{label}.csv"
        if csv.exists():
            lines = [r.split(",") for r in csv.read_text().splitlines()[1:]]
            E = [float(r[2]) for r in lines]
            done = [k for k, r in enumerate(lines) if len(r) < 5 or r[4] == "1"]
            lab = label.replace("_", "/") + " single points"
            if len(done) < len(lines):
                lab += f" ({len(done)} computed, ● ; spline between)"
            curves.append((lab, E, colors[i % 3], done if len(done) < len(lines) else None))
    curves[0] = curves[0] + (None,)
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    for lab, E, c, done in curves:
        E = np.array(E)
        ax.plot(s, E, "-", color=c, lw=2.2, label=lab)
        if done:
            ax.plot(s[done], E[done], "o", color=c, ms=3)
        k = int(np.argmax(E))
        ax.annotate(f"{E[k]:.1f} kJ/mol ({E[k] / 4.184:.1f} kcal/mol)", (s[k], E[k]),
                    textcoords="offset points", xytext=(8, 4), fontsize=7.5, color=c)
    ax.axhline(0, color="#888", lw=1, ls="--")
    ax.text(s[-1], 2, "separated reactants", fontsize=8, color="#555", ha="right", va="bottom")
    ax.set_xlabel("reaction coordinate s (bohr·amu$^{1/2}$), working-level IRC path")
    ax.set_ylabel("energy relative to separated reactants (kJ/mol)")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / 4.184, lambda x: x * 4.184))
    sec.set_ylabel("(kcal/mol)")
    ax.set_title(f"{rx.cfg['equation']}: same path, energies at several levels, SMD({rx.cfg['solvent']})",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(rx.root / "02_irc" / "energy_profile_levels.png")
    plt.close(fig)


STEPS = {"profile_levels": profile_levels, "fragments": fragments, "ts": ts, "irc": irc, "endpoints": endpoints,
         "singlepoints": singlepoints, "barriers": barriers, "profile": profile,
         "diagnostics": diagnostics, "path_singlepoints": path_singlepoints}
