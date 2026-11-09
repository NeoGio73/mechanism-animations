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

    def species(self):
        for role in ("reactants", "products"):
            for s in self.cfg[role]:
                yield role, s

    def key(self, name: str) -> int:
        return self.cfg["key_atoms"][name]

    def bonds(self, g: Geometry) -> tuple[float, float]:
        """(nucleophile–center, center–leaving group) distances in Å."""
        return (g.distance(self.key("nucleophile"), self.key("center")),
                g.distance(self.key("center"), self.key("leaving_group")))


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
        g = from_smiles(s["smiles"])
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
    return {"freq_cm": float(h["freqs"][k]),
            "central_atom_share": float(share[c]),
            "central_atom_alignment_with_axis": float(along_c),
            "largest_atom": geom.symbols[int(np.argmax(share))] + str(int(np.argmax(share)))}


def ts(rx: Reaction) -> None:
    """OptTS + Freq from the TS guess, exact Hessian at start and every 5 steps."""
    wd = rx.d("01_ts/ts")
    g = read_xyz(rx.root / "00_inputs" / rx.cfg["ts_guess"])
    blocks = "%geom Calc_Hess true Recalc_Hess 5 MaxIter 100 end"
    r = orca.run(wd, "ts", rx.inp("OptTS Freq", g, blocks=blocks), rx.log, note="TS search from guess")
    if not r.ok or not r.opt_converged:
        raise orca.OrcaError(f"TS optimization failed; see {r.out}")
    gts = read_xyz(wd / "ts.xyz")
    mode = check_ts_mode(rx, wd / "ts.hess", gts)
    d_nu, d_lg = rx.bonds(gts)
    out = {"E_eh": r.energy_eh, "G_eh": r.gibbs_eh, "imaginary_cm": r.imaginary,
           "n_imaginary": len(r.imaginary), "mode_check": mode,
           "d_nu_center_A": d_nu, "d_center_lg_A": d_lg,
           "angle_nu_center_lg_deg": gts.angle(rx.key("nucleophile"), rx.key("center"),
                                               rx.key("leaving_group"))}
    _write_json(wd.parent / "ts.json", out)
    print(json.dumps(out, indent=2))
    if len(r.imaginary) != 1:
        raise orca.OrcaError(f"TS has {len(r.imaginary)} imaginary frequencies, need exactly 1")


# ---------------------------------------------------------------- IRC
def irc(rx: Reaction) -> None:
    """IRC in both directions from the optimized TS, reading the TS Hessian."""
    wd = rx.d("02_irc")
    shutil.copy(rx.root / "01_ts/ts/ts.hess", wd / "ts.hess")
    g = read_xyz(rx.root / "01_ts/ts/ts.xyz")
    step = rx.cfg.get("irc_scale_displ_sd")       # smaller → more frames along the path
    extra = f"  Scale_Displ_SD {step}\n" if step else ""
    blocks = ('%irc\n  MaxIter 300\n  Direction both\n  InitHess read\n'
              f'  Hess_Filename "ts.hess"\n{extra}  PrintLevel 1\nend')
    note = "IRC both directions" + (f", Scale_Displ_SD {step}" if step else ", default step")
    r = orca.run(wd, "irc", rx.inp("IRC", g, blocks=blocks), rx.log, note=note)
    if not r.ok:
        raise orca.OrcaError(f"IRC failed; see {r.out}")


def endpoints(rx: Reaction) -> None:
    """Optimize the two IRC end geometries into the reactant and product complexes."""
    wd = rx.d("02_irc")
    out = {}
    for end in ("F", "B"):
        g = read_xyz(wd / f"irc_IRC_{end}.xyz")
        d_nu, d_lg = rx.bonds(g)
        label = "reactant_complex" if d_nu > d_lg else "product_complex"
        ew = rx.d(f"02_irc/endpoints/{label}")
        g.write(ew / "start.xyz", f"IRC {end} end")
        r = orca.run(ew, label, rx.inp("Opt Freq", g), rx.log, note=f"from IRC {end} end")
        if not r.ok or not r.opt_converged:
            raise orca.OrcaError(f"{label} optimization failed; see {r.out}")
        name = label
        if r.imaginary:
            # A minimum with a small imaginary mode (flat surface): push the structure along that
            # mode (largest atom moves 0.1 Å) and re-optimize with tighter criteria, once.
            h = orca.read_hess(ew / f"{label}.hess")
            v = h["modes"][:, int(np.argmin(h["freqs"]))].reshape(-1, 3)
            gp = read_xyz(ew / f"{label}.xyz")
            gp = Geometry(gp.symbols, gp.coords + 0.1 * v / np.linalg.norm(v, axis=1).max())
            name = f"{label}_tight"
            r = orca.run(ew, name, rx.inp("TightOpt Freq", gp), rx.log,
                         note=f"rerun: {label} had imaginary mode(s) {[round(f, 1) for f in r.imaginary]} "
                              "cm-1; displaced 0.1 Å along it, TightOpt")
            if not r.ok or not r.opt_converged:
                raise orca.OrcaError(f"{name} optimization failed; see {r.out}")
        go = read_xyz(ew / f"{name}.xyz")
        d_nu2, d_lg2 = rx.bonds(go)
        out[label] = {"irc_end": end, "E_eh": r.energy_eh, "G_eh": r.gibbs_eh,
                      "imaginary_cm": r.imaginary, "xyz": f"02_irc/endpoints/{label}/{name}.xyz",
                      "d_nu_center_A": d_nu2, "d_center_lg_A": d_lg2}
        if r.imaginary:
            print(f"WARNING {label} still has imaginary frequencies {r.imaginary}")
    if set(out) != {"reactant_complex", "product_complex"}:
        raise orca.OrcaError(f"IRC ends do not give one reactant and one product side: {out}")
    _write_json(wd / "endpoints.json", out)
    print(json.dumps(out, indent=2))


# ---------------------------------------------------------------- big-basis single points
def singlepoints(rx: Reaction) -> None:
    """def2-TZVP single points (same functional, same solvent) on all stationary points."""
    out = {}
    for name, (xyz, q, m) in _stationary_points(rx).items():
        r = orca.run(rx.d(f"01_ts/sp_tzvp/{name}"), name,
                     rx.inp("", read_xyz(xyz), q, m, sp=True), rx.log, note="def2-TZVP single point")
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
           "def2-SVP": table(E_svp), "def2-TZVP": table(sp)}
    G = {n: v["G_eh"] for n, v in frag.items()}
    G.update(ts=tsj["G_eh"], reactant_complex=ep["reactant_complex"]["G_eh"],
             product_complex=ep["product_complex"]["G_eh"])
    if all(v is not None for v in G.values()):
        out["def2-SVP_gibbs_298K_1atm"] = table(G)
    _write_json(rx.root / "02_irc" / "barriers.json", out)
    for lvl in ("def2-SVP", "def2-TZVP"):
        t = out[lvl]
        print(f"{lvl}: barrier from complex {fmt_kj(t['barrier_from_reactant_complex'])}, "
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

    Free-energy barriers use the working-level thermal corrections (G - E at def2-SVP)
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

    levels = [{"label": "B3LYP_def2-SVP (working)", "E": {
        **{n: v["E_eh"] for n, v in frag.items()}, "ts": tsj["E_eh"],
        "reactant_complex": ep["reactant_complex"]["E_eh"],
        "product_complex": ep["product_complex"]["E_eh"]}}]
    levels.append({"label": "B3LYP_def2-TZVP (protocol)",
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
    _write_json(rx.root / "01_ts" / "diagnostics.json", {"units": "kJ/mol", "rows": rows})
    md = ["| Level (SMD water, def2-SVP geometries) | TS − separated reactants, ΔE | "
          "TS − reactant complex, ΔE | Reaction energy, ΔE | ΔG‡ (298 K, 1 M) |",
          "|---|---|---|---|---|"]
    for r in rows:
        g = fmt_kj(r["dG_TS_vs_separated_1M_kJ"]) if r["dG_TS_vs_separated_1M_kJ"] is not None else "—"
        md.append(f"| {r['level']} | {fmt_kj(r['dE_TS_vs_separated_kJ'])} | "
                  f"{fmt_kj(r['dE_TS_vs_complex_kJ'])} | {fmt_kj(r['dE_reaction_kJ'])} | {g} |")
    (rx.root / "01_ts" / "diagnostics.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


def path_singlepoints(rx: Reaction) -> None:
    """Energy at every IRC frame at the levels named in `path_levels` (labels from
    [[diagnostic_levels]]), on the working-level geometries. Relative to separated reactants
    at the same level. Writes 02_irc/energies_<label>.csv. Two frames run at a time."""
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

        with ThreadPoolExecutor(max_workers=2) as ex:
            res = list(ex.map(job, range(len(frames))))
        E = [r.energy_eh for r in res]
        total = sum(r.wall_s for r in res if r.wall_s == r.wall_s)
        rows = ["frame,E_eh,E_kJ_per_mol,E_kcal_per_mol"]
        for k, e in enumerate(E):
            rel = kj(e - e_reac)
            rows.append(f"{k},{e:.10f},{rel:.3f},{rel / 4.184:.3f}")
        (rx.root / "02_irc" / f"energies_{label}.csv").write_text("\n".join(rows) + "\n")
        _write_json(rx.root / "02_irc" / f"energies_{label}.json",
                    {"label": label, "method": lv["method"], "E_separated_products_kJ": kj(e_prod - e_reac),
                     "max_frame": int(np.argmax(E)), "E_max_kJ": kj(max(E) - e_reac)})
        rx.log.note(f"Single points at all {len(frames)} IRC frames at {label} "
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
    E = _irc_summary_energies((wd / "irc.out").read_text(errors="replace"))
    if len(E) != len(frames):
        raise orca.OrcaError(f"{len(frames)} IRC frames but {len(E)} energies")
    d0 = rx.bonds(frames[0])
    if d0[0] < d0[1]:                       # first frame is product side: reverse
        frames, E = frames[::-1], E[::-1]
    its = int(np.argmax(E))
    steps = [0.0] + [mass_weighted_step(a, b) for a, b in zip(frames[:-1], frames[1:])]
    s = np.cumsum(steps)
    s -= s[its]

    frag = _read_json(rx.root / "01_ts/fragments.json")
    e_reac = sum(v["E_eh"] for v in frag.values() if v["role"] == "reactants")
    e_prod = sum(v["E_eh"] for v in frag.values() if v["role"] == "products")
    rel = np.array([kj(e - e_reac) for e in E])

    lines = ["frame,s_bohr_amu05,E_eh,E_kJ_per_mol,E_kcal_per_mol,d_nu_center_A,d_center_lg_A"]
    traj = []
    for k, (g, e) in enumerate(zip(frames, E)):
        dn, dl = rx.bonds(g)
        lines.append(f"{k},{s[k]:.4f},{e:.10f},{rel[k]:.3f},{rel[k] / 4.184:.3f},{dn:.4f},{dl:.4f}")
        traj.append(g.to_xyz(f"frame {k} s={s[k]:.4f} bohr*amu^0.5 E={e:.10f} Eh "
                             f"E_rel={rel[k]:.3f} kJ/mol"))
    (wd / "energies.csv").write_text("\n".join(lines) + "\n")
    (wd / "irc_full_trj.xyz").write_text("".join(traj))
    meta = {"n_frames": len(frames), "ts_frame": its, "s_min": float(s[0]), "s_max": float(s[-1]),
            "E_separated_reactants_kJ": 0.0, "E_separated_products_kJ": kj(e_prod - e_reac),
            "units": {"s": "bohr*amu^0.5", "E": "kJ/mol relative to separated reactants"}}
    _write_json(wd / "irc_meta.json", meta)
    plot_profile(rx, s, rel, its, kj(e_prod - e_reac), wd / "energy_profile.png")
    print(json.dumps(meta, indent=2))


def plot_profile(rx: Reaction, s, rel, its, e_prod_kj, path: Path, level="B3LYP-D3(BJ)/def2-SVP") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    ax.plot(s, rel, "-", color="#1f4e79", lw=2.2, zorder=3)
    ax.plot(s, rel, "o", color="#1f4e79", ms=3.2, zorder=4)
    ax.axhline(0, color="#888", lw=1, ls="--")
    ax.axhline(e_prod_kj, color="#888", lw=1, ls=":")
    x1 = s[-1]
    ax.text(x1, 2, "separated reactants  0 kJ/mol", fontsize=8, color="#555", va="bottom", ha="right")
    ax.text(x1, e_prod_kj + 2, f"separated products  {e_prod_kj:.0f} kJ/mol", fontsize=8,
            color="#555", va="bottom", ha="right")
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
    curves = [("B3LYP-D3(BJ)/def2-SVP (working level, IRC)", [float(r.split(",")[3]) for r in rows], "#9aa4b1")]
    colors = ["#1f4e79", "#c0392b", "#2e8b57"]
    for i, label in enumerate(rx.cfg.get("path_levels", [])):
        csv = rx.root / "02_irc" / f"energies_{label}.csv"
        if csv.exists():
            E = [float(r.split(",")[2]) for r in csv.read_text().splitlines()[1:]]
            curves.append((label.replace("_", "/") + " single points", E, colors[i % 3]))
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    for lab, E, c in curves:
        E = np.array(E)
        ax.plot(s, E, "-", color=c, lw=2.2, label=lab)
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
