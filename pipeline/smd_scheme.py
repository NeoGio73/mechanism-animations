"""SMD-scheme check ("HF-density SMD check"), J. Picione's post-E2 order step 4 (2026-11-23).

At the stationary points only, with the same reference states as the paper's barriers:
  Scheme A (current): DLPNO-CCSD(T)/ma-def2-TZVP with SMD (the single points already in 01_ts/diagnostics).
  Scheme B: gas-phase DLPNO-CCSD(T)/ma-def2-TZVP + [E(SMD) - E(gas)] from M06-2X/ma-def2-TZVP, at the same
            geometries, with the same thermal corrections.
ΔG‡ and ΔrG for all six final reactions under both schemes, and the difference. If every difference is
≤ 2 kJ/mol he writes one sentence in Methods; if any is larger, STOP and report before changing any table.

New jobs: the gas-phase twin of every SMD single point at M06-2X/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP that
the barriers use (01_ts/diagnostics/<level>/<species>/, and the zero conformers' sp/<level>/ folders). The gas
input is the SMD input with only the solvent removed (the CPCM(...) keyword and the %cpcm block); geometry,
charge, multiplicity, method, basis, grids and SCF settings are identical. Each source job is checked first: its
.out ends normally and its energy equals the value the paper uses (01_ts/diagnostics.json) to 1e-6 Eh.

Jobs go to <reaction>/01_ts/smd_scheme/gas/<level>/<species>/. Two lanes take jobs from one shared list (a job is
claimed by creating its folder's .claim directory); each job uses 4 cores through the ordinary core queue.

Usage: python -X utf8 -m pipeline.smd_scheme list | run <lane name> [--retry] | report
"""
from __future__ import annotations

import json
import re
import sys
import time
import tomllib
from pathlib import Path

from . import orca
from .steps import Reaction

ROOT = Path(__file__).resolve().parents[1]
REACTIONS = ["sn2_hydroxide_bromomethane_m062x", "e2_ethoxide_2_bromopropane_m062x",
             "nucadd_cyanide_acetone_m062x", "da_butadiene_ethylene_m062x",
             "hshift_3_methyl_2_butyl_m062x", "hshift_2_3_dimethyl_2_butyl_m062x"]
LEVELS = ["DLPNO-CCSD(T)_ma-def2-TZVP", "M06-2X_ma-def2-TZVP"]
CONSOLE = ROOT / "docs/console_smd_scheme_{}.txt"


def say(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}", flush=True)


def _final_energy(out: Path) -> float:
    tx = out.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in tx:
        raise RuntimeError(f"{out}: did not end normally")
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", tx)[-1])


def gas_input(smd: str) -> str:
    """The SMD input with the solvent model removed and nothing else changed."""
    lines = smd.splitlines()
    out = []
    for ln in lines:
        if ln.lstrip().startswith("!"):
            ln = re.sub(r"\s+CPCM\([^)]*\)", "", ln, flags=re.I)
            ln = re.sub(r"\s+SMD\([^)]*\)", "", ln, flags=re.I)
        if re.match(r"\s*%cpcm\b", ln, flags=re.I):
            if not re.search(r"\bend\s*$", ln, flags=re.I):
                raise ValueError("multi-line %cpcm block: not handled")
            continue
        out.append(ln)
    g = "\n".join(out) + "\n"
    assert not re.search(r"cpcm|smd", g, flags=re.I), g[:300]
    return g


def jobs() -> list[dict]:
    """Every gas-phase twin to compute, largest first (DLPNO before M06-2X, more atoms first)."""
    js = []
    for name in REACTIONS:
        R = ROOT / "reactions" / name
        diag = json.loads((R / "01_ts/diagnostics.json").read_text(encoding="utf-8"))["energies_eh"]
        cfg = tomllib.loads((R / "00_inputs/reaction.toml").read_text(encoding="utf-8"))
        for lv in LEVELS:
            for sp, e in diag[lv].items():
                d = R / "01_ts/diagnostics" / lv / sp
                inps = sorted(d.glob("*.inp"))
                assert len(inps) == 1, (d, inps)
                js.append({"reaction": name, "level": lv, "species": sp, "src": inps[0], "E_smd": e})
            for zc in cfg.get("zero_conformers", []):
                res = json.loads((R / zc["result"]).read_text(encoding="utf-8"))
                d = R / Path(res["job"]).parent / "sp" / lv
                inps = sorted(d.glob("*.inp"))
                assert len(inps) == 1, (d, inps)
                js.append({"reaction": name, "level": lv, "species": zc["name"], "src": inps[0],
                           "E_smd": res["sp"][lv]})
    for j in js:
        e = _final_energy(j["src"].with_suffix(".out"))
        if abs(e - j["E_smd"]) > 1e-6:
            raise RuntimeError(f"{j['src']}: energy {e} differs from the value used ({j['E_smd']})")
        j["natoms"] = int(j["src"].read_text().split("* xyz")[1].count("\n")) - 2
        j["wd"] = ROOT / "reactions" / j["reaction"] / "01_ts/smd_scheme/gas" / j["level"] / j["species"]
    js.sort(key=lambda j: (j["level"] != LEVELS[0], -j["natoms"]))
    return js


def unfinished(js: list[dict]) -> list[dict]:
    """Jobs without a normally terminated output (not started, crashed or failed)."""
    out = []
    for j in js:
        o = j["wd"] / f"{j['src'].stem}.out"
        if not (o.exists() and "ORCA TERMINATED NORMALLY" in o.read_text(errors="replace")):
            out.append(j)
    return out


def run(lane: str, retry: bool = False) -> None:
    js = jobs()
    if retry:              # only when no other lane is alive (pipeline/after_smd.py): release unfinished claims
        for j in unfinished(js):
            c = j["wd"] / ".claim"
            if c.exists():
                for f in c.iterdir():
                    f.unlink()
                c.rmdir()
        say(f"lane {lane}: retry of {len(unfinished(js))} unfinished jobs")
    say(f"lane {lane}: {len(js)} gas-phase jobs in the shared list")
    for j in js:
        j["wd"].mkdir(parents=True, exist_ok=True)
        try:
            (j["wd"] / ".claim").mkdir()
        except FileExistsError:
            continue
        (j["wd"] / ".claim" / "lane.txt").write_text(lane, encoding="utf-8")
        rx = Reaction.load(ROOT / "reactions" / j["reaction"])
        t0 = time.time()
        say(f"lane {lane}: {j['reaction']} {j['level']} {j['species']} ({j['natoms']} atoms) start")
        r = orca.run(j["wd"], j["src"].stem, gas_input(j["src"].read_text()), rx.log,
                     note=f"SMD-scheme check: gas-phase twin of {j['src'].relative_to(rx.root).as_posix()}")
        if not r.ok:
            say(f"lane {lane}: {j['reaction']} {j['level']} {j['species']} FAILED; see {r.out}")
            continue
        say(f"lane {lane}: {j['reaction']} {j['level']} {j['species']} done, {(time.time() - t0) / 60:.1f} min, "
            f"E(gas) {r.energy_eh:.9f} Eh")
    say(f"lane {lane}: finished")


def main() -> None:
    cmd = sys.argv[1]
    if cmd == "list":
        js = jobs()
        for j in js:
            print(j["reaction"], j["level"], j["species"], j["natoms"], j["src"].relative_to(ROOT).as_posix())
        print(len(js), "jobs; source energies all match the values used")
    elif cmd == "run":
        run(sys.argv[2], retry="--retry" in sys.argv[3:])
    elif cmd == "report":
        from .smd_scheme_report import report
        report()


if __name__ == "__main__":
    main()
