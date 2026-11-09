"""Write ORCA inputs, run ORCA, parse outputs, keep the per-reaction run log."""
from __future__ import annotations

import datetime as dt
import os
import re
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .geometry import Geometry

ORCA_DIR = Path(os.environ.get("ORCA_DIR", r"F:\ORCA"))
MSMPI_BIN = Path(os.environ.get("MSMPI_BIN", r"C:\Program Files\Microsoft MPI\Bin"))


class OrcaError(RuntimeError):
    pass


@dataclass
class OrcaResult:
    name: str
    workdir: Path
    wall_s: float
    ok: bool
    energy_eh: float | None = None
    gibbs_eh: float | None = None
    opt_converged: bool | None = None
    freqs_cm: list[float] = field(default_factory=list)

    @property
    def out(self) -> Path:
        return self.workdir / f"{self.name}.out"

    @property
    def imaginary(self) -> list[float]:
        return [f for f in self.freqs_cm if f < 0]


def make_input(keywords: str, geom: Geometry, charge: int, mult: int,
               blocks: str = "", nprocs: int = 1, maxcore_mb: int = 3000) -> str:
    pal = f"%pal nprocs {nprocs} end\n" if nprocs > 1 else ""
    return (f"! {keywords}\n{pal}%maxcore {maxcore_mb}\n{blocks.strip()}\n\n"
            f"* xyz {charge} {mult}\n{geom.xyz_block()}\n*\n")


def run(workdir: Path, name: str, inp: str, log: "RunLog", note: str = "",
        reuse: bool = True) -> OrcaResult:
    """Run one ORCA job in workdir/name.inp. Reuses a finished job unless reuse=False."""
    workdir.mkdir(parents=True, exist_ok=True)
    inp_path, out_path = workdir / f"{name}.inp", workdir / f"{name}.out"
    if reuse and out_path.exists() and inp_path.exists() and inp_path.read_text() == inp \
            and "ORCA TERMINATED NORMALLY" in out_path.read_text(errors="replace"):
        res = parse(workdir, name, wall_s=float("nan"))
        print(f"[reuse] {workdir.name}/{name}: finished earlier")
        return res
    inp_path.write_text(inp)
    env = os.environ.copy()
    env["PATH"] = f"{MSMPI_BIN};{ORCA_DIR};" + env.get("PATH", "")
    print(f"[run] {workdir.name}/{name} ...", flush=True)
    t0 = time.time()
    with open(out_path, "w") as fh:
        subprocess.run([str(ORCA_DIR / "orca.exe"), inp_path.name], cwd=workdir,
                       stdout=fh, stderr=subprocess.STDOUT, env=env)
    res = parse(workdir, name, wall_s=time.time() - t0)
    log.add(res, inp.splitlines()[0], note)
    print(f"[done] {workdir.name}/{name}: ok={res.ok} wall={res.wall_s:.0f} s "
          f"E={res.energy_eh}", flush=True)
    return res


def parse(workdir: Path, name: str, wall_s: float) -> OrcaResult:
    text = (workdir / f"{name}.out").read_text(errors="replace")
    ok = "ORCA TERMINATED NORMALLY" in text
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    g = re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", text)
    opt = None
    if "GEOMETRY OPTIMIZATION" in text:
        opt = "THE OPTIMIZATION HAS CONVERGED" in text
    res = OrcaResult(name, workdir, wall_s, ok,
                     float(e[-1]) if e else None, float(g[-1]) if g else None, opt)
    hess = workdir / f"{name}.hess"
    if hess.exists():
        res.freqs_cm = read_hess(hess)["freqs"].tolist()
    return res


def read_hess(path: Path) -> dict:
    """Frequencies (cm-1), normal modes (3N x 3N, columns), symbols, masses from a .hess file."""
    lines = Path(path).read_text().splitlines()
    idx = {ln.strip(): i for i, ln in enumerate(lines) if ln.startswith("$")}
    out: dict = {}
    i = idx["$vibrational_frequencies"]
    n = int(lines[i + 1])
    out["freqs"] = np.array([float(lines[i + 2 + k].split()[1]) for k in range(n)])
    i = idx["$normal_modes"]
    nr, nc = map(int, lines[i + 1].split())
    modes = np.zeros((nr, nc))
    j = i + 2
    while j < len(lines) and not lines[j].startswith("$") and lines[j].strip():
        cols = [int(c) for c in lines[j].split()]
        for r in range(nr):
            vals = lines[j + 1 + r].split()[1:]
            for c, v in zip(cols, vals):
                modes[r, c] = float(v)
        j += 1 + nr
    out["modes"] = modes
    i = idx["$atoms"]
    na = int(lines[i + 1])
    out["symbols"] = [lines[i + 2 + k].split()[0] for k in range(na)]
    out["masses"] = np.array([float(lines[i + 2 + k].split()[1]) for k in range(na)])
    return out


class RunLog:
    """Appends one row per job to run_log.md. Notes record what changed and why."""

    def __init__(self, path: Path):
        self.path = path
        if not path.exists():
            path.write_text(
                "# Run log\n\nEvery ORCA job, its wall time, and what changed when a job "
                "failed. Times are wall-clock on the Omen (Ryzen 7 9800X3D).\n\n"
                "| When | Job | Keywords | Wall time | Result | Note |\n"
                "|---|---|---|---|---|---|\n")

    def add(self, res: OrcaResult, keyword_line: str, note: str = "") -> None:
        when = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
        if not res.ok:
            result = "FAILED"
        else:
            parts = [f"E = {res.energy_eh:.8f} Eh"] if res.energy_eh is not None else []
            if res.opt_converged is not None:
                parts.append("opt converged" if res.opt_converged else "opt NOT converged")
            if res.freqs_cm:
                parts.append(f"{len(res.imaginary)} imaginary freq"
                             + (f" ({', '.join(f'{f:.1f}' for f in res.imaginary)} cm⁻¹)"
                                if res.imaginary else ""))
            result = "; ".join(parts)
        rel = res.workdir.relative_to(self.path.parent).as_posix()
        row = (f"| {when} | `{rel}/{res.name}` | `{keyword_line.lstrip('! ')}` | "
               f"{res.wall_s:.0f} s | {result} | {note} |\n")
        with open(self.path, "a", encoding="utf-8") as fh:
            fh.write(row)

    def note(self, text: str) -> None:
        with open(self.path, "a", encoding="utf-8") as fh:
            fh.write(f"| {dt.datetime.now():%Y-%m-%d %H:%M} | — | — | — | — | {text} |\n")
