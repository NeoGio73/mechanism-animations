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

from . import cores
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
    stopped: str | None = None      # reason, when a watch rule ended the job early

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
        reuse: bool = True, watch=None, poll_s: float = 30.0) -> OrcaResult:
    """Run one ORCA job in workdir/name.inp. Reuses a finished job unless reuse=False.

    watch: optional function called every poll_s seconds while ORCA runs; when it returns a
    reason (string), the job is ended and the reason goes to <name>.stopped.txt and the run log.
    A job ended that way is reused like a finished one."""
    workdir.mkdir(parents=True, exist_ok=True)
    inp_path, out_path = workdir / f"{name}.inp", workdir / f"{name}.out"
    stop_path = workdir / f"{name}.stopped.txt"
    if reuse and out_path.exists() and inp_path.exists() and inp_path.read_text() == inp:
        if "ORCA TERMINATED NORMALLY" in out_path.read_text(errors="replace"):
            res = parse(workdir, name, wall_s=float("nan"))
            print(f"[reuse] {workdir.name}/{name}: finished earlier")
            return res
        if stop_path.exists():
            res = parse(workdir, name, wall_s=float("nan"))
            res.stopped = stop_path.read_text(encoding="utf-8").strip()
            print(f"[reuse] {workdir.name}/{name}: stopped earlier ({res.stopped})")
            return res
    _move_crashed_run(workdir, name, out_path, stop_path)
    inp_path.write_text(inp)
    stop_path.unlink(missing_ok=True)
    env = os.environ.copy()
    env["PATH"] = f"{MSMPI_BIN};{ORCA_DIR};" + env.get("PATH", "")
    # core and memory budget across all lanes (pipeline/cores.py): wait for a free slot first
    try:
        what = str(workdir.relative_to(Path(__file__).resolve().parents[1] / "reactions")).replace("\\", "/")
    except ValueError:
        what = workdir.name
    rid = cores.acquire(*cores.request(inp), what=f"{what}/{name}", folder=str(workdir.resolve()))
    print(f"[run] {workdir.name}/{name} ...", flush=True)
    t0 = time.time()
    reason = None
    try:
        with open(out_path, "w") as fh:
            proc = subprocess.Popen([str(ORCA_DIR / "orca.exe"), inp_path.name], cwd=workdir,
                                    stdout=fh, stderr=subprocess.STDOUT, env=env)
            while True:
                try:
                    proc.wait(timeout=poll_s if watch else None)
                    break
                except subprocess.TimeoutExpired:
                    reason = watch()
                    if reason:
                        # ORCA starts MPI child processes: end the whole process tree
                        subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                                       capture_output=True)
                        proc.wait()
                        break
    finally:
        done = "ORCA TERMINATED NORMALLY" in out_path.read_text(errors="replace") if out_path.exists() else False
        cores.release(rid, "" if done else "  (did not terminate normally)")
    res = parse(workdir, name, wall_s=time.time() - t0)
    if reason:
        res.stopped = reason
        stop_path.write_text(reason + "\n", encoding="utf-8")
        note = f"{note}; stopped early: {reason}" if note else f"stopped early: {reason}"
    log.add(res, inp.splitlines()[0], note)
    print(f"[done] {workdir.name}/{name}: ok={res.ok} wall={res.wall_s:.0f} s "
          f"E={res.energy_eh}", flush=True)
    return res


def _move_crashed_run(workdir: Path, name: str, out_path: Path, stop_path: Path) -> None:
    """Clean-folder rule (J. Picione, 2026-11-11/12): a rerun after a crash must not start among the
    crashed run's files (ORCA would even take its .gbw as the starting guess). Jobs of one production
    chain may share the folder, so the layout stays; the crashed run's files (<name>.*, <name>_*, and
    the fixed-name *.tmp files, which only a crashed job leaves behind) move whole to
    crashed_<time>_<name>/."""
    if not out_path.exists() or stop_path.exists():
        return
    if "ORCA TERMINATED NORMALLY" in out_path.read_text(errors="replace"):
        return
    dest = workdir / f"crashed_{time.strftime('%H%M%S')}_{name}"
    dest.mkdir()
    # other chain jobs with longer names (plus_min beside plus) keep their files
    longer = [q.stem for q in workdir.glob("*.inp") if q.stem != name and q.stem.startswith(name)]
    for p in list(workdir.iterdir()):
        if not p.is_file() or any(p.name.startswith(o + ".") or p.name.startswith(o + "_") for o in longer):
            continue
        if p.name.startswith(f"{name}.") or p.name.startswith(f"{name}_") or p.suffix == ".tmp":
            p.rename(dest / p.name)
    print(f"[crashed earlier] {workdir.name}/{name}: earlier run moved to {dest.name}/", flush=True)


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
        if res.stopped:
            result = "stopped early (see note)"
        elif not res.ok:
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
