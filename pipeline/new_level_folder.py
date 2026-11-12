"""Make a new reaction folder for a full recomputation at the M06-2X working level (J. Picione,
2026-11-11), keeping the old folder untouched.

Copies 00_inputs from the source folder and rewrites reaction.toml:
  method       M06-2X-D3(0)/ma-def2-SVP, DefGrid3 (working level: geometries, TS, IRC, orbitals)
  sp_method    M06-2X-D3(0)/ma-def2-TZVP, DefGrid3
  path levels  M06-2X/ma-def2-TZVP (DefGrid3), B3LYP-D3(BJ)/ma-def2-TZVP and DLPNO-CCSD(T)/ma-def2-TZVP
               at every IRC frame (no sampling)
  TS guess     the M06-2X TS of the functional sensitivity test (01_ts/level_test/M06-2X/ts/ts.xyz)
  IRC          Adapt_Scale_Displ false, first step 0.10 bohr·amu½; Scale_Displ_SD set by calibration
               before the IRC (placeholder: the SN2 value)
Usage: python -m pipeline.new_level_folder <source folder> <new folder>
"""
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "reactions"
M_SVP = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX VeryTightSCF DefGrid3"
M_TZVP = "M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX VeryTightSCF DefGrid3"
SD_PLACEHOLDER = 0.02359


def make(src: str, dst: str) -> None:
    s, d = ROOT / src, ROOT / dst
    if d.exists():
        raise SystemExit(f"{d} exists")
    shutil.copytree(s / "00_inputs", d / "00_inputs")
    shutil.copy(s / "01_ts/level_test/M06-2X/ts/ts.xyz", d / "00_inputs/ts_guess_m062x.xyz")
    t = (d / "00_inputs/reaction.toml").read_text(encoding="utf-8")
    sub = [(r'^irc_init_step_s = .*\n', ''),          # an older first-step line; the new one is written below
           (r'^method = ".*"', f'method = "{M_SVP}"'),
           (r'^sp_method = ".*"', f'sp_method = "{M_TZVP}"'),
           (r'^ts_guess = ".*"', 'ts_guess = "ts_guess_m062x.xyz"   # M06-2X TS of the sensitivity test'),
           (r'^ts_preopt = .*', 'ts_preopt = false'),
           (r'^irc_scale_displ_sd = .*',
            f'irc_scale_displ_sd = {SD_PLACEHOLDER}   # calibrated before the IRC: region median about 0.10\n'
            'irc_adapt_scale_displ = false   # J. Picione, 2026-11-11\n'
            'irc_init_step_s = 0.10          # first step, bohr·amu½'),
           (r'^path_levels = .*',
            'path_levels = ["M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]'),
           (r'^path_sp_sampling = .*\n', ''),
           (r'method = "M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF"', f'method = "{M_TZVP}"')]
    for a, b in sub:
        t, n = re.subn(a, b, t, count=1, flags=re.M)
        if n == 0 and not a.startswith(r"^irc_init_step_s") and not a.startswith(r"^path_sp") \
                and not a.startswith(r"^ts_preopt"):
            raise SystemExit(f"pattern not found in {src}: {a}")
    if "irc_init_step_s = 0.10" not in t:
        raise SystemExit("IRC settings not written")
    head = (f"# Full recomputation at the M06-2X working level (J. Picione, 2026-11-11), from {src}.\n"
            "# Made by pipeline/new_level_folder.py; the source folder stays untouched.\n")
    t = head + t + ('\n[[diagnostic_levels]]\nlabel = "B3LYP_ma-def2-TZVP"\n'
                    'why = "the former working functional, for comparison"\n'
                    'method = "B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF"\n')
    # VeryTightSCF for every IRC gradient and every per-frame single point (J. Picione, 2026-11-12)
    t = re.sub(r'^((?:sp_)?method = ".*?) TightSCF', r'\1 VeryTightSCF', t, flags=re.M)
    (d / "00_inputs/reaction.toml").write_text(t, encoding="utf-8")
    print(f"made {d}")


if __name__ == "__main__":
    make(sys.argv[1], sys.argv[2])
