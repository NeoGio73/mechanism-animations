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
Usage: python -m pipeline.new_level_folder <source folder> <new folder> [--inputs-only]
--inputs-only (2026-11-13): the new folder already exists with its M06-2X stationary points
(pipeline/stationary_m062x.py) but no 00_inputs; write 00_inputs only, TS guess = its own
01_stationary/ts/ts.xyz. The ORCA IRC lines are optional there: the production IRC is the true IRC
(pipeline/true_irc_prod.py), which does not use them.
Nothing is written until the new reaction.toml is complete.
"""
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "reactions"
M_SVP = "M062X D3ZERO ma-def2-SVP def2/J RIJCOSX VeryTightSCF DefGrid3"
M_TZVP = "M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX VeryTightSCF DefGrid3"
SD_PLACEHOLDER = 0.02359
NL = "\n"


def make(src: str, dst: str, inputs_only: bool = False) -> None:
    s, d = ROOT / src, ROOT / dst
    if inputs_only:
        if (d / "00_inputs").exists() or not (d / "01_stationary/ts/ts.xyz").exists():
            raise SystemExit(f"{d}: needs 01_stationary/ts and no 00_inputs")
    elif d.exists():
        raise SystemExit(f"{d} exists")
    t = (s / "00_inputs/reaction.toml").read_text(encoding="utf-8")
    guess_note = "its own 01_stationary/ts" if inputs_only else "M06-2X TS of the sensitivity test"
    sub = [(r"^irc_init_step_s = .*" + NL, ""),        # an older first-step line; the new one is written below
           (r'^method = ".*"', f'method = "{M_SVP}"'),
           (r'^sp_method = ".*"', f'sp_method = "{M_TZVP}"'),
           (r'^ts_guess = ".*"', f'ts_guess = "ts_guess_m062x.xyz"   # {guess_note}'),
           (r"^ts_preopt = .*", "ts_preopt = false"),
           (r"^irc_scale_displ_sd = .*",
            f"irc_scale_displ_sd = {SD_PLACEHOLDER}   # calibrated before the IRC: region median about 0.10" + NL +
            "irc_adapt_scale_displ = false   # J. Picione, 2026-11-11" + NL +
            "irc_init_step_s = 0.10          # first step, bohr·amu½"),
           (r"^path_levels = .*",
            'path_levels = ["M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]'),
           (r"^path_sp_sampling = .*" + NL, ""),
           (r'method = "M062X D3ZERO ma-def2-TZVP def2/J RIJCOSX TightSCF"', f'method = "{M_TZVP}"')]
    optional = [r"^irc_init_step_s", r"^path_sp", r"^ts_preopt"]
    if inputs_only:          # the true IRC uses no ORCA IRC lines; a step-1-only config may lack the rest
        optional += [r"^irc_", r"^ts_guess", r"^path_levels", r'method = "M062X']
    for a, b in sub:
        t, n = re.subn(a, b.replace("\\", "\\\\"), t, count=1, flags=re.M)
        if n == 0 and not any(a.startswith(o) for o in optional):
            raise SystemExit(f"pattern not found in {src}: {a}")
    if inputs_only and not re.search(r"^ts_guess = ", t, flags=re.M):
        t = re.sub(r"^method = ", f'ts_guess = "ts_guess_m062x.xyz"   # {guess_note}' + NL + "method = ",
                   t, count=1, flags=re.M)
    if "irc_init_step_s = 0.10" not in t and not inputs_only:
        raise SystemExit("IRC settings not written")
    head = (f"# Full recomputation at the M06-2X working level (J. Picione, 2026-11-11), from {src}." + NL +
            "# Made by pipeline/new_level_folder.py; the source folder stays untouched." + NL)
    if inputs_only:
        head += ("# 00_inputs written 2026-11-13 (--inputs-only) beside the existing M06-2X stationary points;" + NL +
                 "# TS guess = 01_stationary/ts/ts.xyz. IRC: true IRC (pipeline/true_irc_prod.py)." + NL)
    t = head + t + (NL + "[[diagnostic_levels]]" + NL + 'label = "B3LYP_ma-def2-TZVP"' + NL +
                    'why = "the former working functional, for comparison"' + NL +
                    'method = "B3LYP D3BJ ma-def2-TZVP def2/J RIJCOSX TightSCF"' + NL)
    # VeryTightSCF for every IRC gradient and every per-frame single point (J. Picione, 2026-11-12)
    t = re.sub(r'^((?:sp_)?method = ".*?) TightSCF', r"\1 VeryTightSCF", t, flags=re.M)
    shutil.copytree(s / "00_inputs", d / "00_inputs")
    shutil.copy((d / "01_stationary/ts/ts.xyz") if inputs_only else (s / "01_ts/level_test/M06-2X/ts/ts.xyz"),
                d / "00_inputs/ts_guess_m062x.xyz")
    (d / "00_inputs/reaction.toml").write_text(t, encoding="utf-8")
    print(f"made {d}")


if __name__ == "__main__":
    make(sys.argv[1], sys.argv[2], inputs_only="--inputs-only" in sys.argv)
