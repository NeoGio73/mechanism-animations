"""What runs after the SMD-scheme lanes (J. Picione, 2026-11-25 evening):

1. Wait until no SMD-scheme lane process is alive.
2. If any gas-phase job has no normally terminated output (a lane crashed, or a job failed), run ONE retry lane
   (pipeline.smd_scheme run R --retry: unfinished jobs only; a crashed job is moved aside by orca.run, never deleted).
3. Write the SMD-scheme report (pipeline.smd_scheme report). If jobs are still missing after the retry, write that
   instead (docs/results/smd_scheme_check_FAILED.md): no result from incomplete jobs.
4. Start the cyanide explicit-water test in every case, so the machine is not idle overnight (same 36 h estimate
   stop). When it ends nothing else starts.

Console: docs/console_after_smd.txt. Usage: python -X utf8 -m pipeline.after_smd
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

from .smd_scheme import jobs, unfinished
from .start_after import alive

ROOT = Path(__file__).resolve().parents[1]
PY = str(ROOT / ".venv/Scripts/python.exe")
WATER = ROOT / "reactions/nucadd_cyanide_acetone_m062x/04_explicit_water/console.txt"


def say(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}", flush=True)


def main() -> None:
    say("waiting for the SMD-scheme lanes to end")
    while alive("pipeline.smd_scheme run"):
        time.sleep(60)
    time.sleep(5)
    miss = unfinished(jobs())
    if miss:
        say(f"{len(miss)} gas-phase jobs unfinished; one retry lane")
        with open(ROOT / "docs/console_smd_scheme_R.txt", "a", encoding="utf-8") as fh:
            subprocess.run([PY, "-X", "utf8", "-u", "-m", "pipeline.smd_scheme", "run", "R", "--retry"],
                           cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
        miss = unfinished(jobs())
    if miss:
        lst = "\n".join(f"- {j['reaction']} {j['level']} {j['species']}" for j in miss)
        (ROOT / "docs/results/smd_scheme_check_FAILED.md").write_text(
            "# SMD-scheme check: jobs still unfinished after one retry\n\nNo result is reported from incomplete "
            f"jobs. The water test was started anyway (J. Picione, 2026-11-25).\n\n{lst}\n", encoding="utf-8")
        say(f"STILL UNFINISHED after the retry: {len(miss)} jobs; no SMD-scheme result; see "
            "docs/results/smd_scheme_check_FAILED.md")
    else:
        with open(ROOT / "docs/console_smd_scheme_report.txt", "a", encoding="utf-8") as fh:
            r = subprocess.run([PY, "-X", "utf8", "-u", "-m", "pipeline.smd_scheme", "report"], cwd=ROOT,
                               stdout=fh, stderr=subprocess.STDOUT)
        say(f"SMD-scheme report written (return code {r.returncode}); docs/results/smd_scheme_check.md")
    say("explicit-water test starts")
    with open(WATER, "a", encoding="utf-8") as fh:
        r = subprocess.run([PY, "-X", "utf8", "-u", "-m", "pipeline.explicit_water",
                            "reactions/nucadd_cyanide_acetone_m062x", "run", "--now"], cwd=ROOT,
                           stdout=fh, stderr=subprocess.STDOUT)
    say(f"explicit-water test ended, return code {r.returncode}; nothing else starts")
    sys.exit(0)


if __name__ == "__main__":
    main()
