"""Start a long-running process with no visible window (J. Picione, 2026-11-15).

Why: every process started through WMI as a console program (Git Bash, python.exe, powershell.exe) was handed
to Windows Terminal, the default terminal on this machine, and ran as a visible tab. Closing that window or tab,
or anything that ends Windows Terminal, ends the lane with it. This launcher runs under pythonw.exe (a windowed
program, so it gets no console itself) and starts the target with CREATE_NO_WINDOW: the target gets a console
without a window, which is never handed to Windows Terminal, and its child processes (bash → python → ORCA)
share that console. The launcher exits at once; the target's parent is then gone, as with WMI.

Usage (through WMI, from PowerShell):
  "<base pythonw.exe>" "<repo>\\pipeline\\launch_hidden.py" <base64 of the UTF-8 command line> [<working dir>]
The command line goes base64-encoded so no quoting survives nested shells. One line per start is appended to
docs/hidden_launches.log (time, pid, command).
"""
import base64
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CREATE_NO_WINDOW = 0x08000000


def main() -> None:
    a = sys.argv[1]
    if a.startswith("@"):                     # a long command line (the guard, 27,000 characters) from a file
        cmd = Path(a[1:]).read_text(encoding="utf-8-sig").strip()   # PowerShell 5.1 writes a BOM
    else:
        cmd = base64.b64decode(a).decode("utf-8")
    cwd = sys.argv[2] if len(sys.argv) > 2 else str(ROOT)
    p = subprocess.Popen(cmd, cwd=cwd, creationflags=CREATE_NO_WINDOW, close_fds=True,
                         stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(ROOT / "docs" / "hidden_launches.log", "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  pid {p.pid}  {cmd[:400]}\n")


if __name__ == "__main__":
    main()
