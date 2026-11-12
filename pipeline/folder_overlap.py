"""Check 1a (J. Picione, 2026-11-11): every ORCA job that ran at the same time as another ORCA job in
the same folder.

Times from the outputs themselves: end = last-modified time of <job>.out; start = end − ORCA's
"TOTAL RUN TIME" (jobs that ended normally) or end − the wall time in run_log.md (jobs that did not).
Outputs moved into crashed_*/ or stopped_*/ folders ran in the folder above.
Writes docs/level_test_queue/folder_overlap.json and prints a table.
"""
import json
import re
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RX = ROOT / "reactions"
RT = re.compile(r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec")


def run_log_walls():
    """{(reaction, job path as logged): wall s} from every run_log.md."""
    out = {}
    for log in RX.glob("*/run_log.md"):
        for ln in log.read_text(encoding="utf-8").splitlines():
            m = re.match(r"\| (\d{4}-\d\d-\d\d \d\d:\d\d) \| `([^`]+)` \|[^|]*\| (\d+) s \|", ln)
            if m:
                end = time.mktime(time.strptime(m.group(1), "%Y-%m-%d %H:%M"))
                out.setdefault((log.parent.name, m.group(2)), []).append((end, float(m.group(3))))
    return out


def jobs():
    walls = run_log_walls()
    res = []
    for o in RX.rglob("*.out"):
        if o.name.endswith((".smd.out",)) or "path_sp" in o.parts and False:
            continue
        try:
            head = o.read_text(errors="replace")
        except OSError:
            continue
        if "O   R   C   A" not in head[:20000]:
            continue
        ran_in = o.parent
        moved = None
        while ran_in.name.startswith(("crashed_", "stopped_")):
            moved = ran_in.name
            ran_in = ran_in.parent
        end = o.stat().st_mtime
        ok = "ORCA TERMINATED NORMALLY" in head[-5000:]
        m = RT.search(head[-5000:])
        rx = o.relative_to(RX).parts[0]
        if m:
            d, h, mi, s, ms = map(int, m.groups())
            dur = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000
            src = "ORCA total run time"
        else:
            rel = str(ran_in.relative_to(RX / rx)).replace("\\", "/") + "/" + o.stem
            rows = walls.get((rx, rel), [])      # the row logged when this output ended (same minute)
            near = [w for e, w in rows if abs(e - end) < 120]
            dur = near[0] if near else None
            src = "run_log wall time" if dur else None
        res.append({"out": str(o.relative_to(RX)).replace("\\", "/"), "folder": str(ran_in.relative_to(RX)).replace("\\", "/"),
                    "moved_to": moved, "ok": ok, "end": end, "start": end - dur if dur else None, "dur_s": dur,
                    "time_source": src})
    return res


def main():
    J = jobs()
    by = {}
    for j in J:
        by.setdefault(j["folder"], []).append(j)
    rows, unknown = [], []
    for f, js in by.items():
        if len(js) < 2:
            continue
        for a in js:
            if a["start"] is None:
                unknown.append(a)
                continue
            for b in js:
                if b is a or b["start"] is None:
                    continue
                ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
                if ov > 0:
                    rows.append({"job": a["out"], "ok": a["ok"], "other": b["out"], "other_ok": b["ok"],
                                 "overlap_s": round(ov), "from": time.strftime("%m-%d %H:%M:%S", time.localtime(max(a["start"], b["start"]))),
                                 "to": time.strftime("%m-%d %H:%M:%S", time.localtime(min(a["end"], b["end"])))})
    unknown_shared = [u["out"] for u in unknown]
    out = {"n_orca_outputs": len(J), "overlaps": rows, "no_start_time_in_shared_folders": unknown_shared}
    (ROOT / "docs/level_test_queue/folder_overlap.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"{len(J)} ORCA outputs scanned")
    for r in sorted(rows, key=lambda r: r["job"]):
        print(f"{'OK ' if r['ok'] else 'BAD'} {r['job']}\n     with {r['other']} ({'ok' if r['other_ok'] else 'not ok'}) "
              f"{r['from']} – {r['to']} ({r['overlap_s']} s)")
    print("shared folder, start time unknown:", unknown_shared)


if __name__ == "__main__":
    main()
