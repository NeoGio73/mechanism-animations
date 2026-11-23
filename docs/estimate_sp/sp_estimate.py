"""Remaining path single points per reaction and level; finish-time simulation (scratch, for the 07:00 report)."""
import json, re, statistics, sys, heapq
from pathlib import Path
from datetime import datetime, timedelta
ROOT = Path(r"C:\Users\giova\Mechanism Animations")
sys.path.insert(0, str(ROOT))
from pipeline.steps import Reaction, path_sp_frames
from pipeline.geometry import read_xyz_frames

RX = ["nucadd_cyanide_acetone_m062x", "da_butadiene_ethylene_m062x", "e2_ethoxide_2_bromopropane_m062x",
      "hshift_2_3_dimethyl_2_butyl_m062x", "hshift_3_methyl_2_butyl_m062x"]
LV = ["M06-2X_ma-def2-TZVP", "B3LYP_ma-def2-TZVP", "DLPNO-CCSD(T)_ma-def2-TZVP"]
row = re.compile(r"^\| (\S+ \S+) \| `([^`]*)` \| `[^`]*` \| (\d+) s \|.*\| (.*) \|\s*$")

def done(rx, lv, k):
    o = ROOT / "reactions" / rx / f"02_irc/path_sp/{lv}/f{k:03d}/f.out"
    return o.exists() and "ORCA TERMINATED NORMALLY" in o.read_text(errors="ignore")[-3000:]

def cores_of(note):
    m = re.search(r"on (\d+) cores \(asked (\d+)\)", note)
    return int(m.group(1)) if m else 4

def ts_diag(rx, lv):
    t = c = None
    for ln in (ROOT / "reactions" / rx / "run_log.md").read_text(encoding="utf-8").splitlines():
        m = row.match(ln)
        if m and m.group(2) == f"01_ts/diagnostics/{lv}/ts/ts":
            t, c = int(m.group(3)), cores_of(m.group(4))
    return t * c / 4          # 4-core equivalent, ideal scaling

out, lanes = {}, {}
for rx in RX:
    r = Reaction.load(ROOT / "reactions" / rx)
    n = len(read_xyz_frames(r.root / "02_irc/irc_full_trj.xyz"))
    its = json.loads((r.root / "02_irc/irc_meta.json").read_text())["ts_frame"]
    # measured M06-2X path times (4-core equivalent), last 80 jobs
    m06 = []
    for ln in (r.root / "02_irc/path_sp/run_log_path_sp.md").read_text(encoding="utf-8").splitlines():
        m = row.match(ln)
        if m and m.group(2).startswith("M06-2X"):
            m06.append(int(m.group(3)) * cores_of(m.group(4)) / 4)
    t_m06 = statistics.median(m06[-80:])
    t = {LV[0]: t_m06, LV[1]: t_m06 * 29.5 / 34.5,                     # SN2 path ratio B3LYP/M06-2X
         LV[2]: ts_diag(rx, LV[2])}                                     # TS time (SN2 path frames took ~0.5 of it)
    jobs, info = [], {}
    for lv in LV:
        todo = path_sp_frames(r, n, its, lv)
        left = [k for k in todo if not done(rx, lv, k)]
        info[lv] = {"todo": len(todo), "left": len(left), "t4_s": round(t[lv])}
        jobs += [t[lv]] * len(left)
    out[rx] = {"frames": n, "levels": info, "m06_jobs_measured": len(m06)}
    lanes[rx] = jobs

def simulate(lanes, start, dlpno_scale=1.0, eff=0.88, pause=None):
    """2 workers per lane, 4 cores per job, 8 cores; FIFO; job time / eff (measured overhead)."""
    q = {rx: list(j) for rx, j in lanes.items()}
    t = 0.0; free = 2; ev = []; waiting = []; last = {}
    for rx in q:
        waiting += [rx, rx]
    while waiting or ev:
        while free and waiting:
            rx = waiting.pop(0)
            if q[rx]:
                d = q[rx].pop(0)
                d = d * (dlpno_scale if d > 150 else 1.0) / eff
                heapq.heappush(ev, (t + d, rx)); free -= 1
        if not ev: break
        t, rx = heapq.heappop(ev); free += 1; last[rx] = t
        if q[rx]: waiting.append(rx)
    return {rx: start + timedelta(seconds=s) for rx, s in last.items()}

now = datetime.now()
res = {"now": now.isoformat(timespec="minutes"), "per_reaction": out}
for name, sc in (("late (path DLPNO = TS time)", 1.0), ("early (path DLPNO = 0.5 x TS time, as SN2)", 0.5)):
    fin = simulate(lanes, now, sc)
    res[name] = {rx: f.strftime("%a %Y-%m-%d %H:%M") for rx, f in fin.items()}
print(json.dumps(res, indent=1))
