"""Remaining path single points and finish-time estimate for the reports (copied from the session scratchpad 2026-11-23).
Run: .venv/Scripts/python.exe -X utf8 docs/estimate_sp/sp_estimate_now.py
Prints per-lane DLPNO ratio (path frame time / TS single-point time, all frames and last 24 h), jobs left per level,
and finish times for ratio 1.08 x (0.8 / 1.0 / 1.3 = low / mid / high); hshift lanes use their own ratios."""
import sys, statistics
S=r"C:\Users\giova\Mechanism Animations\docs\estimate_sp"
exec(open(S+r"\sp_estimate_0700b.py",encoding="utf-8").read().split("now=datetime.now()")[0])
exec(open(S+r"\sp_estimate_1700b.py",encoding="utf-8").read().split("import heapq")[1].split("now=datetime.now()")[0])
DA,SEC,TER,E2="da_butadiene_ethylene_m062x","hshift_3_methyl_2_butyl_m062x","hshift_2_3_dimethyl_2_butyl_m062x","e2_ethoxide_2_bromopropane_m062x"
for rx in (SEC,TER,E2):
    p=ROOT/"reactions"/rx/"02_irc/path_sp/run_log_path_sp.md"
    rows=[(m.group(1),int(m.group(3))*cores_of(m.group(4))/4) for ln in p.read_text(encoding="utf-8").splitlines() for m in [row.match(ln)] if m and m.group(2).startswith("DLPNO")]
    ts=ts_diag(rx,LV[2]); x=[v for _,v in rows]; y=[v for d,v in rows if d>=(datetime.now()-timedelta(hours=24)).strftime("%Y-%m-%d %H:%M")]
    print(rx, "n",len(x),"ratio all",round(statistics.median(x)/ts,3) if x else None,"last 24 h n",len(y),"ratio",round(statistics.median(y)/ts,3) if y else None, "TS 4c s", round(ts))
now=datetime.now(); lanes,info=build(True)
for rx in RX: print(rx, info[rx])
for name,f in [("low",0.8),("mid",1.0),("high",1.3)]:
    sc={rx:1.08*f for rx in RX}; sc[SEC]=1.27*f; sc[TER]=1.06*f
    print(name, sim(lanes,now,sc))
sc={rx:1.08 for rx in RX}; sc[SEC]=1.27; sc[TER]=1.06; sc[E2]=0.95
print("E2 at 0.95", sim(lanes,now,sc))
