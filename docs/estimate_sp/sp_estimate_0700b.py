import json, re, statistics, sys, heapq
from pathlib import Path
from datetime import datetime, timedelta
sys.path.insert(0, r"C:\Users\giova\Mechanism Animations\docs\estimate_sp")
import importlib.util
src = open(sys.path[0] + r"\sp_estimate.py", encoding="utf-8").read().split("now = datetime.now()")[0]
g = {}; exec(src, g)
ROOT=g["ROOT"]; RX=g["RX"]; LV=g["LV"]; row=g["row"]; done=g["done"]; cores_of=g["cores_of"]; ts_diag=g["ts_diag"]
from pipeline.steps import Reaction, path_sp_frames, mirror_partners
from pipeline.geometry import read_xyz_frames
# measured DLPNO path / TS ratio (cyanide)
dl=[]
for ln in (ROOT/"reactions/nucadd_cyanide_acetone_m062x/02_irc/path_sp/run_log_path_sp.md").read_text(encoding="utf-8").splitlines():
    m=row.match(ln)
    if m and m.group(2).startswith("DLPNO"): dl.append(int(m.group(3))*cores_of(m.group(4))/4)
tsc=ts_diag("nucadd_cyanide_acetone_m062x", LV[2])
ratio=statistics.median(dl)/tsc
print("cyanide DLPNO path jobs", len(dl), "median 4c s", statistics.median(dl), "TS", tsc, "ratio", round(ratio,3), "min/max", min(dl)/tsc, max(dl)/tsc)
def build(copy):
    lanes={}; info={}
    for rx in RX:
        r=Reaction.load(ROOT/"reactions"/rx); n=len(read_xyz_frames(r.root/"02_irc/irc_full_trj.xyz"))
        its=json.loads((r.root/"02_irc/irc_meta.json").read_text())["ts_frame"]
        mir=mirror_partners(r,n) if copy else {}
        t={}; 
        for lv in LV[:2]:
            x=[]
            for ln in (r.root/"02_irc/path_sp/run_log_path_sp.md").read_text(encoding="utf-8").splitlines():
                m=row.match(ln)
                if m and m.group(2).startswith(lv.split("_")[0]): x.append(int(m.group(3))*cores_of(m.group(4))/4)
            t[lv]=statistics.median(x[-80:]) if x else None
        if t[LV[1]] is None: t[LV[1]]=t[LV[0]]*29.5/34.5
        t[LV[2]]=ts_diag(rx,LV[2])
        jobs=[]; inf={}
        for i,lv in enumerate(LV):
            todo=path_sp_frames(r,n,its,lv)
            if i>0 and mir: todo=[k for k in todo if k not in mir]
            left=[k for k in todo if not done(rx,lv,k)]
            inf[lv]=(len(todo),len(left),round(t[lv]))
            jobs+=[(t[lv], i==2)]*len(left)
        lanes[rx]=jobs; info[rx]=inf
    return lanes, info
def simulate(lanes,start,sc,eff=0.88):
    q={rx:list(j) for rx,j in lanes.items()}; t=0.0; free=2; ev=[]; waiting=[]; last={}
    for rx in q: waiting+=[rx,rx]
    while waiting or ev:
        while free and waiting:
            rx=waiting.pop(0)
            if q[rx]:
                d,isdl=q[rx].pop(0); d=d*(sc if isdl else 1)/eff
                heapq.heappush(ev,(t+d,rx)); free-=1
        if not ev: break
        t,rx=heapq.heappop(ev); free+=1; last[rx]=t
        if q[rx]: waiting.append(rx)
    return {rx:(start+timedelta(seconds=s)).strftime("%a %m-%d %H:%M") for rx,s in last.items()}
now=datetime.now()
for copy in (True,):
    lanes,info=build(copy)
    print("COPY" if copy else "NO COPY"); 
    for rx in RX: print(" ",rx,info[rx])
    for sc in (ratio*0.8, ratio, ratio*1.3):
        print("  scale",round(sc,2), simulate(lanes,now,sc))
for rx in RX:
    p=ROOT/"reactions"/rx/"02_irc/path_sp/run_log_path_sp.md"
    x=[int(m.group(3))*cores_of(m.group(4))/4 for ln in p.read_text(encoding="utf-8").splitlines() for m in [row.match(ln)] if m and m.group(2).startswith("DLPNO")]
    if x: print("DLPNO", rx, len(x), "ratio median", round(statistics.median(x)/ts_diag(rx,LV[2]),3))
