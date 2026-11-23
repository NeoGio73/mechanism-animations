import sys
src=open(sys.argv[1],encoding="utf-8").read().split("now=datetime.now()")[0]
exec(src)
from datetime import datetime, timedelta
import heapq
def sim(lanes,start,scales,eff=0.88):
    q={rx:list(j) for rx,j in lanes.items()}; t=0.0; free=2; ev=[]; waiting=[]; last={}
    for rx in q: waiting+=[rx,rx]
    while waiting or ev:
        while free and waiting:
            rx=waiting.pop(0)
            if q[rx]:
                d,isdl=q[rx].pop(0); d=d*(scales[rx] if isdl else 1)/eff
                heapq.heappush(ev,(t+d,rx)); free-=1
        if not ev: break
        t,rx=heapq.heappop(ev); free+=1; last[rx]=t
        if q[rx]: waiting.append(rx)
    return {rx[:12]:(start+timedelta(seconds=s)).strftime("%a %m-%d %H:%M") for rx,s in last.items()}
now=datetime.now(); lanes,info=build(True)
# DA: its own frames, rising toward the TS like cyanide (x1.2 on the reactant-end time 580 s) = 700 s / 995 s
da=700/995
for name,sc in [("all cyanide 1.07",{rx:ratio for rx in RX}),("DA own 0.70, rest 1.07",{**{rx:ratio for rx in RX},"da_butadiene_ethylene_m062x":da}),("all 0.70",{rx:da for rx in RX})]:
    print(name, sim(lanes,now,sc))
