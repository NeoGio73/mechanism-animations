"""Remaining machine time for the six reactions at the M06-2X working level (estimate, 2026-11-12).

Machine time = sum over jobs of wall time x cores used / 8 (the 8-core budget fully packed), so one
machine-hour is the whole 8-core budget for one hour. Wall times per job are taken at the core
counts they were logged at (IRC sides 2 cores, single points 4 cores). Every input is tagged
M (measured in this project's logs) or X (extrapolated); see the notes at the end.
Run: .venv/Scripts/python.exe -X utf8 docs/estimate_m062x_remaining.py
"""
H = 3600.0
STEP_AVG = {0.10: 0.104, 0.05: 0.052}   # M: mean step along the whole SN2 step-test paths (162 and 323 steps over 16.84)

# per reaction: path length (bohr amu^1/2), M06-2X/ma-def2-SVP gradient at 4 cores (s), per-frame
# single points at 4 cores (s): DLPNO, M06-2X/TZVP (default grid), B3LYP/TZVP; B3LYP orbital SP (s);
# M06-2X TS Opt+Freq wall (s); orbitals shown; tags
R = {
    "SN2": dict(L=16.84, Lt="M", g=25.8, dl=145, dlt="M scaled", m6=8, b3=9, orb=4 * 2, ts=257, n_orb=4,
                note="path length and steps measured on the B3LYP step test; DLPNO per frame 116 s (pilot path) x 154/123 (ma-def2-SVP TS vs pilot TS)"),
    "E2": dict(L=26.97, Lt="M", g=51.3, dl=1262, dlt="M", m6=59, b3=60, orb=35, ts=1523, n_orb=5,
               note="IRC step measured at M06-2X (233 s at 2 cores); B3LYP path length (stopping rule met earlier, at step 98 of 116, so this is an upper bound)"),
    "Diels–Alder": dict(L=11.50, Lt="M", g=42.0, dl=727, dlt="M", m6=42, b3=43, orb=29, ts=686, n_orb=5, note=""),
    "Cyanide + acetone": dict(L=9.30, Lt="M", g=32.5, dl=506, dlt="M", m6=38, b3=33, orb=27, ts=521, n_orb=4, note=""),
    "2,3-dimethyl-2-butyl H shift": dict(L=2.0, Lt="X", g=45.6, dl=640, dlt="M (stationary points)", m6=42, b3=40, orb=35, ts=842,
                                         n_orb=5, note="one IRC open -> TS -> bridged; the second half is its mirror image (no jobs); length = half the B3LYP open-to-open path (3.41) plus the short bridged side"),
    "3-methyl-2-butyl H shift": dict(L=6.8, Lt="X", g=24.3 * 1.18, dl=380, dlt="X", m6=35, b3=35, orb=29, ts=570, n_orb=5,
                                     note="no reaction path yet: length taken as twice the tertiary B3LYP path; DLPNO from the tertiary cation scaled by (16/19)^3; gradient from the gas-phase M06-2X TS x 1.18 (DCM/gas, tertiary scans)"),
}
IRC_PER_GRAD = 233.0 / 51.3   # M: E2 M06-2X IRC step at 2 cores / E2 M06-2X gradient at 4 cores
DEFGRID3 = 1.5                # X: cost of M06-2X/TZVP single points on DefGrid3 vs the default grid (not measured; range 1.2-2)
ORB_M06 = 1.5                 # X: M06-2X/DefGrid3 orbital single point vs B3LYP (not measured)
CUBE_S = 0.6                  # M: cyanide, 480 cubes in 5 min
POST_MIN = 1 + 3 * 3 + 2      # M/X: meshes 1 min (M), three videos ~3 min each (M: 181-180 s), card and figures ~2 min (X)


def estimate(name, r, step):
    steps = round(r["L"] / STEP_AVG[step])
    frames = steps + 1
    irc_step = r["g"] * IRC_PER_GRAD                       # wall per IRC step, 2 cores
    irc = steps * irc_step * 2 / 8                        # machine-s: both sides at 2 cores each
    sp = frames * (r["dl"] + r["m6"] * DEFGRID3 + r["b3"]) * 4 / 8
    orb = frames * r["orb"] * ORB_M06 * 4 / 8
    post = frames * r["n_orb"] * CUBE_S + POST_MIN * 60   # whole machine
    stat = (r["ts"] * 3) * 4 / 8                          # TS re-optimization + 2 end points, Opt+Freq
    return {"steps": steps, "IRC": irc / H, "TS+ends": stat / H, "SP": sp / H, "orbitals": orb / H, "post": post / H}


if __name__ == "__main__":
    for step in (0.10, 0.05):
        print(f"\n== step {step:.2f} bohr·amu½")
        tot = 0.0
        print(f"{'reaction':32s} {'steps':>5s} {'TS+ends':>7s} {'IRC':>6s} {'SP x3':>6s} {'orbit.':>6s} {'cards':>6s} {'total h':>8s}")
        for n, r in R.items():
            e = estimate(n, r, step)
            t = sum(v for k, v in e.items() if k != "steps")
            tot += t
            print(f"{n:32s} {e['steps']:5d} {e['TS+ends']:7.1f} {e['IRC']:6.1f} {e['SP']:6.1f} {e['orbitals']:6.1f} {e['post']:6.1f} {t:8.1f}")
        print(f"{'six reactions':32s} {'':5s} {'':7s} {'':6s} {'':6s} {'':6s} {'':6s} {tot:8.1f}  = {tot / 24:.1f} days")
