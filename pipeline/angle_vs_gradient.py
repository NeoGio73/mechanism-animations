"""Angle check against gradient size, all six production paths (J. Picione, 2026-11-16 evening; for the SI).

For every computed frame of every production side: the mass-weighted angle between the path direction and the
downhill gradient (the angle of the angle check, pipeline.path_angles.side_angles) and the size of the gradient at
that frame (RMS Cartesian gradient, Eh/bohr). Existing data only: each side's path.json (energies, gradients) and
path.xyz; no new calculation. The question: do the large angles occur only where the gradient is nearly zero (near the
TS and in the slow final approach to a minimum)?

Zones: near the TS = |s| <= 1 bohr·amu½ (as the angle check); final approach = the rest of the side from the first frame
whose RMS gradient falls below FLAT and stays below it to the side's end; middle = everything between.
Mirror frames (tertiary shift) have no jobs and are not included. The secondary shift's backward side is still being
extended: its finished first part (true_irc/backward_part1) is used until it ends; rerun this script then.
Writes docs/si_figures/angle_vs_gradient.{png,json,md}.
Usage: .venv/Scripts/python.exe -X utf8 -m pipeline.angle_vs_gradient
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .path_angles import ROOT, _xyz_frames, side_angles

PATHS = [("SN2", "sn2_hydroxide_bromomethane_m062x"), ("Cyanide + acetone", "nucadd_cyanide_acetone_m062x"),
         ("Diels–Alder", "da_butadiene_ethylene_m062x"), ("E2", "e2_ethoxide_2_bromopropane_m062x"),
         ("Tertiary H shift", "hshift_2_3_dimethyl_2_butyl_m062x"), ("Secondary H shift", "hshift_3_methyl_2_butyl_m062x")]
NEAR_S = 1.0          # bohr·amu½, near-TS zone of the angle check
LARGE = 15.0          # degrees: the angle check's 90th-percentile limit
FLAT = 5e-4           # Eh/bohr RMS: "nearly zero" for the final approach (10x the end-minimum test's 5e-5)
OUT = ROOT / "docs/si_figures"


def _side_dir(rx: str, d: str) -> Path:
    base = ROOT / "reactions" / rx / "02_irc/true_irc"
    if not (base / d / "path.json").exists() and (base / f"{d}_part1/path.json").exists():
        return base / f"{d}_part1"
    return base / d


def frames() -> list[dict]:
    rows = []
    for name, rx in PATHS:
        sym, x_ts = _xyz_frames(ROOT / "reactions" / rx / "01_stationary/ts/ts.xyz")[0]
        for d in ("backward", "forward"):
            sd = _side_dir(rx, d)
            js = json.loads((sd / "path.json").read_text(encoding="utf-8"))
            fr = _xyz_frames(sd / "path.xyz")
            keep = [k for k, (_, x) in enumerate(fr) if np.max(np.abs(x - x_ts)) > 1e-6]
            X = [fr[k][1] for k in keep]
            G = [np.asarray(js["gradients"][k]).ravel() for k in keep]
            rms = [float(np.sqrt(np.mean(g ** 2))) for g in G]
            ang = side_angles(X, G, x_ts, sym)          # frames 1 .. n-1 (the last frame has no next step)
            # final approach: from the first frame after which the RMS gradient stays below FLAT to the end
            k0 = len(rms)
            while k0 > 0 and rms[k0 - 1] < FLAT:
                k0 -= 1
            for r in ang:
                i = r["frame_from_ts"] - 1
                zone = "near_ts" if r["s_abs"] <= NEAR_S else ("final_approach" if i >= k0 else "middle")
                rows.append({"path": name, "reaction": rx, "side": d, "folder": sd.relative_to(ROOT).as_posix(),
                             "frame_from_ts": r["frame_from_ts"], "s_abs": r["s_abs"], "angle_deg": r["mw_deg"],
                             "rms_grad_eh_bohr": rms[i], "zone": zone})
    return rows


def summary(rows: list[dict]) -> dict:
    out = {}
    for name, _ in PATHS + [("all", None)]:
        R = [r for r in rows if name == "all" or r["path"] == name]
        big = [r for r in R if r["angle_deg"] > LARGE]
        mid = [r for r in R if r["zone"] == "middle"]
        out[name] = {
            "frames": len(R), "large": len(big),
            "large_by_zone": {z: sum(r["zone"] == z for r in big) for z in ("near_ts", "middle", "final_approach")},
            "large_max_rms_grad": max((r["rms_grad_eh_bohr"] for r in big), default=None),
            "middle_frames": len(mid),
            "middle_median_angle": float(np.median([r["angle_deg"] for r in mid])) if mid else None,
            "middle_max_angle": max((r["angle_deg"] for r in mid), default=None),
            "middle_min_rms_grad": min((r["rms_grad_eh_bohr"] for r in mid), default=None),
            "zones": {z: {"n": sum(r["zone"] == z for r in R),
                          "median_angle": (float(np.median([r["angle_deg"] for r in R if r["zone"] == z]))
                                           if any(r["zone"] == z for r in R) else None),
                          "median_rms_grad": (float(np.median([r["rms_grad_eh_bohr"] for r in R if r["zone"] == z]))
                                              if any(r["zone"] == z for r in R) else None)}
                      for z in ("near_ts", "middle", "final_approach")}}
    return out


def plot(rows: list[dict], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid = "#2b2b29", "#6b6a63", "#e4e3dd"
    col = {"near_ts": "#eb6834", "middle": "#2a78d6", "final_approach": "#1baf7a"}    # categorical slots 2, 1, 3
    mark = {"near_ts": "^", "middle": "o", "final_approach": "s"}                     # shape too, never colour alone
    lab = {"near_ts": "Near the TS (|s| ≤ 1 bohr·amu½)", "middle": "Middle of the side",
           "final_approach": "Final approach (RMS gradient < 5 × 10⁻⁴ Eh/bohr to the end, dotted line)"}
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "font.family": "DejaVu Sans"})
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8), sharex=True, sharey=True, constrained_layout=True)
    g_all = [r["rms_grad_eh_bohr"] for r in rows]
    for ax, (name, _) in zip(axs.ravel(), PATHS):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=grid, linewidth=0.6, which="major")
        ax.axhline(LARGE, color=muted, linewidth=0.8, linestyle=(0, (4, 3)))
        ax.axvline(FLAT, color=muted, linewidth=0.6, linestyle=(0, (1, 2)))
        for z in ("middle", "final_approach", "near_ts"):
            R = [r for r in rows if r["path"] == name and r["zone"] == z]
            if R:
                ax.scatter([r["rms_grad_eh_bohr"] for r in R], [r["angle_deg"] for r in R], s=9, marker=mark[z],
                           c=col[z], edgecolors="#fcfcfb", linewidths=0.4, label=lab[z], zorder=3)
        ax.set_title(name, fontsize=8.5, color=ink, loc="left")
        ax.set_xscale("log")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[0, 0].set_xlim(min(g_all) / 1.5, max(g_all) * 1.5)
    axs[0, 0].set_ylim(0, 120)
    axs[0, 0].set_yticks([0, 15, 30, 60, 90, 120])
    for ax in axs[1]:
        ax.set_xlabel("RMS gradient (Eh/bohr)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Mass-weighted angle (°)")
    h, l = axs.ravel()[4].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=3, frameon=False, fontsize=7.5, markerscale=1.6)
    fig.savefig(path, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = frames()
    s = summary(rows)
    (OUT / "angle_vs_gradient.json").write_text(json.dumps({"summary": s, "frames": rows}, indent=1), encoding="utf-8")
    plot(rows, OUT / "angle_vs_gradient.png")
    return s


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
