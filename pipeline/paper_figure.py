"""Manuscript figure stills, one composite per reaction (J. Picione, 2026-11-27): docs/paper_figures/<id>.png + .md.

Row 1: four IRC frames rendered in Blender with the video's styling (pipeline/render_figure_blender.py, exact IRC
frames, orthographic camera): first frame, the reactant-side frame at half the barrier height (before the TS), the
TS, last frame; the changing orbitals as in the viewer and video (the same meshes); the default-level energy curve
below with every panel marked. Extra rows per reaction (his list):
  E2            Newman view along Cβ–Cα of the first three row-1 frames;
  Diels–Alder   IRC frames 188 and 189 (the step where the computed localization switches), as displayed now (the
                display-only mirror frame at 188) and as computed (meshes built without the mirror frames, in
                05_render/paper_figure/meshes_computed; the viewer data are not touched);
  tertiary      the canonical LUMO at the row-1 frames;
  secondary     row 1 is three panels: the bridged eclipsed form (path start), the C2–C3 rotation TS, the
                tertiary cation (path end), with the three tracked orbitals.
190 mm wide at 300 dpi or better, PNG only (no PDF in the repository). Smallest text 7 pt.
Run: .venv/Scripts/python.exe -X utf8 -m pipeline.paper_figure <id> [<id> ...]   (or "all")
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from .figure import MIN_PT, ORB_ALPHA, _half_height, _irc_no, _line_label
from .steps import Reaction, _read_json
from .video import BLENDER

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "paper_figures"
PANEL_PX = (1000, 800)
W_MM = 190
DPI = 300

LABELS = {
    "sn2_hydroxide_bromomethane": {0: "C", 1: "Br", 2: "O"},
    "e2_ethoxide_2_bromopropane": {0: "Cβ", 1: "Cα", 5: "H", 3: "Br", 13: "O"},
    "da_butadiene_ethylene": {i: f"C{i + 1}" for i in range(6)},
    "nucadd_cyanide_acetone": {0: "C", 1: "O", 4: "C(N)", 5: "N"},
    "hshift_2_3_dimethyl_2_butyl": {0: "C2", 1: "C3", 2: "H"},
    "hshift_3_methyl_2_butyl": {0: "C1", 1: "C2", 3: "C3", 9: "H"},
}


def _row1(rid, M):
    n = len(M["frames"])
    if rid == "hshift_3_methyl_2_butyl":
        return [("bridged eclipsed form", 0), ("C2–C3 rotation TS", M["ts_frame"]), ("tertiary cation", n - 1)]
    return [("first frame", 0), ("before the TS", _half_height(M)[0]),
            ("transition state", M["ts_frame"]), ("last frame", n - 1)]


def panels(rid, M):
    """Rows of panels; each panel names its render group (mesh folder, orbital kind and keys)."""
    loc = {"meshes": "04_meshes", "kind": "localized", "keys": None}
    row1 = [dict(name=f"p{i}", frame=k, view="side", title=t, group="main", **loc) for i, (t, k) in
            enumerate(_row1(rid, M))]
    rows = [row1]
    idx = {f["irc_frame"]: i for i, f in enumerate(M["frames"])}
    if rid == "e2_ethoxide_2_bromopropane":
        rows.append([dict(name=f"q{i}", frame=row1[i]["frame"], view="newman", title=row1[i]["title"], group="main",
                          **loc) for i in range(3)])
        rows[1][0]["row_title"] = "Newman view along Cβ–Cα (Cβ in front)"
    if rid == "da_butadiene_ethylene":
        r = []
        for tag, grp, md in (("as displayed", "main", "04_meshes"),
                             ("as computed", "computed", "05_render/paper_figure/meshes_computed")):
            for fr in (188, 189):
                r.append(dict(name=f"q{len(r)}", frame=idx[fr], view="side", title=f"IRC frame {fr}, {tag}",
                              group=grp, meshes=md, kind="localized", keys=None))
        r[0]["row_title"] = "IRC frames 188 and 189, where the computed localization switches: as displayed and as computed"
        rows.append(r)
    if rid == "hshift_2_3_dimethyl_2_butyl":
        rows.append([dict(name=f"q{i}", frame=p["frame"], view="side", title=p["title"], group="lumo",
                          meshes="04_meshes", kind="canonical", keys=["lumo"]) for i, p in enumerate(row1)])
        rows[1][0]["row_title"] = "Lowest unoccupied orbital (canonical LUMO, not an IBO), the same frames"
    letters = iter("abcdefghijkl")
    for row in rows:
        for p in row:
            p["letter"] = next(letters)
    return rows


def render(rx, M, rows) -> dict:
    base = rx.d("05_render/paper_figure")
    flat = [p for row in rows for p in row]
    groups = []
    for p in flat:
        if p["group"] not in groups:
            groups.append(p["group"])
    proj, setup = {}, None
    for g in groups:
        ps = [p for p in flat if p["group"] == g]
        odir = base / g
        odir.mkdir(parents=True, exist_ok=True)
        st = {"width": PANEL_PX[0], "height": PANEL_PX[1], "samples": 64, "pad_A": 0.35, "orbital_alpha": ORB_ALPHA,
              "dashes_3d": False, "orbital_kind": ps[0]["kind"],
              "panels": [{k: p[k] for k in ("name", "frame", "view")} for p in ps]}
        if ps[0]["keys"]:
            st["orbital_keys"] = ps[0]["keys"]
        if setup:
            st["setup"] = setup                    # same camera scale as the first group
        (odir / "figure_settings.json").write_text(json.dumps(st, indent=1), encoding="utf-8")
        log = odir / "blender_figure.log"
        with open(log, "w") as fh:
            r = subprocess.run([str(BLENDER), "-b", "--factory-startup", "--python",
                                str(Path(__file__).with_name("render_figure_blender.py")), "--",
                                str(rx.root / ps[0]["meshes"]), str(odir)], stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode != 0 or "FIGURE_RENDER_DONE" not in log.read_text(errors="replace"):
            raise RuntimeError(f"Blender figure render failed; see {log}")
        pj = json.loads((odir / "projection.json").read_text(encoding="utf-8"))
        setup = setup or pj.pop("_setup")
        pj.pop("_setup", None)
        for k, v in pj.items():
            v["png"] = str(odir / f"{k}.png")
        proj.update(pj)
    return proj


def compose(rx, M, rows, proj) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from PIL import Image

    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": MIN_PT, "axes.linewidth": 0.6,
                                "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    rid = rx.cfg["id"]
    labels = LABELS.get(rid, {})
    curve = next(c for c in M["curves"] if c["key"] == M["default_curve"])
    s = np.array([f["s"] for f in M["frames"]])
    E = np.array(curve["E_kJ"], dtype=float)
    ok = np.isfinite(E)
    W_in = W_MM / 25.4
    pw = W_in / 4
    imgs, crops = {}, []
    for row in rows:
        lo, hi = PANEL_PX[1], 0
        for p in row:
            a = np.asarray(Image.open(proj[p["name"]]["png"]).convert("RGB"))
            imgs[p["name"]] = a
            used = np.where((a < 248).any(axis=2).any(axis=1))[0]
            lo, hi = min(lo, used.min()), max(hi, used.max())
        crops.append((max(0, lo - 24), min(PANEL_PX[1], hi + 24)))
    phs = [pw * (b - a) / PANEL_PX[0] for a, b in crops]
    used = []                                 # every orbital drawn in any panel, for the legend
    for row in rows:
        for p in row:
            for o in M["orbitals"]:
                if o["kind"] == p["kind"] and (not p["keys"] or o["key"] in p["keys"]) and o not in used:
                    used.append(o)
    title_h, head_h, curve_h, gap = 0.30, 0.20, 2.05, 0.10
    legend_h = 0.26 + 0.16 * len(used)
    n_head = sum(1 for row in rows if row[0].get("row_title"))
    H_in = sum(phs) + len(rows) * title_h + n_head * head_h + legend_h + curve_h + gap
    fig = plt.figure(figsize=(W_in, H_in), dpi=DPI)
    y = H_in
    for r, row in enumerate(rows):
        ph = phs[r]
        top, bot = crops[r]
        if row[0].get("row_title"):
            y -= head_h
            fig.text(0.03 / W_in, (y + head_h - 0.04) / H_in, row[0]["row_title"], ha="left", va="top",
                     fontsize=MIN_PT, style="italic", fontweight="bold", color="#3b4654")
        y -= title_h
        start = (4 - len(row)) / 2 if len(row) < 4 else 0
        for c, p in enumerate(row):
            x0 = (start + c) * pw
            ax = fig.add_axes([x0 / W_in, (y - ph) / H_in, pw / W_in, ph / H_in])
            ax.imshow(imgs[p["name"]][top:bot], extent=[0, PANEL_PX[0], bot, top], interpolation="lanczos")
            ax.set_xlim(0, PANEL_PX[0])
            ax.set_ylim(bot, top)
            ax.axis("off")
            pr = proj[p["name"]]
            for i, j in pr["partial_bonds"]:
                a_, b_ = np.array(pr["atoms_px"][i][:2]), np.array(pr["atoms_px"][j][:2])
                u = (b_ - a_) / np.linalg.norm(b_ - a_)
                a_ = a_ + u * pr["atom_radius_A"][i] * pr["px_per_A"]
                b_ = b_ - u * pr["atom_radius_A"][j] * pr["px_per_A"]
                ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color="#2b2f36", lw=1.1, dashes=(2.4, 1.6),
                        solid_capstyle="butt", zorder=3,
                        path_effects=[pe.Stroke(linewidth=2.3, foreground="white"), pe.Normal()])
            for i, lab in labels.items():
                if pr["view"] == "newman" and M.get("view_axis") and i == M["view_axis"][1]:
                    continue
                px, py, _ = pr["atoms_px"][i]
                ax.annotate(lab, (px, py), xytext=(3.2, 3.2), textcoords="offset points", fontsize=MIN_PT,
                            fontweight="bold", color="#1d2733", ha="left", va="bottom",
                            path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
            k = p["frame"]
            fig.text((x0 + 0.03) / W_in, (y + title_h - 0.05) / H_in, f"({p['letter']}) {p['title']}",
                     ha="left", va="top", fontsize=MIN_PT, fontweight="bold")
            fig.text((x0 + 0.03) / W_in, (y + title_h - 0.17) / H_in,
                     f"IRC frame {_irc_no(M, k)}, s = {s[k]:+.2f}".replace("-", "−"),
                     ha="left", va="top", fontsize=MIN_PT, color="#3b4654")
        y -= ph

    lax =fig.add_axes([0.0, (y - legend_h) / H_in, 1.0, legend_h / H_in])
    y -= legend_h
    lax.set_xlim(0, W_in)
    lax.set_ylim(0, legend_h)
    lax.axis("off")
    sw, line = 0.09, 0.16                     # one orbital per line (the labels are long)
    for n, o in enumerate(used):
        pc, nc = o.get("colors", ["#2f6db5", "#e08a1e"])
        lab = o["label"] if o["kind"] == "localized" else o["label"].replace("canonical LUMO", "LUMO (canonical)")
        yy = legend_h - 0.06 - (n + 1) * line
        lax.add_patch(Rectangle((0.08, yy), sw, sw, color=pc, lw=0))
        lax.add_patch(Rectangle((0.08 + sw, yy), sw, sw, color=nc, lw=0))
        lax.text(0.08 + 2 * sw + 0.05, yy + sw / 2, lab, va="center", fontsize=MIN_PT)
    what = ("Intrinsic bond orbitals that change along the path" if all(o["kind"] == "localized" for o in used) else
            "Intrinsic bond orbitals that change along the path, and the canonical LUMO")
    lax.text(0.08, 0.07, f"{what} (dark = + phase, light = − phase); dashed: bonds forming or breaking.",
             fontsize=MIN_PT, color="#3b4654", va="center")

    y -= gap
    ax = fig.add_axes([0.085, (y - curve_h + 0.38) / H_in, 0.83, (curve_h - 0.55) / H_in])
    ax.plot(s[ok], E[ok], color="#1f4e79", lw=1.2)
    ax.axhline(0, color="#9aa4b1", lw=0.6, ls="--")
    ep = M.get("E_separated_products_kJ")
    ep = curve.get("E_separated_products_kJ", ep)
    if rid == "hshift_2_3_dimethyl_2_butyl":   # one ion, degenerate shift: the zero is the cation (reactant = product)
        ep = None
    if ep is not None:
        ax.axhline(ep, color="#9aa4b1", lw=0.6, ls=":")
    marks = {}
    for row in rows:
        for p in row:
            marks.setdefault(p["frame"], []).append(p["letter"])
    # frames closer than 2 % of the s range share one label (Diels–Alder 188/189: "e–h")
    ks = sorted(marks, key=lambda k: s[k])
    groups = [[ks[0]]]
    for k in ks[1:]:
        if s[k] - s[groups[-1][-1]] < 0.02 * (s[-1] - s[0]):
            groups[-1].append(k)
        else:
            groups.append([k])
    merged = {}
    for g in groups:
        ls = sorted(x for k in g for x in marks[k])
        merged[g[0]] = [f"{ls[0]}–{ls[-1]}"] if len(ls) > 2 else ls
        for k in g[1:]:
            ax.plot([s[k]], [E[k]], "o", ms=3.4, color="#c0392b", mec="white", mew=0.5, zorder=5)
    marks = merged
    span = max(np.nanmax(E), 0) - min(np.nanmin(E), ep if ep is not None else np.nanmin(E))
    for k, ls in marks.items():
        ax.plot([s[k]], [E[k]], "o", ms=3.4, color="#c0392b", mec="white", mew=0.5, zorder=5)
        below = E[k] > np.nanmax(E) - 0.12 * span and k != M["ts_frame"]
        dx, ha = (4, "left") if k == 0 else ((-4, "right") if k == len(s) - 1 else (0, "center"))
        ax.annotate(",".join(ls), (s[k], E[k]), xytext=(dx, -9 if below else 4.5), textcoords="offset points",
                    ha=ha, va="top" if below else "bottom", fontsize=MIN_PT, fontweight="bold", color="#c0392b")
    kmax = int(np.nanargmax(E))
    if kmax not in marks:
        ax.plot([s[kmax]], [E[kmax]], "o", ms=3.6, mfc="white", mec="#1f4e79", mew=0.9, zorder=6)
        ax.annotate(f"maximum, IRC frame {_irc_no(M, kmax)}", (s[kmax], E[kmax]), xytext=(6, 2),
                    textcoords="offset points", ha="left", va="bottom", fontsize=MIN_PT, color="#1f4e79")
    ax.set_xlabel("reaction coordinate s (bohr·amu$^{1/2}$)")
    ax.set_ylabel("ΔE (kJ/mol)")
    sec = ax.secondary_yaxis("right", functions=(lambda v: v / 4.184, lambda v: v * 4.184))
    sec.set_ylabel("ΔE (kcal/mol)")
    ax.set_xlim(s[0], s[-1])
    lo, hi = min(np.nanmin(E), ep if ep is not None else np.nanmin(E)), max(np.nanmax(E), 0)
    ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.16 * (hi - lo))
    ax.spines["top"].set_visible(False)
    zl = M.get("zero_label", "separated reactants")
    if rid == "hshift_2_3_dimethyl_2_butyl":
        zl = "2,3-dimethyl-2-butyl cation (reactant = product)"
    lvl = curve["label"].replace(" single points on the same path", "")
    ax.set_title(f"Electronic energy ΔE relative to {zl}: {lvl}",
                 fontsize=MIN_PT, loc="left", color="#3b4654")
    _line_label(fig, ax, s[ok], E[ok], 0.0, zl, prefer="right")
    if ep is not None:
        pl = "alkoxide (the product, one ion)" if rid == "nucadd_cyanide_acetone" else M.get("product_label", "separated products")
        _line_label(fig, ax, s[ok], E[ok], ep, pl, prefer="left")
    fig.canvas.draw()
    small = [(t.get_text(), t.get_fontsize()) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible() and t.get_fontsize() < MIN_PT]
    if small:
        raise RuntimeError(f"text below {MIN_PT} pt: {small}")
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{rid}.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    return out


def caption(rx, M, rows, out: Path) -> None:
    from PIL import Image
    curve = next(c for c in M["curves"] if c["key"] == M["default_curve"])
    E, s = curve["E_kJ"], [f["s"] for f in M["frames"]]
    rid = rx.cfg["id"]
    fr = "; ".join(f"({p['letter']}) {p['title']}" + ("" if p["title"].startswith("IRC frame") else f", IRC frame {_irc_no(M, p['frame'])}") + f", s = {s[p['frame']]:+.2f} "
                   f"bohr·amu½, ΔE = {E[p['frame']]:+.1f} kJ/mol ({E[p['frame']] / 4.184:+.1f} kcal/mol)"
                   for row in rows for p in row)
    orbs = [o for o in M["orbitals"] if o["kind"] == "localized"]
    w, h = Image.open(out).size
    L = [f"**{M['title']}: {M['equation']}.** Minimum-energy path (IRC) at {M['level']}. Panels: {fr}. Shown: the "
         f"intrinsic bond orbitals (IBOs) that change along the path, tracked frame to frame ({'; '.join(o['label'] for o in orbs)}), "
         f"isovalue ±{M['isovalue']} {M['isovalue_units']}, dark = + phase, light = − phase, drawn at "
         f"{ORB_ALPHA:.0%} opacity (video and viewer 60 %); the same meshes as the viewer and video. Dashed, drawn over "
         "the orbitals: bonds forming or breaking. Bottom: electronic energy ΔE relative to "
         f"{'the cation (reactant = product)' if rid == 'hshift_2_3_dimethyl_2_butyl' else M.get('zero_label', 'separated reactants')}, "
         f"{curve['label']} (the default curve), against the "
         f"mass-weighted reaction coordinate s; red points mark the panels. TS imaginary frequency "
         f"{abs(M['imaginary_frequency_cm']):.0f}i cm⁻¹."]
    if rid != "hshift_3_methyl_2_butyl":
        kh, half, emax, kmax = _half_height(M)
        L.append(f"Panel (b), before the TS: the reactant-side frame whose default-level energy is closest to half the "
                 f"barrier height above the first frame, {half:+.1f} kJ/mol ({half / 4.184:+.1f} kcal/mol), with the "
                 f"highest point of the curve at {emax:+.1f} kJ/mol ({emax / 4.184:+.1f} kcal/mol), IRC frame "
                 f"{_irc_no(M, kmax)}.")
    if rid == "e2_ethoxide_2_bromopropane":
        L.append("Panels (e)–(g): frames (a)–(c) viewed along the Cβ–Cα bond (Newman view, Cβ in front, Cα directly "
                 "behind), framing atoms within 4 Å of the axis.")
    if rid == "da_butadiene_ethylene":
        L.append("Panels (e)–(h): IRC frames 188 and 189, the step where the computed localization switches. (e), (f) as "
                 "displayed in the viewer and video: at frame 188 the display shows the mirror image of the computed "
                 "active orbitals, which continues the previous frame (display-only; the reflected set spans the computed "
                 "orbitals, smallest singular value of the 3×17 overlap ≥ 0.99982 at frames 179–199). (g), (h) as "
                 "computed, without that display substitution.")
    if rid == "hshift_2_3_dimethyl_2_butyl":
        L.append("Panels (e)–(h): the canonical lowest unoccupied orbital (LUMO), not an IBO, at the same frames as "
                 "(a)–(d).")
    if rid == "hshift_3_methyl_2_butyl":
        L.append("Three panels (his choice): the bridged eclipsed form at the path start, the C2–C3 rotation TS, and the "
                 "tertiary cation at the path end; all three tracked orbitals shown. Atom names in the chemical "
                 "numbering CH₃(C1)–C2H⁺–C3H(CH₃)–C4H₃.")
    L.append(f"File: {out.name}, {w} × {h} px = {W_MM} mm wide at {DPI} dpi. Rendered with Blender (EEVEE), "
             "orthographic camera; composite with matplotlib (pipeline/paper_figure.py).")
    import re
    L = [re.sub(r"(?<![\w)\]])-(?=\d)", "−", x) for x in L]          # true minus before digits only
    out.with_suffix(".md").write_text("\n\n".join(L) + "\n", encoding="utf-8")


def figure(rid: str) -> Path:
    rx = Reaction.load(ROOT / "reactions" / f"{rid}_m062x")
    if rid == "da_butadiene_ethylene":
        from .meshes import build_meshes
        if not (rx.root / "05_render/paper_figure/meshes_computed/manifest.json").exists():
            build_meshes(rx, out_sub="05_render/paper_figure/meshes_computed", display_mirror=False)
    M = _read_json(rx.root / "04_meshes" / "manifest.json")
    rows = panels(rid, M)
    proj = render(rx, M, rows)
    out = compose(rx, M, rows, proj)
    caption(rx, M, rows, out)
    print(f"figure: {out}", flush=True)
    return out


def viewer_shot(rid: str = "sn2_hydroxide_bromomethane", port: int = 8766) -> Path:
    """Screenshot of the viewer (production data, ?r=<folder>) at the TS frame, with the energy diagram and the orbital
    switches as the page shows them: headless Chrome, software WebGL, 2244 px wide (190 mm at 300 dpi). A local
    http.server on 127.0.0.1 serves the repository only for the time of the screenshot."""
    import time
    M = _read_json(ROOT / "reactions" / f"{rid}_m062x" / "04_meshes" / "manifest.json")
    out = OUT / f"{rid}_viewer.png"
    srv = subprocess.Popen([str(ROOT / ".venv/Scripts/python.exe"), "-m", "http.server", str(port), "--bind",
                            "127.0.0.1"], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        time.sleep(1.5)
        chrome = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
        url = f"http://127.0.0.1:{port}/web/viewer/index.html?r={rid}_m062x&t={M['ts_frame']}"
        subprocess.run([chrome, "--headless=new", "--hide-scrollbars", "--window-size=1920,1120",
                        f"--force-device-scale-factor={2244 / 1920:.4f}", "--enable-unsafe-swiftshader",
                        "--use-angle=swiftshader", "--virtual-time-budget=20000", f"--screenshot={out}", url],
                       check=True, timeout=180, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    finally:
        srv.terminate()
    out.with_suffix(".md").write_text(
        f"Screenshot of the interactive viewer (web/viewer/index.html?r={rid}_m062x), {M['title']}, opened at the TS "
        f"(IRC frame {M['frames'][M['ts_frame']]['irc_frame']}), default orbitals on, energy diagram and orbital "
        "switches as the page shows them. Headless Chrome with software WebGL, 1920 × 1120 CSS px at "
        f"{2244 / 1920:.3f}× = 2244 px wide (190 mm at 300 dpi).\n", encoding="utf-8")
    print(f"viewer screenshot: {out}", flush=True)
    return out


RX = ["sn2_hydroxide_bromomethane", "e2_ethoxide_2_bromopropane", "nucadd_cyanide_acetone",
      "hshift_2_3_dimethyl_2_butyl", "hshift_3_methyl_2_butyl", "da_butadiene_ethylene"]

if __name__ == "__main__":
    for r in (RX if sys.argv[1:] == ["all"] else sys.argv[1:]):
        figure(r)
