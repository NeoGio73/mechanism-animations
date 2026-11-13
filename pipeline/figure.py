"""Publication composite for one reaction (J. Picione, 2026-11-11): figures/<id>_composite.{png,pdf,md}.

Row 1: four IRC frames — first frame; the reactant-side frame whose energy at the default level
is closest to half the barrier height above the first frame (J. Picione, 2026-11-11); the TS; last
frame — rendered in Blender (pipeline/render_figure_blender.py) with the
tracked orbitals that change (the localized orbitals in the manifest), dashed forming/breaking
bonds and labels on the key atoms. Optional second row (reaction.toml-independent, set in
EXTRA below): E2, the first three frames in the Newman view along Cβ–Cα; Diels–Alder, IRC frames
59 and 60 on either side of the abrupt localization change. Below: the energy curve at the
default level (manifest `default_curve`) with every panel marked.

180 mm wide, 300 dpi PNG; PDF with the panels embedded as images and the curve, labels and
legend as vectors. Caption file states the levels of theory and the isovalue.
Run: .venv/Scripts/python.exe -X utf8 -m pipeline reactions/<id> figure
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np

from .steps import Reaction, _read_json
from .video import BLENDER

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
PANEL_PX = (1000, 800)
ORB_ALPHA = 0.8      # orbital opacity (video and viewer 0.6): phases must stay apart in greyscale
MIN_PT = 7           # smallest text in the figure, at 180 mm width

# atom index → label, per reaction (names as on the cards)
LABELS = {
    "sn2_hydroxide_bromomethane": {0: "C", 1: "Br", 2: "O"},
    "e2_ethoxide_2_bromopropane": {0: "Cβ", 1: "Cα", 5: "H", 3: "Br", 13: "O"},
    "da_butadiene_ethylene": {i: f"C{i + 1}" for i in range(6)},
    "nucadd_cyanide_acetone": {0: "C", 1: "O", 4: "C(N)", 5: "N"},
}
# second row: (kind, frames or None for the row-1 frames, title)
EXTRA = {
    # (view, row-1 panel indices or IRC frame numbers, panel titles or None, row heading)
    "e2_ethoxide_2_bromopropane": ("newman", {"row1": [0, 1, 2]}, None, "Newman view along Cβ–Cα (Cβ in front)"),
    "da_butadiene_ethylene": ("side", {"irc": [59, 60]}, ["before the orbital change", "after the orbital change"],
                              "IRC frames 59 and 60, either side of the abrupt orbital change"),
}


def _half_height(M):
    """Reactant-side frame (first frame … TS) whose default-level energy is closest to
    E(first) + ½ (E(max) − E(first)); E(max) is the highest point of the default-level curve (for
    SN2 one frame after the working-level TS). Refuses if the reactant side crosses that height
    more than once (no clear half-height point)."""
    E = np.array(next(c for c in M["curves"] if c["key"] == M["default_curve"])["E_kJ"], dtype=float)
    ts = M["ts_frame"]
    # computed frames only (a sparse level has no value at the other frames; no interpolation)
    idx = [k for k in range(ts + 1) if np.isfinite(E[k])]
    half = E[0] + (np.nanmax(E) - E[0]) / 2
    d = E[idx] - half
    crossings = [idx[j] for j in range(len(idx) - 1) if np.sign(d[j]) != np.sign(d[j + 1])]
    if len(crossings) != 1:
        raise RuntimeError(f"no clear half-height point on the reactant side (crossings {crossings})")
    return idx[int(np.argmin(np.abs(d)))], half, float(np.nanmax(E)), int(np.nanargmax(E))


def _frames(M):
    n = len(M["frames"])
    return [("first frame", 0), ("half the barrier height", _half_height(M)[0]),
            ("transition state", M["ts_frame"]), ("last frame", n - 1)]


def _irc_no(M, k):
    return M["frames"][k].get("irc_frame", k)


def _panels(rx, M):
    row1 = [{"name": f"p{i}", "frame": k, "view": "side", "title": t} for i, (t, k) in enumerate(_frames(M))]
    rows = [row1]
    ex = EXTRA.get(rx.cfg["id"])
    if ex:
        view, sel, titles, heading = ex
        if "row1" in sel:
            rows.append([{"name": f"q{i}", "frame": row1[j]["frame"], "view": view, "title": row1[j]["title"]}
                         for i, j in enumerate(sel["row1"])])
        else:
            idx = [[f["irc_frame"] for f in M["frames"]].index(fr) for fr in sel["irc"]]
            rows.append([{"name": f"q{i}", "frame": k, "view": view, "title": titles[i]}
                         for i, k in enumerate(idx)])
        rows[1][0]["row_title"] = heading
    letters = iter("abcdefghij")
    for row in rows:
        for p in row:
            p["letter"] = next(letters)
    return rows


def render_panels(rx, M, rows) -> dict:
    odir = rx.d("05_render/figure")
    flat = [p for row in rows for p in row]
    (odir / "figure_settings.json").write_text(json.dumps(
        {"width": PANEL_PX[0], "height": PANEL_PX[1], "samples": 64, "pad_A": 0.35, "orbital_alpha": ORB_ALPHA, "dashes_3d": False,
         "orbital_kind": "localized",
         "panels": [{k: p[k] for k in ("name", "frame", "view")} for p in flat]}, indent=1), encoding="utf-8")
    script = Path(__file__).with_name("render_figure_blender.py")
    log = odir / "blender_figure.log"
    with open(log, "w") as fh:
        r = subprocess.run([str(BLENDER), "-b", "--factory-startup", "--python", str(script), "--",
                            str(rx.root / "04_meshes"), str(odir)], stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0 or "FIGURE_RENDER_DONE" not in log.read_text(errors="replace"):
        raise RuntimeError(f"Blender figure render failed; see {log}")
    return json.loads((odir / "projection.json").read_text(encoding="utf-8"))


def compose(rx, M, rows, proj) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from PIL import Image

    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": MIN_PT, "pdf.fonttype": 42,
                                "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    odir = rx.root / "05_render" / "figure"
    labels = LABELS.get(rx.cfg["id"], {})
    curve = next(c for c in M["curves"] if c["key"] == M["default_curve"])
    s = np.array([f["s"] for f in M["frames"]])
    E = np.array(curve["E_kJ"], dtype=float)
    ok = np.isfinite(E)            # computed frames (all, except on a sparse level)

    W_in = 180 / 25.4
    pw = W_in / 4
    # crop each row to the rows of pixels its panels use (one crop per row keeps one scale)
    imgs, crops = {}, []
    for row in rows:
        lo, hi = PANEL_PX[1], 0
        for p in row:
            a = np.asarray(Image.open(odir / f"{p['name']}.png").convert("RGB"))
            imgs[p["name"]] = a
            used = np.where((a < 248).any(axis=2).any(axis=1))[0]
            lo, hi = min(lo, used.min()), max(hi, used.max())
        crops.append((max(0, lo - 24), min(PANEL_PX[1], hi + 24)))
    phs = [pw * (b - a) / PANEL_PX[0] for a, b in crops]
    title_h, head_h, legend_h, curve_h, gap = 0.30, 0.20, 0.42, 2.05, 0.10
    n_head = sum(1 for row in rows if row[0].get("row_title"))
    H_in = sum(phs) + len(rows) * title_h + n_head * head_h + legend_h + curve_h + gap
    fig = plt.figure(figsize=(W_in, H_in), dpi=300)
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
            # forming/breaking bonds as 2D dashed lines over the lobes (atom centre to centre, trimmed
            # by the drawn atom radius), so they stay visible through the orbitals
            for i, j in pr["partial_bonds"]:
                a_, b_ = np.array(pr["atoms_px"][i][:2]), np.array(pr["atoms_px"][j][:2])
                u = (b_ - a_) / np.linalg.norm(b_ - a_)
                a_ = a_ + u * pr["atom_radius_A"][i] * pr["px_per_A"]
                b_ = b_ - u * pr["atom_radius_A"][j] * pr["px_per_A"]
                ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color="#2b2f36", lw=1.1, dashes=(2.4, 1.6),
                        solid_capstyle="butt", zorder=3,
                        path_effects=[pe.Stroke(linewidth=2.3, foreground="white"), pe.Normal()])
            for i, lab in labels.items():
                if pr["view"] == "newman" and i == M["view_axis"][1]:
                    continue                 # directly behind the front atom
                px, py, _ = pr["atoms_px"][i]
                ax.annotate(lab, (px, py), xytext=(3.2, 3.2), textcoords="offset points", fontsize=MIN_PT,
                            fontweight="bold", color="#1d2733", ha="left", va="bottom",
                            path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
            k = p["frame"]
            fig.text((x0 + 0.03) / W_in, (y + title_h - 0.05) / H_in,
                     f"({p['letter']}) {p['title']}", ha="left", va="top", fontsize=MIN_PT, fontweight="bold")
            fig.text((x0 + 0.03) / W_in, (y + title_h - 0.17) / H_in,
                     f"IRC frame {_irc_no(M, k)}, s = {s[k]:+.2f}".replace("-", "−"),
                     ha="left", va="top", fontsize=MIN_PT, color="#3b4654")
        y -= ph

    # legend: orbitals, phases, dashed bonds
    y -= legend_h
    lax = fig.add_axes([0.0, y / H_in, 1.0, legend_h / H_in])
    lax.set_xlim(0, W_in)
    lax.set_ylim(0, legend_h)
    lax.axis("off")
    orbs = [o for o in M["orbitals"] if o["kind"] == "localized"]
    x = 0.08
    sw = 0.09
    for o in orbs:
        pc, nc = o.get("colors", ["#2f6db5", "#e08a1e"])
        lax.add_patch(Rectangle((x, 0.24), sw, sw, color=pc, lw=0))
        lax.add_patch(Rectangle((x + sw, 0.24), sw, sw, color=nc, lw=0))
        t = lax.text(x + 2 * sw + 0.05, 0.285, o["label"], va="center", fontsize=MIN_PT)
        fig.canvas.draw()
        bb = t.get_window_extent().transformed(lax.transData.inverted())
        x = bb.x1 + 0.22
    lax.text(0.08, 0.07, "Intrinsic bond orbitals that change along the path (dark = + phase, light = − phase); "
             "dashed: bonds forming or breaking.", fontsize=MIN_PT, color="#3b4654", va="center")

    # energy curve
    y -= gap
    ax = fig.add_axes([0.085, (y - curve_h + 0.38) / H_in, 0.83, (curve_h - 0.55) / H_in])
    ax.plot(s[ok], E[ok], color="#1f4e79", lw=1.2)
    ax.axhline(0, color="#9aa4b1", lw=0.6, ls="--")
    ep = curve.get("E_separated_products_kJ")
    if ep is not None:
        ax.axhline(ep, color="#9aa4b1", lw=0.6, ls=":")
    marks = {}
    for row in rows:
        for p in row:
            marks.setdefault(p["frame"], []).append(p["letter"])
    span = max(np.nanmax(E), 0) - min(np.nanmin(E), ep if ep is not None else np.nanmin(E))
    for k, ls in marks.items():
        ax.plot([s[k]], [E[k]], "o", ms=3.4, color="#c0392b", mec="white", mew=0.5, zorder=5)
        below = E[k] > np.nanmax(E) - 0.12 * span and k != M["ts_frame"]
        dx, ha = (4, "left") if k == 0 else ((-4, "right") if k == len(s) - 1 else (0, "center"))
        ax.annotate(",".join(ls), (s[k], E[k]), xytext=(dx, -9 if below else 4.5), textcoords="offset points",
                    ha=ha, va="top" if below else "bottom", fontsize=MIN_PT, fontweight="bold", color="#c0392b")
    kmax = int(np.nanargmax(E))
    if kmax not in marks:        # default-level maximum off the working-level TS frame (IRCMax point)
        ax.plot([s[kmax]], [E[kmax]], "o", ms=3.6, mfc="white", mec="#1f4e79", mew=0.9, zorder=6)
        ax.annotate(f"maximum, IRC frame {_irc_no(M, kmax)}", (s[kmax], E[kmax]), xytext=(6, 2),
                    textcoords="offset points", ha="left", va="bottom", fontsize=MIN_PT, color="#1f4e79")
    ax.set_xlabel("reaction coordinate s (bohr·amu$^{1/2}$)")
    ax.set_ylabel("E (kJ/mol)")
    sec = ax.secondary_yaxis("right", functions=(lambda v: v / 4.184, lambda v: v * 4.184))
    sec.set_ylabel("E (kcal/mol)")
    ax.set_xlim(s[0], s[-1])
    lo, hi = min(np.nanmin(E), ep if ep is not None else np.nanmin(E)), max(np.nanmax(E), 0)
    ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.16 * (hi - lo))
    ax.spines["top"].set_visible(False)
    lvl = curve["label"].replace(" single points on the same path", "")
    ax.set_title(f"Electronic energy relative to separated reactants: {lvl} single points on the "
                 f"{M['level'].split(',')[0]} path", fontsize=MIN_PT, loc="left", color="#3b4654")
    _line_label(fig, ax, s[ok], E[ok], 0.0, "separated reactants", prefer="right")
    if ep is not None:
        _line_label(fig, ax, s[ok], E[ok], ep, "separated products", prefer="left")

    fig.canvas.draw()
    small = [(x.get_text(), x.get_fontsize()) for x in fig.findobj(matplotlib.text.Text)
             if x.get_text().strip() and x.get_visible() and x.get_fontsize() < MIN_PT]
    if small:
        raise RuntimeError(f"text below {MIN_PT} pt: {small}")
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"{rx.cfg['id']}_composite"
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    fig.savefig(f"{out}.pdf", dpi=300, facecolor="white")
    plt.close(fig)
    # black-and-white check copy (luminance only), next to the panel renders
    Image.open(f"{out}.png").convert("L").save(odir / "composite_greyscale_check.png")
    return out


def _line_label(fig, ax, s, E, y, text, prefer="right"):
    """Label a horizontal reference line where its text box does not touch the curve (densely
    interpolated) or a marker label: first at the preferred end, then stepping along the line,
    above it before below it."""
    xs = np.linspace(s[0], s[-1], 2000)
    pts = ax.transData.transform(np.c_[xs, np.interp(xs, s, E)])
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    others = [x.get_window_extent(rend) for x in ax.texts]
    span = s[-1] - s[0]
    fr = np.linspace(0, 1, 41)
    xs_try = [(s[-1] - f * span, "right") for f in fr] if prefer == "right" else [(s[0] + f * span, "left") for f in fr]
    for va in ("bottom", "top"):
        for x, ha in xs_try:
            t = ax.text(x, y, f" {text} ", ha=ha, va=va, fontsize=MIN_PT, color="#5b6878")
            bb = t.get_window_extent(rend).expanded(1.02, 1.15)
            inside = bb.x0 >= ax.bbox.x0 - 1 and bb.x1 <= ax.bbox.x1 + 1
            hit = ((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) & (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)).any()
            if inside and not hit and not any(bb.overlaps(o) for o in others):
                return t
            t.remove()
    return ax.text(s[-1] if prefer == "right" else s[0], y, f" {text} ", ha=prefer, va="bottom",
                   fontsize=MIN_PT, color="#5b6878")


def caption(rx, M, rows, out: Path) -> None:
    curve = next(c for c in M["curves"] if c["key"] == M["default_curve"])
    E = curve["E_kJ"]
    s = [f["s"] for f in M["frames"]]
    orbs = [o for o in M["orbitals"] if o["kind"] == "localized"]
    fr = "; ".join(f"({p['letter']}) {p['title']}, IRC frame {_irc_no(M, p['frame'])}, s = {s[p['frame']]:+.2f} "
                   f"bohr·amu½, E = {E[p['frame']]:.1f} kJ/mol ({E[p['frame']] / 4.184:.1f} kcal/mol)"
                   for row in rows for p in row).replace("-", "−")
    L = [f"**{M['title']}: {M['equation']}.** Minimum-energy path (IRC) at {M['level']}. "
         f"Top (E at the default level, relative to separated reactants): {fr}. Shown: the intrinsic bond orbitals (IBOs) that change along the path, tracked frame "
         f"to frame — {', '.join(o['label'] for o in orbs)} — at isovalue ±{M['isovalue']} "
         f"{M['isovalue_units']}, dark = + phase, light = − phase (drawn at {ORB_ALPHA:.0%} opacity, video and "
         f"viewer 60 %); orbitals and geometries at {M['level']}. "
         "Dashed, drawn over the orbitals: bonds forming or breaking (drawn dashed while longer than 1.15 and shorter than 1.6 times "
         "their own equilibrium length). "
         f"Bottom: electronic energy relative to separated reactants, {curve['label']} (the default level of "
         f"this reaction), against the mass-weighted reaction coordinate s; red points mark the panels. "
         f"TS imaginary frequency {abs(M['imaginary_frequency_cm']):.0f}i cm⁻¹."]
    kh, half, emax, kmax = _half_height(M)
    L.append(f"Panel (b) is the frame on the reactant side whose energy at the default level is closest to half "
             f"the barrier height above the first frame: E(first) + ½[E(max) − E(first)] = {half:.1f} kJ/mol "
             f"({half / 4.184:.1f} kcal/mol), with E(max) = {emax:.1f} kJ/mol ({emax / 4.184:.1f} kcal/mol) the highest "
             f"point of the default-level curve." + ("" if kmax == M["ts_frame"] else
             f" At this level the highest point is IRC frame {_irc_no(M, kmax)}, one step past the working-level "
             f"TS (frame {_irc_no(M, M['ts_frame'])}, {E[M['ts_frame']]:.1f} kJ/mol); panel (c) shows the "
             "working-level TS."))
    rid = rx.cfg["id"]
    if rid == "e2_ethoxide_2_bromopropane":
        L.append("Panels (e)–(g): the same frames as (a)–(c), viewed along the Cβ–Cα bond (Newman view, "
                 "Cβ in front, Cα directly behind it), with the transferring β-H pointing up; the "
                 "anti-periplanar Br points down. The last frame is not repeated in this view: the bromide "
                 "has left and lies outside the Newman framing (atoms within 4 Å of the Cβ–Cα axis).")
        for x in M.get("irc_frames_not_shown") or []:
            L.append(f"IRC frame {x['frame']} ({x['steps_from_TS']} step from the TS, {x['dE_above_TS_kJ']:.1f} "
                     "kJ/mol above it at the working level) is not shown in the animation and is not used here; "
                     "frame numbers are IRC frame numbers.")
    if rid == "da_butadiene_ethylene":
        L.append("Panels (e) and (f): IRC frames 59 and 60, one and two steps after the TS, on either side of "
                 "the abrupt change of the localized orbitals. Lowest frame-to-frame overlaps of the shown "
                 "orbitals, all at this step: π(C1=C2) → σ(C1–C6) 0.835, π(C5=C6) → σ(C4–C5) 0.838, "
                 "π(C3=C4) → π(C2=C3) 0.827 (lowest). An IBO single point at a quarter of the step shows "
                 "the switch happens within its first quarter (overlap 0.785–0.794, then 0.996–0.997).")
        L.append("Display rotation: at the reactant end the ethylene C=C pair, which IBO localizes as two "
                 "bent bonds, is shown as the equivalent σ + π pair (π shown) until it straightens.")
    L.append("Rendered with Blender (EEVEE), orthographic camera, from the same meshes as the viewer and video.")
    Path(f"{out}.md").write_text("\n\n".join(L) + "\n", encoding="utf-8")


def figure(rx: Reaction) -> None:
    M = _read_json(rx.root / "04_meshes" / "manifest.json")
    rows = _panels(rx, M)
    proj = render_panels(rx, M, rows)
    out = compose(rx, M, rows, proj)
    caption(rx, M, rows, out)
    print(f"figure: {out}.png, .pdf, .md")


def viewer_screenshot(rid: str, t: float, out: Path, port: int = 8766) -> None:
    """1920 × 1080 screenshot of the viewer (headless Chrome, software WebGL), opened at frame t
    (viewer URL parameter `t`) with the default orbitals on."""
    import time
    root = Path(__file__).resolve().parents[1]
    srv = subprocess.Popen([str(root / ".venv/Scripts/python.exe"), "-m", "http.server", str(port)],
                           cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        time.sleep(1.5)
        chrome = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
        url = f"http://localhost:{port}/web/viewer/index.html?r={rid}&t={t}"
        subprocess.run([chrome, "--headless=new", "--hide-scrollbars", "--window-size=1920,1080",
                        "--force-device-scale-factor=1", "--enable-unsafe-swiftshader", "--use-angle=swiftshader",
                        "--virtual-time-budget=20000", f"--screenshot={out}", url], check=True, timeout=180)
    finally:
        srv.terminate()
