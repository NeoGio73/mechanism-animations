"""Rendered video: Blender (headless) frames + energy-diagram inset with moving marker → mp4."""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

from .steps import Reaction, _read_json

BLENDER = Path(r"C:\Program Files\Blender Foundation\Blender 5.2\blender.exe")
FPS = 30


def timeline(s: list[float], hold0=1.5, path=16.0, hold1=2.5) -> list[float]:
    """Fractional IRC-frame positions per video frame. Constant speed along the mass-weighted
    path length s (IRC frames bunch up where the path flattens), eased at start and stop."""
    s = np.asarray(s, dtype=float)
    n_frames = len(s)
    T = []
    for i in range(int(round((hold0 + path + hold1) * FPS))):
        tt = i / FPS
        if tt < hold0:
            u = 0.0
        elif tt > hold0 + path:
            u = 1.0
        else:
            x = (tt - hold0) / path
            ease = 0.08                          # fraction of the path spent easing at each end
            if x < ease:
                u = x * x / (2 * ease * (1 - ease))
            elif x > 1 - ease:
                y = 1 - x
                u = 1 - y * y / (2 * ease * (1 - ease))
            else:
                u = (x - ease / 2) / (1 - ease)
        target = s[0] + u * (s[-1] - s[0])
        k = int(np.clip(np.searchsorted(s, target) - 1, 0, n_frames - 2))
        T.append(float(k + (target - s[k]) / (s[k + 1] - s[k])))
    return T


def render(rx: Reaction, preview: bool = False, only: list[int] | None = None) -> Path:
    mdir = rx.root / "04_meshes"
    odir = rx.d("05_render")
    M = _read_json(mdir / "manifest.json")
    (odir / "timeline.json").write_text(json.dumps(timeline([f["s"] for f in M["frames"]])))
    (odir / "render_settings.json").write_text(json.dumps(
        {"width": 1920, "height": 1080, "samples": 48, "orbital_kind": "localized",
         "shift_y": -0.02, "margin_A": 3.4, "only_frames": only}))
    script = Path(__file__).with_name("render_blender.py")
    cmd = [str(BLENDER), "-b", "--factory-startup", "--python", str(script), "--", str(mdir), str(odir)]
    if preview:
        cmd.append("preview")
    log = odir / ("blender_preview.log" if preview else "blender.log")
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0 or "RENDER_DONE" not in log.read_text(errors="replace"):
        raise RuntimeError(f"Blender render failed; see {log}")
    return odir / ("png_preview" if preview else "png")


def _font(size):
    from PIL import ImageFont
    import matplotlib
    p = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"
    return ImageFont.truetype(str(p), size)


def _curve(M, key):
    for c in M.get("curves", []):
        if c["key"] == key:
            return c
    raise KeyError(f"no energy curve '{key}' in manifest (run path_singlepoints, then meshes)")


def _inset_background(M, curve, w_px=700, h_px=300):
    """Static energy diagram; returns RGBA image and a function data→pixel."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    s = np.array([f["s"] for f in M["frames"]])
    E = np.array(curve["E_kJ"])
    ep = curve["E_separated_products_kJ"]
    dpi = 100
    fig = plt.figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    ax = fig.add_axes([0.13, 0.2, 0.84, 0.64])
    ax.plot(s, E, color="#1f4e79", lw=2.2)
    ax.axhline(0, color="#9aa4b1", lw=1, ls="--")
    ax.axhline(ep, color="#9aa4b1", lw=1, ls=":")
    ax.text(s[-1], 3, "separated reactants", ha="right", va="bottom", fontsize=7.5, color="#5b6878")
    ax.text(s[0], ep + 3, "separated products", ha="left", va="bottom", fontsize=7.5, color="#5b6878")
    top = int(np.argmax(E))
    ax.annotate((f"{'TS' if top == M['ts_frame'] else 'highest point'} {E[top]:.1f} kJ/mol "
                 f"({E[top] / 4.184:.1f} kcal/mol)").replace("-", "−"), (s[top], E[top]),
                xytext=(0, 7), textcoords="offset points", ha="center", fontsize=9, color="#1d2733")
    ax.set_xlabel("reaction coordinate s (bohr·amu$^{1/2}$)", fontsize=9)
    ax.set_ylabel("E (kJ/mol)", fontsize=9)
    ax.tick_params(labelsize=8)
    lo, hi = min(E.min(), ep), max(E.max(), 0)
    ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.25 * (hi - lo))
    ax.set_xlim(s[0], s[-1])
    title = curve["label"].replace(" single points on the same path", ", on the same path")
    ax.set_title(f"energy vs separated reactants · {title}", fontsize=8.5, color="#5b6878", loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.canvas.draw()
    img = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).copy()
    trans = ax.transData

    def to_px(x, y):
        px, py = trans.transform((x, y))
        return px, h_px - py

    plt.close(fig)
    return img, to_px


def composite(rx: Reaction, png_dir: Path, out_mp4: Path, curve_key: str = "working",
              stills: dict | None = None) -> None:
    """stills: {video frame index: output PNG path} — composited frames to keep as images."""
    import imageio_ffmpeg
    from PIL import Image, ImageDraw

    M = _read_json(rx.root / "04_meshes" / "manifest.json")
    T = json.loads((rx.root / "05_render" / "timeline.json").read_text())
    s = np.array([f["s"] for f in M["frames"]])
    curve = _curve(M, curve_key)
    E = np.array(curve["E_kJ"])
    xyz = [np.array(f["xyz"]) for f in M["frames"]]
    ka, sy = M["key_atoms"], M["symbols"]
    inset, to_px = _inset_background(M, curve)
    orbs = [o for o in M["orbitals"] if o["kind"] == "localized"]

    frames = sorted(png_dir.glob("f*.png"))
    first = Image.open(frames[0])
    W, H = first.size
    sc = W / 1920
    # macro_block_size=1: keep 1920 × 1080 exactly (H.264 handles 1080 fine)
    writer = imageio_ffmpeg.write_frames(str(out_mp4), (W, H), fps=FPS, codec="libx264",
                                         pix_fmt_out="yuv420p", quality=None, macro_block_size=1,
                                         output_params=["-crf", "18", "-preset", "slow",
                                                        "-movflags", "+faststart"])
    writer.send(None)
    ins = inset.resize((int(inset.width * sc), int(inset.height * sc)))
    ix, iy = W - ins.width - int(30 * sc), H - ins.height - int(24 * sc)
    for n, fp in enumerate(frames):
        t = T[int(fp.stem[1:])]
        if abs(t - round(t)) < 0.01:      # within 1 % of an IRC frame: read out that frame exactly
            t = float(round(t))
        f0 = min(int(math.floor(t)), len(s) - 1)
        f1 = min(f0 + 1, len(s) - 1)
        w = t - f0
        st, Et = s[f0] + (s[f1] - s[f0]) * w, E[f0] + (E[f1] - E[f0]) * w
        p = xyz[f0] + (xyz[f1] - xyz[f0]) * w
        dn = np.linalg.norm(p[ka["nucleophile"]] - p[ka["center"]])
        dl = np.linalg.norm(p[ka["center"]] - p[ka["leaving_group"]])

        im = Image.open(fp).convert("RGB")
        dr = ImageDraw.Draw(im)
        m = int(36 * sc)
        dr.text((m, m), M["equation"], fill="#1d2733", font=_font(int(40 * sc)))
        dr.text((m, m + int(54 * sc)), f"{M['title']} · minimum-energy path (IRC) · path and orbitals: {M['level']}",
                fill="#5b6878", font=_font(int(22 * sc)))
        y = H - int(150 * sc)
        dr.text((m, y), "Intrinsic bond orbitals, tracked along the path:", fill="#1d2733",
                font=_font(int(21 * sc)))
        sw = int(16 * sc)
        for k, o in enumerate(orbs):
            yy = y + int((32 + 30 * k) * sc)
            pc, nc = o.get("colors", ["#2f6db5", "#e08a1e"])
            x0 = m + int(14 * sc)
            dr.rectangle([x0, yy + 3, x0 + sw, yy + 3 + sw], fill=pc)
            dr.rectangle([x0 + sw, yy + 3, x0 + 2 * sw, yy + 3 + sw], fill=nc)
            dr.text((x0 + 2 * sw + int(10 * sc), yy), o["label"], fill="#1d2733", font=_font(int(21 * sc)))
        dr.text((m + int(14 * sc), y + int((36 + 30 * len(orbs)) * sc)),
                "dark shade = + phase, light shade = − phase", fill="#5b6878", font=_font(int(18 * sc)))

        panel = ins.copy()
        pd = ImageDraw.Draw(panel)
        px, py = to_px(st, Et)
        px, py = px * sc, py * sc
        pd.line([(px, 8 * sc), (px, panel.height - 50 * sc)], fill="#c0392b", width=max(1, int(1.5 * sc)))
        r = 8 * sc
        pd.ellipse([px - r, py - r, px + r, py + r], fill="#c0392b", outline="white", width=max(1, int(2 * sc)))
        im.paste(panel, (ix, iy), panel)
        nu = f"{sy[ka['center']]}–{sy[ka['nucleophile']]}"
        lgl = f"{sy[ka['center']]}–{sy[ka['leaving_group']]}"
        txt = (f"E = {Et:.1f} kJ/mol ({Et / 4.184:.1f} kcal/mol)    {nu} {dn:.2f} Å    "
               f"{lgl} {dl:.2f} Å").replace("-", "−")
        dr.text((ix + int(10 * sc), iy - int(30 * sc)), txt, fill="#1d2733", font=_font(int(20 * sc)))
        writer.send(np.asarray(im).tobytes())
        vf = int(fp.stem[1:])
        if stills and vf in stills:
            im.save(stills[vf])
    writer.close()


def key_still_frames(M, T) -> dict:
    """Video frame indices for reactant, halfway to TS, TS, product."""
    ts = M["ts_frame"]
    near = lambda x: int(np.argmin([abs(t - x) for t in T]))
    return {"reactant": 0, "before_TS": near(ts * 0.6), "TS": near(ts), "product": len(T) - 1}


def contact_sheet(paths: list[Path], out: Path) -> None:
    from PIL import Image
    ims = [Image.open(p) for p in paths]
    w, h = ims[0].size
    sheet = Image.new("RGB", (2 * w, 2 * h), "white")
    for i, im in enumerate(ims):
        sheet.paste(im, ((i % 2) * w, (i // 2) * h))
    sheet.resize((w, h)).save(out)


def _suffix(M, key: str) -> str:
    """'' for the default curve (<id>.mp4); otherwise '_B3LYP', '_M06-2X', '_DLPNO-CCSDT'."""
    if key == M.get("default_curve", "working"):
        return ""
    if key == "working":
        return "_" + M["level"].split("-")[0]
    return "_" + key.split("_")[0].replace("(", "").replace(")", "")


def video(rx: Reaction) -> None:
    """Render once in Blender, then composite one mp4 per energy curve."""
    render(rx)
    video_composite(rx)


def video_composite(rx: Reaction) -> None:
    """One mp4 per energy curve in the manifest, from the frames Blender already rendered
    (the molecule is identical; only the inset energies differ). The default curve
    (reaction.toml `default_curve`) is <id>.mp4."""
    png = rx.root / "05_render" / "png"
    M = _read_json(rx.root / "04_meshes" / "manifest.json")
    T = json.loads((rx.root / "05_render" / "timeline.json").read_text())
    if len(list(png.glob("f*.png"))) != len(T):
        raise RuntimeError(f"{png} does not hold all {len(T)} rendered frames; run the 'video' step")
    for c in M.get("curves", [{"key": "working"}]):
        suffix = _suffix(M, c["key"])
        out = rx.d("05_render") / f"{rx.cfg['id']}{suffix}.mp4"
        sd = rx.d(f"05_render/stills{suffix}")
        stills = {vf: sd / f"{i}_{name}.png" for i, (name, vf) in enumerate(key_still_frames(M, T).items())}
        composite(rx, png, out, c["key"], stills)
        contact_sheet(list(stills.values()), sd / "contact_sheet.png")
        print(f"video: {out} ({out.stat().st_size / 1e6:.1f} MB)")


def video_preview(rx: Reaction) -> None:
    """Half-resolution stills of a few key frames, to check framing before the full render."""
    M = _read_json(rx.root / "04_meshes" / "manifest.json")
    T = timeline([f["s"] for f in M["frames"]])
    want = [0, int(np.argmin([abs(t - M["ts_frame"] * 0.6) for t in T])),
            int(np.argmin([abs(t - M["ts_frame"]) for t in T])), len(T) - 1]
    png = render(rx, preview=True, only=want)
    composite(rx, png, rx.d("05_render") / "preview.mp4")
    # also keep the composited stills as PNGs for a quick look
    import imageio_ffmpeg  # noqa: F401
    from PIL import Image
    for i, fp in enumerate(sorted(png.glob("f*.png"))):
        Image.open(fp).save(rx.d("05_render") / f"preview_raw_{i}.png")
