"""Fingerprint of everything that decides what the Blender frames look like.

Standard library only: imported both by the pipeline and by render_blender.py inside Blender.

The renderer writes the fingerprint to <png dir>/inputs.sha256 and deletes every existing
frame when the stored one differs, so it never reuses frames rendered from other orbitals,
geometries or settings. video_composite refuses frames whose fingerprint does not match.
Energy curves are left out: they are drawn later by the compositor, not by Blender.
"""
import hashlib
import json
from pathlib import Path

STAMP = "inputs.sha256"
# manifest keys that only the compositor or the viewer uses (energies); everything else counts
_ENERGY_ONLY = {"curves", "default_curve", "E_separated_products_kJ", "energy_reference"}


def compute(mdir: Path, odir: Path, script: Path) -> str:
    """sha256 over meshes.bin, the render-relevant manifest, timeline.json,
    render_settings.json (without the frame subset) and the render script itself."""
    mdir, odir = Path(mdir), Path(odir)
    h = hashlib.sha256()
    h.update((mdir / "meshes.bin").read_bytes())
    man = json.loads((mdir / "manifest.json").read_text(encoding="utf-8"))
    man = {k: v for k, v in man.items() if k not in _ENERGY_ONLY}
    man["frames"] = [{k: v for k, v in f.items() if k != "E_kJ"} for f in man.get("frames", [])]
    h.update(json.dumps(man, sort_keys=True).encode())
    h.update((odir / "timeline.json").read_bytes())
    st = json.loads((odir / "render_settings.json").read_text())
    st.pop("only_frames", None)
    h.update(json.dumps(st, sort_keys=True).encode())
    h.update(Path(script).read_bytes())
    return h.hexdigest()


def stored(png: Path) -> str | None:
    p = Path(png) / STAMP
    return p.read_text().strip() if p.exists() else None
