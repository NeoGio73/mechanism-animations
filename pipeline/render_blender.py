"""Runs INSIDE Blender (headless):

  blender -b --factory-startup --python pipeline/render_blender.py -- <meshes_dir> <out_dir> [preview]

Reads 04_meshes/manifest.json + meshes.bin and <out_dir>/timeline.json (a list of fractional
IRC-frame positions t, one per video frame). Renders PNG frames <out_dir>/png/f0000.png ...
Ball-and-stick molecule, forming/breaking bonds dashed, tracked IBOs as two-colour
translucent lobes cross-faded between neighbouring IRC frames. White background.
"""
import json
import math
import sys
from pathlib import Path

import bpy
import numpy as np
from mathutils import Matrix, Vector

argv = sys.argv[sys.argv.index("--") + 1:]
MDIR, ODIR = Path(argv[0]), Path(argv[1])
PREVIEW = len(argv) > 2 and argv[2] == "preview"
M = json.loads((MDIR / "manifest.json").read_text(encoding="utf-8"))
BIN = (MDIR / "meshes.bin").read_bytes()
TIMELINE = json.loads((ODIR / "timeline.json").read_text())
SETTINGS = json.loads((ODIR / "render_settings.json").read_text())
N = len(M["frames"])

ELEM = {"H": ((0.86, 0.86, 0.86), 0.24, 0.31), "C": ((0.20, 0.20, 0.20), 0.36, 0.76),
        "N": ((0.18, 0.33, 0.79), 0.34, 0.71), "O": ((0.85, 0.18, 0.10), 0.34, 0.66),
        "Cl": ((0.18, 0.64, 0.29), 0.42, 1.02), "Br": ((0.58, 0.15, 0.15), 0.46, 1.20),
        "I": ((0.42, 0.18, 0.55), 0.50, 1.39), "F": ((0.50, 0.82, 0.42), 0.30, 0.57),
        "S": ((0.85, 0.72, 0.17), 0.42, 1.05)}
ORB_ALPHA = 0.6


def srgb_to_linear(c):
    return tuple(((x + 0.055) / 1.055) ** 2.4 if x > 0.04045 else x / 12.92 for x in c)


# ------------------------------------------------------------------ scene
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x, scene.render.resolution_y = SETTINGS["width"], SETTINGS["height"]
scene.render.resolution_percentage = 50 if PREVIEW else 100
scene.render.image_settings.file_format = "PNG"
scene.render.film_transparent = False
try:
    scene.eevee.taa_render_samples = SETTINGS.get("samples", 48)
except AttributeError:
    pass
try:
    scene.view_settings.view_transform = "Standard"
except TypeError:
    pass
world = bpy.data.worlds.new("white")
scene.world = world
world.use_nodes = True
bg = next(n for n in world.node_tree.nodes if n.type == "BACKGROUND")
bg.inputs[0].default_value = (1, 1, 1, 1)
bg.inputs[1].default_value = 1.0


def light(kind, energy, rot, size=None):
    ld = bpy.data.lights.new(kind, kind)
    ld.energy = energy
    if size is not None and hasattr(ld, "size"):
        ld.size = size
    ob = bpy.data.objects.new(kind, ld)
    ob.rotation_euler = rot
    scene.collection.objects.link(ob)
    return ob


light("SUN", 3.0, (math.radians(40), math.radians(10), math.radians(30)))
light("SUN", 1.0, (math.radians(-60), math.radians(-20), math.radians(200)))


def material(name, color, rough=0.4, alpha=None):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    bsdf = next(n for n in nt.nodes if n.type == "BSDF_PRINCIPLED")
    bsdf.inputs["Base Color"].default_value = (*srgb_to_linear(color), 1)
    bsdf.inputs["Roughness"].default_value = rough
    if alpha is not None:
        # per-object alpha (object.color[3]) times a base alpha → cross-fades between frames
        info = nt.nodes.new("ShaderNodeObjectInfo")
        mul = nt.nodes.new("ShaderNodeMath")
        mul.operation = "MULTIPLY"
        mul.inputs[1].default_value = alpha
        nt.links.new(info.outputs["Alpha"], mul.inputs[0])
        nt.links.new(mul.outputs[0], bsdf.inputs["Alpha"])
        mat.surface_render_method = "BLENDED"
        mat.use_backface_culling = False
    return mat


def mesh_object(name, verts, faces, mat, smooth=True):
    me = bpy.data.meshes.new(name)
    me.from_pydata([tuple(v) for v in verts], [], [tuple(f) for f in faces])
    me.update()
    if smooth:
        me.polygons.foreach_set("use_smooth", [True] * len(me.polygons))
    me.materials.append(mat)
    ob = bpy.data.objects.new(name, me)
    scene.collection.objects.link(ob)
    return ob


def unit_sphere():
    bpy.ops.mesh.primitive_uv_sphere_add(segments=48, ring_count=24, radius=1.0)
    ob = bpy.context.active_object
    me = ob.data
    me.polygons.foreach_set("use_smooth", [True] * len(me.polygons))
    bpy.data.objects.remove(ob)
    return me


def unit_cylinder():
    bpy.ops.mesh.primitive_cylinder_add(vertices=24, radius=1.0, depth=1.0)
    ob = bpy.context.active_object
    me = ob.data
    me.polygons.foreach_set("use_smooth", [len(p.vertices) == 4 for p in me.polygons])
    bpy.data.objects.remove(ob)
    return me


SPHERE, CYL = unit_sphere(), unit_cylinder()
MATS = {s: material(f"atom_{s}", c, 0.35) for s, (c, _, _) in ELEM.items()}
BOND = material("bond", (0.62, 0.64, 0.67), 0.45)
PARTIAL = material("partial", (0.62, 0.64, 0.67), 0.45)
def hex_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


# one material pair per orbital (same colours as the viewer): dark = + phase, light = − phase
LOBE_MATS = {}
for o in M["orbitals"]:
    pc, nc = o.get("colors", ["#2f6db5", "#e08a1e"])
    LOBE_MATS[o["key"]] = (material(f"{o['key']}_pos", hex_rgb(pc), 0.3, alpha=ORB_ALPHA),
                           material(f"{o['key']}_neg", hex_rgb(nc), 0.3, alpha=ORB_ALPHA))

# shared meshes carry one slot; each object links its own material to that slot
SPHERE.materials.append(MATS["C"])
CYL.materials.append(BOND)
atoms = []
for i, s in enumerate(M["symbols"]):
    ob = bpy.data.objects.new(f"atom{i}", SPHERE)
    scene.collection.objects.link(ob)
    ob.material_slots[0].link = "OBJECT"
    ob.material_slots[0].material = MATS[s]
    ob.scale = (ELEM[s][1],) * 3
    atoms.append(ob)

ka = M["key_atoms"]


def _d(f, a, b):
    return float(np.linalg.norm(np.array(M["frames"][f]["xyz"][a]) - np.array(M["frames"][f]["xyz"][b])))


# forming/breaking bonds judged against their own equilibrium lengths (same rule as the viewer)
# (forming bond: product end; breaking bond: reactant end). A bond that only shortens (C–C → C=C)
# stays an ordinary bond.
KB = M.get("key_bonds") or [{"atoms": [ka["center"], ka["nucleophile"]], "role": "forms"},
                            {"atoms": [ka["center"], ka["leaving_group"]], "role": "breaks"}]
key_eq = {tuple(sorted(b["atoms"])): _d(N - 1 if b["role"] == "forms" else 0, *b["atoms"])
          for b in KB if b["role"] in ("forms", "breaks")}
bonds = []
for i in range(len(atoms)):
    for j in range(i + 1, len(atoms)):
        solid = bpy.data.objects.new(f"bond{i}_{j}", CYL)
        scene.collection.objects.link(solid)
        dashes = []
        for k in range(5):
            d = bpy.data.objects.new(f"dash{i}_{j}_{k}", CYL)
            scene.collection.objects.link(d)
            d.material_slots[0].link = "OBJECT"
            d.material_slots[0].material = PARTIAL
            dashes.append(d)
        cov = ELEM[M["symbols"][i]][2] + ELEM[M["symbols"][j]][2]
        bonds.append({"i": i, "j": j, "solid": solid, "dashes": dashes, "cov": cov,
                      "eq": key_eq.get((i, j))})

Z = Vector((0, 0, 1))


def place(ob, a, b, r):
    d = b - a
    ob.location = a + d * 0.5
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = Z.rotation_difference(d.normalized())
    ob.scale = (r, r, d.length)


def positions(t):
    f0 = min(int(math.floor(t)), N - 1)
    f1 = min(f0 + 1, N - 1)
    w = t - f0
    a = np.array(M["frames"][f0]["xyz"])
    b = np.array(M["frames"][f1]["xyz"])
    return [Vector(p) for p in a + (b - a) * w]


def set_molecule(t):
    pos = positions(t)
    for ob, p in zip(atoms, pos):
        ob.location = p
    for b in bonds:
        d = (pos[b["i"]] - pos[b["j"]]).length
        full = d < 1.15 * (b["eq"] or b["cov"])
        partial = (not full) and b["eq"] is not None and d < 1.6 * b["eq"]
        b["solid"].hide_render = not full
        if full:
            place(b["solid"], pos[b["i"]], pos[b["j"]], 0.085)
        for k, dob in enumerate(b["dashes"]):
            dob.hide_render = not partial
            if partial:
                a, c = pos[b["i"]], pos[b["j"]]
                u0, u1 = (k + 0.2) / 5, (k + 0.8) / 5
                place(dob, a + (c - a) * u0, a + (c - a) * u1, 0.06)


# ------------------------------------------------------------------ orbitals
mid = np.array(M["quant"]["mid"])
scale = np.array(M["quant"]["scale"])
ORBS = [o["key"] for o in M["orbitals"] if o["kind"] == SETTINGS.get("orbital_kind", "localized")]
cache = {}


def lobe(entry):
    ov, nv = entry["v"]
    of, nf = entry["f"]
    if not nv:
        return None, None
    q = np.frombuffer(BIN, dtype="<i2", count=nv * 3, offset=ov).reshape(-1, 3)
    f = np.frombuffer(BIN, dtype="<u2", count=nf * 3, offset=of).reshape(-1, 3)
    return mid + q * scale, f


def orbital_objects(key, frame):
    ck = (key, frame)
    if ck not in cache:
        obs = []
        e = M["mesh_index"][key][frame]
        pos_mat, neg_mat = LOBE_MATS[key]
        for part, mat in (("pos", pos_mat), ("neg", neg_mat)):
            v, f = lobe(e[part])
            if v is not None:
                obs.append(mesh_object(f"{key}_{frame}_{part}", v, f, mat))
        cache[ck] = obs
    return cache[ck]


def set_orbitals(t):
    f0 = min(int(math.floor(t)), N - 1)
    f1 = min(f0 + 1, N - 1)
    w = t - f0
    shown = {}
    for key in ORBS:
        for fr, a in ((f0, math.sqrt(1 - w)), (f1, math.sqrt(w))):
            if a > 0.02:
                for ob in orbital_objects(key, fr):
                    shown[ob.name] = max(shown.get(ob.name, 0), a)
    for obs in cache.values():
        for ob in obs:
            a = shown.get(ob.name, 0.0)
            ob.hide_render = a <= 0.0
            ob.color = (1, 1, 1, a)


# ------------------------------------------------------------------ camera (side view)
ts = M["frames"][M["ts_frame"]]["xyz"]
nu, lg, cc = Vector(ts[ka["nucleophile"]]), Vector(ts[ka["leaving_group"]]), Vector(ts[ka["center"]])
axis = (lg - nu).normalized()
side = Vector((0, 0, 1)).cross(axis)
if side.length < 1e-3:
    side = Vector((1, 0, 0)).cross(axis)
side.normalize()
# Several changing bonds (not all on one line): look straight at the plane they lie in.
kat = sorted({a for b in KB if b["role"] in ("forms", "breaks") for a in b["atoms"]})
P = np.array([ts[a] for a in kat])
if len(kat) >= 4:
    sv = np.linalg.svd(P - P.mean(0))
    if sv[1][1] > 0.25 * sv[1][0]:
        n = Vector(sv[2][2].tolist())
        n = n - axis * n.dot(axis)
        if n.length > 1e-3:
            side = n.normalized() if n.dot(side) >= 0 else -n.normalized()
third = axis.cross(side).normalized()
view_dir = side                     # camera sits on +side, looks back at the molecule
up = -third                         # nucleophile appears on the left, leaving group on the right
forward = -view_dir
right = forward.cross(up).normalized()
up2 = right.cross(forward).normalized()
center = Vector(np.mean([np.array(M["frames"][k]["xyz"]).mean(0) for k in (0, N - 1)], axis=0))
center = center + (cc - center) * 0.5
cam_data = bpy.data.cameras.new("cam")
cam_data.lens = 50
cam_data.shift_y = SETTINGS.get("shift_y", -0.10)
cam = bpy.data.objects.new("cam", cam_data)
scene.collection.objects.link(cam)
scene.camera = cam
extent = max(max((Vector(p) - center).dot(axis) for p in M["frames"][k]["xyz"]) -
             min((Vector(p) - center).dot(axis) for p in M["frames"][k]["xyz"]) for k in (0, N - 1))
half_w = extent / 2 + SETTINGS.get("margin_A", 2.4)
fov_h = 2 * math.atan(cam_data.sensor_width / (2 * cam_data.lens))
dist = half_w / math.tan(fov_h / 2)
# Tall molecules (e.g. a diene below a dienophile): if the atoms plus room for orbital lobes
# would reach the title above or the energy inset (bottom right), move back and recentre, in
# whichever free area shows the molecule larger: centred above the inset, or in the left part of
# the frame beside the inset and above the orbital legend.
fr = SETTINGS.get("fit_rows")
if fr:
    Wpx, Hpx = SETTINGS["width"], SETTINGS["height"]
    pad = SETTINGS.get("fit_pad_A", 1.6)
    allp = [Vector(p) - center for f in M["frames"] for p in f["xyz"]]
    u_hi = max(q.dot(up2) for q in allp) + pad
    u_lo = min(q.dot(up2) for q in allp) - pad
    v_hi = max(q.dot(right) for q in allp) + pad
    v_lo = min(q.dot(right) for q in allp) - pad
    tanh = math.tan(fov_h / 2)
    ppa0 = (Wpx / 2) / (dist * tanh)                                   # pixels per Å now
    top, bot_centre, bot_left = fr
    c0, c1 = SETTINGS["fit_cols_left"]
    y_hi = Hpx / 2 - u_hi * ppa0 + cam_data.shift_y * Wpx
    y_lo = Hpx / 2 - u_lo * ppa0 + cam_data.shift_y * Wpx
    if y_hi < top or y_lo > bot_centre:
        side = 2 * max(v_hi, -v_lo)
        opts = {"centre": (min(ppa0, (bot_centre - top) / (u_hi - u_lo), (Wpx - 80) / side), Wpx / 2, (top + bot_centre) / 2),
                "left": (min(ppa0, (bot_left - top) / (u_hi - u_lo), (c1 - c0) / (v_hi - v_lo)), (c0 + c1) / 2, (top + bot_left) / 2)}
        name, (ppa, xc, yc) = max(opts.items(), key=lambda kv: kv[1][0])
        dist = (Wpx / 2) / (ppa * tanh)
        vc = 0.0 if name == "centre" else (v_hi + v_lo) / 2
        cam_data.shift_x = (Wpx / 2 + vc * ppa - xc) / Wpx
        cam_data.shift_y = (yc - Hpx / 2 + (u_hi + u_lo) / 2 * ppa) / Wpx
        print(f"FRAMING: {name}; {ppa:.0f} px/Å (was {ppa0:.0f}); shift_x {cam_data.shift_x:.3f}, "
              f"shift_y {cam_data.shift_y:.3f}")
cam.location = center + view_dir * dist
cam.matrix_world = Matrix.Translation(cam.location) @ Matrix((right, up2, -forward)).transposed().to_4x4()

# ------------------------------------------------------------------ render loop
png = ODIR / ("png_preview" if PREVIEW else "png")
png.mkdir(parents=True, exist_ok=True)
# Never reuse frames rendered from other inputs: compare the stored fingerprint of meshes,
# geometry, timeline, settings and this script; on any mismatch delete every old frame.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_stamp  # noqa: E402

stamp = render_stamp.compute(MDIR, ODIR, Path(__file__).resolve())
old_stamp = render_stamp.stored(png)
old_frames = sorted(png.glob("f*.png"))
if old_frames and old_stamp != stamp:
    for p in old_frames:
        p.unlink()
    print(f"STALE_FRAMES_REMOVED {len(old_frames)} (stored fingerprint "
          f"{(old_stamp or 'none')[:12]}, current {stamp[:12]})", flush=True)
(png / render_stamp.STAMP).write_text(stamp)
print(f"INPUTS_SHA256 {stamp}", flush=True)
only = SETTINGS.get("only_frames")
gpu_logged = False
for n, t in enumerate(TIMELINE):
    if only is not None and n not in only:
        continue
    out = png / f"f{n:04d}.png"
    if out.exists():
        continue
    set_molecule(t)
    set_orbitals(t)
    scene.render.filepath = str(out)
    bpy.ops.render.render(write_still=True)
    if not gpu_logged:
        # the GPU context only exists after the first render; record what EEVEE ran on
        gpu_logged = True
        try:
            import gpu
            print(f"RENDER_DEVICE engine={scene.render.engine} backend={gpu.platform.backend_type_get()} "
                  f"gpu={gpu.platform.renderer_get()}", flush=True)
        except Exception as e:
            print(f"RENDER_DEVICE engine={scene.render.engine} gpu=unknown ({e})", flush=True)
    if n % 25 == 0:
        print(f"RENDERED {n}/{len(TIMELINE)}", flush=True)
print("RENDER_DONE", flush=True)
