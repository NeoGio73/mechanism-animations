"""Runs INSIDE Blender (headless), for publication figures (pipeline/figure.py):

  blender -b --factory-startup --python pipeline/render_figure_blender.py -- <meshes_dir> <out_dir>

Reads 04_meshes/manifest.json + meshes.bin and <out_dir>/figure_settings.json:
  {"width": px, "height": px, "samples": n, "pad_A": Å,
   "panels": [{"name": str, "frame": int, "view": "side" | "newman"}, ...]}
Renders one PNG per panel (<out_dir>/<name>.png) at an exact IRC frame (no cross-fade), and
<out_dir>/projection.json with the pixel position of every atom in every panel and the forming/
breaking bonds drawn dashed (for labels and, with "dashes_3d": false, 2D dashed overlays that stay
visible over the orbital lobes).

Same styling as the video (pipeline/render_blender.py: atom colours and radii, bond rule, dashed
forming/breaking bonds, two-colour translucent lobes, white background), but an orthographic
camera framed on the panels themselves, with no room left for video overlays. Kept as a separate
script so the video's render fingerprint (pipeline/render_stamp.py) is not touched.
"""
import json
import math
import sys
from pathlib import Path

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Matrix, Vector

argv = sys.argv[sys.argv.index("--") + 1:]
MDIR, ODIR = Path(argv[0]), Path(argv[1])
M = json.loads((MDIR / "manifest.json").read_text(encoding="utf-8"))
BIN = (MDIR / "meshes.bin").read_bytes()
S = json.loads((ODIR / "figure_settings.json").read_text(encoding="utf-8"))
N = len(M["frames"])
ORB_ALPHA = S.get("orbital_alpha", 0.6)

ELEM = {"H": ((0.86, 0.86, 0.86), 0.24, 0.31), "C": ((0.20, 0.20, 0.20), 0.36, 0.76),
        "N": ((0.18, 0.33, 0.79), 0.34, 0.71), "O": ((0.85, 0.18, 0.10), 0.34, 0.66),
        "Cl": ((0.18, 0.64, 0.29), 0.42, 1.02), "Br": ((0.58, 0.15, 0.15), 0.46, 1.20),
        "I": ((0.42, 0.18, 0.55), 0.50, 1.39), "F": ((0.50, 0.82, 0.42), 0.30, 0.57),
        "S": ((0.85, 0.72, 0.17), 0.42, 1.05)}


def srgb_to_linear(c):
    return tuple(((x + 0.055) / 1.055) ** 2.4 if x > 0.04045 else x / 12.92 for x in c)


bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x, scene.render.resolution_y = S["width"], S["height"]
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.film_transparent = False
try:
    scene.eevee.taa_render_samples = S.get("samples", 64)
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

LIGHTS = []


def light(kind, energy, rot):
    ld = bpy.data.lights.new(kind, kind)
    ld.energy = energy
    ob = bpy.data.objects.new(kind, ld)
    ob.rotation_euler = rot
    scene.collection.objects.link(ob)
    LIGHTS.append(ob)
    return ob


light("SUN", 3.0, (math.radians(40), math.radians(10), math.radians(30)))
light("SUN", 1.0, (math.radians(-60), math.radians(-20), math.radians(200)))
LIGHT_ROT0 = [ob.rotation_euler.to_matrix() for ob in LIGHTS]


def material(name, color, rough=0.4, alpha=None):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    bsdf = next(n for n in nt.nodes if n.type == "BSDF_PRINCIPLED")
    bsdf.inputs["Base Color"].default_value = (*srgb_to_linear(color), 1)
    bsdf.inputs["Roughness"].default_value = rough
    if alpha is not None:
        bsdf.inputs["Alpha"].default_value = alpha
        mat.surface_render_method = "BLENDED"
        mat.use_backface_culling = False
    return mat


def mesh_object(name, verts, faces, mat):
    me = bpy.data.meshes.new(name)
    me.from_pydata([tuple(v) for v in verts], [], [tuple(f) for f in faces])
    me.update()
    me.polygons.foreach_set("use_smooth", [True] * len(me.polygons))
    me.materials.append(mat)
    ob = bpy.data.objects.new(name, me)
    scene.collection.objects.link(ob)
    return ob


def unit_mesh(op, **kw):
    op(**kw)
    ob = bpy.context.active_object
    me = ob.data
    bpy.data.objects.remove(ob)
    return me


SPHERE = unit_mesh(bpy.ops.mesh.primitive_uv_sphere_add, segments=64, ring_count=32, radius=1.0)
SPHERE.polygons.foreach_set("use_smooth", [True] * len(SPHERE.polygons))
CYL = unit_mesh(bpy.ops.mesh.primitive_cylinder_add, vertices=32, radius=1.0, depth=1.0)
CYL.polygons.foreach_set("use_smooth", [len(p.vertices) == 4 for p in CYL.polygons])
MATS = {s: material(f"atom_{s}", c, 0.35) for s, (c, _, _) in ELEM.items()}
BOND = material("bond", (0.62, 0.64, 0.67), 0.45)


def hex_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


ORBS = [o for o in M["orbitals"] if o["kind"] == S.get("orbital_kind", "localized")]
LOBE_MATS = {o["key"]: (material(f"{o['key']}_pos", hex_rgb(o.get("colors", ["#2f6db5", "#e08a1e"])[0]), 0.3, ORB_ALPHA),
                        material(f"{o['key']}_neg", hex_rgb(o.get("colors", ["#2f6db5", "#e08a1e"])[1]), 0.3, ORB_ALPHA))
             for o in ORBS}

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
XYZ = [np.array(f["xyz"]) for f in M["frames"]]
KB = M.get("key_bonds") or [{"atoms": [ka["center"], ka["nucleophile"]], "role": "forms"},
                            {"atoms": [ka["center"], ka["leaving_group"]], "role": "breaks"}]
key_eq = {tuple(sorted(b["atoms"])): float(np.linalg.norm(XYZ[N - 1 if b["role"] == "forms" else 0][b["atoms"][0]]
                                                         - XYZ[N - 1 if b["role"] == "forms" else 0][b["atoms"][1]]))
          for b in KB if b["role"] in ("forms", "breaks")}
bonds = []
for i in range(len(atoms)):
    for j in range(i + 1, len(atoms)):
        solid = bpy.data.objects.new(f"bond{i}_{j}", CYL)
        scene.collection.objects.link(solid)
        dashes = [bpy.data.objects.new(f"dash{i}_{j}_{k}", CYL) for k in range(5)]
        for d in dashes:
            scene.collection.objects.link(d)
        bonds.append({"i": i, "j": j, "solid": solid, "dashes": dashes,
                      "cov": ELEM[M["symbols"][i]][2] + ELEM[M["symbols"][j]][2], "eq": key_eq.get((i, j))})
Z = Vector((0, 0, 1))


def place(ob, a, b, r):
    d = b - a
    ob.location = a + d * 0.5
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = Z.rotation_difference(d.normalized())
    ob.scale = (r, r, d.length)


def set_molecule(f):
    """Same bond rule as the video and the viewer: solid below 1.15 × the reference length;
    a forming/breaking bond dashed up to 1.6 × its own equilibrium length."""
    pos = [Vector(p) for p in XYZ[f]]
    for ob, p in zip(atoms, pos):
        ob.location = p
    for b in bonds:
        d = (pos[b["i"]] - pos[b["j"]]).length
        full = d < 1.15 * (b["eq"] or b["cov"])
        partial = (not full) and b["eq"] is not None and d < 1.6 * b["eq"]
        b["solid"].hide_render = not full
        b["partial"] = partial
        if full:
            place(b["solid"], pos[b["i"]], pos[b["j"]], 0.085)
        for k, dob in enumerate(b["dashes"]):
            dob.hide_render = not (partial and S.get("dashes_3d", True))
            if partial and S.get("dashes_3d", True):
                a, c = pos[b["i"]], pos[b["j"]]
                place(dob, a + (c - a) * (k + 0.2) / 5, a + (c - a) * (k + 0.8) / 5, 0.06)


mid = np.array(M["quant"]["mid"])
qscale = np.array(M["quant"]["scale"])


def lobe(entry):
    ov, nv = entry["v"]
    of, nf = entry["f"]
    if not nv:
        return None, None
    q = np.frombuffer(BIN, dtype="<i2", count=nv * 3, offset=ov).reshape(-1, 3)
    f = np.frombuffer(BIN, dtype="<u2", count=nf * 3, offset=of).reshape(-1, 3)
    return mid + q * qscale, f


orb_obs = {}


def set_orbitals(f):
    for obs in orb_obs.values():
        for ob in obs:
            ob.hide_render = True
    for o in ORBS:
        key = o["key"]
        if (key, f) not in orb_obs:
            obs = []
            for part, mat in zip(("pos", "neg"), LOBE_MATS[key]):
                v, fc = lobe(M["mesh_index"][key][f][part])
                if v is not None:
                    obs.append(mesh_object(f"{key}_{f}_{part}", v, fc, mat))
            orb_obs[(key, f)] = obs
        for ob in orb_obs[(key, f)]:
            ob.hide_render = False


def lobe_points(f):
    pts = [XYZ[f]]
    for o in ORBS:
        for part in ("pos", "neg"):
            v, _ = lobe(M["mesh_index"][o["key"]][f][part])
            if v is not None:
                pts.append(v)
    return np.vstack(pts)


# ------------------------------------------------------------------ camera bases
def side_basis():
    """The video's side view (pipeline/render_blender.py), from the TS geometry."""
    ts = XYZ[M["ts_frame"]]
    nu, lg = Vector(ts[ka["nucleophile"]]), Vector(ts[ka["leaving_group"]])
    axis = (lg - nu).normalized()
    side = Vector((0, 0, 1)).cross(axis)
    if side.length < 1e-3:
        side = Vector((1, 0, 0)).cross(axis)
    side.normalize()
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
    forward, up = -side, -third
    right = forward.cross(up).normalized()
    return right, right.cross(forward).normalized(), forward


def newman_basis(f):
    """Looking along view_axis [front, back] (E2: Cβ in front of Cα); the up direction is the
    front atom's transferring H (key bond 'breaks' at the front atom), so the anti-periplanar
    leaving group points down."""
    a, b = M["view_axis"]
    x = XYZ[f]
    forward = Vector((x[b] - x[a]).tolist()).normalized()
    hb = next((bb["atoms"] for bb in KB if bb["role"] == "breaks" and a in bb["atoms"]), None)
    h = hb[1] if hb and hb[0] == a else (hb[0] if hb else None)
    upv = Vector((x[h] - x[a]).tolist()) if h is not None else Vector((0, 0, 1))
    up = (upv - forward * upv.dot(forward)).normalized()
    right = forward.cross(up).normalized()
    return right, right.cross(forward).normalized(), forward


cam_data = bpy.data.cameras.new("cam")
cam_data.type = "ORTHO"
cam = bpy.data.objects.new("cam", cam_data)
scene.collection.objects.link(cam)
scene.camera = cam
W, H = S["width"], S["height"]
pad = S.get("pad_A", 0.35)


def frame_extent(basis, frames, centre=None, radius=None):
    """Bounding box of atoms and lobes in the view plane; with `radius`, only points within that
    distance of the view axis count (Newman view: frame the Cβ–Cα unit, not distant fragments)."""
    right, up, _ = basis
    R, U = np.array(right), np.array(up)
    xs, ys = [], []
    for f in frames:
        P = lobe_points(f)
        c = np.zeros(3) if centre is None else centre(f)
        x, y = (P - c) @ R, (P - c) @ U
        if radius is not None:
            keep = np.hypot(x, y) <= radius
            x, y = x[keep], y[keep]
        xs.append(x)
        ys.append(y)
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    return xs.min() - pad, xs.max() + pad, ys.min() - pad, ys.max() + pad


panels = S["panels"]
by_view = {}
for p in panels:
    by_view.setdefault(p["view"], []).append(p["frame"])
SIDE = side_basis()
setup = {}
if "side" in by_view:
    x0, x1, y0, y1 = frame_extent(SIDE, by_view["side"])
    setup["side"] = {"scale": max(x1 - x0, (y1 - y0) * W / H), "cx": (x0 + x1) / 2, "cy": (y0 + y1) / 2}
if "newman" in by_view:
    a = M["view_axis"][0]
    cen = lambda f: XYZ[f][a]
    sc = 0
    for f in by_view["newman"]:
        x0, x1, y0, y1 = frame_extent(newman_basis(f), [f], cen, S.get("newman_radius_A", 4.0))
        sc = max(sc, 2 * max(-x0, x1), 2 * max(-y0, y1) * W / H)
    setup["newman"] = {"scale": sc}


def aim(p):
    if p["view"] == "side":
        right, up, forward = SIDE
        st = setup["side"]
        c = Vector(np.array(right) * st["cx"] + np.array(up) * st["cy"])
    else:
        right, up, forward = newman_basis(p["frame"])
        st = setup["newman"]
        c = Vector(XYZ[p["frame"]][M["view_axis"][0]].tolist())
    cam_data.ortho_scale = st["scale"]
    cam_data.clip_end = 1000
    loc = c - forward * 60
    cam.matrix_world = Matrix.Translation(loc) @ Matrix((right, up, -forward)).transposed().to_4x4()
    # lights keep the same direction relative to the camera in every view
    rot = Matrix((right, up, -forward)).transposed()
    side_rot = Matrix((SIDE[0], SIDE[1], -SIDE[2])).transposed()
    for ob, r0 in zip(LIGHTS, LIGHT_ROT0):
        ob.rotation_euler = (rot @ side_rot.inverted() @ r0).to_euler()


proj = {}
for p in panels:
    set_molecule(p["frame"])
    set_orbitals(p["frame"])
    aim(p)
    bpy.context.view_layer.update()
    out = ODIR / f"{p['name']}.png"
    scene.render.filepath = str(out)
    bpy.ops.render.render(write_still=True)
    pts = []
    for ob in atoms:
        v = world_to_camera_view(scene, cam, ob.location)
        pts.append([v.x * W, (1 - v.y) * H, v.z])
    proj[p["name"]] = {"frame": p["frame"], "view": p["view"], "atoms_px": pts,
                       "px_per_A": W / cam_data.ortho_scale,
                       "atom_radius_A": [ELEM[s][1] for s in M["symbols"]],
                       "partial_bonds": [[b["i"], b["j"]] for b in bonds if b["partial"]]}
    print(f"RENDERED {p['name']} frame {p['frame']} ({p['view']})", flush=True)
(ODIR / "projection.json").write_text(json.dumps(proj, indent=1), encoding="utf-8")
print("FIGURE_RENDER_DONE", flush=True)
