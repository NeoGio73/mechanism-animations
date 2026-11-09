// Mechanism viewer: molecule + tracked orbitals follow the energy-diagram scrubber.
// Data: ../../reactions/<id>/04_meshes/{manifest.json, meshes.bin}  (see pipeline/meshes.py)
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

const rid = new URLSearchParams(location.search).get('r') || 'sn2_hydroxide_bromomethane';
const base = `../../reactions/${rid}/04_meshes/`;
const KCAL = 4.184;

const ELEM = {
  H: { c: 0xd8d8d8, r: 0.24, cov: 0.31 }, C: { c: 0x3a3a3a, r: 0.36, cov: 0.76 },
  N: { c: 0x2f55c9, r: 0.34, cov: 0.71 }, O: { c: 0xd8321e, r: 0.34, cov: 0.66 },
  F: { c: 0x7fd06a, r: 0.30, cov: 0.57 }, S: { c: 0xd9b72b, r: 0.42, cov: 1.05 },
  Cl: { c: 0x2fa34a, r: 0.42, cov: 1.02 }, Br: { c: 0x9b2a2a, r: 0.46, cov: 1.20 },
  I: { c: 0x6b2d8c, r: 0.50, cov: 1.39 },
};
const el = (s) => ELEM[s] || { c: 0xff66cc, r: 0.35, cov: 0.8 };

const $ = (id) => document.getElementById(id);
const state = { t: 0, kind: 'localized', orbOn: true, labels: false, dash: true, playing: false,
                hidden: new Set() };

let M, BIN, N;
let scene, camera, renderer, labelRenderer, controls;
const atomMeshes = [], labelObjs = [], bonds = [];
const orbCache = new Map();
let orbGroup;

init().catch((e) => { $('loading').textContent = 'Could not load data: ' + e.message; console.error(e); });

async function init() {
  const [m, b] = await Promise.all([
    fetch(base + 'manifest.json').then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); }),
    fetch(base + 'meshes.bin').then((r) => { if (!r.ok) throw new Error(r.status); return r.arrayBuffer(); }),
  ]);
  M = m; BIN = b; N = M.frames.length;
  document.title = `${M.title} · Mechanism Viewer`;
  $('title').textContent = M.title;
  $('equation').textContent = M.equation;
  $('cardLink').href = `../../reactions/${rid}/card.md`;
  $('footnote').textContent = `${M.level}. Orbital isosurfaces at ±${M.isovalue} ${M.isovalue_units}. ` +
    `TS imaginary frequency ${Math.abs(M.imaginary_frequency_cm).toFixed(0)}i cm⁻¹. ` +
    `Energies are electronic, in kJ/mol relative to ${M.energy_reference}.`;
  setupScene();
  buildMolecule();
  buildOrbitalList();
  setupDiagram();
  setupControls();
  state.t = 0;
  setCamera('side');
  update();
  $('loading').classList.add('hidden');
  animate();
}

// ------------------------------------------------------------------ scene
function setupScene() {
  const host = $('three');
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0xffffff);
  camera = new THREE.PerspectiveCamera(32, 1, 0.1, 200);
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  host.appendChild(renderer.domElement);
  labelRenderer = new CSS2DRenderer();
  labelRenderer.domElement.style.position = 'absolute';
  labelRenderer.domElement.style.inset = '0';
  labelRenderer.domElement.style.pointerEvents = 'none';
  host.appendChild(labelRenderer.domElement);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xffffff, 0x8899aa, 1.6));
  const d1 = new THREE.DirectionalLight(0xffffff, 1.6); d1.position.set(5, 8, 10); scene.add(d1);
  const d2 = new THREE.DirectionalLight(0xffffff, 0.6); d2.position.set(-6, -3, -8); scene.add(d2);
  orbGroup = new THREE.Group(); scene.add(orbGroup);
  const ro = new ResizeObserver(() => {
    const w = host.clientWidth, h = host.clientHeight;
    renderer.setSize(w, h); labelRenderer.setSize(w, h);
    camera.aspect = w / h; camera.updateProjectionMatrix();
  });
  ro.observe(host);
}

function animate(now) {
  requestAnimationFrame(animate);
  if (state.playing) stepPlay(now);
  controls.update();
  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
}

// ------------------------------------------------------------------ molecule
const sphereGeo = new THREE.SphereGeometry(1, 32, 20);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 16, 1, true);
const bondMat = new THREE.MeshStandardMaterial({ color: 0x9aa0a6, roughness: 0.5 });
const partialMat = new THREE.MeshStandardMaterial({ color: 0x9aa0a6, roughness: 0.5, transparent: true, opacity: 0.55 });

function buildMolecule() {
  M.symbols.forEach((s, i) => {
    const e = el(s);
    const mesh = new THREE.Mesh(sphereGeo, new THREE.MeshStandardMaterial({ color: e.c, roughness: 0.4 }));
    mesh.scale.setScalar(e.r);
    scene.add(mesh); atomMeshes.push(mesh);
    const div = document.createElement('div');
    div.className = 'atom-label'; div.textContent = `${s}${i}`;
    const lab = new CSS2DObject(div); lab.center.set(-0.35, 1.1); lab.visible = false;
    mesh.add(lab); labelObjs.push(lab);
  });
  // The forming and breaking bonds are judged against their own equilibrium lengths, taken from
  // the path ends (forming bond: product end; breaking bond: reactant end).
  const ka = M.key_atoms;
  const dist = (f, a, b) => Math.hypot(...[0, 1, 2].map((k) => M.frames[f].xyz[a][k] - M.frames[f].xyz[b][k]));
  const eq = new Map([[pairKey(ka.nucleophile, ka.center), dist(N - 1, ka.nucleophile, ka.center)],
                      [pairKey(ka.center, ka.leaving_group), dist(0, ka.center, ka.leaving_group)]]);
  for (let i = 0; i < M.symbols.length; i++) {
    for (let j = i + 1; j < M.symbols.length; j++) {
      const solid = new THREE.Mesh(cylGeo, bondMat);
      const dashed = new THREE.Group();
      for (let k = 0; k < 5; k++) {
        const seg = new THREE.Mesh(cylGeo, partialMat);
        seg.scale.set(1, 0.11, 1); seg.position.y = -0.5 + (k + 0.5) / 5;
        dashed.add(seg);
      }
      scene.add(solid); scene.add(dashed);
      bonds.push({ i, j, solid, dashed, eq: eq.get(pairKey(i, j)),
                   cov: el(M.symbols[i]).cov + el(M.symbols[j]).cov });
    }
  }
}
const pairKey = (a, b) => (a < b ? `${a}-${b}` : `${b}-${a}`);

function framePositions(t) {
  const f0 = Math.floor(t), f1 = Math.min(N - 1, f0 + 1), w = t - f0;
  const a = M.frames[f0].xyz, b = M.frames[f1].xyz;
  return a.map((p, i) => new THREE.Vector3(
    p[0] + (b[i][0] - p[0]) * w, p[1] + (b[i][1] - p[1]) * w, p[2] + (b[i][2] - p[2]) * w));
}

const UP = new THREE.Vector3(0, 1, 0);
function placeBond(obj, a, b, radius) {
  const d = new THREE.Vector3().subVectors(b, a);
  const len = d.length();
  obj.position.copy(a).addScaledVector(d, 0.5);
  obj.quaternion.setFromUnitVectors(UP, d.clone().normalize());
  obj.scale.set(radius, len, radius);
}

function updateMolecule(pos) {
  pos.forEach((p, i) => atomMeshes[i].position.copy(p));
  for (const b of bonds) {
    const d = pos[b.i].distanceTo(pos[b.j]);
    // key bonds: full below 1.15 × equilibrium length, dashed (partial) up to 1.6 ×
    const full = b.eq ? d < 1.15 * b.eq : d < 1.15 * b.cov;
    const partial = !full && !!b.eq && d < 1.6 * b.eq;
    b.solid.visible = full || (partial && !state.dash);
    b.dashed.visible = partial && state.dash;
    if (b.solid.visible) {
      b.solid.material = full ? bondMat : partialMat;
      placeBond(b.solid, pos[b.i], pos[b.j], full ? 0.085 : 0.05);
    }
    if (b.dashed.visible) placeBond(b.dashed, pos[b.i], pos[b.j], 0.06);
  }
  labelObjs.forEach((l) => { l.visible = state.labels; });
}

// ------------------------------------------------------------------ orbitals
// one material pair per orbital: dark shade = + phase, light shade = − phase
const lobeMats = new Map();
function materialsFor(o) {
  if (!lobeMats.has(o.key)) {
    const mk = (c) => new THREE.MeshPhongMaterial({ color: new THREE.Color(c), transparent: true,
      opacity: 0.6, side: THREE.DoubleSide, depthWrite: false, shininess: 60 });
    const [p, n] = o.colors || ['#2f6db5', '#e08a1e'];
    lobeMats.set(o.key, [mk(p), mk(n)]);
  }
  return lobeMats.get(o.key);
}

function lobeGeometry(entry) {
  const [ov, nv] = entry.v, [of, nf] = entry.f;
  if (!nv) return null;
  const q = new Int16Array(BIN, ov, nv * 3);
  const p = new Float32Array(nv * 3);
  const { mid, scale } = M.quant;
  for (let i = 0; i < nv; i++) for (let k = 0; k < 3; k++) p[3 * i + k] = mid[k] + q[3 * i + k] * scale[k];
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(p, 3));
  g.setIndex(new THREE.BufferAttribute(new Uint16Array(BIN, of, nf * 3).slice(), 1));
  g.computeVertexNormals();
  return g;
}

function orbitalObject(o, f) {
  const ck = `${o.key}:${f}`;
  if (orbCache.has(ck)) return orbCache.get(ck);
  const e = M.mesh_index[o.key][f];
  const grp = new THREE.Group();
  const [posMat, negMat] = materialsFor(o);
  const gp = lobeGeometry(e.pos), gn = lobeGeometry(e.neg);
  if (gp) grp.add(new THREE.Mesh(gp, posMat));
  if (gn) grp.add(new THREE.Mesh(gn, negMat));
  orbCache.set(ck, grp);
  return grp;
}

function updateOrbitals(f) {
  orbGroup.clear();
  if (!state.orbOn) return;
  for (const o of M.orbitals) {
    if (o.kind !== state.kind || state.hidden.has(o.key)) continue;
    orbGroup.add(orbitalObject(o, f));
  }
}

function buildOrbitalList() {
  const box = $('orbList');
  box.innerHTML = '';
  for (const o of M.orbitals.filter((x) => x.kind === state.kind)) {
    const lab = document.createElement('label');
    const cb = document.createElement('input');
    cb.type = 'checkbox'; cb.checked = !state.hidden.has(o.key);
    cb.addEventListener('change', () => { cb.checked ? state.hidden.delete(o.key) : state.hidden.add(o.key); update(); });
    const [p, n] = o.colors || ['#2f6db5', '#e08a1e'];
    const sw = document.createElement('span');
    sw.className = 'pair';
    sw.title = '+ phase / − phase';
    sw.innerHTML = `<i style="background:${p}"></i><i style="background:${n}"></i>`;
    lab.append(cb, ' ', sw, ' ' + o.label);
    box.appendChild(lab);
  }
}

// ------------------------------------------------------------------ energy diagram
const svg = $('energy');
const NS = 'http://www.w3.org/2000/svg';
const D = { ml: 58, mr: 64, mt: 30, mb: 36 };
let xs, ys, sArr, eArr, marker, mline;

let curve;   // the energy curve currently shown (working level or single points on the same path)
function setupDiagram() {
  sArr = M.frames.map((f) => f.s);
  const curves = M.curves || [{ key: 'working', label: M.level, E_kJ: M.frames.map((f) => f.E_kJ),
                                E_separated_products_kJ: M.E_separated_products_kJ }];
  const sel = $('curveSel');
  for (const c of curves) sel.add(new Option(c.label, c.key));
  const pick = (key) => {
    curve = curves.find((c) => c.key === key) || curves[0];
    eArr = curve.E_kJ;
    drawDiagram(); update();
  };
  sel.addEventListener('change', () => pick(sel.value));
  curve = curves.find((c) => c.key === M.default_curve) || curves[0];
  eArr = curve.E_kJ;
  sel.value = curve.key;
  drawDiagram();
  new ResizeObserver(drawDiagram).observe(svg);
}

function niceStep(range, n = 6) {
  const raw = range / n, p = 10 ** Math.floor(Math.log10(raw));
  return [1, 2, 2.5, 5, 10].map((m) => m * p).find((v) => v >= raw);
}

function drawDiagram() {
  const W = svg.clientWidth || 800, H = svg.clientHeight || 230;
  const narrow = W < 560;           // phones: shorter labels, tighter margins
  D.ml = narrow ? 44 : 58; D.mr = narrow ? 40 : 64;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = '';
  const s0 = sArr[0], s1 = sArr[N - 1];
  const ep = curve.E_separated_products_kJ;
  const eMin = Math.min(...eArr, ep), eMax = Math.max(...eArr, 0);
  const pad = (eMax - eMin) * 0.12;
  const y0 = eMin - pad, y1 = eMax + pad * 1.8;
  xs = (s) => D.ml + (s - s0) / (s1 - s0) * (W - D.ml - D.mr);
  ys = (e) => H - D.mb - (e - y0) / (y1 - y0) * (H - D.mt - D.mb);
  const add = (tag, attrs, text) => {
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (text !== undefined) n.textContent = text;
    svg.appendChild(n); return n;
  };
  const step = niceStep(y1 - y0);
  for (let e = Math.ceil(y0 / step) * step; e <= y1; e += step) {
    add('line', { x1: D.ml, x2: W - D.mr, y1: ys(e), y2: ys(e), class: 'grid' });
    add('text', { x: D.ml - 6, y: ys(e) + 4, 'text-anchor': 'end' }, `${+e.toFixed(1)}`.replace('-', '−'));
  }
  // kcal/mol axis on the right with its own round numbers
  const kstep = niceStep((y1 - y0) / KCAL);
  for (let k = Math.ceil(y0 / KCAL / kstep) * kstep; k * KCAL <= y1; k += kstep) {
    add('line', { x1: W - D.mr, x2: W - D.mr + 4, y1: ys(k * KCAL), y2: ys(k * KCAL), class: 'tick' });
    add('text', { x: W - D.mr + 7, y: ys(k * KCAL) + 4 }, `${+k.toFixed(1)}`.replace('-', '−'));
  }
  add('line', { x1: W - D.mr, x2: W - D.mr, y1: D.mt, y2: H - D.mb, class: 'axis' });
  add('text', { x: 12, y: D.mt - 12 }, 'kJ/mol');
  add('text', { x: W - 4, y: D.mt - 12, 'text-anchor': 'end' }, narrow ? 'kcal' : '(kcal/mol)');
  add('line', { x1: D.ml, x2: W - D.mr, y1: H - D.mb, y2: H - D.mb, class: 'axis' });
  const sStep = niceStep(s1 - s0);
  for (let s = Math.ceil(s0 / sStep) * sStep; s <= s1; s += sStep) {
    add('line', { x1: xs(s), x2: xs(s), y1: H - D.mb, y2: H - D.mb + 4, class: 'tick' });
    add('text', { x: xs(s), y: H - D.mb + 16, 'text-anchor': 'middle' }, `${+s.toFixed(1)}`);
  }
  add('text', { x: (D.ml + W - D.mr) / 2, y: H - 4, 'text-anchor': 'middle' },
      'reaction coordinate s (bohr·amu½)');
  // reference levels; labels sit at the end of the plot away from the curve
  const awayFrom = (e) => (Math.abs(eArr[0] - e) > Math.abs(eArr[N - 1] - e)
    ? { x: D.ml + 4, a: 'start' } : { x: W - D.mr - 4, a: 'end' });
  add('line', { x1: D.ml, x2: W - D.mr, y1: ys(0), y2: ys(0), class: 'ref' });
  let p = awayFrom(0);
  add('text', { x: p.x, y: ys(0) - 4, 'text-anchor': p.a }, narrow ? 'reactants' : 'separated reactants, 0');
  add('line', { x1: D.ml, x2: W - D.mr, y1: ys(ep), y2: ys(ep), class: 'ref' });
  p = awayFrom(ep);
  add('text', { x: p.x, y: ys(ep) - 4, 'text-anchor': p.a },
      narrow ? 'products' : `separated products, ${ep.toFixed(0)}`.replace('-', '−'));
  add('path', { d: sArr.map((s, k) => `${k ? 'L' : 'M'}${xs(s).toFixed(1)},${ys(eArr[k]).toFixed(1)}`).join(''), class: 'curve' });
  sArr.forEach((s, k) => add('circle', { cx: xs(s), cy: ys(eArr[k]), r: 1.8, class: 'pt' }));
  // label the highest point of the curve shown (single points on the working-level path can
  // peak one or two frames away from the working-level TS)
  const ts = eArr.indexOf(Math.max(...eArr));
  const lab = (k, text, dy, anchor) => add('text', { x: xs(sArr[k]), y: ys(eArr[k]) + dy, 'text-anchor': anchor, class: 'lab' },
    (narrow ? `${text} ${eArr[k].toFixed(0)}` : `${text} ${eArr[k].toFixed(1)} (${(eArr[k] / KCAL).toFixed(1)})`).replaceAll('-', '−'));
  lab(ts, ts === M.ts_frame ? 'TS' : narrow ? 'max' : 'highest point', -10, 'middle');
  // the path ends are minima only at the working level; at other levels call them path ends
  const working = curve.key === 'working';
  lab(0, narrow ? 'RC' : working ? 'reactant complex' : 'reactant end of path', 18, 'start');
  lab(N - 1, narrow ? 'PC' : working ? 'product complex' : 'product end of path', -8, 'end');
  mline = add('line', { class: 'mline' });
  marker = add('circle', { r: 6.5, class: 'marker' });
  placeMarker();
}

function interp(arr, t) {
  const f0 = Math.floor(t), f1 = Math.min(N - 1, f0 + 1), w = t - f0;
  return arr[f0] + (arr[f1] - arr[f0]) * w;
}

function placeMarker() {
  if (!marker) return;
  const x = xs(interp(sArr, state.t)), y = ys(interp(eArr, state.t));
  marker.setAttribute('cx', x); marker.setAttribute('cy', y);
  const H = svg.viewBox.baseVal.height;
  mline.setAttribute('x1', x); mline.setAttribute('x2', x);
  mline.setAttribute('y1', D.mt - 6); mline.setAttribute('y2', H - D.mb);
}

function tFromClientX(cx) {
  const r = svg.getBoundingClientRect();
  const W = svg.viewBox.baseVal.width;
  const x = (cx - r.left) * W / r.width;
  return tFromS(sArr[0] + (x - D.ml) / (W - D.ml - D.mr) * (sArr[N - 1] - sArr[0]));
}

// ------------------------------------------------------------------ controls
function setT(t) { state.t = Math.max(0, Math.min(N - 1, t)); update(); }

function update() {
  const pos = framePositions(state.t);
  updateMolecule(pos);
  const f = Math.round(state.t);
  updateOrbitals(f);
  placeMarker();
  const ka = M.key_atoms, sy = M.symbols;
  const dn = pos[ka.nucleophile].distanceTo(pos[ka.center]);
  const dl = pos[ka.center].distanceTo(pos[ka.leaving_group]);
  const E = interp(eArr, state.t);
  const nu = `${sy[ka.center]}–${sy[ka.nucleophile]}`, lg = `${sy[ka.center]}–${sy[ka.leaving_group]}`;
  $('readout').textContent =
    `frame ${f} of ${N - 1}   s = ${interp(sArr, state.t).toFixed(2)} bohr·amu½\n` +
    `E = ${E.toFixed(1)} kJ/mol (${(E / KCAL).toFixed(1)} kcal/mol)\n` +
    `${nu} ${dn.toFixed(2)} Å   ${lg} ${dl.toFixed(2)} Å`;
  svg.setAttribute('aria-valuemax', N - 1);
  svg.setAttribute('aria-valuenow', f);
  svg.setAttribute('aria-valuetext',
    `frame ${f} of ${N - 1}, energy ${E.toFixed(1)} kilojoules per mole, ${nu} ${dn.toFixed(2)} angstrom, ${lg} ${dl.toFixed(2)} angstrom`);
}

let lastNow = null;
function stepPlay(now) {
  if (lastNow === null) { lastNow = now; return; }
  const dt = (now - lastNow) / 1000; lastNow = now;
  // constant speed along the path length s (IRC frames bunch up where the path flattens)
  const s = interp(sArr, state.t) + dt * (sArr[N - 1] - sArr[0]) / 12;   // whole path in about 12 s
  if (s >= sArr[N - 1]) { setT(N - 1); togglePlay(false); return; }
  setT(tFromS(s));
}

function tFromS(s) {
  if (s <= sArr[0]) return 0;
  if (s >= sArr[N - 1]) return N - 1;
  let k = 0; while (k < N - 2 && sArr[k + 1] < s) k++;
  return k + (s - sArr[k]) / (sArr[k + 1] - sArr[k]);
}

function togglePlay(on = !state.playing) {
  state.playing = on; lastNow = null;
  if (on && state.t >= N - 1) state.t = 0;
  $('play').textContent = on ? 'Pause' : 'Play';
  $('play').setAttribute('aria-pressed', String(on));
}

function setupControls() {
  $('orbOn').addEventListener('change', (e) => { state.orbOn = e.target.checked; update(); });
  $('labelsOn').addEventListener('change', (e) => { state.labels = e.target.checked; update(); });
  $('dashOn').addEventListener('change', (e) => { state.dash = e.target.checked; update(); });
  document.querySelectorAll('input[name=kind]').forEach((r) => r.addEventListener('change', (e) => {
    state.kind = e.target.value; buildOrbitalList(); update();
  }));
  document.querySelectorAll('[data-cam]').forEach((b) => b.addEventListener('click', () => setCamera(b.dataset.cam)));
  $('play').addEventListener('click', () => togglePlay());
  $('toR').addEventListener('click', () => { togglePlay(false); setT(0); });
  $('toTS').addEventListener('click', () => { togglePlay(false); setT(M.ts_frame); });
  $('toP').addEventListener('click', () => { togglePlay(false); setT(N - 1); });
  let dragging = false;
  svg.addEventListener('pointerdown', (e) => { dragging = true; svg.setPointerCapture(e.pointerId); togglePlay(false); setT(tFromClientX(e.clientX)); svg.focus(); });
  svg.addEventListener('pointermove', (e) => { if (dragging) setT(tFromClientX(e.clientX)); });
  svg.addEventListener('pointerup', () => { dragging = false; });
  svg.addEventListener('keydown', (e) => {
    const big = e.shiftKey ? 5 : 1;
    const k = e.key;
    if (k === 'ArrowRight' || k === 'ArrowUp') setT(Math.floor(state.t + 1e-6) + big);
    else if (k === 'ArrowLeft' || k === 'ArrowDown') setT(Math.ceil(state.t - 1e-6) - big);
    else if (k === 'Home') setT(0);
    else if (k === 'End') setT(N - 1);
    else if (k === ' ') togglePlay();
    else return;
    if (k !== ' ') togglePlay(false);
    e.preventDefault();
  });
}

function setCamera(mode) {
  const ka = M.key_atoms, ts = M.frames[M.ts_frame].xyz;
  const v = (p) => new THREE.Vector3(...p);
  const nu = v(ts[ka.nucleophile]), lg = v(ts[ka.leaving_group]), c = v(ts[ka.center]);
  const axis = new THREE.Vector3().subVectors(lg, nu).normalize();
  let side = new THREE.Vector3(0, 0, 1).cross(axis);
  if (side.lengthSq() < 1e-3) side = new THREE.Vector3(1, 0, 0).cross(axis);
  side.normalize();
  const third = new THREE.Vector3().crossVectors(axis, side).normalize();
  let r = 0;
  for (const f of [M.frames[0], M.frames[N - 1], M.frames[M.ts_frame]])
    for (const p of f.xyz) r = Math.max(r, v(p).distanceTo(c));
  const dist = (r + 1.6) / Math.sin(THREE.MathUtils.degToRad(camera.fov / 2)) * 0.9;
  let dir, up;
  if (mode === 'axis') { dir = axis.clone().negate(); up = third; }
  else if (mode === 'top') { dir = third.clone(); up = axis.clone(); }
  else { dir = side.clone(); up = third.clone().negate(); }   // nucleophile left, leaving group right
  camera.position.copy(c).addScaledVector(dir, dist);
  camera.up.copy(up);
  controls.target.copy(c);
  camera.lookAt(c);
  controls.update();
}
