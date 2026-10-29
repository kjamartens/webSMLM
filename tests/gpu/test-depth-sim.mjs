#!/usr/bin/env node
// Correctness gate for the depth-dependent simulation (Li et al., Biomed. Opt. Express 10, 2708
// (2019): beads on the coverslip and in a gel, clathrin-coated pits, Nup107 pores at depth).
// Every number the documentation quotes for this feature comes from a key written here through
// writeResults('depth-sim', …), so it traces to a JSON value (SMLMpipelinelessons §9).
//
//  (0) byte identity: seeded movies with every new setting at its default hash to the values
//      pinned below — the build before the depth work (2026-09-28u). A changed hash means a
//      default movie changed; never update the constants without knowing why.
//  (1) pupil: supercritical annulus and Fresnel apodisation reach the PSF (on/off kernels and
//      z-CRLBs differ), the annulus power fraction on the grid matches a 1-D quadrature and the
//      published numbers, the mismatch phase projects onto the report's spherical Zernikes.
//  (2) focus map: f0(z_p), its slope and the model's own RIMF.
//  (3) physical depth mode: emitters at the sample depth reproduce the objective-mode movie;
//      CPU workers = GPU.
//  (4) gel bead stack: truth, skipped frames, glass-bead slope 1.
//  (5) structures: ccp and nup107 geometry.
//  (6) per-object scores move in the published direction.
//  (7) presets: every one simulates, survives a settings round trip, and an edit flips it.
//
// Run: node tests/gpu/test-depth-sim.mjs   (Playwright via tools/package.json, see tests/README.md)
import { launchPage } from '../lib/launch.mjs';
import { writeResults } from '../lib/report.mjs';

// Pinned on build 2026-09-28u (before any depth work), CPU path, harness settings in SIM_BASE.
const BASELINE = { uniform3D: 2860914707, nup: 555640748 };

const { browser, page } = await launchPage({ headless: false });
let failed = false;
const results = {};
const check = (label, ok, detail) => {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}${detail ? '  — ' + detail : ''}`);
  if (!ok) failed = true;
};

// In-page harness: set a field (or its paramOverrides entry when it has no control), simulate,
// hash. The same body is what the interactive verification ran in the Browser pane.
await page.evaluate(() => {
  window.__set = (id, v) => {
    const el = document.getElementById(id);
    if (!el) { paramOverrides[id] = v; return; }
    if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
    el.dispatchEvent(new Event('change'));
  };
  window.__simBase = () => {
    __set('useGpu', false); paramOverrides.simulation_3d = true;
    __set('simulation_seed', 31); __set('simulation_fov', 64); __set('frames', 50); __set('simulation_zRange', 300);
    __set('simulation_psfZRange', 400); __set('simulation_psfZStep', 20);
  };
  window.__simHash = async (cfg) => {
    __simBase();
    for (const [k, v] of Object.entries(cfg || {})) __set(k, v);
    const st = await generateSynthetic();
    return { hash: await stackPixelHash(st), n: st.n, events: groundTruthEvents.length, path: lastSimTimings.path };
  };
});

try {
  // ---- (0) byte identity ---------------------------------------------------------------------
  const base = await page.evaluate(async () => ({
    uniform3D: await __simHash({ simulation_structureType: 'uniform3D' }),
    nup: await __simHash({ simulation_structureType: 'nup' }),
  }));
  results.baseline = base;
  for (const k of Object.keys(BASELINE))
    check(`(0) seeded ${k} movie byte-identical to 2026-09-28u`, base[k].hash === BASELINE[k], `hash ${base[k].hash} vs ${BASELINE[k]}`);

  // ---- (1) pupil -------------------------------------------------------------------------------
  const p1 = await page.evaluate(async () => {
    const out = {};
    const li = { NA: 1.43, lambda: 670e-9, ns: 1.3345, ni: 1.518 };
    // 1-D quadrature reference, uniform amplitude and with Fresnel apodisation
    for (const apod of [false, true]) {
      const p = { ...li, saf: true, apod };
      const P0 = psfPupilPower(p, 0).total;
      out[apod ? 'quadApod' : 'quad'] = [0, 50, 100, 340, 2000].map(z => {
        const r = psfPupilPower(p, z * 1e-9);
        return { z, safPct: 100 * r.saf / r.total, collected: r.total / P0 };
      });
    }
    // the same fraction on the Cartesian pupil grid the PSF is built from
    out.grid = [0, 50, 100, 340, 2000].map(z => {
      const cfg = { ...li, ti0: 150e-6, particleAxialPosition: z * 1e-9, zStep: 10e-9, zernikeCoeffs: new Array(PSF_NZERNIKE).fill(0),
                    maskType: 'none', maskModes: 5, maskWaist: 1, saf: true, apod: false };
      const M = psfFftM(cfg), dk = psfFftDkForCfg(cfg);
      const plan = psfPupilPlan(cfg, M, dk), k0 = 2 * Math.PI / cfg.lambda, c0 = Math.floor(M / 2), rc = psfBMax(cfg);
      let tot = 0, saf = 0;
      for (let q = 0; q < plan.idx.length; q++) {
        const i = plan.idx[q], kx = ((i % M) - c0) * dk, ky = (Math.floor(i / M) - c0) * dk;
        const rho = Math.hypot(kx, ky) / (k0 * cfg.NA), a2 = plan.amp[q] * plan.amp[q];
        tot += a2; if (rho > rc) saf += a2;
      }
      return { z, safPct: 100 * saf / tot, M };
    });
    // penetration depth of the intensity at the pupil edge
    out.penetrationNm = 1e9 / (2 * (2 * Math.PI / li.lambda) * Math.sqrt(li.NA * li.NA - li.ns * li.ns));
    // Projection of the index-mismatch phase onto ORTHONORMAL radial Zernikes over the
    // subcritical disc r = rho/rho_c, at the focus where the defocus term vanishes (report §1).
    {
      const k0 = 2 * Math.PI / li.lambda, rc = li.ns / li.NA, N = 4000, zp = 1e-6;
      const R = {
        2: r => Math.sqrt(3) * (2 * r * r - 1),
        4: r => Math.sqrt(5) * (6 * r ** 4 - 6 * r * r + 1),
        6: r => Math.sqrt(7) * (20 * r ** 6 - 30 * r ** 4 + 12 * r * r - 1),
        8: r => 3 * (70 * r ** 8 - 140 * r ** 6 + 90 * r ** 4 - 20 * r * r + 1),
      };
      const proj = (f, Z) => { let s = 0; for (let i = 0; i < N; i++) { const r = (i + 0.5) / N; s += f(r) * Z(r) * r; } return 2 * s / N; };
      const sample = r => Math.sqrt(li.ns ** 2 - (li.NA * rc * r) ** 2), imm = r => Math.sqrt(li.ni ** 2 - (li.NA * rc * r) ** 2);
      const c2s = proj(sample, R[2]), c2i = proj(imm, R[2]);
      const f0 = zp * c2s / c2i;                      // objective travel that cancels defocus
      const phi = r => k0 * (zp * sample(r) - f0 * imm(r));
      out.zernikeRadPerUm = { Z4: proj(phi, R[4]), Z6: proj(phi, R[6]), Z8: proj(phi, R[8]), rimfFullPupil: zp / f0 };
    }
    return out;
  });
  results.pupil = p1;
  const pub = [12.87, 9.72, 7.37, 2.22, 0.05], quadTarget = [12.91, 9.77, 7.42, 2.27, 0.08];
  p1.quad.forEach((r, i) => check(`(1) SAF fraction ${r.z} nm (quadrature)`, Math.abs(r.safPct - quadTarget[i]) < 0.05,
    `${r.safPct.toFixed(2)} % (published ${pub[i]} %)`));
  p1.grid.forEach((r, i) => check(`(1) SAF fraction ${r.z} nm (M=${r.M} grid vs quadrature)`, Math.abs(r.safPct - p1.quad[i].safPct) < 0.3,
    `${r.safPct.toFixed(2)} % vs ${p1.quad[i].safPct.toFixed(2)} %`));
  [36.3, 29.7, 24.1, 9.0, 0.4].forEach((t, i) => check(`(1) SAF fraction with apodisation ${p1.quadApod[i].z} nm`,
    Math.abs(p1.quadApod[i].safPct - t) < 0.1, `${p1.quadApod[i].safPct.toFixed(2)} %`));
  [[50, 0.965], [340, 0.891], [2000, 0.872]].forEach(([z, t]) => { const r = p1.quad.find(q => q.z === z);
    check(`(1) collected fraction ${z} nm`, Math.abs(r.collected - t) < 0.002, r.collected.toFixed(3)); });
  check('(1) evanescent 1/e penetration at the pupil edge', Math.abs(p1.penetrationNm - 103.8) < 0.1, `${p1.penetrationNm.toFixed(1)} nm`);
  const zr = p1.zernikeRadPerUm;
  check('(1) mismatch phase → Z4⁰ (report −0.2765 rad/µm)', Math.abs(zr.Z4 + 0.2765) < 0.01, zr.Z4.toFixed(4));
  check('(1) mismatch phase → Z6⁰ (report −0.1651 rad/µm)', Math.abs(zr.Z6 + 0.1651) < 0.01, zr.Z6.toFixed(4));
  check('(1) mismatch phase → Z8⁰ (report −0.0986 rad/µm)', Math.abs(zr.Z8 + 0.0986) < 0.01, zr.Z8.toFixed(4));
  check('(1) full-subcritical-pupil RIMF (report 0.7255)', Math.abs(zr.rimfFullPupil - 0.7255) < 0.002, zr.rimfFullPupil.toFixed(4));

  // The switches reach the PSF on every path (the 2026-09-20g lesson: identical output for
  // different settings is the alarm). Same small kernel built with each flag, CPU pool and GPU.
  const reach = await page.evaluate(async () => {
    __simBase(); __set('simulation_psfNA', 1.43); __set('simulation_psfNs', 1.3345); __set('simulation_psfDepth', 100);
    const sig = async () => { const cfg = readPsfConfigFromUI(); const k = await getOrBuildPsfKernelStack(cfg);
      const s = k.slices[k.focusIndex]; let h = 0; for (let i = 0; i < s.length; i += 97) h = (h * 31 + Math.round(s[i] * 1e9)) % 2147483647;
      return { h, zCrlb: k.zCrlb ? k.zCrlb.atFocusZSig : null, path: lastPsfBuildTimings.path }; };
    const out = {};
    for (const gpu of [false, true]) {
      __set('useGpu', gpu);
      __set('simulation_psfSaf', false); __set('simulation_psfApodization', false); const off = await sig();
      __set('simulation_psfSaf', true); const saf = await sig();
      __set('simulation_psfApodization', true); const both = await sig();
      out[gpu ? 'gpu' : 'cpu'] = { off, saf, both };
    }
    __set('simulation_psfSaf', false); __set('simulation_psfApodization', false); __set('useGpu', false);
    __set('simulation_psfNA', 1.4); __set('simulation_psfNs', 1.33); __set('simulation_psfDepth', 0);
    return out;
  });
  results.reach = reach;
  for (const k of ['cpu', 'gpu']) {
    const r = reach[k];
    check(`(1) SAF reaches the kernel (${k}, ${r.saf.path})`, r.off.h !== r.saf.h && r.off.zCrlb !== r.saf.zCrlb,
      `z-CRLB at focus off ${r.off.zCrlb?.toFixed(2)} / on ${r.saf.zCrlb?.toFixed(2)} nm`);
    check(`(1) apodisation reaches the kernel (${k})`, r.saf.h !== r.both.h, `z-CRLB ${r.both.zCrlb?.toFixed(2)} nm`);
  }

  // Grid convergence (SMLMpipelinelessons §1.5): one plane at the emitter's paraxial focus, the
  // Li 2019 optics with the annulus on, on pupil grids of 64…512. Relative L2 distance of the
  // photon-normalised plane from the next-finer grid, and the share of the plane's light in a
  // 10-pixel border of the 6 µm kernel (light the kernel would cut). The production grid is
  // psfFftM(): 256 with the annulus on.
  const conv = await page.evaluate(async () => {
    const out = [];
    for (const zNm of [0, 2000, 5000]) {
      const cfg = { NA: 1.43, lambda: 670e-9, ns: 1.3345, ni: 1.518, ti0: 150e-6, particleAxialPosition: zNm * 1e-9, zStep: 10e-9,
                    zernikeCoeffs: new Array(PSF_NZERNIKE).fill(0), maskType: 'none', maskModes: 5, maskWaist: 1, saf: true, apod: false };
      cfg.zernikeCoeffs[5] = 0.183;
      const res = 25e-9, nx = 241, planes = {};
      for (const M of [64, 128, 256, 512]) {
        const kMax = 2 * Math.PI / cfg.lambda * cfg.NA * psfApertureMax(cfg), dk = 2 * kMax / (M - 4);
        const cart = computePsfPupilCartesianForZPlane(0, 1, cfg, M, dk);
        const s = computePsfIntensityPlaneFFT(nx, nx, res, cart, M, dk);
        let sum = 0; for (const v of s) sum += v; for (let i = 0; i < s.length; i++) s[i] /= sum;
        planes[M] = s;
      }
      const l2 = (a, b) => { let d = 0, n = 0; for (let i = 0; i < a.length; i++) { d += (a[i] - b[i]) ** 2; n += b[i] ** 2; } return Math.sqrt(d / n); };
      let border = 0; const s = planes[256];
      for (let y = 0; y < nx; y++) for (let x = 0; x < nx; x++) if (x < 10 || y < 10 || x >= nx - 10 || y >= nx - 10) border += s[y * nx + x];
      out.push({ zNm, l2_64_128: l2(planes[64], planes[128]), l2_128_256: l2(planes[128], planes[256]), l2_256_512: l2(planes[256], planes[512]), borderShare256: border });
    }
    return out;
  });
  results.convergence = conv;
  for (const c of conv)
    check(`(1) pupil grid converged at ${c.zNm} nm (256 vs 512 < 1 %)`, c.l2_256_512 < 0.01,
      `L2 64→128 ${(100 * c.l2_64_128).toFixed(2)} %, 128→256 ${(100 * c.l2_128_256).toFixed(2)} %, 256→512 ${(100 * c.l2_256_512).toFixed(3)} %, kernel border ${(100 * c.borderShare256).toFixed(2)} %`);

  // ---- (2) focus map: Li 2019 optics, astigmatism 183 mλ at OSA j=5, annulus on ----------------
  const fmap = await page.evaluate(async () => {
    const z = new Array(PSF_NZERNIKE).fill(0); z[5] = 183;
    __set('simulation_psfNA', 1.43); __set('simulation_psfWavelength', 670); __set('simulation_psfNs', 1.3345);
    __set('simulation_psfNi', 1.518); __set('simulation_psfSaf', true); __set('simulation_psfZernikeCustom', z.join(','));
    const fm = await psfFocusMap(readPsfConfigFromUI(), 0, 5000);
    return { j5: psfIndexToNM(5), anchors: fm.anchors.map(a => ({ d: a.dNm, f0: a.f0Nm, how: a.how, atEdge: a.atEdge })),
             at: [340, 2000, 3500, 5000].map(d => ({ d, f0: fm.f0(d), rimf: fm.rimf(d), slope: fm.f0prime(d) })) };
  });
  results.focusMap = fmap;
  check('(2) OSA j=5 is vertical astigmatism (n=2, m=+2)', fmap.j5[0] === 2 && fmap.j5[1] === 2, JSON.stringify(fmap.j5));
  check('(2) every focus anchor found by the width zero, none at the sweep edge',
    fmap.anchors.every(a => a.how === 'width' && !a.atEdge), `${fmap.anchors.length} anchors`);
  check('(2) f0 at the coverslip is 0', Math.abs(fmap.anchors[0].f0) < 2, `${fmap.anchors[0].f0.toFixed(1)} nm`);
  const at = d => fmap.at.find(a => a.d === d);
  check('(2) f0(2000) within [paraxial 2283, full-pupil 2760] nm (report MLE 2619)', at(2000).f0 > 2283 && at(2000).f0 < 2760, `${at(2000).f0.toFixed(0)} nm`);
  check('(2) f0(5000) within [5688, 6892] nm', at(5000).f0 > 5688 && at(5000).f0 < 6892, `${at(5000).f0.toFixed(0)} nm`);
  check('(2) model RIMF at 340/2000/3500 nm near report 0.70/0.76/0.79',
    [[340, 0.70], [2000, 0.76], [3500, 0.79]].every(([d, t]) => Math.abs(at(d).rimf - t) < 0.02),
    [340, 2000, 3500].map(d => at(d).rimf.toFixed(3)).join(' / '));

  // ---- (3) physical depth mode -----------------------------------------------------------------
  // (a) the default movies are untouched by the new code paths
  const base3 = await page.evaluate(async () => {
    __set('simulation_psfSaf', false); __set('simulation_psfZernikeCustom', ''); __set('simulation_psfNA', 1.4);
    __set('simulation_psfWavelength', 660); __set('simulation_psfNs', 1.33); __set('simulation_psfNi', 1.518);
    paramOverrides.simulation_depthMode = 'objective';
    return { uniform3D: await __simHash({ simulation_structureType: 'uniform3D' }), nup: await __simHash({ simulation_structureType: 'nup' }) };
  });
  for (const k of Object.keys(BASELINE))
    check(`(3) seeded ${k} movie still byte-identical after the depth code`, base3[k].hash === BASELINE[k], `hash ${base3[k].hash}`);
  // (b) at the coverslip, a flat sample in physical mode is the objective-mode movie, bit for bit
  const same = await page.evaluate(async () => {
    const cfg = { simulation_structureType: 'tiltedPlane', simulation_zRange: 0, simulation_psfSaf: true, simulation_psfDepth: 0 };
    paramOverrides.simulation_depthMode = 'objective'; const a = await __simHash(cfg);
    paramOverrides.simulation_depthMode = 'physical'; const b = await __simHash(cfg);
    paramOverrides.simulation_depthMode = 'objective';
    return { a, b };
  });
  results.physicalAtCoverslip = same;
  check('(3) physical = objective for a flat sample at the coverslip', same.a.hash === same.b.hash, `${same.a.hash} / ${same.b.hash}`);
  // (c) residual mapping: the plane an emitter at depth D+delta is drawn from against the PSF built
  // directly at that depth and the same objective position — and against the opposite sign.
  const resid = await page.evaluate(async () => {
    const z = new Array(PSF_NZERNIKE).fill(0); z[5] = 183;
    __set('simulation_psfNA', 1.43); __set('simulation_psfWavelength', 670); __set('simulation_psfNs', 1.3345);
    __set('simulation_psfNi', 1.518); __set('simulation_psfSaf', true); __set('simulation_psfZernikeCustom', z.join(','));
    __set('simulation_psfZRange', 400); __set('simulation_psfZStep', 10);
    paramOverrides.simulation_depthMode = 'physical';
    const cfg = readPsfConfigFromUI(), D = 2000;
    const fObj = (await psfFocusMap(cfg, D, D)).f0(D);
    const deltas = [-120, -60, 60, 120], zps = deltas.map(d => D + d);
    const ds = await buildPhysicalDepthStack(cfg, D, fObj, zps);
    const corr = (a, b) => { let sa = 0, sb = 0, sab = 0, saa = 0, sbb = 0, n = a.length;
      for (let i = 0; i < n; i++) { sa += a[i]; sb += b[i]; } const ma = sa / n, mb = sb / n;
      for (let i = 0; i < n; i++) { const x = a[i] - ma, y = b[i] - mb; sab += x * y; saa += x * x; sbb += y * y; } return sab / Math.sqrt(saa * sbb); };
    const out = [];
    for (const zp of zps) {
      const direct = await buildPsfPlanesQuiet({ ...cfg, particleAxialPosition: zp * 1e-9, nz: 1, sweepCenterM: cfg.ti0 - fObj * 1e-9, fftDk: psfFftDkForCfg(cfg) }, ds.nx, ds.nx);
      const used = ds.slices[ds.planeOf(zp)];
      // the opposite sign: mirror the plane index about the bin's own centre plane
      const B = ds.bins.find(b => ds.planeOf(zp) >= b.offset && ds.planeOf(zp) < b.offset + b.n);
      const mirror = ds.slices[B.offset + B.half - (ds.planeOf(zp) - B.offset - B.half)];
      const centre = ds.slices[B.offset + B.half];
      out.push({ zp, corrUsed: corr(direct.slices[0], used), corrMirror: corr(direct.slices[0], mirror), corrNoResidual: corr(direct.slices[0], centre) });
    }
    paramOverrides.simulation_depthMode = 'objective';
    return { bins: ds.bins.map(b => ({ zb: b.zb, n: b.n })), out };
  });
  results.residualMapping = resid;
  for (const r of resid.out)
    check(`(3) residual mapping at ${r.zp} nm reproduces the direct PSF, and beats the opposite sign`,
      r.corrUsed > 0.995 && r.corrUsed > r.corrMirror && r.corrUsed >= r.corrNoResidual,
      `corr used ${r.corrUsed.toFixed(4)}, opposite sign ${r.corrMirror.toFixed(4)}, no residual ${r.corrNoResidual.toFixed(4)}`);
  // (d) a real physical-mode movie, 2 um deep, ±500 nm of structure, on the CPU workers and the GPU
  const phys = await page.evaluate(async () => {
    const cfg = { simulation_structureType: 'uniform3D', simulation_zRange: 500, simulation_psfDepth: 2000, simulation_psfZRange: 1000 };
    paramOverrides.simulation_depthMode = 'physical';
    const movies = [];
    for (const g of [false, true]) {
      __simBase(); for (const [k, v] of Object.entries(cfg)) __set(k, v); __set('useGpu', g);
      const st = await generateSynthetic();
      movies.push({ path: lastSimTimings.path, frames: await st.getFrames(0, st.n),
        depths: groundTruthEvents.map(e => e.depth), zs: groundTruthEvents.map(e => e.z), rates: groundTruthEvents.map(e => e.rate) });
    }
    paramOverrides.simulation_depthMode = 'objective'; __set('useGpu', false);
    let close = 0, tot = 0, maxAbs = 0;
    for (let fi = 0; fi < movies[0].frames.length; fi++) for (let i = 0; i < movies[0].frames[fi].length; i++) {
      const d = Math.abs(movies[0].frames[fi][i] - movies[1].frames[fi][i]); tot++; if (d <= 1e-3) close++; if (d > maxAbs) maxAbs = d; }
    const m = movies[0];
    const depthOk = m.depths.every((d, i) => Math.abs(d - Math.max(0, 2000 + m.zs[i])) < 1e-9);
    return { paths: movies.map(x => x.path).join('/'), fracClose: close / tot, maxAbs, depthOk, nEvents: m.zs.length,
             depthRange: [Math.min(...m.depths), Math.max(...m.depths)] };
  });
  results.physicalMovie = phys;
  check('(3) physical movie: every event carries depth = 2000 + z', phys.depthOk, `${phys.nEvents} events, depth ${phys.depthRange.map(v => v.toFixed(0)).join('…')} nm`);
  check('(3) physical movie: CPU workers = GPU', phys.paths === 'cpu/gpu' && phys.fracClose >= 0.999,
    `paths ${phys.paths}, ${(100 * phys.fracClose).toFixed(3)}% of pixels within 1e-3 ADU, max|Δ| ${phys.maxAbs.toExponential(2)}`);
} finally {
  writeResults('depth-sim', results);
  await browser.close();
}
if (failed) { console.log('\nSOME CHECKS FAILED'); process.exitCode = 1; } else console.log('\nALL CHECKS PASSED');
