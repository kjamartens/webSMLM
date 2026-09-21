#!/usr/bin/env node
// Correctness gate for the simulator's summed-kernel splat, its counter-based camera-noise RNG, and
// (when WebGPU is available) the GPU frame stage and GPU PSF build that mirror them.
// Deliberately SMALL (64 px frames, a 41-plane PSF, a few dozen emitters) so it runs in seconds on a
// laptop; bench-simulation.mjs is where the sizes that matter for speed live.
//
//  (a) pcg4d: the WGSL port produces the same u32 stream as pcg4dInto(), bit for bit.
//  (b) splatZernikeEmitter() (block sums, one interpolation per camera pixel) against the OLD
//      algorithm (interpolate every oversample^2 sub-cell, then sum), which is injected below as a
//      reference — it no longer exists in webSMLM.html. All four placement modes, emitters inside,
//      straddling and outside the frame edge, in-focus and defocused planes.
//  (c) GPU splat against CPU splat, noise-free.
//  (d) GPU camera noise against CPU camera noise on the same expectation image (sCMOS and EMCCD),
//      plus the noise statistics themselves.
//  (e) GPU PSF planes (direct evaluator) against CPU direct planes.
import assert from 'node:assert/strict';
import { launchPage, checkGpu } from '../lib/launch.mjs';

const { browser, page } = await launchPage({ headless: true });
let failed = false;
const check = (label, ok, detail) => {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}${detail ? '  — ' + detail : ''}`);
  if (!ok) failed = true;
};
try {
  // Small, fast PSF: 41 planes (±400 nm at 20 nm), default kernel width and oversample.
  await page.evaluate(() => {
    const set = (id, v) => { const el = document.getElementById(id); el.value = v; el.dispatchEvent(new Event('change')); };
    set('simulation_psfZRange', 400); set('simulation_psfZStep', 20); set('simulation_psfZernikePreset', 'astigModerate');
  });

  // ---- (b) summed-kernel splat vs the old per-sub-cell algorithm --------------------------------
  const b = await page.evaluate(async () => {
    // ---- reference: the pre-2026-09-21b splat, verbatim apart from the ref_ prefixes ----
    function ref_psfKernelAt(kernel, nx, ny, ix, iy){ if(ix<0||iy<0||ix>=nx||iy>=ny) return 0; return kernel[iy*nx+ix]; }
    function ref_nearest(kernel, nx, ny, x, y){ return ref_psfKernelAt(kernel, nx, ny, Math.round(x), Math.round(y)); }
    function ref_bilinear(kernel, nx, ny, x, y){
      const x0=Math.floor(x), y0=Math.floor(y), fx=x-x0, fy=y-y0;
      const v00=ref_psfKernelAt(kernel,nx,ny,x0,y0),   v10=ref_psfKernelAt(kernel,nx,ny,x0+1,y0);
      const v01=ref_psfKernelAt(kernel,nx,ny,x0,y0+1), v11=ref_psfKernelAt(kernel,nx,ny,x0+1,y0+1);
      return (v00*(1-fx)+v10*fx)*(1-fy) + (v01*(1-fx)+v11*fx)*fy;
    }
    function ref_bicubic(kernel, nx, ny, x, y){
      const x0=Math.floor(x), y0=Math.floor(y), wx=catmullRomWeights(x-x0), wy=catmullRomWeights(y-y0);
      let acc=0;
      for(let j=-1;j<=2;j++){ let row=0; for(let i=-1;i<=2;i++) row+=wx[i+1]*ref_psfKernelAt(kernel,nx,ny,x0+i,y0+j); acc+=wy[j+1]*row; }
      return acc;
    }
    function ref_sampleKernelAt(mode, kernel, nx, ny, x, y){
      if(mode==='nearest') return ref_nearest(kernel,nx,ny,x,y);
      if(mode==='cubic')   return ref_bicubic(kernel,nx,ny,x,y);
      return ref_bilinear(kernel,nx,ny,x,y);
    }
    function ref_splatFFT(img, w, h, mx, my, photons, kernelObj, camRad){
      const {kernel, nx, ny, oversample}=kernelObj;
      const kcx=(nx-1)/2, kcy=(ny-1)/2, X0=Math.round(mx), Y0=Math.round(my);
      const refX=kcx+(X0-mx-0.5)*oversample+0.5, refY=kcy+(Y0-my-0.5)*oversample+0.5;
      const shiftX=refX-Math.round(refX), shiftY=refY-Math.round(refY);
      const shifted=fftShiftKernelTile(kernel, nx, ny, shiftX, shiftY);
      for(let dy=-camRad; dy<=camRad; dy++){ const Y=Y0+dy; if(Y<0||Y>=h) continue; const baseY=kcy+(Y-my-0.5)*oversample;
        for(let dx=-camRad; dx<=camRad; dx++){ const X=X0+dx; if(X<0||X>=w) continue; const baseX=kcx+(X-mx-0.5)*oversample;
          let sum=0;
          for(let sy=0; sy<oversample; sy++){ const iy=Math.round(baseY+sy+0.5-shiftY);
            for(let sx=0; sx<oversample; sx++){ const ix=Math.round(baseX+sx+0.5-shiftX); sum+=ref_psfKernelAt(shifted, nx, ny, ix, iy); } }
          img[Y*w+X]+=photons*sum; } }
    }
    function ref_splat(img, w, h, mx, my, photons, kernelObj, camRad, interpMode){
      if(interpMode==='fft') return ref_splatFFT(img, w, h, mx, my, photons, kernelObj, camRad);
      const {kernel, nx, ny, oversample}=kernelObj;
      const kcx=(nx-1)/2, kcy=(ny-1)/2, X0=Math.round(mx), Y0=Math.round(my);
      for(let dy=-camRad; dy<=camRad; dy++){ const Y=Y0+dy; if(Y<0||Y>=h) continue; const baseY=kcy+(Y-my-0.5)*oversample;
        for(let dx=-camRad; dx<=camRad; dx++){ const X=X0+dx; if(X<0||X>=w) continue; const baseX=kcx+(X-mx-0.5)*oversample;
          let sum=0;
          for(let sy=0; sy<oversample; sy++){ const ky=baseY+sy+0.5;
            for(let sx=0; sx<oversample; sx++) sum+=ref_sampleKernelAt(interpMode, kernel, nx, ny, baseX+sx+0.5, ky); }
          img[Y*w+X]+=photons*sum; } }
    }
    // ---- fixture ----
    const cfg = readPsfConfigFromUI();
    const kstack = await getOrBuildPsfKernelStack(cfg);
    const W = 64, H = 64, os = cfg.oversample, camRad = Math.max(1, Math.ceil(cfg.halfWidthPx / os));
    const rng = mulberry32(99);
    const inside = Array.from({ length: 6 }, () => [4 + rng() * 56, 4 + rng() * 56]);
    const edges = [[-2.3, 30.1], [0.49, 0.51], [63.5, 12.25], [66.2, 70.7], [31.999, -0.5], [20.5, 63.5]];
    const out = [];
    for (const zi of [kstack.focusIndex, 3]) {
      const norm = normalizePsfPlane(kstack.slices[zi]).norm;
      for (const mode of ['nearest', 'linear', 'cubic', 'fft']) {
        const pts = mode === 'fft' ? [...inside.slice(0, 2), edges[0], edges[2]] : [...inside, ...edges];
        const plane = simKernelPlane(norm, kstack.nx, kstack.ny, os, mode);
        const objNew = { ...plane, nx: kstack.nx, ny: kstack.ny, oversample: os };
        const objOld = { kernel: norm, nx: kstack.nx, ny: kstack.ny, oversample: os };
        const a = new Float64Array(W * H), c = new Float64Array(W * H);
        let tNew = 0, tOld = 0;
        for (const [x, y] of pts) {
          let t = performance.now(); splatZernikeEmitter(a, W, H, x, y, 1000, objNew, camRad, mode); tNew += performance.now() - t;
          t = performance.now(); ref_splat(c, W, H, x, y, 1000, objOld, camRad, mode); tOld += performance.now() - t;
        }
        let peak = 0, maxd = 0, sa = 0, sc = 0;
        for (let i = 0; i < a.length; i++) { peak = Math.max(peak, Math.abs(c[i])); maxd = Math.max(maxd, Math.abs(a[i] - c[i])); sa += a[i]; sc += c[i]; }
        out.push({ zi, mode, n: pts.length, relMax: maxd / peak, relSum: Math.abs(sa - sc) / sc, tNew, tOld });
      }
    }
    return { camRad, nx: kstack.nx, os, out };
  });
  console.log(`(b) kernel ${b.nx}² oversample ${b.os}, camRad ${b.camRad}`);
  for (const r of b.out)
    check(`(b) summed splat = old splat  plane ${String(r.zi).padStart(2)} ${r.mode.padEnd(7)} (${r.n} emitters)`,
      r.relMax < 1e-5 && r.relSum < 1e-6,
      `max|Δ|/peak ${r.relMax.toExponential(2)}, Δsum ${r.relSum.toExponential(2)}, ${r.tOld.toFixed(0)} → ${r.tNew.toFixed(0)} ms`);

  // ---- noise statistics on the CPU (no GPU needed) ----------------------------------------------
  const st = await page.evaluate(() => {
    const res = {};
    for (const [type, lam] of [['scmos', 5], ['scmos', 200], ['emccd', 5], ['emccd', 200]]) {
      const n = 64 * 64 * 8, img = new Float32Array(n).fill(lam);
      const cam = { type, qe: 1, emGain: 1e6, cic: 0, bitDepth: 30 };
      // gain 1 ADU/e, zero read noise, zero offset: the ADU values ARE the counts
      applySimCameraNoise(img, new Float32Array(n), 1, 0, 777, 3, cam);
      let m = 0, v = 0; for (const x of img) m += x; m /= n; for (const x of img) v += (x - m) ** 2; v /= n - 1;
      res[`${type} λ=${lam}`] = { mean: m, fano: v / m };
    }
    // reproducible and address-only: the same (seed, frame) twice gives the same frame
    const a = new Float32Array(4096).fill(12), c = new Float32Array(4096).fill(12);
    applySimCameraNoise(a, new Float32Array(4096), 1, 2, 5, 9, { type: 'scmos' });
    applySimCameraNoise(c, new Float32Array(4096), 1, 2, 5, 9, { type: 'scmos' });
    res.repeat = a.every((x, i) => x === c[i]);
    return res;
  });
  for (const [k, r] of Object.entries(st)) if (k !== 'repeat') {
    const lam = Number(k.split('=')[1]), fanoWant = k.startsWith('emccd') ? 2 : 1;
    check(`noise statistics ${k}`, Math.abs(r.mean - lam) / lam < 0.02 && Math.abs(r.fano - fanoWant) < 0.08,
      `mean ${r.mean.toFixed(3)}, var/mean ${r.fano.toFixed(3)} (want ${fanoWant})`);
  }
  check('noise reproducible from (seed, frame)', st.repeat);

  const gpu = await checkGpu(page);
  if (!gpu.available) console.log('GPU checks: SKIP (WebGPU unavailable)');
  else if (typeof (await page.evaluate(() => typeof simulateFramesGpu)) === 'string' &&
           (await page.evaluate(() => typeof simulateFramesGpu)) === 'undefined') console.log('GPU checks: SKIP (no GPU simulator in this build)');
  else {
    // ---- (a) pcg4d bit-exact ----
    const a = await page.evaluate(async () => simGpuSelfTestRng(await getGpuEngine(), 1 << 18));
    check('(a) pcg4d WGSL = JS, bit-exact', a.mismatches === 0, `${a.n.toLocaleString()} u32 compared, ${a.mismatches} mismatches`);

    // ---- (c)+(d) GPU frame stage vs CPU frame stage ----
    const cd = await page.evaluate(async () => simGpuSelfTestFrames(await getGpuEngine()));
    for (const r of cd.splat)
      check(`(c) GPU splat = CPU splat  ${r.label}`, r.relMax < 2e-5 && r.relSum < 1e-5,
        `max|Δ|/peak ${r.relMax.toExponential(2)}, Δsum ${r.relSum.toExponential(2)}`);
    for (const r of cd.noise)
      check(`(d) GPU noise = CPU noise  ${r.label}`, r.fracClose >= 0.999,
        `${(100 * r.fracClose).toFixed(3)}% of pixels within 1e-3 ADU, max|Δ| ${r.maxAbs.toExponential(2)}`);

    // ---- (e) GPU PSF (direct) vs CPU direct ----
    const e = await page.evaluate(async () => simGpuSelfTestPsf(await getGpuEngine()));
    check('(e) GPU direct PSF = CPU direct PSF', e.relMax < 1e-4, `max|Δ|/peak ${e.relMax.toExponential(2)} over ${e.nz} planes`);
  }
} finally {
  await browser.close();
}
if (failed) { console.log('\nSOME CHECKS FAILED'); process.exitCode = 1; } else console.log('\nALL CHECKS PASSED');
