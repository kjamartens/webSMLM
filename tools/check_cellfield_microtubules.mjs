#!/usr/bin/env node
// Containment + slope-following regression checker for cell_field_sim's
// microtubule generator (cell_field_sim/microtubules.js). This targets the
// separate `cell_field_sim/` prototype, NOT the main webSMLM.html app --
// reuses tools/node_modules (the same Playwright dependency webSMLM-cli.mjs
// already brings in) rather than adding a second install step.
//
// What it checks, and why:
//  1. CONTAINMENT -- every generated point must sit (a) within the cell's own
//     blobby footprint (cellRadiusAt), (b) outside the nucleus ellipsoid, and
//     (c) at or below sampleCytoMeshHeight(cell,p,x,y) -- the ACTUAL rendered
//     (smoothed) cytoplasm height, not the raw analytic cytoHeightAt(). This
//     is not circular even though mtClampIntoCytoplasm() now clamps against
//     the very same sampleCytoMeshHeight(): it verifies the clamp is actually
//     being applied end-to-end, catching a future regression where a code
//     path skips the clamp or reverts to the raw cytoHeightAt().
//  2. SLOPE-FOLLOWING -- histograms z as a FRACTION of the local ceiling
//     (frac = pt.z / sampleCytoMeshHeight(cell,p,pt.x,pt.y)) across every
//     point of every microtubule. The bug this catches is bimodality: a
//     population pinned almost entirely near frac~1 (still at nucleus
//     height) and frac~0 (pinned to the floor), with a near-empty middle
//     band, is exactly the "hard snap, no slope" artifact a user reported and
//     the fraction-of-local-ceiling z model (see microtubules.js's own
//     mtGenerateOne comment) exists to fix.
//  3. STEP SLOPE -- for every consecutive pair of points on a path, the ratio
//     |dz| / lateral distance. The bug this catches is a real, reported
//     "incredibly sudden, fully vertical drop": an earlier version of the
//     over/under-the-nucleus crossing rule (see mtGenerateOne's own comment)
//     blended toward a FIXED absolute height regardless of what the
//     surrounding fraction-based height happened to be, so leaving the
//     nucleus's blend zone could snap z from that fixed value back to a very
//     different natural one over a short lateral distance -- a near-vertical
//     segment. A handful of genuinely steep single steps is expected (e.g.
//     the very first step off a near-polar nucleus-surface start point), so
//     this only fails on a non-trivial SHARE of steps being that steep, not
//     any single one.
//  4. TURN RADIUS -- for every interior point, the implied local radius of
//     curvature (avg segment length / turn angle between the incoming and
//     outgoing segment) against `p.mtMinTurnRadius`. The bug this catches is
//     a real, reported "very tight turn radii": the raw random walk's own
//     per-step heading cap alone left ~16% of turns violating the intended
//     radius (measured directly, some down to ~0.001 um) because the
//     Brownian-bridge drift correction added on top of that walk has no
//     curvature bound of its own -- mtEnforceMinTurnRadius() (see its own
//     comment in microtubules.js) fixes this on the FINAL rendered geometry.
//     A small residual share is expected (the independent containment clamp,
//     which runs after this and moves a point without any awareness of its
//     neighbours, can reintroduce a tight local turn) -- this only fails on a
//     non-trivial share, not any single violation.
//
// Usage:
//   cd tools && npm install   (once -- installs Playwright + Chromium)
//   node check_cellfield_microtubules.mjs
//
// Exits non-zero on any containment violation or failed slope check, so this
// is usable as a real regression gate after touching microtubule or
// cytoplasm-height-field code.

import { chromium } from 'playwright';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const htmlPath = path.resolve(__dirname, '..', 'cell_field_sim', 'index.html');
const htmlUrl = 'file://' + htmlPath.replace(/\\/g, '/');

// A tiny centre-region view -- enough chunks to see several cells (some
// present, some not, per p.density) without spending too long on packing.
const CHUNK_RANGE = 3; // -CHUNK_RANGE..CHUNK_RANGE in both cx,cy
const SEED = 1234;

// The tool's own current defaults, plus the exact bug-triggering
// high-density/large-offset scenario from the report.
const SCENARIOS = [
  { name: 'defaults', overrides: {} },
  {
    name: 'bug-report (high density, large offset)',
    overrides: { mtDensity: 0.825, mtStartFracMax: 0.3, mtStartOffsetXY: 2.0, mtEndFracMin: 0.3, mtMinSeparation: 0.5 },
  },
  {
    // Wide End direction jitter + Wobble turn strength are what actually
    // exercise the over/under-the-nucleus crossing rule (see
    // mtNucleusFootprintBlend's own comment) -- neither scenario above
    // touches it, so the step-slope check below would otherwise never run
    // against the one code path it exists to catch a regression in.
    name: 'nucleus-crossing (wide jitter)',
    overrides: { mtDensity: 0.3, mtEndJitterDeg: 180, mtWobbleTurn: 0.6, mtStartFracMax: 0.15, mtEndFracMax: 0.6 },
  },
];

const CONTAIN_TOL = 1e-3; // um -- tiny float-rounding slack, not a real violation budget
const MID_BAND_MIN_SHARE = 0.12; // well under the ~1/3 a uniform distribution would give
const STEEP_SLOPE_RATIO = 6; // |dz| more than 6x the lateral step distance (>~80.5 deg from horizontal)
const STEEP_STEP_MAX_SHARE = 0.01; // at most 1% of steps allowed to be that steep (a handful near start points is expected)
// Below this, |dz| is small enough on its own (a few tenths of a um, against
// typical cell/nucleus dimensions of tens of um) not to read as a visible
// "sudden drop" regardless of the RATIO -- the ratio alone spikes whenever
// the free walk's own realized lateral step happens to land near zero (a
// real, expected occurrence, not a defect -- see microtubules.js's own
// per-step z cap, deliberately sized against the NOMINAL step length rather
// than the realized lateral distance for exactly this reason), and without
// this floor the check would flag plenty of genuinely tiny, invisible
// height wiggles as if they were the artifact it exists to catch.
const STEEP_MIN_DZ = 0.3; // um
const TURN_RADIUS_SLACK = 0.9; // flag only when noticeably tighter than the configured minimum (float/geometry slack)
const TIGHT_TURN_MAX_SHARE = 0.02; // at most 2% of interior turns allowed to violate

async function main() {
  console.log(`Launching Chromium (headless)...`);
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage();
  page.on('pageerror', err => console.error('[page error]', err));

  await page.goto(htmlUrl);
  await page.waitForFunction(() => typeof buildMicrotubulesForCell === 'function' && typeof sampleCytoMeshHeight === 'function');

  let anyFail = false;

  for (const scenario of SCENARIOS) {
    const result = await page.evaluate(({ overrides, seed, chunkRange, containTol, midBandMinShare, steepSlopeRatio, steepMinDz, turnRadiusSlack }) => {
      for (const [id, val] of Object.entries(overrides)) {
        const el = document.getElementById(id);
        if (!el) throw new Error(`Unknown control id: ${id}`);
        el.value = String(val);
      }

      const p = params();
      ensureCytoCacheFresh(seed, p);

      const map = buildCandidateMap(seed, -chunkRange, -chunkRange, chunkRange, chunkRange, p);
      if (p.enablePacking) {
        relax(map, p, p.relaxIters);
        for (let round = 0; round < 6; round++) {
          const r = prune(map, p);
          if (r === 0) break;
          relax(map, p, Math.min(8, p.relaxIters));
        }
      }

      let nCells = 0, nMt = 0, nPts = 0;
      let containViolations = 0;
      let worstContainViolation = 0;
      const containExamples = [];
      // Coarse 3-bin histogram of frac = z / localCeiling across every point.
      const bins = [0, 0, 0]; // [0,1/3), [1/3,2/3), [2/3,1]
      let totalFracPts = 0;
      let totalSteps = 0, steepSteps = 0, worstSlope = 0;
      const steepExamples = [];
      let totalTurns = 0, tightTurns = 0, worstRadius = Infinity;
      const tightExamples = [];
      const minTurnRadius = p.mtMinTurnRadius;

      for (const cell of map.values()) {
        if (!cell.present) continue;
        nCells++;
        const mts = buildMicrotubulesForCell(seed, cell.cx, cell.cy, cell, p);
        for (const pts of mts) {
          nMt++;
          for (const pt of pts) {
            nPts++;

            // (a) radial footprint containment
            const ang = Math.atan2(pt.y, pt.x);
            const rc = cellRadiusAt(cell, ang);
            const dist = Math.hypot(pt.x, pt.y);
            const radialViol = dist - rc;

            // (b) nucleus ellipsoid exclusion (must be OUTSIDE, norm >= 1)
            const dxN = pt.x - cell.nucOffX, dyN = pt.y - cell.nucOffY;
            const cr = Math.cos(-cell.nucRot), sr = Math.sin(-cell.nucRot);
            const lx = dxN * cr - dyN * sr, ly = dxN * sr + dyN * cr;
            const a = cell.nucLong / 2, b = cell.nucShort / 2, rz = cell.nucHeight / 2;
            const nx = lx / Math.max(1e-6, a), ny = ly / Math.max(1e-6, b), nz = (pt.z - cell.nucZ) / Math.max(1e-6, rz);
            const ellNorm = Math.sqrt(nx * nx + ny * ny + nz * nz);
            const nucViol = 1 - ellNorm; // positive means INSIDE (a violation)

            // (c) height ceiling -- the SAME function the clamp now calls
            const ceilH = Math.max(0, sampleCytoMeshHeight(cell, p, pt.x, pt.y));
            const heightViol = pt.z - ceilH; // positive means ABOVE the surface (a violation)
            const floorViol = -pt.z; // positive means below 0

            const worstHere = Math.max(radialViol, nucViol, heightViol, floorViol);
            if (worstHere > containTol) {
              containViolations++;
              worstContainViolation = Math.max(worstContainViolation, worstHere);
              if (containExamples.length < 5) {
                containExamples.push({ cx: cell.cx, cy: cell.cy, radialViol, nucViol, heightViol, floorViol });
              }
            }

            if (ceilH > 1e-6) {
              const frac = Math.min(1, Math.max(0, pt.z / ceilH));
              totalFracPts++;
              if (frac < 1 / 3) bins[0]++;
              else if (frac < 2 / 3) bins[1]++;
              else bins[2]++;
            }
          }

          // (3) step slope: |dz| vs. lateral distance between consecutive
          // points on THIS path -- see this script's own header comment.
          for (let k = 1; k < pts.length; k++) {
            const a = pts[k - 1], b = pts[k];
            const lateral = Math.hypot(b.x - a.x, b.y - a.y);
            const dz = Math.abs(b.z - a.z);
            const slope = dz / Math.max(lateral, 1e-6);
            totalSteps++;
            if (slope > worstSlope) worstSlope = slope;
            if (slope > steepSlopeRatio && dz > steepMinDz) {
              steepSteps++;
              if (steepExamples.length < 5) {
                steepExamples.push({ cx: cell.cx, cy: cell.cy, lateral, dz, slope });
              }
            }
          }

          // (4) turn radius: implied local radius of curvature at every
          // interior point -- see this script's own header comment.
          if (minTurnRadius > 0) {
            for (let i = 1; i < pts.length - 1; i++) {
              const a = pts[i - 1], b = pts[i], c = pts[i + 1];
              const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
              const uLen = Math.hypot(ux, uy), vLen = Math.hypot(vx, vy);
              if (uLen < 1e-9 || vLen < 1e-9) continue;
              const cosAng = Math.min(1, Math.max(-1, (ux * vx + uy * vy) / (uLen * vLen)));
              const turnAngle = Math.acos(cosAng);
              const avgLen = (uLen + vLen) / 2;
              const radius = turnAngle > 1e-6 ? avgLen / turnAngle : Infinity;
              totalTurns++;
              if (radius < worstRadius) worstRadius = radius;
              if (radius < minTurnRadius * turnRadiusSlack) {
                tightTurns++;
                if (tightExamples.length < 5) {
                  tightExamples.push({ cx: cell.cx, cy: cell.cy, radius, turnAngleDeg: turnAngle * 180 / Math.PI });
                }
              }
            }
          }
        }
      }

      const midShare = totalFracPts > 0 ? bins[1] / totalFracPts : 0;
      const steepShare = totalSteps > 0 ? steepSteps / totalSteps : 0;
      const tightShare = totalTurns > 0 ? tightTurns / totalTurns : 0;
      return {
        nCells, nMt, nPts, containViolations, worstContainViolation, containExamples,
        bins, totalFracPts, midShare,
        slopeFail: totalFracPts > 0 && midShare < midBandMinShare,
        totalSteps, steepSteps, steepShare, worstSlope, steepExamples,
        minTurnRadius, totalTurns, tightTurns, tightShare, worstRadius, tightExamples,
      };
    }, { overrides: scenario.overrides, seed: SEED, chunkRange: CHUNK_RANGE, containTol: CONTAIN_TOL, midBandMinShare: MID_BAND_MIN_SHARE, steepSlopeRatio: STEEP_SLOPE_RATIO, steepMinDz: STEEP_MIN_DZ, turnRadiusSlack: TURN_RADIUS_SLACK });

    console.log(`\n=== Scenario: ${scenario.name} ===`);
    console.log(`  cells=${result.nCells}  microtubules=${result.nMt}  points=${result.nPts}`);
    console.log(`  containment violations: ${result.containViolations} (worst=${result.worstContainViolation.toFixed(4)} um, tol=${CONTAIN_TOL})`);
    if (result.containExamples.length) {
      for (const ex of result.containExamples) {
        console.log(`    e.g. cell (${ex.cx},${ex.cy}): radial=${ex.radialViol.toFixed(4)} nuc=${ex.nucViol.toFixed(4)} height=${ex.heightViol.toFixed(4)} floor=${ex.floorViol.toFixed(4)}`);
      }
    }
    const [lo, mid, hi] = result.bins;
    console.log(`  frac histogram: [0,1/3)=${lo}  [1/3,2/3)=${mid}  [2/3,1]=${hi}  (mid share=${(result.midShare * 100).toFixed(1)}%, need >= ${(MID_BAND_MIN_SHARE * 100).toFixed(0)}%)`);

    const stepSlopeFail = result.steepShare > STEEP_STEP_MAX_SHARE;
    console.log(`  step slope: worst=${result.worstSlope.toFixed(2)}x  steep steps (>${STEEP_SLOPE_RATIO}x AND dz>${STEEP_MIN_DZ}um)=${result.steepSteps}/${result.totalSteps} (${(result.steepShare * 100).toFixed(2)}%, need <= ${(STEEP_STEP_MAX_SHARE * 100).toFixed(0)}%)`);
    if (result.steepExamples.length) {
      for (const ex of result.steepExamples) {
        console.log(`    e.g. cell (${ex.cx},${ex.cy}): lateral=${ex.lateral.toFixed(3)}um dz=${ex.dz.toFixed(3)}um slope=${ex.slope.toFixed(1)}x`);
      }
    }

    let tightTurnFail = false;
    if (result.minTurnRadius > 0) {
      tightTurnFail = result.tightShare > TIGHT_TURN_MAX_SHARE;
      const worstStr = Number.isFinite(result.worstRadius) ? result.worstRadius.toFixed(4) : 'n/a';
      console.log(`  turn radius: worst=${worstStr}um (min=${result.minTurnRadius}um)  tight turns (<${(result.minTurnRadius * TURN_RADIUS_SLACK).toFixed(3)}um)=${result.tightTurns}/${result.totalTurns} (${(result.tightShare * 100).toFixed(2)}%, need <= ${(TIGHT_TURN_MAX_SHARE * 100).toFixed(0)}%)`);
      if (result.tightExamples.length) {
        for (const ex of result.tightExamples) {
          console.log(`    e.g. cell (${ex.cx},${ex.cy}): radius=${ex.radius.toFixed(4)}um turnAngle=${ex.turnAngleDeg.toFixed(1)}deg`);
        }
      }
    }

    const failed = result.containViolations > 0 || result.slopeFail || stepSlopeFail || tightTurnFail;
    if (failed) anyFail = true;
    console.log(`  ${failed ? 'FAIL' : 'PASS'}${result.containViolations > 0 ? ' -- containment violations found' : ''}${result.slopeFail ? ' -- middle band under-populated (hard-snap artifact)' : ''}${stepSlopeFail ? ' -- too many steep single-step z jumps (sudden-drop artifact)' : ''}${tightTurnFail ? ' -- too many tight turns (min-turn-radius regression)' : ''}`);
  }

  await browser.close();

  if (anyFail) {
    console.error('\nOne or more scenarios FAILED.');
    process.exit(1);
  }
  console.log('\nAll scenarios PASSED.');
}

main().catch(err => { console.error(err); process.exit(1); });
