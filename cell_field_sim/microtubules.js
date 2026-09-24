// Wrapped in a function ONLY so index.html can read this file's own source text
// (Function.prototype.toString) and hand it to the generator Web Workers: a page
// opened straight from disk (file://) cannot fetch() or importScripts() a sibling
// file, but it can always read a function's source. index.html re-runs the body
// as a plain global script, so everything below behaves as if it were top-level.
window.__MT_SRC = function () {
'use strict';

// ---- Microtubules ------------------------------------------------------------
// Kept in its own file, deliberately independent of the cell/nucleus/cytoplasm
// generator in index.html's own inline <script> -- so this whole module can be
// swapped, deleted, or (eventually) ported on its own. It is a PURE
// computational core: every function here takes plain numbers/plain cell-
// geometry fields (plus a (seed,cx,cy) chunk address) and returns plain
// {x,y,z} point arrays. Nothing in this file touches canvas/DOM -- only
// index.html's own draw() loop (which calls buildMicrotubulesForCell below)
// projects/strokes the result.
//
// Units: micrometres (um) for every length, radians for every angle, unless
// named "Deg". Coordinate frame: the cell's cached LOCAL frame -- centre-
// relative (origin at the cell's own centre, i.e. c.x/c.y already subtracted),
// pre-packRot -- the exact same frame cellRadiusAt()/nucleusSignedDistLocal()/
// cellOutlineLocal()/sampleCytoMeshHeight() already use in index.html. The
// caller (draw()) rotates/translates each returned point into world space via
// the existing localToWorld(cell, x, y) before projecting, exactly like the
// cytoplasm mesh and nucleus already do -- this file never needs to know
// about c.x/c.y/c.packRot at all.
//
// Depends on a small, fixed set of helpers/fields defined in index.html's own
// script (shared global scope, plain <script> tags -- see that file's own
// comment on load order): hashUnit(seed,cx,cy,k), hashStream(seed,cx,cy,base),
// lerp(a,b,t), smoothstep(t), cellRadiusAt(cell,thetaLocal), cellOutlineLocal(cell,n),
// nucleusSignedDistLocal(cell,lx,ly),
// sampleCytoMeshHeight(cell,p,x,y) -- the ACTUAL rendered (smoothed) height
// field, not the raw analytic cytoHeightAt(), so containment matches what's
// drawn; and reads cell.semiMajor/semiMinor/rot/
// harmAmp/harmPh/modFloor/rOuter/nucOffX/nucOffY/nucLong/nucShort/nucRot/
// nucZ/nucHeight (all already resolved by rawCandidate()/envelopNucleus()).
//
// RNG: every random draw goes through the same address-based hashUnit/
// hashStream(seed,cx,cy,channel) convention the rest of the file uses (ported
// from demoCam_SMLM_MM's sim::Pcg4d / SMLMCounterRng.h -- see index.html's own
// top-of-file comment) -- never a mutable-state PRNG. A draw depends only on
// its own address, so a given seed always reproduces byte-identical
// microtubules regardless of pan/zoom/draw order. Channel numbering below is
// namespaced well above the cell generator's own CH.* (which tops out at 44)
// and its hashStream bases (none currently used elsewhere) so nothing collides.
//
// Ground-up portability note: this file is written so the SAME algorithm could
// be reimplemented in demoCam_SMLM_MM's own C++ simulation engine (which
// already shares this project's pcg4d hash) without a redesign -- plain
// scalar math, no JS-only idioms, named constants for every channel/bin-count/
// iteration-cap, and this header documenting units/frame/RNG convention
// explicitly since a C++ port won't have this repo's plan file to refer back
// to.

// Direct hashUnit channels this file uses (never multiplied by hashStream's
// own *4096 addressing -- see hashStream's own comment in index.html) --
// picked well above the cell generator's CH.* (max 44).
const MT_CH_COUNT = 500;

// hashStream base for a single microtubule's OWN generation draws (direction,
// start/end offsets, per-step wobble, priority) -- every draw for microtubule
// `mtIndex` on `resampleRound` attempt comes from one continuous stream, same
// pattern the removed buildMicrotubule() used (hashStream(seed,cx,cy,100+i)),
// just re-spaced (mtIndex directly, resampleRound*100000) so a resample
// attempt draws genuinely different numbers rather than repeating the first
// attempt's own sequence.
const MT_STREAM_BASE = 1000;
const MT_RESAMPLE_SPACING = 100000;

// Direction weighting: how many azimuth bins to score around the nucleus
// centre (see buildMtDirectionTable) -- 64 matches the default cytoTheta's own
// resolution choice in index.html for the same reason (smooth enough at any
// reasonable zoom without the ray-march cost blowing up).
const MT_N_DIR = 64;
// Steps for the coarse linear march (+ a bisection refine) that finds where a
// ray from the nucleus centre exits the cell's own blobby footprint -- see
// rayCellBoundaryFromNucleus.
const MT_MARCH_STEPS = 24;
const MT_BISECT_ITERS = 24;

// Minimum-separation pass: bounded like the cell generator's own relax()/
// prune() (see index.html's own comments on PRUNE_ROUNDS) -- best-effort, not
// a hard guarantee, so a pathological configuration can't hang the tab.
// Kept deliberately SMALL (measured, not guessed): near the nucleus, many
// microtubules' own FIXED start points are packed tighter than a typical
// mtMinSeparation by construction once density gets even moderately high
// (mtStartFracMax confines every start to a thin band right at the nucleus
// edge) -- those specific conflicts can never actually resolve (neither
// point is allowed to move), so a round involving them never converges to
// "no violation found" and always burns its FULL round budget. Measured
// directly: at mtDensity=0.2 (well under half of the new 2/um^2 ceiling) the
// original MT_NUDGE_ROUNDS=15 + a 3-round resample fallback froze the tab's
// very first draw() call for well over 30s. A handful of rounds is enough to
// fix the realistic case this pass targets (a few stray, genuinely
// resolvable crossings), not the structurally-infeasible one.
const MT_NUDGE_ROUNDS = 5;
const MT_RESAMPLE_ROUNDS = 2;
const MT_POST_RESAMPLE_NUDGE_ROUNDS = 3;
// Per-round pairwise-check ceiling for the grid-based nudge pass (see
// mtNudgeRoundGrid's own comment) -- bounds a single round's cost even when
// many points end up crammed into one crowded neighbourhood (routine near
// the nucleus at high mtDensity), same idiom as CHUNK_CAP/RENDER_CAP
// elsewhere in index.html. This bounds the CANDIDATE-COMPARISON cost once a
// bucket is found; it does NOT bound how many buckets get looked up in the
// first place (every point checks its own 3x3x3 neighbourhood regardless) --
// that cost is what MT_COLLISION_MAX_TOTAL_POINTS below is really for.
const MT_COLLISION_OP_BUDGET = 300000;
// Skips collision resolution entirely for a cell whose microtubules add up to
// more total points than this. The grid (mtBuildSpatialIndex) turns the
// per-round cost from O(M^2*S^2) into O(N), but N*27 bucket lookups per
// round, repeated over several rounds, is still real work -- measured
// directly: 40000 was NOT low enough on its own (see MT_NUDGE_ROUNDS's own
// comment on the actual hang this produced) until the round counts above
// were also cut down; the two fixes together are what keeps this fast.
// Every microtubule is still individually contained inside the cell volume
// (mtClampIntoCytoplasm runs regardless of this cap), just not necessarily
// separated from its neighbours once it's hit -- a documented best-effort
// limit, not silent data loss.
const MT_COLLISION_MAX_TOTAL_POINTS = 8000;
// Hard ceiling on a SINGLE microtubule's own step count (see mtGenerateOne) --
// independent of the per-cell collision caps above, so a small Step length
// combined with a large Wobble path length multiplier can't alone balloon
// one path's point count without bound.
const MT_MAX_STEPS_PER_MT = 1500;
// Hard ceiling on microtubule COUNT per cell -- mtDensity's own 2/um^2
// ceiling times a large blobby cell's real footprint area can otherwise ask
// for many thousands of paths; this bounds generation cost directly rather
// than relying only on the collision-pass caps above (which apply to
// SEPARATION, not to whether the paths get built/drawn at all). Same
// "just cap it" idiom as index.html's own CHUNK_CAP/RENDER_CAP.
const MT_MAX_PER_CELL = 5000;

// Containment clamp keeps a nudged/wobbled point just inside the cell's own
// blobby footprint/nucleus surface rather than exactly on it (a point sitting
// EXACTLY on a boundary can re-trigger the same clamp on the next pass due to
// floating-point rounding).
const MT_CONTAIN_MARGIN = 0.98;

// How far beyond the nucleus's own lateral (equatorial) footprint the
// "must clear the nucleus vertically" rule (see mtNucleusFootprintBlend,
// used by mtGenerateOne below) ramps back down to 0 -- expressed as a
// fraction of the nucleus's own equatorial radius in that direction, so it
// scales with nucleus size rather than being a fixed um value. Wide enough
// for a smooth, visually continuous rise/fall as a path enters/leaves the
// footprint (no visible kink at the exact ellipse boundary), narrow enough
// that it only affects paths genuinely passing near the nucleus, not the
// whole cell.
const MT_NUCLEUS_CLEAR_BLEND = 0.4;

// How much vertical clearance a microtubule keeps from the nucleus surface
// while riding over/under it (see mtGenerateOne) -- half of the cell's own
// `nucMargin` (the SAME margin envelopNucleus()/cytoHeightAt() already
// guarantee the cytoplasm keeps above the nucleus everywhere within that
// margin band), not an independent hardcoded value, so raising/lowering
// nucMargin in the sidebar scales this too rather than the two drifting out
// of relation to each other. Floored at 0.05 um so a user-set nucMargin of
// exactly 0 (its slider minimum) still leaves a non-zero, visible gap.
// Minimum radius of curvature (um, `p.mtMinTurnRadius`, "Min turn radius" --
// 0 disables the constraint) a microtubule's own path is allowed to bend at,
// regardless of Wobble turn strength -- see mtEnforceMinTurnRadius's own
// comment below for why this has to run on the FINAL rendered geometry, not
// just the raw random walk that feeds it.

// Slope limiter ceiling -- now the settable `p.mtMaxZSlope` ("Max height
// slope (xxy)", slider 1-20 step 0.5, default 5, both wired in index.html): z may
// change at most this many times the lateral distance moved between two
// consecutive points. Used by TWO separate passes -- see each one's own
// comment for why one alone isn't enough:
//  1. The post-generation smoothing pass inside mtGenerateOne (below, `const
//     maxDz = ... * stepLen`), which bounds dz against the path's own
//     NOMINAL step length, not the REALIZED lateral distance actually moved
//     -- deliberately, since a persistent random walk can double back to
//     near-zero net lateral movement between two consecutive points, and
//     bounding purely by realized distance there would mask a genuine height
//     change rather than spread it out (see that pass's own comment).
//  2. `mtLimitZSlopeRealized` (below), a FINAL safety net run once per path
//     in buildMicrotubulesForCell, after collision resolution -- on the
//     REALIZED lateral distance between the truly final points. Added
//     because #1's own nominal-stepLen basis, plus points collision
//     resolution (`mtNudgeRoundGrid`) moves AFTER #1 already ran with no
//     slope recheck of its own, could still leave a real, reported "massive
//     change in z in one or a few steps" on the rendered geometry -- exactly
//     the ratio a viewer's own eye reads as "how steep is this segment", and
//     the same ratio tools/check_cellfield_microtubules.mjs's own step-slope
//     check measures.

function mtNucleusClearance(p) {
  return Math.max(0.05, 0.5 * p.nucMargin);
}

// 1 when (x,y) sits within the nucleus's own LATERAL (equatorial) elliptical
// footprint -- ignoring z entirely, unlike mtClampIntoCytoplasm's full 3D
// ellipsoid check -- ramping smoothly down to 0 by MT_NUCLEUS_CLEAR_BLEND
// beyond its edge. Used to decide how strongly a path point's z should be
// pulled toward clearing the nucleus vertically (mtGenerateOne) rather than
// following the ordinary fraction-of-local-ceiling height model, which on
// its own has no idea the nucleus sits in the way and would happily
// interpolate straight through it.
function mtNucleusFootprintBlend(cell, x, y) {
  const dxN = x - cell.nucOffX, dyN = y - cell.nucOffY;
  const cr = Math.cos(-cell.nucRot), sr = Math.sin(-cell.nucRot);
  const lx = dxN * cr - dyN * sr, ly = dxN * sr + dyN * cr;
  const a = cell.nucLong / 2, b = cell.nucShort / 2;
  const norm = Math.hypot(lx / Math.max(1e-6, a), ly / Math.max(1e-6, b));
  const t = Math.min(1, Math.max(0, (norm - 1) / MT_NUCLEUS_CLEAR_BLEND));
  return smoothstep(1 - t);
}

// Post-process curvature limiter, run on the FINAL rendered (x,y) of a
// microtubule's own path -- see mtGenerateOne's own comment on why the
// per-step heading-turn cap in the raw random walk isn't enough on its own
// (the Brownian-bridge drift correction added afterward is a separate,
// additive, fixed-direction pull with no curvature bound of its own, and the
// REALIZED point-to-point direction can still bend sharply even when the raw
// walk's own curvature is well within bounds).
//
// Works in SEGMENT-VECTOR space, not by nudging point positions directly: a
// first version tried pulling an over-curved point toward the midpoint of
// its two neighbours, a soft correction needing many repeated passes to
// converge on anything but a mild violation -- measured directly, it left
// ~3% of turns still violating the cap after 20 passes, some barely reduced
// from a near-total reversal at all. This version reads off each segment's
// own direction+length, clamps the ANGLE step between consecutive segments
// to the cap EXACTLY (one shot, not an iterative pull -- same "rotate to an
// explicit clamped angle" idiom mtGenerateOne's own raw-walk heading cap
// uses, just applied to the realized geometry instead of the walk that
// produced it), then rebuilds every point from the fixed start by summing
// the (now-clamped) segment vectors. That reconstruction generally no longer
// lands exactly on the true end point, so a fresh Brownian-bridge-style
// drift correction (same smoothstep taper as mtGenerateOne's own) pins it
// back -- which can reintroduce a little curvature of its own, so this whole
// clamp+redrift cycle repeats a few rounds; each round's own drift shortfall
// is much smaller than the last (most of the curvature was already fixed),
// so it converges quickly. Endpoints (i=0/i=steps) are never moved.
// `minRadius<=0` (the slider's own "off" position) skips this entirely.
function mtEnforceMinTurnRadius(pts, minRadius) {
  const n = pts.length;
  if (!(minRadius > 0) || n < 3) return;
  const endX = pts[n - 1].x, endY = pts[n - 1].y;

  for (let round = 0; round < 8; round++) {
    const segX = new Array(n - 1), segY = new Array(n - 1), segLen = new Array(n - 1);
    for (let i = 1; i < n; i++) {
      segX[i - 1] = pts[i].x - pts[i - 1].x;
      segY[i - 1] = pts[i].y - pts[i - 1].y;
      segLen[i - 1] = Math.hypot(segX[i - 1], segY[i - 1]);
    }

    let changed = false;
    for (let i = 1; i < segX.length; i++) {
      const prevLen = segLen[i - 1], curLen = segLen[i];
      if (prevLen < 1e-9 || curLen < 1e-9) continue;
      const prevAng = Math.atan2(segY[i - 1], segX[i - 1]);
      const curAng = Math.atan2(segY[i], segX[i]);
      let dAng = curAng - prevAng;
      while (dAng > Math.PI) dAng -= 2 * Math.PI;
      while (dAng < -Math.PI) dAng += 2 * Math.PI;
      const avgLen = (prevLen + curLen) / 2;
      const cap = 2 * Math.asin(Math.min(1, avgLen / (2 * minRadius)));
      if (Math.abs(dAng) <= cap) continue;
      const clampedAng = prevAng + Math.sign(dAng) * cap;
      segX[i] = Math.cos(clampedAng) * curLen;
      segY[i] = Math.sin(clampedAng) * curLen;
      changed = true;
    }

    for (let i = 1; i < n; i++) {
      pts[i].x = pts[i - 1].x + segX[i - 1];
      pts[i].y = pts[i - 1].y + segY[i - 1];
    }

    const dxErr = endX - pts[n - 1].x, dyErr = endY - pts[n - 1].y;
    if (Math.abs(dxErr) > 1e-9 || Math.abs(dyErr) > 1e-9) {
      for (let i = 1; i < n - 1; i++) {
        const t = i / (n - 1);
        const s = t * t * (3 - 2 * t); // smoothstep
        pts[i].x += dxErr * s;
        pts[i].y += dyErr * s;
      }
      pts[n - 1].x = endX; pts[n - 1].y = endY;
    }

    if (!changed) break;
  }
}

// A steep FIRST/LAST segment is a genuinely different case from a steep
// INTERIOR one, and `mtLimitZSlopeRealized` below (which only ever nudges
// i=1..n-2) structurally cannot fix it: i=0 and i=n-1 are never moved (their
// z is a deliberate, meaningful value -- a real point on the nucleus
// ellipsoid's surface, or an independently-drawn end fraction), so when the
// offending segment IS the endpoint one, the interior relax pass is bounded
// by that fixed z no matter how many rounds it runs -- it can only ever pull
// the NEIGHBOUR closer, and the neighbour's own budget is itself limited by
// ITS other neighbour. Measured directly: this is exactly why a single-digit
// mtMaxZSlope could still leave triple-digit realized slopes after the
// interior-only relax pass -- every worst offender traced back to the very
// first or last segment of its path. Fixing an impossible endpoint jump by
// deleting the endpoint (promoting its former neighbour to be the new,
// effective start/end) reads as the microtubule simply not growing that one
// extra bit -- a small, usually invisible truncation -- rather than forcing
// a still-too-steep segment to exist regardless, or fighting to nudge a
// value that structurally can't move far enough to matter. Run BEFORE
// mtLimitZSlopeRealized so that pass's own interior relax works from
// already-valid endpoints. Capped at MT_MAX_END_TRIM points per end so a
// pathological path (e.g. its whole first quarter genuinely near-vertical)
// can't be trimmed down to nothing.
const MT_MAX_END_TRIM = 50;
function mtTrimSteepEnds(pts, maxSlope) {
  if (!(maxSlope > 0)) return;
  for (let n = 0; n < MT_MAX_END_TRIM && pts.length > 2; n++) {
    const a = pts[0], b = pts[1];
    const lateral = Math.hypot(b.x - a.x, b.y - a.y);
    if (Math.abs(b.z - a.z) <= maxSlope * Math.max(lateral, 1e-6)) break;
    pts.shift();
  }
  for (let n = 0; n < MT_MAX_END_TRIM && pts.length > 2; n++) {
    const last = pts.length - 1;
    const a = pts[last], b = pts[last - 1];
    const lateral = Math.hypot(a.x - b.x, a.y - b.y);
    if (Math.abs(a.z - b.z) <= maxSlope * Math.max(lateral, 1e-6)) break;
    pts.pop();
  }
}

// Final safety-net slope limiter -- see MT_MAX_Z_SLOPE's own comment above
// for why this is a SEPARATE pass from mtGenerateOne's own nominal-stepLen
// version, not a duplicate: this one runs ONCE PER PATH in
// buildMicrotubulesForCell, AFTER collision resolution has finished moving
// points around, and measures the ratio against the REALIZED lateral
// distance between the truly final pair of points -- catching whatever a
// post-generation nudge (or a doubled-back realized step the nominal check's
// own stepLen basis was never meant to police) leaves behind. Same
// alternating-sweep relaxation idiom as the nominal-stepLen version (a
// violation can span several points; each pass propagates the correction one
// step further). Endpoints (i=0/i=n-1) are never moved -- same reasoning as
// every other postprocess pass in this file: their z is a deliberate,
// meaningful value (nucleus-surface start, independently-drawn end). See
// `mtTrimSteepEnds` above for the complementary fix when the endpoint ITSELF
// is the unfixable offender, which must run first.
// `maxSlope<=0` (not reachable via the slider, whose minimum is 1, but kept
// as a safe no-op for a directly-scripted config) skips this entirely.
function mtLimitZSlopeRealized(pts, maxSlope) {
  const n = pts.length;
  if (!(maxSlope > 0) || n < 3) return;
  for (let pass = 0; pass < 20; pass++) {
    let changed = false;
    const forward = pass % 2 === 0;
    for (let k = 1; k < n - 1; k++) {
      const i = forward ? k : n - 1 - k;
      const prev = pts[i - 1], cur = pts[i], nxt = pts[i + 1];
      const dPrev = Math.hypot(cur.x - prev.x, cur.y - prev.y) * maxSlope;
      const dNxt = Math.hypot(nxt.x - cur.x, nxt.y - cur.y) * maxSlope;
      const lo = Math.max(prev.z - dPrev, nxt.z - dNxt);
      const hi = Math.min(prev.z + dPrev, nxt.z + dNxt);
      // Same "leave it, containment clamp still guards it" reasoning as the
      // nominal-stepLen pass's own lo>hi case: the two neighbours are
      // themselves too far apart in z for any single value to satisfy both
      // slope budgets against their OWN realized distances at once.
      if (lo > hi) continue;
      if (cur.z < lo) { cur.z = lo; changed = true; }
      else if (cur.z > hi) { cur.z = hi; changed = true; }
    }
    if (!changed) break;
  }
}

function mtNormalize3(v) {
  const n = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / n, v[1] / n, v[2] / n];
}
// Exact polygon area (shoelace) over cellOutlineLocal's own sampled points --
// reused directly rather than a circle-equivalent approximation, since the
// footprint is already sampled at fixed angular steps for rendering.
function mtShoelaceArea(pts) {
  let area = 0;
  for (let i = 0; i < pts.length; i++) {
    const [x1, y1] = pts[i], [x2, y2] = pts[(i + 1) % pts.length];
    area += x1 * y2 - x2 * y1;
  }
  return Math.abs(area) / 2;
}

// Nucleus radius at LOCAL azimuth theta (measured from the nucleus's own
// centre, world/local-frame convention -- same "theta minus own rotation"
// form as cellRadiusAt's base-ellipse term).
function mtNucleusRadiusAt(cell, theta) {
  const a = cell.nucLong / 2, b = cell.nucShort / 2;
  const phi = theta - cell.nucRot;
  return (a * b) / Math.hypot(b * Math.cos(phi), a * Math.sin(phi));
}

// Distance from the NUCLEUS centre to the cell's own (blobby, possibly
// off-centre-relative) footprint boundary along azimuth theta -- marches a
// ray p(t) = nucCenter + t*(cos theta, sin theta) outward from the nucleus's
// own edge (rNucStart) and tests hypot(p) > cellRadiusAt(cell, angle-from-
// CELL-centre) (cellRadiusAt's own convention measures from the local origin,
// i.e. the cell centre, not the nucleus centre -- p is already in that same
// local frame since the nucleus offset is centre-relative), refined by
// bisection once the boundary is bracketed. Falls back to the last sampled t
// in the (never expected in practice, envelopNucleus guarantees clearance)
// case the march never exits within the capped search radius.
function mtRayCellBoundaryFromNucleus(cell, theta, rNucStart) {
  const ux = Math.cos(theta), uy = Math.sin(theta);
  const maxT = (cell.rOuter + Math.hypot(cell.nucOffX, cell.nucOffY)) * 1.3 + 1;
  const dt = Math.max(1e-3, (maxT - rNucStart) / MT_MARCH_STEPS);
  let prevT = rNucStart, tLo = null, tHi = null;
  for (let s = 0; s <= MT_MARCH_STEPS; s++) {
    const t = rNucStart + s * dt;
    const px = cell.nucOffX + ux * t, py = cell.nucOffY + uy * t;
    const ang = Math.atan2(py, px);
    const rc = cellRadiusAt(cell, ang);
    const dist = Math.hypot(px, py);
    if (dist > rc) { tHi = t; tLo = s === 0 ? t : prevT; break; }
    prevT = t;
  }
  if (tHi === null) return prevT;
  for (let it = 0; it < MT_BISECT_ITERS; it++) {
    const mid = (tLo + tHi) / 2;
    const px = cell.nucOffX + ux * mid, py = cell.nucOffY + uy * mid;
    const ang = Math.atan2(py, px);
    const rc = cellRadiusAt(cell, ang);
    const dist = Math.hypot(px, py);
    if (dist > rc) tHi = mid; else tLo = mid;
  }
  return (tLo + tHi) / 2;
}

// Once per cell: a discrete weight table over MT_N_DIR azimuth bins around the
// nucleus centre, weight(theta) = max(0, rCellFromNucleus(theta) -
// rNuc(theta)) -- the cytoplasm "thickness" in that direction, i.e. how much
// room a microtubule has to grow into. Directions with more cytoplasm are
// proportionally more likely to be drawn (see mtSampleDirection). Returns a
// cumulative-weight array for fast inverse-CDF sampling.
function buildMtDirectionTable(cell) {
  const bins = new Array(MT_N_DIR);
  const cum = new Array(MT_N_DIR + 1);
  cum[0] = 0;
  for (let i = 0; i < MT_N_DIR; i++) {
    const theta = (i / MT_N_DIR) * Math.PI * 2;
    const rNuc = mtNucleusRadiusAt(cell, theta);
    const rCell = mtRayCellBoundaryFromNucleus(cell, theta, rNuc);
    const w = Math.max(0, rCell - rNuc);
    bins[i] = { theta, rNuc, rCell, w };
    cum[i + 1] = cum[i] + w;
  }
  return { bins, cum, total: cum[MT_N_DIR] };
}

// Draws one azimuth from the direction table: u1 picks a weighted bin
// (inverse-CDF, binary search), u2 jitters continuously within that bin's own
// angular width so microtubules aren't confined to MT_N_DIR fixed headings.
// Falls back to a uniform draw if every bin has zero weight (a degenerate
// cell where the nucleus fills essentially the whole footprint).
function mtSampleDirection(table, u1, u2) {
  if (table.total <= 1e-9) return u1 * Math.PI * 2;
  const target = u1 * table.total;
  let lo = 0, hi = MT_N_DIR;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (table.cum[mid + 1] < target) lo = mid + 1; else hi = mid;
  }
  const i = Math.min(MT_N_DIR - 1, lo);
  const binWidth = (Math.PI * 2) / MT_N_DIR;
  return table.bins[i].theta + (u2 - 0.5) * binWidth;
}

// Per-cell geometry cache (direction table, local outline, footprint area) --
// same reasoning as index.html's own cytoCache: these fields are a pure
// function of the cell's own already-resolved shape (semiMajor/semiMinor/rot/
// harmAmp/harmPh/modFloor/nucleus fields), not of its current (packing-
// relaxed) x,y, so recomputing this ray-marched table from scratch on every
// draw() call (including every frame of a plain pan/zoom/tilt drag where
// nothing shape-related has changed) would reintroduce the exact
// "incredibly slow when moving around" problem cytoCache itself was built to
// fix. Keyed on chunk id, invalidated by a signature over the cell's own
// resolved numeric fields (not on index.html's `p` directly, keeping this
// file self-contained).
const MT_GEOM_CACHE_MAX = 6000;
const mtGeomCache = new Map();
function mtCellShapeSig(cell) {
  return [cell.semiMajor, cell.semiMinor, cell.rot, cell.harmAmp.join(','), cell.harmPh.join(','),
    cell.modFloor, cell.nucOffX, cell.nucOffY, cell.nucLong, cell.nucShort, cell.nucRot,
    cell.nucZ, cell.nucHeight, cell.rOuter].join('|');
}
function getMtCellGeometry(cell) {
  const key = cell.cx + ',' + cell.cy;
  const sig = mtCellShapeSig(cell);
  const entry = mtGeomCache.get(key);
  if (entry && entry.sig === sig) return entry;
  if (mtGeomCache.size > MT_GEOM_CACHE_MAX) mtGeomCache.clear();
  const localOutline = cellOutlineLocal(cell, 48);
  const fresh = {
    sig,
    dirTable: buildMtDirectionTable(cell),
    localOutline,
    areaUm2: mtShoelaceArea(localOutline),
  };
  mtGeomCache.set(key, fresh);
  return fresh;
}

// Density is microtubules per um^2 of the cell's own footprint area (not per
// cell, and not tied to semiMajor the way the old placeholder was) --
// fractional expected counts are resolved via a per-cell hash draw
// (probabilistic rounding) so low densities don't band every cell to the same
// integer count.
function mtCountForCell(seed, cx, cy, cell, p, geom) {
  const g = geom || getMtCellGeometry(cell);
  const expected = Math.max(0, g.areaUm2 * p.mtDensity);
  const base = Math.floor(expected);
  const frac = expected - base;
  const roll = hashUnit(seed, cx, cy, MT_CH_COUNT);
  return base + (roll < frac ? 1 : 0);
}

// Pushes `pt` (mutated in place) back inside the cell's own cytoplasm volume
// if the free/wobbled walk carried it out -- the whole microtubule has to
// stay inside the cell, not just its (already-valid-by-construction)
// endpoints. Three checks, in order:
//  1. Out of the NUCLEUS (a real 3D ellipsoid, not just its 2D footprint --
//     rotated into the nucleus's own local frame, then pushed back along the
//     radial direction in that normalized ellipsoid space if inside).
//  2. Inside the cell's own blobby FOOTPRINT (radial clamp against
//     cellRadiusAt at this point's own azimuth).
//  3. Between 0 and the cytoplasm HEIGHT-FIELD at this (possibly
//     footprint-clamped) xy position.
// A hard clamp, not a physically-motivated bounce/reflection -- simple,
// always exactly satisfies containment, and reads as the microtubule
// bending along the nuclear envelope/cell cortex when wobble would have
// carried it through, a reasonable stand-in for a real prototype.
function mtClampIntoCytoplasm(cell, p, geom, pt) {
  const dxN = pt.x - cell.nucOffX, dyN = pt.y - cell.nucOffY;
  const cr = Math.cos(-cell.nucRot), sr = Math.sin(-cell.nucRot);
  const lx = dxN * cr - dyN * sr, ly = dxN * sr + dyN * cr;
  const a = cell.nucLong / 2, b = cell.nucShort / 2, rz = cell.nucHeight / 2;
  const nz = (pt.z - cell.nucZ) / Math.max(1e-6, rz);
  const nx = lx / Math.max(1e-6, a), ny = ly / Math.max(1e-6, b);
  const ellNorm = Math.sqrt(nx * nx + ny * ny + nz * nz);
  if (ellNorm < 1) {
    // Pushing OUT of the nucleus needs the normalized radius to end up
    // GREATER than 1 (outside), so the margin divides rather than multiplies
    // here -- the opposite direction from the footprint/height clamps below,
    // which pull a point back INSIDE their own outer bound.
    const scale = (ellNorm > 1e-9 ? 1 / ellNorm : 1) / MT_CONTAIN_MARGIN;
    const lx2 = lx * scale, ly2 = ly * scale;
    pt.z = cell.nucZ + nz * rz * scale;
    const cr2 = Math.cos(cell.nucRot), sr2 = Math.sin(cell.nucRot);
    pt.x = cell.nucOffX + lx2 * cr2 - ly2 * sr2;
    pt.y = cell.nucOffY + lx2 * sr2 + ly2 * cr2;
  }

  const ang = Math.atan2(pt.y, pt.x);
  const rc = cellRadiusAt(cell, ang);
  const dist = Math.hypot(pt.x, pt.y);
  if (dist > rc) {
    const scale = (rc * MT_CONTAIN_MARGIN) / Math.max(1e-9, dist);
    pt.x *= scale; pt.y *= scale;
  }

  // Against the ACTUAL rendered (smoothed) mesh height, not the raw analytic
  // cytoHeightAt() -- see sampleCytoMeshHeight's own comment in index.html:
  // buildCytoMesh()'s Laplacian smoothing lowers the sharp nucleus-adjacent
  // dome peak below what cytoHeightAt() alone would return, so clamping
  // against the raw function let a point sit above the surface actually
  // drawn (a real, reported "pokes out of the dome" bug).
  const topH = Math.max(0, sampleCytoMeshHeight(cell, p, pt.x, pt.y));
  if (pt.z < 0) pt.z = 0;
  else if (pt.z > topH) pt.z = topH * MT_CONTAIN_MARGIN;
}

// Builds one microtubule's full geometry (steps 2-6 of the design): a weighted
// starting direction from the nucleus centre, a start position near the
// nucleus edge, an end direction (jittered -- see End direction jitter,
// widened up to 180 deg so an end can genuinely land on the far side of the
// nucleus, not just a nearby azimuth) and end position near the cell edge,
// and a free correlated-random-walk path in XY (real loops/U-bends allowed)
// between them, corrected to land exactly on both endpoints. z comes from
// the fraction-of-local-ceiling model (see the per-point loop below) with an
// over/under-the-nucleus override applied wherever a point's (x,y) passes
// near or over the nucleus's own lateral footprint -- see
// mtNucleusFootprintBlend's own comment -- and the whole path is then
// clamped to stay inside the cell's own cytoplasm volume throughout.
// `resampleRound`
// (0 = first attempt) reseeds the whole draw sequence when the collision pass
// below needs to regenerate this microtubule from scratch.
function mtGenerateOne(seed, cx, cy, mtIndex, resampleRound, cell, p, geom) {
  const next = hashStream(seed, cx, cy, MT_STREAM_BASE + mtIndex + resampleRound * MT_RESAMPLE_SPACING);

  const theta0 = mtSampleDirection(geom.dirTable, next(), next());
  // The START point is a genuine 3D point on the nucleus ELLIPSOID's surface
  // (standard parametrization: azimuth theta0 -- still the weighted-by-
  // cytoplasm direction from buildMtDirectionTable, unchanged -- plus a polar
  // angle psi drawn uniform-on-sphere via cosPsi=1-2u), not azimuth-with-an-
  // independently-random-Z. That older version put every start point on the
  // nucleus's own EQUATORIAL ring regardless of its (unrelated) Z draw --
  // geometrically inconsistent (a ring-radius point paired with a near-pole Z
  // sits outside the true ellipsoid, or just reads as "microtubules only
  // start at the equator" once the containment clamp corrects it) and, more
  // visibly, a real cause of the whole population sitting flat near one
  // height: with every start already pinned near mid-height, a path only
  // reached the cytoplasm dome's own mid-height slope (see cytoHeightAt) if
  // its own random wobble happened to climb there, so the slope's own
  // mid-region routinely went uncovered (confirmed against a real screenshot
  // showing exactly that gap). Starting genuinely above/below the nucleus
  // too means a path travelling from a raised start down to its low, near-
  // the-edge end has to CROSS that slope by construction.
  const phi0 = theta0 - cell.nucRot;
  const cosPsi = 1 - 2 * next();
  const sinPsi = Math.sqrt(Math.max(0, 1 - cosPsi * cosPsi));
  const aN = cell.nucLong / 2, bN = cell.nucShort / 2, rzN = cell.nucHeight / 2;
  const lx0 = aN * sinPsi * Math.cos(phi0), ly0 = bN * sinPsi * Math.sin(phi0), lz0 = rzN * cosPsi;
  const crN = Math.cos(cell.nucRot), srN = Math.sin(cell.nucRot);
  const surfX = cell.nucOffX + lx0 * crN - ly0 * srN;
  const surfY = cell.nucOffY + lx0 * srN + ly0 * crN;
  const surfZ = cell.nucZ + lz0;
  const outward = mtNormalize3([surfX - cell.nucOffX, surfY - cell.nucOffY, lz0]);

  const rNuc0 = mtNucleusRadiusAt(cell, theta0); // 2D equatorial reference, only used below to scale how far "out" means
  const rCell0 = mtRayCellBoundaryFromNucleus(cell, theta0, rNuc0);
  const startFrac = lerp(p.mtStartFracMin, p.mtStartFracMax, next());
  const startDist = startFrac * rCell0;
  let startX = surfX + outward[0] * startDist;
  let startY = surfY + outward[1] * startDist;
  let startZ = surfZ + outward[2] * startDist;

  // Small random XY offset (settable, up to mtStartOffsetXY) on top of the
  // surface-derived position -- otherwise every microtubule starting near
  // the same (theta0,psi) still emerges from EXACTLY the same point, so nothing
  // can ever cross above/below another one right at the nucleus. Applied in
  // the WORLD xy plane (not the nucleus's own local surface tangent) since
  // it only needs to scatter starts apart, not stay tangent to the ellipsoid.
  const offR = next() * Math.max(0, p.mtStartOffsetXY);
  const offAng = next() * Math.PI * 2;
  startX += Math.cos(offAng) * offR;
  startY += Math.sin(offAng) * offR;

  const jitterRad = (next() * 2 - 1) * p.mtEndJitterDeg * Math.PI / 180;
  const thetaEnd = theta0 + jitterRad;
  const rNucEnd = mtNucleusRadiusAt(cell, thetaEnd);
  const rCell1 = mtRayCellBoundaryFromNucleus(cell, thetaEnd, rNucEnd);
  const endFrac = lerp(p.mtEndFracMin, p.mtEndFracMax, next());
  const endR = Math.max(0, rCell1 - endFrac * rCell1);
  const endX = cell.nucOffX + Math.cos(thetaEnd) * endR;
  const endY = cell.nucOffY + Math.sin(thetaEnd) * endR;

  // z is tracked as a FRACTION of the LOCAL cytoplasm ceiling (0 = floor,
  // 1 = the actual rendered height right at that (x,y)), not an absolute
  // height driven by a fixed-rate linear walk -- the old approach decreased z
  // at a constant rate per path-STEP while the true ceiling (sampleCytoMeshHeight)
  // stays near full nucleus height for a while, then drops over a short band,
  // then goes flat at the rim for the remaining, typically much longer,
  // distance -- so a path's own z was usually still near its starting value
  // once (x,y) had already crossed that drop-off band, forcing
  // mtClampIntoCytoplasm to snap it down hard and hold it near the floor for
  // the rest of the path (a real, reported "hard snap, not a slope"
  // artifact). Interpolating a FRACTION and re-deriving z from the LOCAL
  // ceiling at every point instead means z always rides whatever slope is
  // actually there, by construction.
  const startTopH = Math.max(0, sampleCytoMeshHeight(cell, p, startX, startY));
  const fracStart = startTopH > 1e-9 ? Math.min(1, Math.max(0, startZ / startTopH)) : 0;
  const fracEnd = next(); // endZ = fracEnd * (local ceiling at endX,endY) -- see the per-point loop below

  const dx = endX - startX, dy = endY - startY;
  const straightLen = Math.hypot(dx, dy);
  const turnMag = Math.max(0, p.mtWobbleTurn);
  // Path length BUDGET -- how far the walk is allowed to wander before it
  // has to arrive -- deliberately larger than the straight-line distance
  // (settable multiplier) so there's real room for loops/U-bends, not just a
  // perturbed straight line. Skipped entirely at turnMag=0 (no wobble asked
  // for): with a constant heading the raw walk below is already dead
  // straight, and budgeting extra length there would just make it overshoot
  // in the initial direction and bend backwards to correct, a worse-looking
  // result than simply not padding the length at all.
  const pathBudget = turnMag > 0 ? Math.max(straightLen, straightLen * Math.max(1, p.mtWobbleFactor)) : straightLen;
  // Capped so a small mtStepLen combined with a large mtWobbleFactor can't
  // blow a single microtubule up to an unbounded point count on its own
  // (independent of mtDensity/mtMinSeparation, which the collision-pass caps
  // above already cover) -- stepLen widens past what the slider asked for
  // once the cap binds, a graceful degradation rather than a hang.
  const steps = Math.min(MT_MAX_STEPS_PER_MT, Math.max(4, Math.round(pathBudget / Math.max(0.02, p.mtStepLen))));
  const stepLen = pathBudget / steps;

  // Free 2D correlated random walk in XY only (a persistent random walk: each
  // step's heading is the OLD heading nudged by a random perpendicular kick
  // of magnitude turnMag, then renormalized), started pointing at the target
  // but otherwise completely untethered from it -- this is what actually
  // allows real U-bends/loops in the lateral wandering that's actually
  // visible in the oblique view (a heading-based perpendicular-offset
  // wobble, tried first, can only ever deviate sideways from a straight line
  // and can never double back). z no longer rides along with this walk --
  // see the fraction-of-local-ceiling model above -- so there's no need for
  // the 3D cross-product basis (mtCross3/mtNormalize3) here, just the single
  // perpendicular direction a 2D heading has. Tracked relative to the start
  // (raw[0] = origin) so the drift-correction below can be a simple, exact
  // vector subtraction.
  //
  // `smoothLp` is a PERSISTENCE LENGTH (um, `p.mtSmoothLen`, "Path
  // smoothing") shared by this walk's own heading noise and z's fraction
  // noise below -- floored at stepLen (the slider's own minimum is now 1 um,
  // comfortably above any realistic Step length, so this floor is mostly a
  // safety net rather than something the slider's own minimum normally hits)
  // so smoothLp never drops below stepLen, `headingKickScale` never exceeds
  // 1, i.e. never an ever-shrinking-stepLen blow-up. Applying it here (not
  // just to z) is the fix for a real,
  // reported asymmetry: at a SMALLER Step length (more, smaller steps over
  // the same real distance), a flat per-step kick magnitude -- independent
  // of how physically long each step is -- makes the path accumulate MORE
  // total wander over a fixed real distance (each step contributes a fresh,
  // same-size kick; more steps per um means more kicks per um). Scaling the
  // kick's own magnitude by sqrt(stepLen/smoothLp) is the standard
  // worm-like-chain relation (tangent-angle variance grows linearly with arc
  // length, i.e. per-step variance must scale with stepLen for the total to
  // come out step-length-INVARIANT over a fixed real distance) -- this is
  // what actually keeps "how wiggly per um" constant regardless of Step
  // length, rather than just capping the worst case (below).
  const smoothLp = Math.max(stepLen, p.mtSmoothLen);
  const headingKickScale = Math.sqrt(stepLen / smoothLp);
  // Minimum radius of curvature (`p.mtMinTurnRadius`, "Min turn radius", um
  // -- 0 disables it): heading_new = normalize(heading_old + kick*perp)
  // rotates heading by exactly atan2(kick, 1) -- unbounded as the kick
  // grows, which at a high Wobble turn strength let a single step bend by
  // tens of degrees, several such steps in a row (a real possibility, not
  // just a tail-risk edge case, since a fresh kick is drawn every step)
  // producing an implausibly tight hairpin/near-180 kink (a real, reported
  // artifact) with no relation to how far the step actually moved. Recast as
  // an explicit ROTATION by that same angle (mathematically identical to the
  // old vector-add-then-renormalize for any UNCLAMPED angle -- heading/perp
  // are orthonormal, so the new vector's angle off heading is exactly
  // atan2(kick,1) either way) makes the angle a first-class value that CAN
  // be clamped: `stepTurnCap`, this path's own per-step ceiling, is derived
  // from the minimum radius via the chord-angle relation for a step of
  // length stepLen -- expressed as a radius (um) rather than a flat
  // degrees-per-step limit specifically so it stays a genuine geometric
  // floor regardless of Step length. Deliberately independent of
  // turnMag/mtWobbleTurn, per the report asking for a curvature floor
  // "without changing wobble turn strength" -- turnMag still controls how
  // OFTEN/how close to the cap a path turns; this only ever pulls in the
  // tail where an uncapped kick would have exceeded it.
  //
  // THIS CAP ON ITS OWN IS NOT ENOUGH, measured directly: ~16% of a real
  // run's interior turns still violated the intended radius, some down to
  // ~0.001 um, i.e. a near-total reversal -- because it only bounds the RAW
  // walk's own curvature, and the Brownian-bridge drift correction added
  // below is a SEPARATE, additive, fixed-direction pull with no curvature
  // bound of its own; the REALIZED point-to-point direction (raw step +
  // drift increment, which grows from 0 at the start to its full value at
  // the end) can still bend sharply wherever the two don't point the same
  // way, however gently the raw walk itself curves. `mtEnforceMinTurnRadius`
  // (below, run on the FINAL rendered x,y once drift is already added) is
  // the fix that actually catches this -- this raw-walk cap is kept anyway
  // as a cheap first line of defence that reduces how much work that pass
  // has to do.
  const stepTurnCap = p.mtMinTurnRadius > 0 ? 2 * Math.asin(Math.min(1, stepLen / (2 * p.mtMinTurnRadius))) : Math.PI;
  let heading = straightLen > 1e-9 ? [dx / straightLen, dy / straightLen] : [1, 0];
  const raw = [[0, 0]];
  for (let i = 1; i <= steps; i++) {
    if (turnMag > 0) {
      const n1 = next() * 2 - 1;
      let turnAngle = Math.atan2(turnMag * headingKickScale * n1, 1);
      if (turnAngle > stepTurnCap) turnAngle = stepTurnCap;
      else if (turnAngle < -stepTurnCap) turnAngle = -stepTurnCap;
      const cosT = Math.cos(turnAngle), sinT = Math.sin(turnAngle);
      heading = [heading[0] * cosT - heading[1] * sinT, heading[0] * sinT + heading[1] * cosT];
    }
    const prev = raw[i - 1];
    raw.push([prev[0] + heading[0] * stepLen, prev[1] + heading[1] * stepLen]);
  }

  // Brownian-bridge-style correction: the free walk above generally does NOT
  // land exactly on (dx,dy), so distribute the shortfall across every point
  // via a smoothstep taper (0 at the start, 1 at the end) -- a smooth,
  // low-frequency correction that pins both endpoints EXACTLY without
  // fighting the walk's own local wiggles/loops, the same "de-trend a free
  // random walk" construction a Brownian bridge uses.
  const driftX = dx - raw[steps][0], driftY = dy - raw[steps][1];

  // Finalize xy BEFORE computing z (which needs to sample the ceiling/
  // nucleus footprint at the path's own ACTUAL, post-correction position,
  // not a soon-to-be-adjusted one) -- build the plain x,y first, run the
  // minimum-turn-radius corrector on them, then compute z from the result.
  const ptsXY = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const s = t * t * (3 - 2 * t); // smoothstep
    ptsXY[i] = { x: startX + raw[i][0] + driftX * s, y: startY + raw[i][1] + driftY * s };
  }
  mtEnforceMinTurnRadius(ptsXY, p.mtMinTurnRadius);

  // Fraction noise: a small persistent-random-walk scalar process (same idea
  // as the heading kick above, just applied to a scalar instead of a
  // direction), tapered to 0 at both ends (sin(pi*t), 0 at t=0 and t=1) so it
  // never disturbs the exact fracStart/fracEnd endpoints -- gives the
  // interior a bit of natural roughness without losing the slope-following
  // behaviour the interpolation itself provides. Amplitude scaled off the
  // existing mtWobbleTurn slider (down-weighted) rather than a new one.
  //
  // Its CORRELATION LENGTH is the SAME `smoothLp` persistence length the
  // heading kick above now uses (both ultimately driven by the one Path
  // smoothing slider) -- fixing a real, reported asymmetry: z used to be a
  // flat AR(1) with a fixed 0.8-per-STEP coefficient, decorrelating after a
  // roughly constant NUMBER of steps regardless of how physically long each
  // step was, so at a small Step length (more steps per um) the same
  // few-step correlation window covered much less real distance and z
  // visibly got MORE jagged per um exactly where a smaller Step length was
  // chosen to look more realistic -- while xy's heading walk, despite having
  // the same "redraws every step" shape, looked far less affected, because
  // (before this round) its minimum-curvature clamp already tied a good
  // share of ordinary steps to stepLen via the chord-angle relation; z had
  // no equivalent tie at all. Reparametrized as a proper discretized
  // Ornstein-Uhlenbeck process in ARC LENGTH: decay = exp(-stepLen/corrLen)
  // is the exact relation that keeps the process's correlation length fixed
  // in real um regardless of how many steps a given physical distance is
  // chopped into, and the innovation's own sqrt(1-decay^2) scaling is the
  // standard OU identity that keeps the STATIONARY variance constant
  // regardless of corrLen/stepLen too -- so Path smoothing controls only how
  // quickly the noise wanders (frequency), never how far it wanders
  // (amplitude, still fracNoiseAmp). If corrLen ever drops to stepLen (the
  // floor this shares with headingKickScale above -- not reachable via the
  // slider's own current minimum of 1 um for any realistic Step length, but
  // still a safe floor for a directly-scripted config), decay = exp(-1) =~
  // 0.37 per step -- a deliberately modest, not zero, floor (a flat
  // white-noise reset every single step reads as static, not roughness) that
  // still redraws fast enough to look like "smoothing off".
  const fracNoiseAmp = 0.15 * turnMag;
  const fracNoiseDecay = Math.exp(-stepLen / smoothLp);
  const fracNoiseInnovScale = Math.sqrt(Math.max(0, 1 - fracNoiseDecay * fracNoiseDecay));
  let fracNoise = 0;

  // Over/under-the-nucleus crossings: the fraction-of-local-ceiling model
  // above has no idea the nucleus sits in the way -- left alone, a path
  // whose (x,y) happens to swing over the nucleus footprint (more likely the
  // wider End direction jitter/Wobble settings are pushed, since those are
  // what let an end azimuth land on a genuinely different side of the
  // nucleus than the start) would just interpolate straight through empty
  // fraction-space with no awareness it's passing over solid nucleus volume
  // -- a real, reported gap ("basically no microtubules cross over/under the
  // nucleus"). `goOverNucleus` commits to ONE side for the whole path (not
  // drawn per point) so it doesn't flip-flop -- but it is tied to `cosPsi`
  // (the START point's own hemisphere on the nucleus ellipsoid, drawn
  // above), NOT an independent coin flip. An independent flip was a real,
  // reported bug: the start point sits ON the nucleus surface, so it's
  // still deep inside the blend zone at the very next (i=1) point -- if the
  // flip disagreed with which hemisphere the start actually landed on (e.g.
  // start drawn near the nucleus's underside but the flip said "go over"),
  // the override forced an enormous, instantaneous jump right at the
  // beginning of the path (an "incredibly sudden, fully vertical drop/rise")
  // to reconcile the two. Matching the flip to cosPsi's own sign means the
  // override's very first application already agrees with where the path
  // actually starts. `canGoUnder` requires real clearance beneath the
  // nucleus (basal side) before that option is even offered, overriding a
  // below-the-equator start if there's nowhere to go.
  const nucTopZ = cell.nucZ + cell.nucHeight / 2;
  const nucBottomZ = cell.nucZ - cell.nucHeight / 2;
  const nucClearance = mtNucleusClearance(p);
  const canGoUnder = nucBottomZ - nucClearance > 0;
  const goOverNucleus = !canGoUnder || cosPsi >= 0;

  // The target is a FLOOR (over) / CEILING (under) applied to zNormal, not an
  // independent absolute height -- a real, reported bug in an earlier
  // version blended toward a FIXED value (nucTopZ+clearance regardless of
  // what zNormal happened to be nearby), so once a path left the blend zone
  // it had to snap from that fixed value back to whatever zNormal
  // independently was there -- often very different, since zNormal only
  // depends on the frac interpolation, not on distance from the nucleus, and
  // can be small even where the ceiling is still tall right next to the
  // nucleus. Squeezed into the blend zone's short lateral width, that
  // mismatch showed up as an "incredibly sudden, fully vertical drop."
  // Anchoring to zNormal instead means the target EQUALS zNormal whenever
  // zNormal already clears the nucleus (no intervention, no mismatch to snap
  // back from) and only pulls it toward the clearance boundary when it
  // doesn't -- continuous by construction at the point intervention
  // starts/stops.
  function applyNucleusOverride(x, y, zNormal, ceilH) {
    const blend = mtNucleusFootprintBlend(cell, x, y);
    if (blend <= 0) return zNormal;
    const target = goOverNucleus
      ? Math.max(zNormal, Math.min(ceilH, nucTopZ + nucClearance))
      : Math.min(zNormal, Math.max(0, nucBottomZ - nucClearance));
    return lerp(zNormal, target, blend);
  }

  const pts = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const s = t * t * (3 - 2 * t); // smoothstep
    const x = ptsXY[i].x, y = ptsXY[i].y;
    const ceilH = Math.max(0, sampleCytoMeshHeight(cell, p, x, y));
    let z;
    if (i === 0) {
      // The START point is deliberately excluded from the override -- its z
      // already comes from a specific, meaningful construction (a genuine
      // point on the nucleus ellipsoid's own surface, uniform-on-sphere) that
      // this shouldn't second-guess; applying a floor/ceiling here would
      // pull every "over" start up toward the pole regardless of where psi
      // actually placed it, reintroducing the equatorial-clustering bug that
      // construction exists to avoid (see its own comment above).
      z = fracStart * ceilH;
    } else if (i === steps) {
      // The END point, unlike the start, has NO structural tie to the
      // nucleus -- fracEnd is a free, independent draw -- so there is
      // nothing meaningful about it this override could second-guess. A
      // real, reported bug: leaving it out entirely meant that whenever the
      // path was STILL inside the blend zone right up to the last point
      // (the second-to-last point correctly pulled up/down to clear the
      // nucleus), the true end simply reverted to its own unrelated
      // fracEnd*ceilH with no continuity at all -- the same class of "sudden
      // vertical drop" as the start-point bug, just at the finish instead of
      // the beginning. Routing it through the same override fixes that.
      z = applyNucleusOverride(x, y, fracEnd * ceilH, ceilH);
    } else {
      fracNoise = fracNoiseDecay * fracNoise + fracNoiseAmp * fracNoiseInnovScale * (next() * 2 - 1);
      const frac = Math.min(1, Math.max(0, lerp(fracStart, fracEnd, s) + fracNoise * Math.sin(Math.PI * t)));
      z = applyNucleusOverride(x, y, frac * ceilH, ceilH);
    }
    pts[i] = { x, y, z };
  }

  // Explicit z smoothing (postprocess, on top of the OU reparametrization
  // above) -- a real, separate finding: the OU process keeps fracNoise's
  // STATIONARY amplitude constant regardless of corrLen (Path smoothing),
  // shaping only its FREQUENCY, and its own per-step innovation scale
  // (`fracNoiseInnovScale = sqrt(1-decay^2)`, decay=exp(-stepLen/corrLen))
  // only shrinks as sqrt(stepLen/corrLen) as corrLen grows past stepLen --
  // a SLOW rate. At a small Step length (0.1 um by default) this stays
  // substantial (~26% of the full noise amplitude still injected FRESH
  // every single step) even at the slider's own maximum, since reducing it
  // to, say, 5% would need corrLen roughly 800x stepLen -- far outside any
  // practical slider range. A real, reported "even with very high Path
  // smoothing there is a lot of small jitter" confirmed this directly
  // (measured: z's own local curvature only dropped from a comparable
  // frac-noise-driven contribution, not toward zero, as Path smoothing was
  // pushed to its max). This is NOT the same bug the OU reparametrization
  // fixed (that one was an ASYMMETRY -- z decorrelating faster than xy at a
  // smaller Step length -- not present-but-too-weak smoothing at any Step
  // length) so it needed a separate fix: a plain box-filter moving average
  // directly over z, window size in POINTS derived from the same Path
  // smoothing slider and this path's own stepLen (so it stays a real-um
  // window regardless of Step length, same reasoning as everywhere else in
  // this file), endpoints excluded. This directly suppresses whatever
  // high-frequency content survives the OU process's own frequency-shaping,
  // regardless of how small Step length is set -- complementary to, not a
  // replacement for, the OU fix (which is still what keeps z from getting
  // WORSE at a smaller Step length in the first place).
  const smoothWinPts = Math.max(0, Math.round(p.mtSmoothLen / stepLen));
  if (smoothWinPts > 0) {
    const zOrig = pts.map(pt => pt.z);
    const half = Math.max(1, Math.round(smoothWinPts / 2));
    for (let i = 1; i < pts.length - 1; i++) {
      const lo = Math.max(0, i - half), hi = Math.min(pts.length - 1, i + half);
      let sum = 0, n = 0;
      for (let j = lo; j <= hi; j++) { sum += zOrig[j]; n++; }
      pts[i].z = sum / n;
    }
  }

  // Slope limiter (postprocess, on top of the anchored floor/ceiling above):
  // even anchored to zNormal, a microtubule whose free wobble happens to
  // swing back near the nucleus at a point where the surrounding
  // fraction-based height is naturally low still needs a genuinely large z
  // change concentrated in the short lateral distance the blend zone spans
  // to clear the nucleus in time -- still a real, visible "sudden, near-
  // vertical drop" even though the anchor fix above removed the WORSE
  // fixed-target mismatch case. A few iterative smoothing passes (same
  // "relax toward neighbours" idiom as mtNudgeRoundGrid's own collision
  // resolution elsewhere in this file) cap z's rate of change per step,
  // pulling an over-steep interior point toward whatever range both its
  // neighbours' own slope budgets allow, rather than leaving one segment to
  // absorb an entire height change alone. Endpoints (i=0/i=steps) are
  // excluded, same reasoning as the override above -- their z is a
  // deliberate, meaningful value this shouldn't second-guess.
  //
  // Uses the microtubule's own NOMINAL per-step distance (stepLen), not the
  // REALIZED lateral distance between each specific pair of points, as the
  // reference -- a persistent random walk with real turning can occasionally
  // double back enough that two consecutive points land almost on top of
  // each other in xy (realized lateral distance near 0), and bounding z
  // purely by THAT distance would then force z arbitrarily close to its
  // neighbours too, masking a genuine height change rather than smoothing it
  // out over a few more steps. stepLen is fixed for the whole path, so this
  // is a per-step budget, not a per-realized-distance one. Sweeps back and
  // forth (alternating direction each pass, the standard way to propagate a
  // 1D range constraint faster than always scanning the same direction) so a
  // violation spanning several points converges in a handful of passes
  // rather than needing one pass per point.
  const maxDz = Math.max(1, p.mtMaxZSlope) * stepLen;
  for (let pass = 0; pass < 20; pass++) {
    let changed = false;
    const forward = pass % 2 === 0;
    for (let k = 1; k < pts.length - 1; k++) {
      const i = forward ? k : pts.length - 1 - k;
      const prev = pts[i - 1], cur = pts[i], nxt = pts[i + 1];
      const lo = Math.max(prev.z - maxDz, nxt.z - maxDz);
      const hi = Math.min(prev.z + maxDz, nxt.z + maxDz);
      // lo>hi means the two neighbours are themselves too far apart in z for
      // ANY single value to satisfy both slope budgets at once -- left alone
      // rather than force an arbitrary compromise; the containment clamp
      // below still keeps whatever value remains valid regardless.
      if (lo > hi) continue;
      if (cur.z < lo) { cur.z = lo; changed = true; }
      else if (cur.z > hi) { cur.z = hi; changed = true; }
    }
    if (!changed) break;
  }

  // Truncate at the cell's own OUTER EDGE rather than clamp-and-continue past
  // it. The free walk above has genuine excursions beyond the footprint
  // (more of them the larger Wobble path length/turn strength are set), and
  // the old behaviour -- radially clamping every such point back onto the
  // boundary via mtClampIntoCytoplasm, then letting the walk continue -- let
  // a single microtubule cross the boundary many times, each crossing pulled
  // onto almost exactly the same rim radius: a real, reported "microtubules
  // hugging the outline" artifact, a dense traced ring right at the cell
  // edge that isn't a real cytoskeletal structure. A growing microtubule
  // reaching the cortex is more faithfully represented as simply ending
  // there than as bouncing along it. Scans forward from the start (built to
  // be interior, on/near the nucleus surface) so the FIRST exit ends the
  // path outright -- once truncated, nothing further is generated, so there
  // is no way for a second or third crossing to build up a hugging run. Does
  // NOT apply to the nucleus boundary (still a push-back clamp, below) or the
  // height ceiling (still a clamp, since the fraction-of-local-ceiling z
  // model already keeps height close to right by construction) -- this is
  // specifically about the outer footprint, the one boundary the free walk
  // can wander back and forth across on its own.
  let cutLen = pts.length;
  for (let i = 0; i < pts.length; i++) {
    const ang = Math.atan2(pts[i].y, pts[i].x);
    if (Math.hypot(pts[i].x, pts[i].y) > cellRadiusAt(cell, ang)) { cutLen = i; break; }
  }
  pts.length = cutLen;

  // The SURVIVING points -- now guaranteed within the footprint by the
  // truncation above -- must still stay outside the nucleus and below the
  // cytoplasm height-field; the fraction-of-local-ceiling z already keeps
  // height close to right by construction, but this (and the collision-nudge
  // pass elsewhere in this file, which moves points after generation) still
  // need it as a safety net. mtClampIntoCytoplasm's own footprint clamp is a
  // no-op here in the ordinary case (nothing left in `pts` violates it) --
  // kept as-is rather than split out, since the collision-nudge pass still
  // needs the full three-way clamp for points it moves after this point.
  for (let i = 0; i < pts.length; i++) mtClampIntoCytoplasm(cell, p, geom, pts[i]);

  return { pts, priority: next() };
}

// Spatial-hash grid over every point of every microtubule in the cell, bucket
// size = minSep. Two points further apart than minSep in every axis can never
// be within minSep of each other, so checking a point against only its own
// bucket + the 26 neighbouring buckets (searchNeighbors below) is guaranteed
// to find every real conflict -- turning what used to be an O(M^2*S^2)
// all-pairs scan into an O(N) one (N = total points in the cell). This is the
// difference between "fine for a handful of microtubules" and actually
// working at the up-to-2/um^2 density ceiling (mtDensity), which can put
// thousands of points in a single cell -- a real, measured hang (Chromium's
// tab froze mid-drag) at that combination with the original all-pairs
// version is why this exists.
// Plain integer spatial-hash key (Teschner et al.'s constants) rather than a
// string ("ix,iy,iz") -- a numeric key lets the grid be a Map<number,array>
// instead of Map<string,array>, which V8 hashes/compares meaningfully faster
// than a freshly-concatenated string, and this function is called N times
// per round, N being the whole point of this file's own performance pass. A
// hash COLLISION (two different buckets landing on the same key) is
// harmless here -- it only costs a few extra real-distance checks against
// points that turn out not to be nearby, never an incorrect result.
function mtBucketKey(ix, iy, iz) {
  return (Math.imul(ix, 73856093) ^ Math.imul(iy, 19349663) ^ Math.imul(iz, 83492791)) >>> 0;
}
function mtBuildSpatialIndex(mts, cellSize) {
  const grid = new Map();
  for (let i = 0; i < mts.length; i++) {
    const pts = mts[i].pts;
    for (let j = 0; j < pts.length; j++) {
      const pt = pts[j];
      const key = mtBucketKey(Math.floor(pt.x / cellSize), Math.floor(pt.y / cellSize), Math.floor(pt.z / cellSize));
      let arr = grid.get(key);
      if (!arr) { arr = []; grid.set(key, arr); }
      arr.push(i, j); // endpoints (j=0 or pts.length-1) are included as obstacles too, just never the point that gets moved (see the caller's own loop bounds)
    }
  }
  return grid;
}

// One Jacobi-style round (same "every point updates from a single shared
// snapshot" idea as index.html's own relax() for cell packing): scans every
// point once via the grid, accumulates a push-apart vector for each
// NON-endpoint point from every other-microtubule point within minSep found
// in its own 3x3x3 bucket neighbourhood, then applies every push at once and
// re-clamps into the cell volume. Returns false once no point had a
// violation (converged).
//
// The grid turns the TOTAL cost into O(N), but a single bucket can still
// hold many points if a lot of microtubules happen to pass close together
// (routine right near the nucleus, where every path's start end is confined
// to a narrow band by mtStartFracMax) -- checking all of THOSE against each
// other is locally O(bucket^2), the same blow-up the grid exists to avoid,
// just scoped to one crowded neighbourhood instead of the whole cell.
// MT_COLLISION_OP_BUDGET bails out of the CURRENT round once total pairwise
// checks cross a fixed ceiling, same "just cap it" idiom as CHUNK_CAP/
// RENDER_CAP elsewhere in index.html -- a real, measured hang (Chromium
// froze) at mtDensity's own new 2/um^2 ceiling combined with a small
// mtMinSeparation is why this exists, not a hypothetical.
function mtNudgeRoundGrid(cell, p, geom, mts, minSep) {
  const cellSize = Math.max(minSep, 1e-6);
  const grid = mtBuildSpatialIndex(mts, cellSize);
  const pushes = mts.map(m => m.pts.map(() => ({ x: 0, y: 0, z: 0 })));
  let any = false, ops = 0;
  scan: for (let i = 0; i < mts.length; i++) {
    const pts = mts[i].pts;
    for (let j = 1; j < pts.length - 1; j++) { // endpoints never move, so never need a push accumulated
      const pt = pts[j];
      const ix = Math.floor(pt.x / cellSize), iy = Math.floor(pt.y / cellSize), iz = Math.floor(pt.z / cellSize);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const arr = grid.get(mtBucketKey(ix + dx, iy + dy, iz + dz));
        if (!arr) continue;
        for (let k = 0; k < arr.length; k += 2) {
          const i2 = arr[k], j2 = arr[k + 1];
          if (i2 === i) continue; // same microtubule is never a conflict target
          if (++ops > MT_COLLISION_OP_BUDGET) break scan;
          const p2 = mts[i2].pts[j2];
          const ddx = pt.x - p2.x, ddy = pt.y - p2.y, ddz = pt.z - p2.z;
          const d = Math.hypot(ddx, ddy, ddz);
          if (d >= minSep) continue;
          any = true;
          let nx, ny, nz;
          if (d < 1e-6) { nx = 1; ny = 0; nz = 0; } else { nx = ddx / d; ny = ddy / d; nz = ddz / d; }
          const mag = (minSep - d) / 2 + 1e-4;
          pushes[i][j].x += nx * mag; pushes[i][j].y += ny * mag; pushes[i][j].z += nz * mag;
        }
      }
    }
  }
  if (!any) return false;
  // Apply whatever was found even if the budget cut the scan short -- a
  // partial round still makes progress; returning `true` below (never
  // `false` on a budget-cut round) means the caller keeps calling this until
  // MT_NUDGE_ROUNDS runs out rather than mistaking a cut-short round for
  // convergence.
  for (let i = 0; i < mts.length; i++) {
    const pts = mts[i].pts;
    for (let j = 1; j < pts.length - 1; j++) {
      const push = pushes[i][j];
      if (push.x === 0 && push.y === 0 && push.z === 0) continue;
      pts[j].x += push.x; pts[j].y += push.y; pts[j].z += push.z;
      mtClampIntoCytoplasm(cell, p, geom, pts[j]);
    }
  }
  return true;
}

// Step 7: within ONE cell's own microtubule set only (each starts near its
// own nucleus and ends near its own cell edge, so it never leaves that cell's
// territory -- no cross-chunk/spatial-hash checking needed, just within-cell).
// Iteratively nudges apart (grid-accelerated, see mtNudgeRoundGrid) any pair
// closer than mtMinSeparation in 3D (a purely-XY pair separated in Z is not a
// violation -- "ok to overlap in 2D" per the design). After MT_NUDGE_ROUNDS,
// any microtubule that STILL has a violating point gets its lower-priority
// member regenerated from scratch (a fresh hash-stream attempt) -- bounded to
// MT_RESAMPLE_ROUNDS so a pathological configuration can't loop indefinitely;
// a residual violation past that is a documented best-effort limit, same
// spirit as index.html's own prune() PRUNE_ROUNDS cap.
function mtResolveCollisions(seed, cx, cy, cell, p, mts, geom) {
  const minSep = Math.max(0, p.mtMinSeparation);
  if (minSep <= 0 || mts.length < 2) return;
  let totalPoints = 0;
  for (const m of mts) totalPoints += m.pts.length;
  if (totalPoints > MT_COLLISION_MAX_TOTAL_POINTS) return; // see MT_COLLISION_MAX_TOTAL_POINTS's own comment

  function nudgePass(rounds) {
    for (let round = 0; round < rounds; round++) {
      if (!mtNudgeRoundGrid(cell, p, geom, mts, minSep)) return true;
    }
    return false;
  }

  if (nudgePass(MT_NUDGE_ROUNDS)) return;

  for (let resampleRound = 0; resampleRound < MT_RESAMPLE_ROUNDS; resampleRound++) {
    // One more grid pass, used only to find which microtubules still have a
    // violating point (not to move anything yet) -- same O(N) scan mtNudgeRoundGrid
    // already does, just reporting per-microtubule instead of applying a push.
    const cellSize = Math.max(minSep, 1e-6);
    const grid = mtBuildSpatialIndex(mts, cellSize);
    const violating = new Set();
    for (let i = 0; i < mts.length; i++) {
      const pts = mts[i].pts;
      outer: for (let j = 0; j < pts.length; j++) {
        const pt = pts[j];
        const ix = Math.floor(pt.x / cellSize), iy = Math.floor(pt.y / cellSize), iz = Math.floor(pt.z / cellSize);
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
          const arr = grid.get(mtBucketKey(ix + dx, iy + dy, iz + dz));
          if (!arr) continue;
          for (let k = 0; k < arr.length; k += 2) {
            const i2 = arr[k], j2 = arr[k + 1];
            if (i2 === i) continue;
            const p2 = mts[i2].pts[j2];
            if (Math.hypot(pt.x - p2.x, pt.y - p2.y, pt.z - p2.z) < minSep) { violating.add(i); break outer; }
          }
        }
      }
    }
    if (violating.size === 0) return;
    // Resample the lower half (by priority) of the still-violating set --
    // regenerating every violator at once risks two violators swapping
    // places forever; only ever replacing the losing half converges.
    const sorted = [...violating].sort((a, b) => mts[a].priority - mts[b].priority);
    const loserCount = Math.max(1, Math.ceil(sorted.length / 2));
    for (let k = 0; k < loserCount; k++) {
      const idx = sorted[k];
      mts[idx] = mtGenerateOne(seed, cx, cy, idx, resampleRound + 1, cell, p, geom);
    }
    if (nudgePass(MT_POST_RESAMPLE_NUDGE_ROUNDS)) return;
  }
  // Still-violating microtubules after this point are left as-is --
  // best-effort, see this function's own header comment.
}

// Final-result cache: the same "recomputed from scratch every draw() call
// made panning/dragging incredibly slow" lesson index.html's own cytoCache
// comment documents applies even more here -- generation plus the O(M^2*S^2)
// collision-resolution pass is real work, and draw() calls buildMicrotubulesForCell
// for every visible cell on every single pointermove/wheel frame, not just
// once per parameter change. Keyed on chunk id, invalidated by a signature
// over the seed, the cell's own resolved shape (mtCellShapeSig) and every
// mt*-prefixed control -- a pan/zoom/tilt drag that doesn't touch any of
// those hits this cache for every already-seen chunk, exactly like
// getCytoGeometry does for the cell body/nucleus mesh.
const MT_RESULT_CACHE_MAX = 6000;
const mtResultCache = new Map();
function mtResultSig(seed, cell, p) {
  // p.nucMargin isn't an mt*-prefixed control, but it feeds the cytoplasm
  // height field (via cytoHeightAt, MODULE: index.html) that
  // sampleCytoMeshHeight reads AND the over/under-nucleus clearance
  // (mtNucleusClearance) directly -- a real, previously-latent staleness gap:
  // without it here, changing nucMargin alone would silently keep serving a
  // cached microtubule set built against the old clearance/height field.
  return [seed, mtCellShapeSig(cell), p.mtDensity, p.mtStartFracMin, p.mtStartFracMax,
    p.mtStartOffsetXY, p.mtEndFracMin, p.mtEndFracMax, p.mtEndJitterDeg, p.mtWobbleTurn,
    p.mtWobbleFactor, p.mtStepLen, p.mtSmoothLen, p.mtMinTurnRadius, p.mtMinSeparation, p.mtMaxZSlope, p.nucMargin, p.cytoMaxSlope, p.cytoDomeSlope, p.cytoSmoothPasses, p.cytoRings, p.cytoTheta].join('|');
}

// The one entry point index.html's draw() calls: builds every microtubule for
// one cell (already collision-resolved) and returns their point arrays
// ({x,y,z}, cell-local frame) ready for localToWorld()+project().
function buildMicrotubulesForCell(seed, cx, cy, cell, p) {
  if (!cell || !cell.present) return [];
  const key = cx + ',' + cy;
  const sig = mtResultSig(seed, cell, p);
  const cached = mtResultCache.get(key);
  if (cached && cached.sig === sig) return cached.paths;
  if (mtResultCache.size > MT_RESULT_CACHE_MAX) mtResultCache.clear();

  const geom = getMtCellGeometry(cell);
  // Capped independently of MT_COLLISION_MAX_TOTAL_POINTS -- that one only
  // protects the collision PASS; without this, a very large/blobby cell at
  // mtDensity's own 2/um^2 ceiling would still GENERATE thousands of
  // multi-hundred-point paths (and index.html's draw() would still have to
  // stroke all of them) even with collision-checking skipped.
  const count = Math.min(MT_MAX_PER_CELL, mtCountForCell(seed, cx, cy, cell, p, geom));
  let paths = [];
  if (count > 0) {
    const mts = [];
    for (let i = 0; i < count; i++) mts.push(mtGenerateOne(seed, cx, cy, i, 0, cell, p, geom));
    mtResolveCollisions(seed, cx, cy, cell, p, mts, geom);
    paths = mts.map(m => m.pts);
    // Final realized-slope safety net, run AFTER collision resolution -- see
    // mtLimitZSlopeRealized's own comment for why this has to happen here,
    // not just inside mtGenerateOne. Re-clamped afterward (same three-way
    // containment check every other z-adjusting pass in this file re-applies)
    // since pulling z toward a neighbour's own value is not guaranteed to
    // land inside THIS point's own footprint/nucleus/height bounds.
    for (const pts of paths) {
      mtTrimSteepEnds(pts, p.mtMaxZSlope);
      mtLimitZSlopeRealized(pts, p.mtMaxZSlope);
      for (let i = 0; i < pts.length; i++) mtClampIntoCytoplasm(cell, p, geom, pts[i]);
    }
  }
  mtResultCache.set(key, { sig, paths });
  return paths;
}

// ---- Surface labels (nanobody/antibody + dye) --------------------------------
// Each centerline is treated as a MT_RADIUS_NM cylinder carrying the 13_3
// protofilament lattice. Per lattice site: attachment point on the surface ->
// binder tip (stalk of MT_BINDER_NM, radially outward) -> dye, displaced from
// the tip by displaceByLinker (uniform in volume, ported from webSMLM.html's
// NUP model). Full-network output can reach millions of points, so this is
// only ever called on a short window (see the debug preview in index.html).
const MT_RADIUS_NM = 12.5;
const MT_N_PROTOFILAMENTS = 13;
const MT_DIMER_NM = 8;
const MT_LATTICE_START = 3;       // 13_3 lattice: neighbouring protofilaments stagger by 3/13 dimer
const MT_BINDER_NM = 12;
const MT_LINKER_MIN_NM = 2;
const MT_LINKER_MAX_NM = 5;
const MT_LABEL_STREAM_BASE = 7777777;

function mtDisplaceByLinker(p, minNm, maxNm, rng) {
  const u = rng() * 2 - 1, phi = rng() * 2 * Math.PI, sn = Math.sqrt(1 - u * u);
  const r = Math.cbrt(Math.pow(minNm, 3) + (Math.pow(maxNm, 3) - Math.pow(minNm, 3)) * rng());
  return { a: p.a + r * sn * Math.cos(phi), b: p.b + r * sn * Math.sin(phi), c: p.c + r * u };
}

// pts: centerline [{x,y,z}] (um). opts: {startUm, lenUm, efficiency, binderNm,
// linkerMinNm, linkerMaxNm}. rng: () => [0,1). Returns [{x,y,z (dye, um),
// ax..az (attachment), bx..bz (binder tip), k (protofilament), s (nm along the
// window), att/tip/dye: {u,v,s} offsets in nm in the axis frame}].
function buildMicrotubuleLabelPoints(pts, opts, rng) {
  const o = opts || {};
  const binderNm = o.binderNm != null ? o.binderNm : MT_BINDER_NM;
  const lMin = o.linkerMinNm != null ? o.linkerMinNm : MT_LINKER_MIN_NM;
  const lMax = o.linkerMaxNm != null ? o.linkerMaxNm : MT_LINKER_MAX_NM;
  const eff = o.efficiency != null ? o.efficiency : 1;
  const n = pts.length;
  if (n < 2) return [];

  // Per-segment tangent + parallel-transported normal frame, cumulative length (um).
  const cum = [0], T = [], U = [], V = [];
  let prevU = null;
  for (let i = 0; i < n - 1; i++) {
    const dx = pts[i + 1].x - pts[i].x, dy = pts[i + 1].y - pts[i].y, dz = pts[i + 1].z - pts[i].z;
    const len = Math.hypot(dx, dy, dz);
    cum.push(cum[i] + len);
    let t = len > 1e-12 ? [dx / len, dy / len, dz / len] : (T.length ? T[i - 1] : [1, 0, 0]);
    let u;
    if (!prevU) {
      u = Math.abs(t[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    } else u = prevU;
    const d = u[0] * t[0] + u[1] * t[1] + u[2] * t[2];
    u = [u[0] - d * t[0], u[1] - d * t[1], u[2] - d * t[2]];
    const ul = Math.hypot(u[0], u[1], u[2]) || 1;
    u = [u[0] / ul, u[1] / ul, u[2] / ul];
    T.push(t); U.push(u);
    V.push([t[1] * u[2] - t[2] * u[1], t[2] * u[0] - t[0] * u[2], t[0] * u[1] - t[1] * u[0]]);
    prevU = u;
  }
  const total = cum[n - 1];
  const lenUm = Math.min(o.lenUm != null ? o.lenUm : 1, total);
  const startUm = o.startUm != null ? Math.min(Math.max(0, o.startUm), total - lenUm) : (total - lenUm) / 2;
  const phase = o.phase != null ? o.phase : rng() * 2 * Math.PI;
  const nm = 1e-3;
  const out = [];
  for (let k = 0; k < MT_N_PROTOFILAMENTS; k++) {
    const th = phase + k * 2 * Math.PI / MT_N_PROTOFILAMENTS;
    const ct = Math.cos(th), st = Math.sin(th);
    const off = (k * MT_LATTICE_START * MT_DIMER_NM / MT_N_PROTOFILAMENTS) % MT_DIMER_NM;
    let seg = 0;
    for (let sNm = off; sNm < lenUm * 1000 - 1e-9; sNm += MT_DIMER_NM) {
      if (rng() >= eff) continue;
      const S = startUm + sNm * nm;
      while (seg < n - 2 && cum[seg + 1] < S) seg++;
      const t = T[seg], u = U[seg], v = V[seg], f = S - cum[seg];
      const cx = pts[seg].x + t[0] * f, cy = pts[seg].y + t[1] * f, cz = pts[seg].z + t[2] * f;
      const r = [ct * u[0] + st * v[0], ct * u[1] + st * v[1], ct * u[2] + st * v[2]];
      const R = MT_RADIUS_NM * nm, B = (MT_RADIUS_NM + binderNm) * nm;
      const a = [cx + r[0] * R, cy + r[1] * R, cz + r[2] * R];
      const b = [cx + r[0] * B, cy + r[1] * B, cz + r[2] * B];
      const dye = mtDisplaceByLinker({ a: b[0], b: b[1], c: b[2] }, lMin * nm, lMax * nm, rng);
      const dv = [dye.a - cx, dye.b - cy, dye.c - cz];
      out.push({
        x: dye.a, y: dye.b, z: dye.c,
        ax: a[0], ay: a[1], az: a[2],
        bx: b[0], by: b[1], bz: b[2],
        k, s: sNm,
        att: { u: R / nm * ct, v: R / nm * st, s: sNm },
        tip: { u: B / nm * ct, v: B / nm * st, s: sNm },
        dye: {
          u: (dv[0] * u[0] + dv[1] * u[1] + dv[2] * u[2]) / nm,
          v: (dv[0] * v[0] + dv[1] * v[1] + dv[2] * v[2]) / nm,
          s: sNm + (dv[0] * t[0] + dv[1] * t[1] + dv[2] * t[2]) / nm,
        },
      });
    }
  }
  return out;
}
};
