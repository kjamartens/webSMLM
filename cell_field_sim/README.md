# Cell field simulation — exploration

Standalone prototype, deliberately kept out of `webSMLM.html` while the approach is
still being explored. Everything for this lives in this folder; nothing here is
wired into the main app.

Open `index.html` directly in a browser — no build step, no server.

## The idea being tested

A continuous, panable field of cells (later: populated with microtubules, lamins,
mitochondria, NUPs, ...) where **a structure's position is a pure function of a
world coordinate, not a stored array**. That means:

- No upper bound on field size — panning to `(1_000_000, -400_000)` µm costs the
  same as panning to the origin.
- The same seed always regenerates the same structure at the same coordinate, so
  revisiting a location after panning away reproduces it exactly (needed for both
  interactive exploration and, later, drift-corrected re-analysis of a
  demoCam_SMLM_MM-style simulated acquisition, where the sample is fixed and the
  stage/FOV moves under it).
- Nothing needs to be precomputed for the whole field — only the chunks currently
  intersecting the viewport (+1 chunk margin) are ever evaluated.

## Method

**Hash**: `pcg4d`, ported line-for-line from webSMLM's `makeSimNoiseRng` /
demoCam_SMLM_MM's `sim::Pcg4d` (`SMLMCounterRng.h`) — the same counter/address-based
hash both codebases already use for camera noise. A draw depends only on its own
`(seed, a, b, c)` address, never on how many draws came before it, which is exactly
what "seed → coordinate" placement needs. Reusing it (rather than a sequential PRNG
like `mulberry32`) keeps this exploration on the same primitive as the rest of the
simulation stack, in case it's worth porting back later.

**Cell placement**: the world is divided into a square chunk grid (`chunkSize`, µm).
Each chunk gets **at most one** candidate cell centre, placed by jittering within a
`jitter` fraction of the chunk (`rawCandidate()` in `index.html`) — the standard
"stratified / jittered grid" trick for cheap, roughly blue-noise point fields used in
procedural generation (one candidate per grid cell, then jitter it). `jitter < 1`
guarantees a minimum spacing of `chunkSize*(1-jitter)` between any two
adjacent-chunk centres; `jitter = 1` allows them to touch. This is weaker than a true
Poisson-disk process, but it's O(1) per chunk with no global state, which is the
actual requirement here — a true Bridson/Poisson-disk sampler needs the whole
previously-placed point set to reject a new candidate against, which doesn't compose
with "any chunk, any order, any thread."

An `occupancy`/`density` draw independently decides whether a chunk's candidate is
even realized as a cell, so the field can be sparser than "one cell per chunk"
without changing the spacing statistics of the cells that *are* present.

**Cell footprint**: an ellipse (`semiMajor`/`semiMinor`, random short/long ratio,
random rotation) whose *geometric mean* diameter — not either single axis — is drawn
from `cellDiamMin`/`cellDiamMax`, so "diameter ~25-30 µm" stays true regardless of how
elongated the cell ends up. That ellipse is then modulated by a few angular harmonics
(`cellRadiusAt()`, amplitude set by `cellBlob`, slider goes to 4) to get a wobbly,
non-convex outline closer to a real confluent-culture cell than a plain ellipse — see
the reference phase-contrast images this round was built against.

**Keeping "diameter" meaning the same size at any blobbiness took three separate
fixes, not one** — worth knowing before touching this code again:

1. The modulation clamp (`CELL_MOD_MIN`/`CELL_MOD_MAX`) is **symmetric around 1**: an
   earlier asymmetric range visibly enlarged the average cell as blobbiness went up
   (more clamp headroom above 1 than below biases the mean radius upward once the
   clamp saturates).
2. Even a perfectly symmetric clamp doesn't fix it on its own — adding *any* angular
   wobble increases enclosed area on average (Jensen's inequality: area is a sum of
   r², and distinct-frequency cosines are exactly orthogonal over a full period, so
   the pre-clamp E[mod²] = 1 + 0.5·Σ(harmAmp²) grows with blobbiness regardless of
   clamp symmetry). `rawCandidate()` corrects for this **numerically per cell**, not
   analytically: an analytic pre-correction based on the *unclamped* variance was
   tried first and failed once the clamp started saturating (routine by blob≈1.5),
   because a saturated signal's true mean-square looks nothing like the unclamped
   formula predicts. The actual fix samples the cell's own already-clamped outline
   (48 points) and rescales `semiMajor`/`semiMinor` by the exact ratio needed so the
   sampled mean-square area matches the unmodulated ellipse's — verified numerically
   flat (equivalent diameter within 0.1 µm) across blob 0–4 in isolation.
3. `envelopNucleus()`'s own lateral correction (below) used to scale the *whole*
   ellipse uniformly to fix any single problem angle — which re-inflated cell size
   exactly in proportion to how deep a lobed outline's dips got, undoing fix #2 by a
   different path (confirmed: with fix #2 alone but the old scale-based envelopment,
   equivalent diameter still grew 27→48 µm from blob 0→3). Fixed by raising a
   *per-cell modulation floor* (`c.modFloor`) instead — see its own comment.

`rOuter` (the radius packing/overlap checks actually use) is pinned to the true
clamped ceiling (`semiMajor * CELL_MOD_MAX`), not the unclamped worst-case sum of
harmonic amplitudes — the unclamped sum routinely overshoots what `cellRadiusAt()`
ever actually draws once the clamp saturates, which was quietly forcing packing to
space cells much further apart than their drawn outlines needed.

**Nucleus**: a true 3D ellipsoid per cell — long axis, short/long ratio, and height
(as a fraction of the long axis) are tunable in the sidebar (defaults: long axis
8–12 µm, ratio 0.6–1, height 0.3–0.6× long axis). **Nucleus size is correlated with
cell size, not an independent draw**: both reuse the exact same underlying hash draw
(`sizeT` in `rawCandidate()`), so a cell at the 80th percentile of `cellDiamMin`–
`cellDiamMax` gets a nucleus at the 80th percentile of `nucLongMin`–`nucLongMax` too —
same percentile, not just "correlated noise." Offset from the cell's own centroid by a
random amount up to `nucOffsetFrac × cell radius` (default 0.1 — subtle, "somewhere in
the middle" rather than dead centre).

`envelopNucleus()` guarantees the full ellipsoid stays inside the cell with at least
`nucMargin` µm of clearance (floored at 1 µm) on every side. Vertically this grows
`height` and re-clamps `nucZ` (cheap, no interaction with blobbiness). Laterally it
samples the nucleus's own boundary (32 points) and, for any point the cell's outline
would otherwise fall short of, raises a **per-cell modulation floor**
(`c.modFloor`, read by `cellRadiusAt()` in place of the global `CELL_MOD_MIN`) just
high enough to cover it — filling in a specific dip rather than scaling
`semiMajor`/`semiMinor` (see point 3 above for why that distinction matters). Only
falls back to whole-cell scaling in the rare case even maximum modulation can't cover
the nucleus. Rendered as a stack of horizontal ellipse slices (`drawNucleus()`), sized
by the true ellipsoid equation at each slice height, coloured by absolute z with a
small blue→red depth ramp (`depthColor()`) — a deliberate echo of the "colour by
depth" convention in this project's own reference imagery, so the 3D-ness reads at a
glance instead of looking like a flat blob.

**3D view**: a simple oblique projection (`project()`, one `tilt` parameter, no
azimuth) — x is untouched, y is foreshortened by `cos(tilt)` and z lifts the point on
screen by `sin(tilt)`. `tilt = 0` is exactly the original top-down view (verified: the
same seed/view produces pixel-identical cell footprints at tilt 0 as the pre-3D pass
did). The cell body itself is drawn as a simple vertical extrusion of its 2D outline
(flat top/bottom, connecting edges at a subsample of vertices) — not a dome, a
known simplification, see below.

**Packing — move, don't shrink** (this round's main ask): cells are never resized to
avoid overlap. Instead, `relax()` runs a fixed number of **Jacobi** iterations: every
candidate looks at its neighbours within a fixed chunk radius (`interactionChunks()`)
and, if overlapping, moves away by half the overlap. Because each candidate only ever
updates *its own* position from a shared same-iteration snapshot, a pair pushes apart
symmetrically with no double-counting and no dependency on Map iteration order —
see the comment above `relax()` for why that's true. If a handful of cells are still
badly stuck after relaxation (e.g. three candidates landing in the same tight
corner), `prune()` drops the lower-priority one of any pair still overlapping beyond a
stricter threshold — "moved apart, or removed entirely," never shrunk. A fixed
per-candidate priority hash means both members of a stuck pair agree on which one
gets dropped regardless of evaluation order. The HUD's "removed (stuck)" count makes
this visible; tune `Density`/`Chunk size`/`Relaxation iterations` until it's ~0 for a
demo, or crank density up to see it kick in.

Two safety nets, both real fixes found by stress-testing, not just precautions:
`CHUNK_CAP` skips packing (falls back to raw, unresolved positions) once the
relaxation window gets too large; `RENDER_CAP` skips drawing entirely once the
*visible* chunk count would mean tens of thousands of detailed polygons per frame —
without it, zooming out far enough genuinely hangs the tab (found by an automated
stress test, not a theoretical concern). Zoom in if either message appears.

**Microtubules** (`microtubules.js`, kept independent of the cell/nucleus/cytoplasm
generator in `index.html` on purpose — see that file's own header comment): rebuilt
from scratch, genuinely 3D and anchored to the real nucleus/cytoplasm geometry rather
than a free-wandering flat random walk. Per cell: a discrete azimuth-weight table is
built once around the nucleus centre (`buildMtDirectionTable()`), weighted by how much
cytoplasm sits in each direction (a ray-marched distance from the nucleus edge to the
cell's own blobby footprint edge) — directions with more room are proportionally more
likely. Each microtubule draws a weighted azimuth from that table, then a genuine 3D
start point on the nucleus ELLIPSOID's own surface (azimuth + a polar angle drawn
uniform-on-sphere, so starts land anywhere over the nucleus — including straight above
or below it — not only on its equatorial ring, which is what an earlier version did and
which turned out to be a real cause of the whole population sitting flat: with every
start pinned near one height, a path only reached the cytoplasm dome's own mid-height
slope if its wobble happened to climb there, routinely leaving that slope's mid-region
uncovered), pushed outward by a settable random distance plus a settable small XY
jitter (so multiple microtubules starting near the same nucleus region don't all emerge
from exactly one point and can cross above/below each other), an end direction (the
start direction plus a settable jitter), and an end position a settable random offset
inside the cell edge. Between them, a free correlated random walk in XY only (heading
nudged by a random perpendicular kick each step, magnitude "Wobble turn strength") is
generated with a settable path-length budget ("Wobble path length ×", larger than the
straight-line distance so there's real room to loop/double back — genuine U-bends, not
just a perturbed straight line), then corrected via a Brownian-bridge-style drift
subtraction so it lands exactly on both endpoints regardless of how far it wandered.

**z is not part of that walk.** It's tracked as a FRACTION of the LOCAL cytoplasm
ceiling (0 = floor, 1 = the actual rendered height right at that point's own x,y),
interpolated smoothly between the start and end fraction (plus a little tapered
AR(1)-ish noise for roughness) and only converted to an absolute height using the local
ceiling at *that specific point*. An earlier version derived z from the same free walk
as x/y — a roughly constant rate of descent per path-step — while the true ceiling
(`sampleCytoMeshHeight()`) stays near full nucleus height for a while, then drops over a
short band near the nucleus, then goes flat at the rim for the remaining, usually much
longer, distance to the cell edge. That mismatch meant a path's z was typically still
near its starting value once its x/y had already crossed the drop-off band, so the
containment clamp had to snap it down hard and hold it near the floor for the rest of
the path — a real, reported "hard transition, not a slope" artifact. Re-deriving z from
the local ceiling at every point instead makes it ride whatever slope is actually there,
by construction, rather than needing a late external clamp to fix it up.

**Over/under-the-nucleus crossings**: the fraction-of-local-ceiling z model has no idea
the nucleus sits in the way, so on its own a path whose x/y happens to swing across the
nucleus's footprint would just interpolate straight through solid nucleus volume — a
real, reported gap ("basically no microtubules cross over/under the nucleus"), since
every microtubule's end azimuth used to be a small jitter away from its own start azimuth
(so paths stayed radial, never actually needing to cross the nucleus at all). Fixed two
ways together: **End direction jitter now ranges up to 180°** (from a 90° ceiling), so an
end can genuinely land on the far side of the nucleus from where a microtubule started,
not just a nearby azimuth on the same side; and any INTERIOR path point (start/end
themselves are left alone) whose x/y falls within, or close to, the nucleus's own lateral
elliptical footprint (`mtNucleusFootprintBlend()`, ignoring z, smoothly blending to 0 by
40% of the nucleus's own radius beyond its edge) has its z pulled toward clearing the
nucleus vertically — above its top or below its bottom, whichever side that specific
microtubule committed to (drawn once per path, not per point, so it doesn't flip-flop
mid-crossing; "below" is only offered when there's real clearance beneath the nucleus to
begin with). The clearance itself is half of the existing **Cell margin around nucleus**
(`nucMargin`) slider — the same margin the cytoplasm height field already guarantees
above the nucleus everywhere nearby — rather than an unrelated hardcoded value, so
raising/lowering that slider scales this too. This needed no change to the containment
clamp below: a point cleared this way already sits outside the nucleus ellipsoid by
construction (z alone puts it there), so the ordinary nucleus-push clamp simply never
triggers for it.

**This crossing rule shipped with its own real bug — a "sudden, fully vertical drop"**
right where a crossing began or ended, found by the user against a real screenshot. Three
compounding causes, all fixed:

1. The clearance target was blended toward as a FIXED absolute height (`nucTopZ +
   clearance`) regardless of what the ordinary fraction-based height happened to be right
   there — so leaving the blend zone meant snapping from that fixed value back to a
   (often very different) natural one, squeezed into the blend zone's short lateral
   width. Fixed by making the clearance a FLOOR (over) / CEILING (under) applied ON TOP
   of the normal fraction-based height instead of an independent target — it now equals
   the normal height whenever that height already clears the nucleus (no snap at all),
   and only pulls upward/downward when it doesn't.
2. Which side to go — `goOverNucleus` — used to be an INDEPENDENT coin flip, disconnected
   from where the microtubule's own start point actually sits on the nucleus ellipsoid. A
   start drawn on the nucleus's underside combined with a flip that said "go over" forced
   an enormous jump right at the path's very first step, while it was still essentially
   sitting on the nucleus surface. Fixed by tying the choice to `cosPsi` — the same draw
   that already placed the start point above or below the nucleus's equator — so the
   override's first application already agrees with where the path actually begins.
3. The END point was excluded from the override entirely (reasoned, at the time, that its
   `fracEnd` draw was "a specific, meaningful construction" not to be second-guessed —
   true for the START, which is genuinely anchored to the nucleus surface, but NOT true
   for the end, which has no structural tie to the nucleus at all). So whenever the path
   was still inside the blend zone right up to its last interior point (correctly pulled
   to clear the nucleus), the true end simply reverted to its own unrelated
   `fracEnd*ceilH` with no continuity — the same class of bug as #2, just at the finish
   instead of the start. Fixed by routing the end point through the same override too.
A final safety net (`MT_MAX_Z_SLOPE`, `microtubules.js`) also caps z's own rate of change
per nominal step (not per REALIZED lateral distance — a persistent random walk can
occasionally double back to near-zero net lateral movement between two consecutive
points, and bounding purely by that would mask a genuine height change rather than
spread it out) via a few back-and-forth smoothing passes, the same "relax toward
neighbours" idiom `mtNudgeRoundGrid()`'s own collision resolution already uses elsewhere
in this file — catches whatever the three fixes above don't, rather than relying on them
being exhaustive. `tools/check_cellfield_microtubules.mjs`'s own step-slope check (see
below) exists specifically to catch a regression in this; its own ratio metric additionally
needs an absolute `|dz|` floor (`STEEP_MIN_DZ`) before flagging a step, since the RATIO
alone spikes on plenty of genuinely tiny, invisible height wiggles whenever the walk's
own realized lateral step happens to land near zero — an expected occurrence, not a
defect, given the z cap above is deliberately sized against the nominal step length
rather than the realized one.

Every point of every microtubule — including its own start/end — is then clamped to
stay inside the cell's actual cytoplasm volume (inside the blobby footprint, below the
cytoplasm height-field, and outside the nucleus ellipsoid), not just checked at the
endpoints; with the fraction-of-local-ceiling z model this is mostly a safety net for
the collision-nudge pass (which moves points after generation) rather than the primary
height mechanism. **The clamp's own height ceiling is the ACTUAL rendered (smoothed)
mesh height (`sampleCytoMeshHeight()`, bilinearly sampled from the same mesh
`buildCytoMesh()` draws), not the raw analytic `cytoHeightAt()`** — the mesh's Laplacian
smoothing pass systematically lowers the sharp nucleus-adjacent dome peak below its raw
analytic value, so clamping against the raw function used to let a point sit above the
surface actually drawn (a real, reported "pokes out of the dome" bug, most visible right
near that peak). **The cell's own outer edge is a cutoff, not a clamp**: the free walk
has genuine excursions past the footprint (more of them the larger Wobble path
length/turn strength are set), and radially clamping every such point back onto the
boundary — the original behaviour — let a single microtubule cross the edge many times,
each crossing pulled onto almost exactly the same rim radius, tracing a dense ring right
along the outline (a real, reported "microtubules hugging the edge" artifact, not a real
cytoskeletal structure). A microtubule is now simply truncated at its first exit from
the footprint, scanning forward from its start — a growing microtubule reaching the
cortex reads as ending there, not bouncing along it, and there's no way for a second or
third crossing to build up a hugging run once the first one already cut the path short.
A per-cell minimum-separation pass (`mtResolveCollisions()`,
spatial-hash-grid accelerated — see its own comment for why) then nudges — and, if still
stuck after a capped number of rounds, regenerates — any microtubule within a settable
minimum 3D distance of another (2D overlap is allowed; only true 3D closeness counts).
Density is expressed as microtubules per µm² of the cell's own footprint area (up to
2/µm²), not a per-cell count.

**Minimum turn radius** (`p.mtMinTurnRadius`, **Min turn radius (µm)**, default 0.35 µm, 0
disables it): the free walk's own per-step heading kick (`heading_new =
normalize(heading_old + kick*perp)`) rotates the heading by `atan2(kick, 1)` — unbounded
as `turnMag` (Wobble turn strength) rises, so a high-wobble path could bend tens of
degrees in a single step, several such steps in a row producing an implausibly tight
hairpin/near-180° kink (a real, reported artifact — "very tight turn radii", not the
intended loopy-but-plausible wobble). Fixed WITHOUT touching Wobble turn strength's own
meaning, as asked: the heading update is recast as an explicit rotation by that same
angle (mathematically identical to the old vector-add-then-renormalize for any unclamped
value), which makes the angle a first-class value that can be clamped — `stepTurnCap` (per
path) derives from the slider via the chord-angle relation for that path's own Step
length, expressed as a radius rather than a flat degrees-per-step number specifically so
it stays a genuine geometric floor regardless of Step length.

**That cap on the raw walk alone was NOT enough — a real, reported "this hasn't happened,
I still see tight turns" follow-up, measured and confirmed**: ~16% of a real run's
interior turns still violated the intended radius, some down to a radius of ~0.001 µm
(essentially a full reversal). The raw-walk cap only bounds the RAW random walk's own
curvature; the Brownian-bridge drift correction added afterward (a separate, additive,
fixed-direction pull growing smoothly from 0 at the start to its full value at the end, no
curvature bound of its own) can still make the REALIZED point-to-point direction bend
sharply wherever the two don't point the same way, however gently the raw walk itself
curves. Fixed by `mtEnforceMinTurnRadius()`, a postprocess that runs on the FINAL rendered
(x,y) — after drift is already added — directly in segment-direction space: it reads off
each segment's own direction and length, clamps the ANGLE step between consecutive
segments to the cap EXACTLY (an explicit rotation, not an iterative pull — a first version
tried softly nudging an over-curved point toward the midpoint of its neighbours, which
needed far too many repeated passes to unwind anything but a mild violation), rebuilds
every point from the fixed start by summing the clamped segment vectors, then re-pins the
true end point with a fresh Brownian-bridge-style drift correction (which can reintroduce
a little curvature of its own, so the whole clamp+redrift cycle repeats a few rounds — each
round's own shortfall is much smaller than the last, so it converges quickly). Reduced the
measured violation rate from 16% to ~0.05% (residual attributed to the independent
containment clamp below, which moves a point's position without any awareness of its
neighbours and runs after this pass).

**Line thickness**: a settable per-cell draw thickness (**Line thickness (µm)**,
`mtLineWidth`, default 0.1 µm — the stroke width sent to Canvas2D is `mtLineWidth *
pxPerUm`, so it's a genuine physical µm thickness that scales with zoom rather than a
fixed screen-pixel width) replaces what used to be a hardcoded `0.1 * pxPerUm` in
`draw()`.

**Path smoothing** (`p.mtSmoothLen`, **Path smoothing (µm)**, default 0.3 µm, up to 3 µm)
controls a PERSISTENCE LENGTH shared by both the xy heading walk and z's fraction noise —
a real, reported follow-up ("it seems to be only in z — is the logic different here?").
It was: z's fraction noise was a flat AR(1) with a fixed 0.8-per-STEP decay coefficient,
decorrelating after a roughly constant NUMBER of steps regardless of how physically long
each step was, so at a smaller Step length (more, smaller steps over the same real
distance) the same few-step correlation window covered much less real distance and z got
visibly MORE jagged per µm exactly where a smaller Step length was chosen to look more
realistic. The xy heading walk has the identical "redraws every step regardless of step
length" shape underneath, but looked far less affected for a specific reason: its
minimum-turn-radius cap (above) already ties a large share of ordinary steps to Step
length via the chord-angle relation, incidentally giving xy a physically-anchored scale
z never had. Both are now reparametrized properly: z as a discretized Ornstein-Uhlenbeck
process in ARC LENGTH (`decay = exp(-stepLen/corrLen)`, the exact relation keeping its
correlation length fixed in real µm regardless of how many steps a given distance is
chopped into, with the innovation's own `sqrt(1-decay²)` scaling — the standard OU
identity — keeping the STATIONARY variance constant so the slider controls only frequency,
never amplitude); xy's heading kick magnitude is scaled by `sqrt(stepLen/corrLen)`, the
worm-like-chain relation for a persistent random walk's tangent angle (variance grows
linearly with arc length, so per-step variance must scale with stepLen for the total over
a fixed real distance to come out step-length-invariant). Both share ONE persistence
length (`corrLen`, floored at Step length so the slider's own minimum, 0, maps to
`corrLen=stepLen` — z's decay ≈ 0.37/step, xy's kick scale = 1 — a deliberately modest,
not zero, "off" floor, not an ever-shrinking blow-up).

Path smoothing and Min turn radius are two DIFFERENT, complementary constraints, not the
same knob under two names: smoothing shapes the TYPICAL, everyday wiggliness (how far the
walk wanders per µm, all cases, including the ones far under the turn-radius cap), while
minimum turn radius is a hard FLOOR only the rare tail ever hits (the worst-case single
step/segment). Verified together (`tools/check_cellfield_microtubules.mjs`'s own scenarios
and a direct single-cell, low-density Playwright screenshot comparison): off (both 0) shows
visibly sharp kinks and zigzags; defaults (0.3 µm / 0.35 µm) already read as gently curving
strands with no visible sharp turns; pushed further (1.5 µm / 0.6 µm) reads as long, sweeping
curves.

**z still had a lot of small-scale jitter even at Path smoothing's own maximum — a real,
reported, separately-investigated follow-up, and a genuinely different bug from the
step-length asymmetry above.** The OU reparametrization keeps fracNoise's STATIONARY
amplitude constant regardless of `corrLen`, shaping only its frequency — but its own
per-step innovation scale (`sqrt(1-decay²)`, `decay=exp(-stepLen/corrLen)`) only shrinks
as `sqrt(stepLen/corrLen)` once `corrLen` exceeds `stepLen` — a SLOW rate. At the default
0.1 µm Step length, measured directly: even at the slider's own maximum (3 µm), ~26% of
the full noise amplitude was still injected fresh every single step (getting that down to
a visually-smooth ~5% would need `corrLen` roughly 800× Step length — tens of µm, far past
any practical slider range). Fixed with a genuinely separate, complementary mechanism: a
plain box-filter MOVING AVERAGE directly over z (window size in points derived from the
same Path smoothing slider and this path's own Step length, so it stays a real-µm window
regardless of Step length; endpoints excluded), applied after the frac/override
computation and before the slope limiter. This directly suppresses whatever
high-frequency content survives the OU process's own frequency-shaping, independent of how
small Step length is set — it does not replace the OU fix (which is still what keeps z
from getting WORSE at a smaller Step length in the first place), it catches what that fix
structurally cannot on its own. Measured: average local z curvature (`|z[i+1]-2z[i]+z[i-1]|`
across a real run) dropped from 0.18 (Path smoothing off) to 0.009 at the slider's maximum
(3 µm) — a ~95% reduction — while the worst single case (the deliberately steep
over/under-nucleus crossing transitions) stayed a comparable, expected magnitude, i.e. this
smooths the everyday wiggle without flattening genuine large crossings.

**Verification**: `tools/check_cellfield_microtubules.mjs` (Playwright) checks three
things independently after generation, across several parameter scenarios including one
with widened End direction jitter/Wobble turn strength specifically to exercise the
over/under-nucleus crossing rule: every point against the same containment bounds the
clamp is supposed to enforce (footprint, nucleus, and `sampleCytoMeshHeight`); a coarse
histogram of z-as-fraction-of-local-ceiling across every point, which should have real
density in the middle third rather than being bimodal (pinned near 0 and 1) — the
signature of the hard-snap artifact above; and, for every consecutive pair of points on a
path, the ratio of |dz| to the lateral distance moved, which should only rarely exceed a
generous multiple — the signature of the "sudden, fully vertical drop" artifact above.
Run it after touching microtubule or cytoplasm-height-field code (`cd tools && node
check_cellfield_microtubules.mjs`; needs `tools/node_modules` — `cd tools && npm
install` once).

**Scale bar**: bottom-right, fixed in screen space (untouched by pan/tilt), picks a
round 1-2-5 length that renders near a target on-screen size at the current zoom
(`drawScaleBar()`). Uses the projection's x-axis scale specifically, which stays
`view.scale * devicePixelRatio` at any tilt (only y is foreshortened) — metrically
correct regardless of the tilt slider.

**Dragging while packing is on**: `relax()` always runs its full configured iteration
count now, including on every `pointermove` during an active drag. It used to drop to
6 iterations while dragging for perceived responsiveness, but that meant packing
looked visibly broken (still-overlapping, under-relaxed cells) for the whole drag and
only "snapped" correct on release — worse than just paying the full cost, which turned
out cheap enough at the chunk counts `CHUNK_CAP` allows through anyway.

## Known limitations / not yet done

- **Wobble path length ×**'s own slider now goes down to 0.9 (from a 1.0 floor), but
  `params()` still floors the EFFECTIVE value at 1 (`Math.max(1, p.mtWobbleFactor)`,
  where the path-length budget is computed) — a value in `[0.9, 1.0)` currently behaves
  identically to 1.0, since the budget can't meaningfully go below the straight-line
  distance without a larger redesign (removing that floor would need pathBudget<straightLen
  to be handled deliberately, not just left to fall out of the math). Not fixed here since
  it wasn't asked for — flagging it so the 0.9–1.0 range doesn't read as broken.
- Microtubule collision avoidance checks point-to-point distance along each path
  (sampled ~one step length apart) rather than true segment-to-segment distance — an
  adequate proxy at the step lengths this is tuned for, not an exact geometric test.
- The minimum-separation pass is deliberately bounded to a small number of rounds
  (`MT_NUDGE_ROUNDS`/`MT_RESAMPLE_ROUNDS`, `microtubules.js`), not run to convergence —
  measured directly: many microtubules' own FIXED start points end up packed tighter
  than a typical minimum separation once density gets even moderately high (every start
  is confined to a thin band right at the nucleus edge), which can never actually
  resolve since neither point is allowed to move; running more rounds against that just
  burns time for no gain. A residual violation is expected background, not a bug to
  chase, at anything beyond a modest density/separation combination.
- Collision resolution is skipped entirely (microtubules stay individually contained,
  just not mutually separated) once a cell's total point count crosses
  `MT_COLLISION_MAX_TOTAL_POINTS` — real headroom is needed for `mtDensity`'s own
  2/µm² ceiling, which can put hundreds of thousands of points in one large cell.
- Direction weighting only accounts for lateral room (cytoplasm thickness in a given
  azimuth), not the cytoplasm's own height at that azimuth — a direction that's wide
  but very thin (close to the tapered cell edge) is weighted the same as an equally
  wide but tall one.
- Over/under-the-nucleus crossings (see above) need **End direction jitter** pushed well
  past its old 90° range to become common — at the default 20° almost every microtubule
  is still structurally radial (start and end azimuths close together) and rarely needs
  to cross the nucleus's footprint at all, so the new clearance rule mostly sits idle
  until a user actually widens the jitter (or pushes Wobble turn strength/path length
  high enough for the free walk to swing across on its own).
- Cell body extrusion is a flat-topped prism, not a dome — real adherent cells are
  much thinner at the edges than the centre. Would need a height *field* over the
  footprint rather than one scalar per cell.
- Packing resolves circle-equivalent (`rOuter`) overlap, not true polygon-vs-polygon
  intersection, so two wobbly outlines can still graze each other slightly even when
  their bounding circles don't overlap — and conversely the relaxation is sometimes
  more conservative than the true outlines would need. Good enough for "close but
  mostly non-overlapping"; not an exact tessellation.
- Relaxation uses a fixed local neighbourhood radius rebuilt from the current
  viewport + a margin, not a globally fixed window — so a cell right at the edge of
  that margin could in principle resolve to a very slightly different position
  depending on how far the margin extends, unlike the raw jittered position (exactly
  reproducible) or the microtubules (exactly reproducible). Not visible in practice
  at the margins currently used (checked: pan-away-and-return reproduces identical
  HUD/positions), but worth knowing if this is pushed further.
- Only one structure type (microtubules) plus now the cell body + nucleus. Lamins,
  mitochondria, NUPs etc. are not started — the intent is that each becomes another
  per-cell generator keyed off the same `(seed, chunkX, chunkY, purpose-channel)`
  address, same pattern as everything else here.
- Not validated against demoCam_SMLM_MM's own structure builders
  (`SMLMStructures.cpp`) or webSMLM's `buildStructure()` in any quantitative way —
  this is a placement-and-navigation prototype, not a density/parameter match yet.

## Why a separate folder / branch

This lives on the `cell-field-simulation` branch, split out from `webSMLM.html`
on purpose — it's an open-ended exploration of a different generation paradigm
(infinite, coordinate-addressed field vs. webSMLM's current finite, whole-FOV
`buildStructure()`), not yet a decided feature. If/when a direction is picked, the
useful pieces (the address-based placement scheme, in particular) can be ported into
either webSMLM's simulator or demoCam_SMLM_MM's `SMLMStructures`, whichever ends up
wanting it.

## Microtubule surface labels (debug preview)

`buildMicrotubuleLabelPoints()` (microtubules.js) treats a centerline as a 25 nm-diameter cylinder
carrying the 13_3 protofilament lattice (13 sites per 8 nm dimer ring, neighbouring protofilaments
staggered by 3/13 dimer). Each site becomes surface attachment point -> binder tip (12 nm radial stalk,
nanobody/antibody) -> dye, displaced from the tip by a uniform-in-volume 2-5 nm linker (ported from
`displaceByLinker()` in webSMLM.html). Constants (`MT_RADIUS_NM`, `MT_BINDER_NM`, `MT_LINKER_*`) sit at
the top of that section. Full-network labeling can reach millions of points, so it is not done in
`draw()`; **Debug: view MT labels** labels only 1 µm of the longest microtubule in view.

## demoCam_SMLM_MM port spec (keep in sync)

[`DEMOCAM_PORT.md`](DEMOCAM_PORT.md) is the hand-off spec for reimplementing this in
`C:\GitHub\demoCam_SMLM_MM` (cells, microtubules, dyes/blinks, dummy XY stage). **While that file
exists, update it in the same commit as any change here that it describes** (hash channels, defaults,
geometry, dye model, function names).
