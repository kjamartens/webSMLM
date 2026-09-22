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
the reference phase-contrast images this round was built against. The modulation
clamp (`CELL_MOD_MIN`/`CELL_MOD_MAX`) is deliberately **symmetric around 1**: an
earlier asymmetric range visibly enlarged the average cell as blobbiness went up (more
clamp headroom above 1 than below biases the mean radius upward once the clamp
saturates, plus adding any symmetric wobble increases enclosed area regardless, by
Jensen's inequality on the r² area integral) — see the comment above `CELL_MOD_MIN`
for the full reasoning. `rOuter` (the radius packing/overlap checks actually use) is
likewise pinned to the true clamped ceiling (`semiMajor * CELL_MOD_MAX`), not the
unclamped worst-case sum of harmonic amplitudes — the unclamped sum routinely
overshoots what `cellRadiusAt()` ever actually draws once the clamp saturates, which
was quietly forcing packing to space cells much further apart than their drawn
outlines needed.

**Nucleus**: a true 3D ellipsoid per cell — long axis, short/long ratio, and height
(as a fraction of the long axis) are tunable in the sidebar (defaults: long axis
8–12 µm, ratio 0.6–1, height 0.3–0.6× long axis). **Nucleus size is correlated with
cell size, not an independent draw**: both reuse the exact same underlying hash draw
(`sizeT` in `rawCandidate()`), so a cell at the 80th percentile of `cellDiamMin`–
`cellDiamMax` gets a nucleus at the 80th percentile of `nucLongMin`–`nucLongMax` too —
same percentile, not just "correlated noise." Offset from the cell's own centroid by a
random amount up to `nucOffsetFrac × cell radius` (default 0.1 — subtle, "somewhere in
the middle" rather than dead centre). `envelopNucleus()` then guarantees the full
ellipsoid stays inside the cell with at least `nucMargin` µm of clearance (floored at
1 µm) on every side, both laterally (grows `semiMajor`/`semiMinor` uniformly, solving
for the exact scale factor needed rather than a worst-case guess — see its own
comment) and vertically (grows `height` and re-clamps `nucZ`) — a big nucleus forces a
bigger cell around it, never the other way round. Rendered as a stack of horizontal
ellipse slices (`drawNucleus()`), sized by the true ellipsoid equation at each slice
height, coloured by absolute z with a small blue→red depth ramp (`depthColor()`) — a
deliberate echo of the "colour by depth" convention in this project's own reference
imagery, so the 3D-ness reads at a glance instead of looking like a flat blob.

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

**Microtubules**: unchanged from the first pass — still the rough, 2D-only
random-walk placeholder (`buildMicrotubule()`), deliberately not touched this round.

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

- Microtubules are still 2D-only and unaware of the new cell footprint shape (they
  use the cell's semi-axes for extent but don't clip to the actual wobbly outline).
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
