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
`jitter` fraction of the chunk (`chunkCandidate()` in `index.html`) — the standard
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

**Overlap avoidance** (`getCell()`): each realized cell looks at its 8 immediate
neighbouring chunks' own (raw) candidates and shrinks its own radius so it doesn't
overlap a present neighbour, clamped to a small minimum. This is one level deep, not
recursive — a neighbour's raw candidate is cheap to recompute and never itself
depends on this cell's shrink decision, so there's no circularity. Gives "close but
mostly non-overlapping" without a hard geometric guarantee, which is what was asked
for; toggle off (`Shrink to avoid overlap`) to see the raw jittered field with real
overlaps.

**Microtubules**: per cell, a small number of gently-curved random walks
(`buildMicrotubule()`), seeded off `(seed, chunkX, chunkY, 100+filamentIndex)` so a
given chunk's filaments are identical every time it's revisited. Currently rough
positioning only, 2D (z is drawn but discarded — kept as a real field rather than
bolted on later, so extending to 3D is a smaller change than it would otherwise be).

## Known limitations / not yet done

- 2D only. z is present in the data shape but unused.
- Only one structure type (microtubules). Lamins, mitochondria, NUPs etc. are not
  started — the intent is that each becomes another per-cell generator keyed off the
  same `(seed, chunkX, chunkY, purpose-channel)` address, same pattern as
  microtubules.
- No cell shape beyond a circle (no nucleus, no membrane irregularity).
- The overlap-avoidance radius shrink doesn't account for microtubules extending
  past the cell radius into a neighbour — filaments can still visually cross a
  neighbouring cell.
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
