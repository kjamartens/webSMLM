# Porting the cell field simulation into demoCam_SMLM_MM

> **Keep this file up to date while it exists.** It is the hand-off spec between the two repos. Whenever
> the cell-field code in `C:\GitHub\websmlm\cell_field_sim\` changes in a way that affects anything
> described here (hash channels, defaults, geometry, dye model, function names), update this file in
> the same commit. Whoever implements the port in `C:\GitHub\demoCam_SMLM_MM` should correct this file
> when reality diverges from it (a decision changed, a section turned out wrong), and mark finished
> milestones in section 11. Delete the file only once the port is complete and `cell_field_sim/README.md`
> plus demoCam's `CLAUDE.md` carry everything that is still true.

**Audience:** a coding agent working in `C:\GitHub\demoCam_SMLM_MM` (C++ Micro-Manager device adapter).
**Goal:** demoCam simulates a *field of cells* (cell body, nucleus, microtubules, fluorophores on the
microtubules) that is effectively infinite, plus a **dummy XY stage** so the field-of-view can be moved
over it, with **blinks generated from the actual dye positions** and rendered by the existing PSF/noise
pipeline.

Read demoCam's `CLAUDE.md` first (property naming convention, the "two draws in one expression" RNG
gotcha, GPU/CPU render paths). Then read the JS source listed in section 1.

---

## 1. Source material (the JS prototype is the reference implementation)

All in `C:\GitHub\websmlm\cell_field_sim\`. Line numbers drift; grep the function names.

| File | What to port |
|---|---|
| `index.html` `pcg4d`, `hashUnit`, `hashStream` | Address-based RNG. Same primitive as demoCam's `Pcg4d` in `SMLMCounterRng.h`. |
| `index.html` `CH`, `rawCandidate`, `cellRadiusAt`, `envelopNucleus`, `cellOutlineLocal`, `nucleusSignedDistLocal` | Per-chunk cell candidate: position, ellipse + angular-harmonic outline, nucleus ellipsoid. |
| `index.html` `cytoHeightAt`, `buildCytoMesh`, `getCytoGeometry`, `sampleCytoMeshHeight`, `smoothCytoGrid` | Cytoplasm height field. Microtubules are clamped against the *smoothed mesh*, not the analytic function, so the mesh must be ported too. |
| `index.html` `buildCandidateMap`, `relax`, `prune`, `interactionChunks` | Packing (cells are moved apart, never shrunk). **See the viewport-dependence trap in 4.3.** |
| `microtubules.js` `buildMicrotubulesForCell` and everything it calls (`mtGenerateOne`, `mtResolveCollisions`, `mtEnforceMinTurnRadius`, `mtClampIntoCytoplasm`, `buildMtDirectionTable`, ...) | Per-cell 3D microtubule centrelines (`{x,y,z}` µm, cell-local frame). |
| `microtubules.js` `buildMicrotubuleLabelPoints`, `MT_*` constants | Lattice site -> binder tip -> dye. Currently only a windowed debug preview in JS. |
| `README.md` in the same folder | The *why* behind every non-obvious decision (blobbiness area correction, z-as-fraction-of-local-ceiling, over/under-nucleus crossing, turn-radius enforcement). Read it before touching an algorithm; do not "simplify" what it says was fixed on purpose. |

The JS header of `microtubules.js` already states it is written to be reimplemented in C++ (plain scalar
math, named constants). The cell generator in `index.html` is the same in spirit.

---

## 2. Architecture

```
                    SharedStageState (process-wide, no MMDevice dependency)
   SMLMDemoXYStage --writes--> x,y target + motion model      zPositionUm (existing)
        (MM::XYStage)                    |                         |
                                         v                         v
  CSMLMDemoCamera live loop / stack worker: per frame ->  stage(x,y) at frame time, z
                                         |
                                         v
      CellFieldSource::EventsForFrame(f, stageXY, fovSize, frameDur, params)
        |  cells (hash by chunk) -> microtubules (per cell, cached) -> dye blocks (per 1 um of MT, cached)
        |  block schedule = every dye's full blink lifetime, sorted by time, pure function of address
        v
   std::vector<sim::BlinkEvent>   (FOV-relative um, zNm, tStart/tEnd in frames)   <-- existing type
        |
        v
   RenderPhotonImage / CollectGpuEmitters / ApplyNoiseChain   (UNCHANGED)
```

Design rules that make this work:

1. **Everything is a pure function of `(seed, address)`.** Cells, microtubules, dyes and their blink
   schedules are never stored in a global array and never depend on draw order, viewport, or
   what was visited before. Caches exist for speed only and can be dropped at any time.
2. **The renderer stays untouched.** The new source emits ordinary `BlinkEvent`s already translated
   into FOV-relative coordinates. `RenderPhotonImage`, the GPU path, camera noise, drift, illumination
   all keep working.
3. **Geometry is C++ in `Simulation/`, no MMDevice includes**, like the rest of the engine, so it can be
   tested standalone.

New files (suggested): `Simulation/SMLMCellField.h/.cpp` (hash, cells, packing, cytoplasm mesh),
`Simulation/SMLMMicrotubules.h/.cpp`, `Simulation/SMLMCellFieldDyes.h/.cpp` (labels, schedules, event
query), `SMLMDemoXYStage.h/.cpp`. Add them to `SMLMDemoCam.vcxproj` and `.filters`.

---

## 3. RNG and addressing

* `hashUnit(seed, cx, cy, k)` = `((Pcg4d(seed, cx, cy, k).a >> 9) + 0.5) * 2^-23`. Identical to
  `CounterRng::Uniform()` in `SMLMCounterRng.h` (same constants). Add a free function
  `HashUnit(uint32_t, uint32_t, uint32_t, uint32_t)`.
* **Integer semantics.** JS does `seed|0`, `cx|0`, `k|0` (int32 wrap) then `Math.imul` (uint32 wrap).
  In C++ take each argument as `int64_t`/`double`, reduce mod 2^32, cast to `uint32_t`. `hashStream`
  computes `k = base * 4096 + ctr`; with `base` up to ~9.5e6 that exceeds 2^32 and **wraps**, so the
  C++ must wrap identically (do the multiply in `uint64_t`, then truncate to `uint32_t`). Negative chunk
  coordinates are normal.
* **Channels.** Copy `CH` (0..44) and the `MT_*` stream bases (`MT_STREAM_BASE = 1000`,
  `MT_RESAMPLE_SPACING = 100000`, label bases) exactly, or the geometry stops matching the JS reference.
  `NUC_LONG (30)` is deliberately retired; leave it unused.
* **Never write two draws from the same sequential stream in one expression.** C++ leaves operand
  evaluation order unspecified (demoCam already hit this with `CombinedShotAndReadNoise`). JS is
  left-to-right, so a line-for-line port can silently diverge. Draw into named locals in the JS order.
* Use `double` everywhere the JS does. `Float32` is only for what demoCam already holds in float.

**Golden vectors.** Before anything else, add `tools/cellfield_parity_check/` (mirror
`tools/psf_parity_check/`): a Node script that `require`s the JS (or evals `index.html`'s script) and
dumps, for a fixed seed, `hashUnit` values including negative `cx,cy` and large `k`; `rawCandidate` for
~20 chunks; and one full `buildMicrotubulesForCell` result. Compare against the C++ output. Targets:
hash bit-exact, cell fields within 1e-9, MT points within 1e-6 µm (they pass through `sqrt`/`atan2`/
iterative clamps, so demand tolerance, not bit-exactness).

---

## 4. Cells, nucleus, cytoplasm

### 4.1 Units and frames
Everything is µm and radians. Cell fields (`semiMajor`, `nucLong`, `height`, ...) are in the cell's
cached **local** frame (origin at the cell centre, before `packRot`). `localToWorld(cell, lx, ly)`
rotates by `packRot` and translates by `cell.x, cell.y`. Microtubules and dyes are generated in the local
frame and mapped to world at the end; `z` is height above the coverslip (0 = floor).

### 4.2 Default parameters (from `index.html`)
Put these in a `CellFieldParams` struct with these defaults; expose only a curated subset as MM
properties (section 9) and keep the rest as constants until asked.

Field: `chunkSize 26`, `jitter 0.8`, `density (occupancy) 0.25`.
Cell: diameter 25-35, elongation (short/long) 0.2-0.8, `cellBlob 1.5`, height 3-6.
Nucleus: long axis 8-12, short/long 0.6-1, height 0.3-0.5 x long axis, offset 0.1, margin 1.5.
Cytoplasm: rim height 0.1-0.3, edge rise 0.1-0.5, mid height 1-2, mid distance 0.1-0.3 x cell radius, `nucMargin` 0.6. `cytoSmoothPasses` (default 12, slider 0-12) is the number of 3x3 binomial smoothing passes `smoothCytoGrid` applies to the mesh; `cytoRings` (default 60, slider 4-60) is the mesh's radial ring count; `cytoTheta` (default 128, slider 32-128) is its angular sample count.
Packing: enabled, min gap 1.0, relax iterations 12, step (damping) 0.2, rotation allowed.
Microtubules: density 0.45 /µm², start offset 0-0.3, start XY jitter 0, end offset 0.01-0.4,
end direction jitter 145 deg, wobble turn 0.8, wobble path x1.05, step length 0.05, path smoothing 1.5,
min turn radius 0.15, min separation 0.05, max z slope 5, line width (draw only, skip).

### 4.3 The viewport-dependence trap (do not port this part as-is)
`draw()` builds the candidate map over "visible chunks + margin", runs `relax`/`prune` over *that
window*, and the README admits a cell near the window edge can resolve to a slightly different position
depending on how far the margin extends. That is fine for a viewer and **wrong for a moving stage**:
panning the FOV away and back must reproduce identical cells.

Do this instead: pack on **fixed, absolutely aligned blocks** (e.g. 8x8 chunks). For a block, build
candidates over the block plus a margin of `interactionChunks(p) + 2`, run `relax`/`prune` (same
iteration counts), and *keep only the cells whose home chunk lies in the block core*. The result then
depends only on the block address. Cells near block borders can differ slightly from what the neighbour
block would have decided for the same pair; accept that (it is deterministic), and log it. Cache per
block. Provide a switch `packing=off` that uses raw jittered positions (exactly reproducible, some
overlaps) for a first milestone.

Keep the safety caps (`CHUNK_CAP`, `PRUNE_ROUNDS`); they exist because a pathological setting hangs.

### 4.4 Microtubules
Port `buildMicrotubulesForCell` faithfully, including the caches keyed by a signature of the cell shape
and parameters (`mtResultSig`, `mtCellShapeSig`). Key constants: `MT_MAX_PER_CELL 5000`,
`MT_COLLISION_MAX_TOTAL_POINTS 8000` (collision resolution is *skipped* above this; documented, not a
bug), `MT_NUDGE_ROUNDS 5`, `MT_RESAMPLE_ROUNDS 2`. Density is per µm² of the cell footprint.

Additionally store, per microtubule, what the dye stage needs and JS recomputes each time: cumulative
arc length and the **parallel-transport frame** `(T, U, V)` per segment (see
`buildMicrotubuleLabelPoints`). Parallel transport is sequential along the path, so this must be a
per-microtubule table built once, not evaluated per dye.

---

## 5. Fluorophore labels (dyes)

### 5.1 Geometry per lattice site (from `buildMicrotubuleLabelPoints`)
13 protofilaments (`MT_N_PROTOFILAMENTS`), 13_3 lattice: protofilament `k` sits at angle
`phase + k*2*pi/13` and axial offset `(k * 3 * 8 / 13) mod 8` nm, then every 8 nm (`MT_DIMER_NM`).
Surface attachment at radius `MT_RADIUS_NM = 12.5` nm from the axis; binder tip a further
`MT_BINDER_NM = 12` nm radially outward; dye displaced from the tip by a uniform-in-volume distance in
`[MT_LINKER_MIN_NM, MT_LINKER_MAX_NM] = [2, 5]` nm in a uniform direction (`mtDisplaceByLinker`,
originally `displaceByLinker` in webSMLM). Per-microtubule seam phase = `2*pi*hashUnit(seed,cx,cy,8000000+i)`.
That is about **1625 sites per µm of microtubule** at 100% labelling.

### 5.2 Do not materialise the sites
A 30 µm cell at the default density holds millions of sites. A FOV of 13x13 µm covers on the order of
800 µm of microtubule = ~1.3 M sites at 100% labelling. Rules:

* **Hash first, geometry later.** A dye's identity is
  `H1 = Pcg4d(seed ^ DYE_SALT, cx, cy, mtIndex)`, `H2 = Pcg4d(H1.a, k, n, purpose)` with `k` the
  protofilament and `n` the dimer index along the microtubule. Every per-dye draw (labelled?, first
  activation time, blink durations, brightness, linker displacement) is `HashUnit`-style off `H2` with a
  different `purpose`. Decide "labelled" and "ever activates before the horizon" from hashes only, and
  compute the 3D position (frame lookup, linker) **only for dyes that will emit**.
* **Blocks.** The unit of generation and caching is a *block*: one microtubule, 1 µm of arc length
  (~1625 sites). Blocks are addressed `(cx, cy, mtIndex, blockIndex)`.
* **Labelling efficiency** is one hash draw per site. Default the new property to something sparse
  (start at 5-10%; the webSMLM structures use 70% but their site counts are tiny).
* The JS preview uses a sequential `hashStream` for labels, so **JS and C++ dye-for-dye positions will
  not match**. The C++ hashing above is normative; only the *statistics* (ring at 12.5 nm, tip at 24.5 nm,
  linker in range, lattice angles/stagger) must agree. If exact match is ever wanted, change the JS to
  the same `H1/H2` scheme and update this file.
* Dyes whose `|z - focus|` lies beyond the PSF kernel's z range must be **culled**, not clamped (the
  existing renderer clamps to the end plane, which would draw a cell's out-of-range dyes as bright
  in-focus-looking blobs). Cull in the query, count them, and log once.

---

## 6. Blinking from dye positions

### 6.1 Per-dye schedule = pure function of the dye's hash
For a labelled dye: first activation `tAct = -ln(U) * activationMeanSec`; then repeat:
ON for `Exp(onLifetimeSec)`, then bleach with probability `blinkBleachProb`, else dark for
`Exp(offLifetimeSec)` and blink again. Per-blink brightness log-normal with CV `photonCV`, mean 1.
Cap blinks per dye (e.g. 1000) so a tiny `blinkBleachProb` cannot loop forever. Reuse the existing
properties `FluoParam_OnLifetimeSec`, `FluoParam_OffLifetimeSec`, `FluoParam_BlinkBleachProb`,
`FluoParam_PhotonCV`, `FluoParam_PhotonsPerSecond`. This is the same three-state model as
`EmitterModel`'s rich path, but with **persistent identity**: the same dye is the same dye every time
the stage returns.

The whole lifetime is finite (bleaching), so a block's schedule can be generated **once**, for all
time, as a list of `{tOnSec, tOffSec, brightness, dyePosLocal}` sorted by `tOn`, plus `maxOnSec` for the
block. No time windows, no replay.

### 6.2 Clock (two modes, one code path)
The schedule is in "dye clock" seconds `tau`. A block maps real simulated time to `tau`:

* **`Global`** (implement first, purely functional): `tau = t`. Blocks you are not looking at bleach
  anyway. Simple, reproducible, parallelisable, but a sample you leave for ten minutes comes back
  depleted.
* **`Illuminated`** (milestone 5, realistic default): `tau = E_block`, the accumulated seconds the block
  has been inside the illuminated FOV. Keep `E_block` in a sparse map `blockId -> double` that lives
  for the acquisition; advance it by `frameDurationSec` for every block that intersects the FOV
  rectangle (+ PSF margin) in a rendered frame; reset on acquisition start (`StartLiveProducer`,
  same place the drift origin resets). Block granularity (1 µm) is a deliberate approximation.
  Reproducible for a given stage trajectory.

Simulated time is `frameIndex * frameDurationSec` (as drift already does), **not** wall-clock time.

### 6.3 Activation rate (keep `density` meaning ON-density)
webSMLM's rule, which demoCam's `EmitterModel` also follows: the user-facing density means *mean
emitters ON per µm² per frame*. Per-dye rate `k = rho * A_fov / (N_dyes_in_fov * meanBlinks * tauOn)`
with `meanBlinks = 1/blinkBleachProb`. `N_dyes_in_fov` is estimated once (labelled-site count over the
FOV at the initial stage position, via the cheap hash-only pass). Expose `General_EmitterDensityPerSec`
as the target; log the resulting `activationMeanSec = 1/k`. Rough magnitude at 5% labelling, 0.2 /µm²,
5 blinks, 90 ms ON: `activationMeanSec` of order 100-300 s. Provide an override property to set the
mean directly.

### 6.4 Query, per frame
```
std::vector<BlinkEvent> CellFieldSource::EventsForFrame(long f, StagePose pose, double fovWUm, double fovHUm,
                                                        double frameDurSec, const ...& params);
```
1. World rectangle = FOV centred on `pose` (x,y) expanded by a PSF margin (`>= 3 * kernel half-width`
   is plenty; ~2 µm).
2. Enumerate cells/blocks intersecting it (cell footprint bounding circle `rOuter` first, then per
   microtubule bounding box, then per 1 µm block centre). Generate/fetch cached blocks.
3. Per block, get `tauNow` (section 6.2). Binary-search the sorted schedule for `tOn` in
   `(tauNow - maxOnSec, tauNow + frameDurSec)`; keep events with `tOff > tauNow`.
4. For each hit emit a `BlinkEvent`:
   * `xUm = xWorld - (pose.x - fovW/2)`, `yUm = yWorld - (pose.y - fovH/2)`  (FOV top-left origin, as
     every existing pattern uses; the renderer then adds drift),
   * `zNm = (zWorld - focusHeightUm) * 1000`  (the Z stage still adds `globalZOffsetUm` in the renderer),
   * `tStart = f + (tOn - tauNow)/frameDurSec`, `tEnd = f + (tOff - tauNow)/frameDurSec`,
   * `brightness` from the schedule.
   Because the FOV moves between frames, **emit a fresh event per frame with that frame's translation**;
   never keep one event across frames. The renderer's frame-overlap weighting then still works because
   `tStart/tEnd` keep their original meaning.
5. Cull events outside the FOV + margin and beyond the kernel z range (section 5.2).

Cost budget: a rendered frame should touch a few thousand blocks and emit hundreds of events. Profile
the first-visit block generation separately from the steady state; the steady state must not
regenerate anything. Cache blocks in an LRU (bounded, like `MT_RESULT_CACHE_MAX`), and generate the
blocks of the next FOV-edge ring ahead of need if first-visit latency shows up in live mode.

---

## 7. Dummy XY stage

### 7.1 Shared state (`Simulation/SharedStageState.h`)
Keep `std::atomic<double> zPositionUm` untouched. Add an XY motion model behind a `std::mutex`:
committed position, target, `moveStart` (`std::chrono::steady_clock`), `speedUmPerSec`, `settleSec`.
Provide, all MMDevice-free:

* `void SetXyTarget(double x, double y)` (starts a move from the *current interpolated* position),
* `void PositionXyAt(steady_clock::time_point now, double& x, double& y)` (constant-velocity linear
  interpolation, clamped at the target),
* `bool XyBusy(now)` (true until arrival + settle),
* `void SetXyOrigin`/limits as plain fields.

The camera samples `PositionXyAt(now)` **once per produced frame** (live) and uses the pose at that
frame. Motion blur during the exposure is ignored on purpose; note it in the docs.

### 7.2 MM device `SMLMDemoXYStage`
Derive from `CXYStageBase<SMLMDemoXYStage>` (`DeviceBase.h`). It is in `third_party/mmCoreAndDevices`.
**Implement every pure virtual it declares**: an abstract class fails at the `new SMLMDemoXYStage()` in
`SMLMDemoCameraModule.cpp`, not in the stage's own files (same trap already recorded for
`IsStageSequenceable`). At minimum: `Initialize/Shutdown/GetName`, `Busy` (return
`XyBusy(now)` so MM waits for the move), `SetPositionSteps/GetPositionSteps`, `SetPositionUm/GetPositionUm`
(the base provides the step-based defaults; overriding the Um versions directly and using a step size
of 0.1 µm is simplest), `SetRelativePositionUm`, `Home` (0,0), `Stop`, `SetOrigin`, `GetLimitsUm`,
`GetStepLimits`, `GetStepSizeXUm/YUm`, `IsXYStageSequenceable` (false). Register it in
`InitializeModuleData`/`CreateDevice` as `MM::XYStageDevice`, name `"SMLMDemoXYStage"`, and add to the
vcxproj/filters.

Properties (group prefix rule from `CLAUDE.md`: no standard keyword applies, so use `General_`):
`General_StageSpeedUmPerSec` (default 5000, i.e. a few mm/s like a real stage), `General_StageSettleMs`
(default 20), `General_StageLimitUm` (default 1e6, the field is effectively unbounded). Add both stage
and camera through the Hardware Configuration Wizard; they communicate only through the singleton, as
the Z stage already does.

### 7.3 Coordinates
The stage position is the **world coordinate of the FOV centre** in µm. World `(0,0)` is the origin of the
cell-field chunk grid. Increasing X moves the FOV to the right over the sample, so features move *left* in
the image; same for Y downward. Real setups vary (camera mirroring/rotation), so add
`General_StageInvertX/Y` (default off) rather than hard-coding a convention. **Drift** composes as an extra
offset on the FOV centre: `effectiveCentre = stage + drift`, using the existing `ComputeDriftOffsetPx`
result converted to µm (do not shift dye positions separately or you will double-apply it).

---

## 8. Camera integration

* Add `PATTERN_CELL_FIELD = 15` to `SMLMPatternType` and a `SimType_Pattern` value `CellField`. It is
  **not** an `IPatternGenerator` (it does not sample sites per blink); dispatch on
  `CurrentPatternType()` in the two places that produce events:
  `StackGenerationWorker` (instead of `model.GenerateAllEvents`) and `LiveProducerLoop` (instead of
  `liveEmitterModel_.AdvanceOneFrame`). The live loop already reads the Z stage per frame; read the XY
  pose the same way and call `EventsForFrame`. The out-of-focus population (`Background_OutOfFocusRatio`)
  is site-list based; for `CellField` disable it with a one-line log, or later feed it from culled dyes.
* **Live/MDA is the intended mode.** A multi-position MDA moves the XY stage between snaps, which is the
  use case. **Precomputed-stack mode** generates all frames up front, so it can only use one stage pose:
  snapshot the pose at generation start, apply it to every frame, and document that moving the stage
  afterwards does not change the stack (the Z stage behaves the same way today). In `Illuminated` mode
  use `E = f * frameDurationSec` for blocks in the FOV.
* `BuildShapingFields` calls `EmitterModel::SampleSitesForHaze` for the cell-contrast/haze background.
  For `CellField` either return an empty site list (background falls back to flat) or sample the dyes in
  the initial FOV. Check that function before deciding; do not let it dereference a null pattern.
* Any property that changes the field must call `InvalidateStack()` (existing convention) so live mode
  rebuilds; the block/cell caches must key on a signature of the parameters, like `mtResultSig`.
* Set the focus plane: `CellField_FocusHeightUm`-style property (default ~1.5 µm above the coverslip,
  inside the nucleus/cytoplasm range). Cells are 3-6 µm tall, so `PSFParam_PsfZRangeUm` needs to be
  >= ~7 µm to cover them; warn (like the existing structure z-extent warning) if it is not.

---

## 9. MM properties

Follow the group-prefix rule in `CLAUDE.md` and update that file's bullet list. Suggested curated set,
all with the defaults in 4.2:

* `SimType_Pattern` gains `CellField`.
* `SimType_CellFieldChunkSizeUm`, `SimType_CellFieldOccupancy`, `SimType_CellFieldPacking` (On/Off),
  `SimType_CellFieldCellDiameterMinUm/MaxUm`, `SimType_CellFieldMicrotubuleDensityPerUm2`,
  `SimType_CellFieldFocusHeightUm`.
* `General_LabelingEfficiencyPct` (existing; reuse for the dyes, default for this pattern 5-10) and
  `General_EmitterDensityPerSec` (target ON-density, section 6.3).
* `FluoParam_...` (existing kinetics) plus `SimType_CellFieldActivationMode` (`Global`/`Illuminated`) and an
  optional `SimType_CellFieldActivationMeanSec` override.
* Stage: `General_StageSpeedUmPerSec`, `General_StageSettleMs`, `General_StageInvertX/Y`.

Seed: reuse `SimType_RandomSeed`; derive the cell-field seed as its own XOR-constant stream (like
`structureSeed`), never the arrival/noise stream. Do not add UI beyond properties.

---

## 10. Verification

1. **Parity** (section 3): hashes, cells, one microtubule set against the JS dump.
2. **Determinism**: same seed, stage moved away 1 mm and back, frames identical (Global mode). Also the
   same after dropping every cache.
3. **Dye statistics** on a straight synthetic microtubule: attachment radius 12.5 nm, tip 24.5 nm, dye
   within `[2,5]` nm of the tip, radial CDF of the linker displacement ~ r^3, per-protofilament angles
   are multiples of 2*pi/13 + phase, axial stagger 3*8/13 nm, count ~ `1625 * length_um * efficiency`.
4. **Kinetics**: measured ON-density in a FOV matches the target within ~10% (this needs the rate
   formula of 6.3), blink-count distribution geometric with mean `1/blinkBleachProb`.
5. **Stage**: extend `tools/test_smlmcam.py` (pymmcore-plus): add the XY stage, set positions, assert
   `Busy` transitions, that a `speed` move takes ~`dist/speed`, that a known feature shifts by the
   expected pixels between two snaps, and that position-in-image flips with `StageInvertX`.
6. **Performance**: report first-visit block generation time, steady-state frame time at the default FOV
   (128x128 @ 100 nm), and memory after panning across a 1 mm path.
7. **Visual**: someone must look at it in Micro-Manager Studio; headless checks do not replace that
   (the repo's other features are all marked "not visually verified" for this reason).

---

## 11. Suggested order of work

Mark each done here.

1. [ ] Hash + golden-vector tooling; port `rawCandidate`/outline/nucleus; parity test on cells.
2. [ ] Cytoplasm mesh + microtubules (unpacked cells first); parity test on one microtubule set.
3. [ ] Dyes: block generation, schedule, `EventsForFrame`; static FOV, `Global` clock; renders through the
       existing pipeline. First visual check.
4. [ ] Shared XY state + `SMLMDemoXYStage` + camera wiring; stage test in `test_smlmcam.py`.
5. [ ] Fixed-block packing (4.3); `Illuminated` clock; activation-rate auto-calibration (6.3).
6. [ ] Properties, `CLAUDE.md` update, performance pass, docs.

## 12. Known gaps to keep in mind (not for the first pass)

* Motion blur during an exposure while the stage moves; per-frame stage jitter.
* Dyes beyond the kernel z range are culled, not added as diffuse background haze.
* `Illuminated` clock is per 1 µm block, not per dye.
* Packing near block borders can differ slightly from what the neighbouring block would decide.
* Only microtubules carry labels; nucleus/cytoplasm labels (lamin, mitochondria, NUP) do not exist yet
  in the JS either.
* The JS prototype and this port are not validated quantitatively against real SMLM data.
