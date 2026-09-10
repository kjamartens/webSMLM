# Simulation parity: webSMLM vs. demoCam_SMLM_MM

This is a point-in-time snapshot of the parameter/feature correspondence
between this project's simulator (`webSMLM.html`) and the Micro-Manager
device adapter at `C:\GitHub\demoCam_SMLM_MM`
(`DeviceAdapter/SMLMDemoCam/Simulation/`). **It is not auto-updated.**

- webSMLM snapshot: commit `74654775f2182` (2026-09-08)
- demoCam_SMLM_MM snapshot: commit `35bcbda7c19c2` (2026-08-31) **plus**
  uncommitted working-tree changes made in the session that wrote this
  file (per-emitter 3D + 3D/NPC structures + labeling efficiency,
  `PsfInterp`, `PsfEvalMethod`/chirp-Z) -- see that repo's `CLAUDE.md`
  "webSMLM parity feature" section for the authoritative up-to-date
  status once those changes are committed.

**Staleness warning:** this project has other collaborators and no
mechanism here notifies demoCam_SMLM_MM of new commits. Before starting
*any* simulation-engine work in either repo that assumes parity (or a gap)
described below, diff this table against both projects' current state and
refresh it.

Direction so far: **demoCam_SMLM_MM catches up to webSMLM.** This project
has not been modified for parity (this file is the only addition here).

## PSF models

| Concept | webSMLM (`simulation_psfModel`) | demoCam_SMLM_MM (`PsfModel`) |
|---|---|---|
| Fixed-sigma Gaussian | `gaussian`, sigma hardcoded 1.3px, not NA/lambda-derived | `Gaussian`, sigma = 0.21*lambda/NA (physics-derived) -- **demoCam ahead** |
| Richards-Wolf (scalar Debye-Kirchhoff) | absent | `RichardsWolf` (PSFGenerator) -- **demoCam ahead** |
| Gibson-Lanni | absent as a standalone model (subsumed into the Zernike model at zero coefficients) | `GibsonLanni` (PSFGenerator) -- **demoCam ahead** |
| Gibson-Lanni + Zernike (scalar, full 2D pupil) | `zernike` (default), N_RHO=20 x N_PHI=40 polar quadrature or chirp-Z (`simulation_psfEvalMethod`: `direct`\|`fft`, default `fft`) | `GibsonLanniZernike`, same N_RHO/N_PHI, `PsfEvalMethod`: `Direct`\|`ChirpZ`, default `Direct` -- ported this session, ~3.7x measured speedup, 0.22-0.29% relative L2 agreement with Direct |
| Zernike coefficients | 15 (OSA 0-14), milliwaves in UI, presets + custom | 15 (OSA 0-14), waves, `PsfZernikeCoefficients` + `PsfZernikePreset` (10 presets) -- same convention/index mapping |
| Sub-pixel kernel placement | `simulation_psfInterp`: `nearest`\|`linear`\|`cubic`\|`fft` (Fourier-shift), default `cubic` | `PsfInterp`: `Nearest`\|`Linear`\|`Cubic`, default `Nearest` (unchanged legacy behavior) -- ported this session (Linear/Cubic); FFT-shift mode not ported (deferred, see below) |
| Z-stack / defocus | per-emitter Z (site's own z), nearest-plane lookup, no blend | per-emitter Z (site depth + Z-stage global offset add) -- ported this session, nearest-plane only (two-plane blend deliberately not ported, matching webSMLM's own removed-and-not-reintroduced decision) |
| Measured/experimental PSF, cubic-spline PSF, DH/tetrapod/biplane | absent in both | absent in both |

## Emitter placement / structures

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Sampling model | finite pre-generated site list (2000-50000 sites, area-scaled), arrival picks a random index | 9 patterns sample continuously per blink (no site list) **plus** (this session) `SiteListPattern` for the 4 new 3D structures below -- see that project's `SMLMPatterns.h` header comment for why the distinction is deliberate |
| 2D patterns | filaments+ring, NPC-only-as-3D (no dedicated 2D resolution-target pattern) | Circle, Lines, Grid, Random, CustomPoints, Spiral, Star, Heart, ResolutionTarget -- **demoCam ahead** on 2D pattern variety |
| 3D structures | `filaments_ring`, `nup`, `tiltedPlane`, `uniform3D`, `shell` | `TiltedPlane`, `Uniform3D`, `Shell`, `NUP` -- ported this session (formulas/constants ported from webSMLM's `buildStructure()`/`buildNupStructure()`); filaments+ring's 3D variant not ported (demoCam's 9 existing 2D patterns instead get z via the `ZSpreadPattern` decorator, a uniform spread rather than webSMLM's sine-correlated z) |
| NUP/NPC geometry | Thevathasan 2019 Nup96, Wanninger 2023 CIR4MICS params, 64 sites/pore, rejection-sampled spacing, bowl curvature, topdown/sideways | same -- ported verbatim this session, with one added corelog warning (rejection-sampling shortfall) webSMLM doesn't have |
| Labeling efficiency | `simulation_labelEfficiency` (%, default 70), applied once to the site list before arrivals | `LabelingEfficiencyPct` (%, default 100), applied once to the 4 site-list structures -- ported this session; **not** applied to the 9 continuous patterns or CustomPoints (scope limit, documented) |
| Emitter Z spread on 2D patterns | filaments' z correlates with the same sine driving y (documented as an accuracy limitation in webSMLM's own docs) | `ZSpreadPattern`: uniform z spread, uncorrelated with the 2D geometry -- a deliberate simplification, not a port of webSMLM's correlated approach |

## Blinking kinetics / photophysics

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Arrival process | areal-density Poisson, sites reused across the movie | areal-density Poisson (per-frame or per-tick), continuous patterns draw a fresh position -- same core math, different site-reuse semantics per the sampling-model row above |
| ON duration | single exponential (mean = `simlifetime`), exactly one blink per arrival | single exponential (mean = `OnLifetimeSec`/`OnLifetimeFrames`), exactly one blink per arrival -- same model |
| Photon count | constant per emitter (no distribution), scaled by frame-overlap fraction | constant per emitter (no distribution), scaled by frame-overlap fraction -- same |
| Bleaching, dark/triplet states, multi-blink kinetics | absent (single finite ON period *is* the "bleach") | absent, same |
| Sub-pixel position | continuous float positions both paths | continuous float positions both paths |

## Camera / noise model

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Chain | bg -> Poisson(signal) -> +Gaussian read noise (post-gain conversion... actually pre-gain, see webSMLM inventory) -> /gain -> +per-pixel offset -> clamp at 0 | photons -> *QE +dark current -> Poisson+Gaussian read noise (combined) -> /gain -> +offset -> clamp [0,65535] -- **demoCam ahead**: QE, dark current, 16-bit clamp, per-pixel sCMOS-style gain AND read-noise maps (not just offset) |
| QE / dark current | absent | `QuantumEfficiency` (0.85), `DarkCurrentElectronsPerSec` (1.03) -- **demoCam ahead** |
| Per-pixel maps | offset only (`simulation_offset_std`) | offset, gain (`PixelGainStdPct`), read noise (`PixelReadNoiseStdPct`) -- **demoCam ahead** |
| EMCCD / excess noise factor | absent | absent |
| Quantization / bit depth | none (Float32 output, no rounding) | 16-bit uint clamp -- **demoCam ahead** (more camera-realistic) |
| Gain units | photons/ADU | photons/ADU -- same convention |

## Background / drift / focus

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Background | uniform scalar, Poisson-fluctuating | uniform scalar (`BackgroundPhotonsPerSec`), Poisson-fluctuating -- same |
| Structured background, vignetting, autofluorescence | absent | absent |
| Drift | linear, single random direction (X/Y), stored ground truth | linear, fixed diagonal (`DriftNmPerSecX`, Y = half rate) -- not random-direction; no ground-truth export (see below) |
| Focus / Z stage | no separate stage concept -- structure Z is absolute | real `SMLMDemoZStage` MM::Stage device, global focus offset, live-drivable -- **demoCam ahead** (no equivalent in webSMLM) |

## Outputs / ground truth

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Raw movie export | **absent** (in-memory only) | uint16 stack via standard MM camera API -- **demoCam ahead** |
| Ground-truth emitter/localization export | in-memory only (`groundTruthEvents`/`groundTruthLocs`), consumed by "View GT" and "Score vs truth"; no file export | **absent entirely** -- `BlinkEvent` positions are discarded after rendering; this is demoCam's single biggest ground-truth gap |
| GT scoring vs. localization output | yes (`scoreTruthCore`, recall/precision/Jaccard/lateral+axial error stats) | absent |

## Not tracked here

webSMLM's ThunderSTORM-compatible CSV export, Settings-JSON round-trip, and
reconstruction/rendering pipeline have no demoCam equivalent (demoCam is a
camera device, not an analysis tool) and are out of scope for this
comparison by design, not a gap.

## Numeric cross-check

Two separate checks, both run this session:

1. **Within demoCam_SMLM_MM**: `GibsonLanniZernikePSF`'s ported chirp-Z
   evaluator vs. its original direct-quadrature evaluator (both Java, same
   process) -- 0.22-0.29% relative L2 (sum-normalized) at NA 1.4/660nm,
   65x65px, zero-Zernike and 0.15-wave astigmatism cases. Recorded in that
   repo's `CLAUDE.md`.
2. **Cross-project** (`demoCam_SMLM_MM/tools/psf_parity_check/`): webSMLM's
   own JS `computePsfPupilForZPlane`/`computePsfIntensityPlane` (the direct
   polar-quadrature path) vs. demoCam's Java port of the same functions
   (`GibsonLanniZernikePSF.computeSliceDirect`), same physical parameters
   (NA 1.4, lambda 660nm, ns=ni=1.518, ti0 150um, in-focus, 0.15-wave
   vertical astigmatism, 65x65px @ 25nm/px oversampled resolution).
   **Result: 0.0000% relative L2 difference** (bit-for-bit agreement,
   modulo ordinary floating-point noise) -- expected and confirms the
   direct-quadrature port from webSMLM into Java is faithful, since it's
   the same algorithm/formulas/constants in both languages. Re-run this
   check (`tools/psf_parity_check/README.md`) after any future change to
   either engine's direct-quadrature PSF math.
