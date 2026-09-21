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

**Partial refresh, 2026-09-19 (webSMLM side only).** webSMLM builds
`2026-09-19a`–`d` added multi-blink photophysics, a structured background,
blinking out-of-focus emitters and scoring criteria. The rows those builds
invalidate are updated below and marked *(webSMLM 2026-09-19)*; **the
demoCam_SMLM_MM column was NOT re-checked** — that repo was not available to
the session making this edit — so wherever such a row says "same" or
"absent" for demoCam it describes the 2026-09-08 snapshot, and every one of
them now reads **webSMLM ahead** until demoCam is diffed again.

**Partial refresh, 2026-09-20 (webSMLM side only).** webSMLM builds
`2026-09-20a`-`h` added an EMCCD sensor path (QE, Gamma gain register,
clock-induced charge, integer ADU at a settable bit depth) with its
analysis-side excess-noise correction, a non-uniform illumination profile,
engineered PSFs (higher-order astigmatism presets and a Gauss-Laguerre
double-helix phase mask) with a per-PSF Cramer-Rao bound, and a PSF-model
fitter (`psfmle`). Rows those builds invalidate or add are marked
*(webSMLM 2026-09-20)*; **the demoCam_SMLM_MM column was again NOT
re-checked** for the same reason as above.

**Partial refresh, 2026-09-21 (webSMLM side only).** webSMLM builds
`2026-09-21b`-`d` changed how the simulator computes, not what it models:
a block-summed splat (same result to 3e-8, ~10x faster), counter-based
camera noise (pcg4d, so seeded movies changed once), and WebGPU paths for
splat + noise and the `direct` PSF build. They also found that the `direct`
evaluator is **wrong on wide kernels** — this matters to demoCam, whose
default is `Direct` (see the PSF-models row and "Numeric cross-check" §3).
Rows marked *(webSMLM 2026-09-21)*; the demoCam column was again NOT
re-checked.

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
| Gibson-Lanni + Zernike (scalar, full 2D pupil) | `zernike` (default), N_RHO=20 x N_PHI=40 polar quadrature or chirp-Z (`simulation_psfEvalMethod`: `direct`\|`fft`, default `fft`); `direct` runs on the GPU with Use GPU acceleration, and logs a warning on kernels wider than its valid radius *(webSMLM 2026-09-21)* | `GibsonLanniZernike`, same N_RHO/N_PHI, `PsfEvalMethod`: `Direct`\|`ChirpZ`, default `Direct` -- ported this session, ~3.7x measured speedup, 0.22-0.29% relative L2 agreement with Direct **on a 1.6 µm kernel only**. **demoCam's default `Direct` inherits the wide-kernel aliasing error** (see "Numeric cross-check" §3) -- switch its default to `ChirpZ`, or keep its kernel under ~2.2 µm *(webSMLM 2026-09-21)* |
| Zernike coefficients | 28 (OSA 0-27, n<=6), milliwaves in UI, presets + custom; a 15-value custom string is zero-padded, so older settings files are unchanged *(webSMLM 2026-09-20)* | 15 (OSA 0-14), waves, `PsfZernikeCoefficients` + `PsfZernikePreset` (10 presets) -- same convention/index mapping |
| Sub-pixel kernel placement | `simulation_psfInterp`: `nearest`\|`linear`\|`cubic`\|`fft` (Fourier-shift), default `cubic` | `PsfInterp`: `Nearest`\|`Linear`\|`Cubic`, default `Nearest` (unchanged legacy behavior) -- ported this session (Linear/Cubic); FFT-shift mode not ported (deferred, see below) |
| Z-stack / defocus | per-emitter Z (site's own z), nearest-plane lookup, no blend | per-emitter Z (site depth + Z-stage global offset add) -- ported this session, nearest-plane only (two-plane blend deliberately not ported, matching webSMLM's own removed-and-not-reintroduced decision) |
| Engineered PSFs | higher-order astigmatism presets (`saddlePoint`/`extendedRange`/`extendedRangeStrong`, astigmatism stacked at OSA j=5/13/25) and a real double-helix phase mask (`simulation_psfMaskType`, Gauss-Laguerre superposition along l=2p+1, measured ~60 deg rotation over +/-800 nm) *(webSMLM 2026-09-20)* | absent |
| PSF figure of merit | `psfZCramerRao()`: x/y/z Cramer-Rao bound from the kernel stack, PSF-shape agnostic (the astigmatism-specific `zUsableNm` is kept alongside it) *(webSMLM 2026-09-20)* | absent |
| PSF-model ("vector") fitting | `psfmle`: fits the camera-pixel-integrated modelled PSF itself, theta=[x,y,N,bg,z], tricubic over samples, coarse z scan then Newton, joint z CRLB, no width calibration needed; single-threaded for now *(webSMLM 2026-09-20)* | absent (fitting is out of scope for a device adapter) |
| Measured/experimental PSF, bead-calibrated spline PSF, biplane | absent in both | absent in both |

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
| ON duration | single exponential (mean = `simlifetime`); one blink per molecule by default, a geometric number (mean 1/`simulation_blinkBleachProb`) otherwise *(webSMLM 2026-09-19)* | single exponential (mean = `OnLifetimeSec`/`OnLifetimeFrames`), exactly one blink per arrival -- same model |
| Photon count | constant by default; per-blink log-normal rate with CV `simulation_photCV`, mean preserved; scaled by frame-overlap fraction *(webSMLM 2026-09-19)* | constant per emitter (no distribution), scaled by frame-overlap fraction -- same |
| Bleaching, dark/triplet states, multi-blink kinetics | three-state ON ⇄ dark → bleached per molecule (`simulation_blinkBleachProb`, `simulation_offLifetime`); `dens` keeps meaning ON-density; no separate triplet state; default = original single-blink model, byte-identical for a seed *(webSMLM 2026-09-19)* | absent, same |
| Sub-pixel position | continuous float positions both paths | continuous float positions both paths |

## Camera / noise model

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Chain | bg -> Poisson(signal) -> +Gaussian read noise (post-gain conversion... actually pre-gain, see webSMLM inventory) -> /gain -> +per-pixel offset -> clamp at 0 | photons -> *QE +dark current -> Poisson+Gaussian read noise (combined) -> /gain -> +offset -> clamp [0,65535] -- **demoCam ahead**: QE, dark current, 16-bit clamp, per-pixel sCMOS-style gain AND read-noise maps (not just offset) |
| QE / dark current | QE on the EMCCD path (`simulation_qe`, 0.9); no dark current *(webSMLM 2026-09-20)* | `QuantumEfficiency` (0.85), `DarkCurrentElectronsPerSec` (1.03) -- **demoCam ahead on dark current** |
| Per-pixel maps | offset only (`simulation_offset_std`) | offset, gain (`PixelGainStdPct`), read noise (`PixelReadNoiseStdPct`) -- **demoCam ahead** |
| EMCCD / excess noise factor | `simulation_cameraType='emccd'`: Poisson -> Gamma gain register, so variance = 2x mean (measured 2.00); clock-induced charge (`simulation_cic`); analysis side corrects for it via `cameraExcessNoise` *(webSMLM 2026-09-20)* | absent -- **webSMLM ahead** |
| Quantization / bit depth | integer ADU clipped at 2^`simulation_bitDepth`-1 on the EMCCD path; sCMOS path still Float32 (kept for byte-identity with earlier builds) *(webSMLM 2026-09-20)* | 16-bit uint clamp |
| Gain units | photons/ADU | photons/ADU -- same convention |
| Illumination profile | `simulation_illumProfile`: `flat` (default) \| `gaussian` \| `sigmoid`, width `simulation_illumFwhmPct`; attenuation normalised to peak 1, applied to emitters AND background *(webSMLM 2026-09-20)* | not checked (this column is the 2026-09-08 snapshot) |

## Background / drift / focus

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Background | `simbg` = FOV-mean photons/px, Poisson-fluctuating; flat by default *(webSMLM 2026-09-19)* | uniform scalar (`BackgroundPhotonsPerSec`), Poisson-fluctuating -- same |
| Structured background, vignetting, autofluorescence | cell-shaped autofluorescence field (`simulation_bgCellContrast`), static out-of-focus haze from the structure's projected density (`simulation_bgHazeWeight`/`Width`), fade to a 30% floor (`simulation_bgDecayFrames`), and blinking out-of-focus emitters through the real defocused PSF (`simulation_hazeRatio`/`Depth`); no vignetting *(webSMLM 2026-09-19)* | absent |
| Drift | linear, single random direction (X/Y), stored ground truth | linear, fixed diagonal (`DriftNmPerSecX`, Y = half rate) -- not random-direction; no ground-truth export (see below) |
| Focus / Z stage | no separate stage concept -- structure Z is absolute | real `SMLMDemoZStage` MM::Stage device, global focus offset, live-drivable -- **demoCam ahead** (no equivalent in webSMLM) |

## Outputs / ground truth

| Concept | webSMLM | demoCam_SMLM_MM |
|---|---|---|
| Raw movie export | **absent** (in-memory only) | uint16 stack via standard MM camera API -- **demoCam ahead** |
| Ground-truth emitter/localization export | in-memory only (`groundTruthEvents`/`groundTruthLocs`), consumed by "View GT" and "Score vs truth"; no file export | **absent entirely** -- `BlinkEvent` positions are discarded after rendering; this is demoCam's single biggest ground-truth gap |
| GT scoring vs. localization output | yes (`scoreTruthCore`): recall/precision/Jaccard, median+percentile lateral/axial error, photon-threshold + edge "don't care" classes, isolated/crowded split, recall-vs-photons curve with 50% point, challenge efficiency *(webSMLM 2026-09-19)*; plus per-molecule metrics and an effective z range, and a `validation_preset` that reproduces the published SMLM-Challenge-2016 rules (3D cylinder matching, quantile photon threshold, border excluded before matching) *(webSMLM 2026-09-20)* | absent |

## Presets (webSMLM only, 2026-09-19)

`simulation_realism` (`min`/`med`/`max`) and `simulation_densityPreset`
(`low`/`med`/`high`) are UI shorthands that only WRITE the ordinary parameters
above; they carry no physics of their own, so a port needs the underlying
parameters, not the presets.

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
3. **Direct vs. chirp-Z vs. an exact Airy disk, on a wide kernel**
   *(webSMLM 2026-09-21)*. Checks 1 and 2 used a 65x65 px @ 25 nm
   (1.6 µm) kernel. At webSMLM's default 6 µm kernel (241x241 @ 25 nm,
   NA 1.4, 660 nm, ns 1.33 / ni 1.518, unaberrated, in focus) the two
   evaluators disagree by **21% relative L2** after sum-normalization. The
   cores agree: FWHM 255.7 nm chirp-Z vs 255.6 nm direct vs 255.1 nm Airy,
   and the first zero is within 1 nm of 0.61·λ/NA_eff. The disagreement is
   in the tails:

   | Light beyond | chirp-Z | direct | exact Airy, same grid |
   |---|---|---|---|
   | 1 µm | 3.44% | 20.07% | 3.39% |
   | 3 µm (corners) | 0.17% | 16.98% | 0.157% |

   Chirp-Z is the physical one. Direct samples the pupil at N_PHI=40
   angles, which resolves exp(i·k·r·cos φ) only while k·NA_eff·r < ~N_PHI/2,
   i.e. out to ~1.6 µm here; beyond that the sum aliases into spurious
   light. After normalization every emitter's core is ~17-20% too dim.
   This is a property of the shared algorithm (same N_RHO/N_PHI), so it
   applies to demoCam's Java `Direct` port too, which is demoCam's
   default. Fix options: use chirp-Z, keep the kernel under ~2.2 µm wide,
   or scale N_PHI with the kernel radius (~128 at 6 µm; not done in either
   project). webSMLM guards its default with a regression check
   (`tests/gpu/test-sim-gpu.mjs`: chirp-Z's light beyond 3 µm < 1%).
