# OREM — Optimal Regularized re-Entry estimation Method

Optimal re-entry time prediction for resident space objects from highly elliptical orbits, using Response Surface Methodology (RSM) and Genetic Algorithm (GA) optimization with the KSROP regularized orbit propagator.

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com

---

## 1. Overview

OREM predicts re-entry times of HEO debris (GTO, Molniya, SSTO upper stages) by:

1. Processing TLE history for a target NORAD ID
2. Selecting TLE zones with a clean linear apogee-decay trend (up to 8 zones, distributed across the object's decay by R² ranking)
3. Generating mean apogee surfaces via RSM (varying eccentricity and ballistic number) for each zone
4. Optimizing (e, BN) per zone with a GA that minimizes RMS error between propagated and observed TLE apogee trajectory — with a physics-based BN floor (G2) and a trust-gated BN-range carryover between zones
5. Propagating each zone's fit with KSROP until re-entry (altitude < 80 km)
6. Reporting the **latest-zone prediction as the primary estimate** plus the all-zone ensemble mean ± spread (`output/OREM_<NORAD>_<DATE>.txt`)

Accuracy (v1.21, 7-object validation campaign, full force model): latest-zone RPE **median 2.4%, mean 4.1%, worst object 10.4%** — see `scratch_rpe/`.

---

## 2. Project Structure

```
OREM/
├── fpm.toml                         fpm package manifest; depends on KSROP (git, pinned tag)
├── src/
│   ├── propagate_ks.F               KS propagator as callable subroutine (OREM's own fork of
│   │                                KSROP's driver_KS.F core, not vendored KSROP code)
│   ├── tle_evolution.F              Batch TLE → orbital evolution
│   ├── zone_select.F                Zone selection — linear apogee decay, top-R² candidates
│   ├── ga.F                         Binary-coded GA optimizer (pop=20; trajectory-matching fitness)
│   ├── rsm.F                        RSM surface generation (9 surfaces per zone)
│   ├── orem.F                       OREM driver: pipeline + G2 BN floor + zone diagnostics
│   │                                + trust-gated BN carryover + compute_rpe
│   ├── report.F                     Prediction report writer (latest-zone primary + ensemble)
│   ├── tle_filter.F                 TLE outlier/maneuver/gap filtering (issue #10)
│   └── swx.F                        Space weather / ATM2D lookup (issue #26)
├── app/
│   └── main_orem.F                  Standalone runner (reads orem.cfg, writes the report)
│
├── input/
│   ├── const_new.dat               Physical constants
│   ├── ATM.DAT                     Atmosphere density table (60-500 km)
│   ├── example_27526.tle.txt       Ariane 5 R/B (i=17.7°, e=0.59, re-entry ~2012-05-09)
│   ├── example_32007.tle.txt       GSLV R/B (i=25.9°, e=0.29, re-entry ~2010-06-06)
│   ├── example_35497.tle.txt       Ariane 5 ESC-A (i=5.7°, e=0.63, re-entry ~2016-10-31)
│   ├── example_35497_zone2.tle.txt Ariane 5 ESC-A zone-2 (12 TLEs, e=0.60, epoch 2015-06-06)
│   ├── example_37151.tle.txt       Long March 3B (i=24.9°, e=0.56, re-entry ~2015-12-03)
│   ├── example_37819.tle.txt       Proton-M R/B (i=63.4°, e=0.47, re-entry ~2013-09-12)
│   ├── example_39615.tle.txt       Proton-M Briz-M (i=48.5°, e=0.68, re-entry ~2017-09-15)
│   ├── example_39615_zone1.tle.txt Proton-M Briz-M zone-1 (10 TLEs, e=0.68, epoch 2015-07-16)
│   ├── example_42928.tle.txt       PSLV-C39 R/B (i=19.2°, e=0.33, re-entry ~2019-02-28)
│   ├── example_42928_zone0.tle.txt PSLV-C39 zone-0 (14 TLEs, e=0.32, epoch 2017-09-22)
│   ├── example_42928_zone12.tle.txt PSLV-C39 zone-12 (12 TLEs, e=0.28, epoch 2018-01-21)
│   ├── orem_42928.cfg              Example config (IDRAG=0, fast)
│   └── orem_42928_drag.cfg         Example config (IDRAG=1, full prediction)
│
├── output/                         Prediction reports (OREM_<NORAD>_<DATE>.txt)
│
├── scratch_gmat/                   GMAT cross-validation artifacts (issues #11/#12/#25):
│                                   xval grid script, re-entry arc scripts + results,
│                                   density probe, drag_ref.py exact-integration reference
├── scratch_legacy_validation/      Ground-truth harness vs the 2017-2021 GA runs (issue #12)
├── scratch_rpe/                    7-object RPE campaigns (4-zone / 8-zone / 8-zone gated)
│                                   + ensemble_eval.py estimator comparison
│
├── test/
│   ├── test_propagate_ks.F          Propagator tests (10)
│   ├── test_tle_evolution.F         TLE evolution tests (56)
│   ├── test_zone_select.F           Zone selection tests (71)
│   ├── test_ga.F                    GA optimizer tests (74)
│   ├── test_ga_sensitivity.F        GA parameter sensitivity study (not in test suite)
│   ├── test_rsm.F                   RSM integration tests (39)
│   ├── test_sw.F                    Space weather / ATM2D tests (18)
│   ├── test_tle_filter.F            TLE filter tests (14)
│   ├── test_orem.F                  Driver + diagnostics + G2 + report tests (31)
│   ├── test_reentry.F               7-object re-entry validation (35)
│   ├── test_e2e.F                   End-to-end integration, IDRAG=1 + full force (20)
│   └── test_gmat.F                  GMAT cross-validation + exact-model drag reference (14)
└── README.md                        (382 tests total)
```

---

## 3. Propagator Interface

The core propagator is `propagate_ks` (`src/propagate_ks.F`) — a callable subroutine forked from
KSROP's `driver_KS.F` core, with OREM-specific re-entry-detection logic layered on top. It is
**not** vendored KSROP code — KSROP itself has no file by this name. Everything it calls into
(`Subrouts.F`, `Legendre.F`, `TLEread.F`) comes from KSROP as a real fpm dependency (see §4), not
from copied-and-drifted local files as in earlier revisions of this repo.

```fortran
call propagate_ks(
     &   x0, xd0, cal0,              ! Initial state + epoch
     &   nrev, istep, tole,           ! Propagation config
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   BN, IDRAG, WE_rot, EPS_f, FR_rot,  ! Drag params
     &   CR_srp, AM_srp, IPSR, ISHAD,       ! SRP params
     &   PSR_srp, amuS, amuM,               ! Constants
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,  ! Atmosphere
     &   max_pts, idump,              ! Output buffer
     &   traj_jd, traj_x, traj_xd,   ! Trajectory output
     &   exit_code)                   ! Status
```

**Exit codes:** 0 = normal completion, 1 = re-entry (alt < 80 km), 2 = divergence (NaN)

**Prerequisites:** Caller must call `init_constants()` before `propagate_ks` to populate the `/xy/` common block.

---

## 4. Building

Requires **Intel oneAPI Fortran** (`ifx`), on Windows or Linux — the only toolchain this project is actually validated against. `propagate_ks.F` currently does **not** compile under plain `gfortran` (implicit-real array dimensions, a function/array name collision on `R`, non-standard function-result assignment syntax that gfortran 13 rejects outright — see issue #28). CI (`.github/workflows/ci.yml`, issue #22) installs `ifx` on the Linux runner rather than `gfortran` for this reason.

### fpm (Fortran Package Manager) — recommended

KSROP (`Subrouts.F`, `Legendre.F`, `TLEread.F`) is no longer vendored locally — it's declared as a
real git dependency in `fpm.toml`, pinned to a KSROP release tag:

```bash
fpm build --compiler ifx --flag "-heap-arrays"
fpm test  --compiler ifx --flag "-heap-arrays"
fpm run orem --compiler ifx --flag "-heap-arrays"
```

`-heap-arrays` (16 MB stack) is required for every target that links `rsm.F`/`ga.F` — the
`surfaces(5000,3,3)` arrays overflow the default stack without it. `fpm.toml` sets
`[fortran] source-form = "fixed"` (fpm defaults `.F` to free-form). On first build, fpm clones
KSROP into `build/dependencies/ksrop/` at the pinned tag automatically (needs network access once,
then it's cached).

### Manual (without fpm)

Only needed if you're not using fpm — e.g. testing against an unpublished local KSROP checkout.
Since KSROP's library files are no longer present under this repo, point at wherever KSROP is
checked out (`../KSROP` if it's a sibling directory, matching `fpm.toml`'s dependency path during
local development):

```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx /heap-arrays /F:16777216 test\test_propagate_ks.F src\propagate_ks.F ..\KSROP\src\Subrouts.F ..\KSROP\src\Legendre.F /exe:test_propagate_ks.exe
```

`/heap-arrays /F:16777216` on Windows, `-heap-arrays` + `ulimit -s unlimited` on Linux (stack size is a shell/OS setting there, not a linker flag). See `test_all.sh` for the complete manual source-list-per-executable table (local/manual fallback only — needs a sibling `../KSROP` checkout, which a CI runner doesn't have; CI uses `fpm` instead, see `.github/workflows/ci.yml`) and `fpm.toml` for the fpm target list.

---

## 5. How to Run (Quick Start)

### Step 1: Compile

```bash
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

fpm build --compiler ifx --flag "-heap-arrays"
```

See §4 for the manual (non-fpm) build path and the Linux/`ulimit` equivalent.

### Step 2: Create a config file

Copy an example and edit:
```
cp input/orem_42928.cfg input/my_run.cfg
```

Config file format (`input/orem_42928.cfg`):
```
input/example_42928.tle.txt          <- TLE file path
42928                                <- NORAD ID
2019 3 3 0 0 0.0                    <- Observed re-entry (yr mo dy hr mn sc). Use 0 0 0 0 0 0.0 if unknown
8                                    <- Max number of zones (8 recommended: later zones sharpen the primary estimate, v1.21)
8 10.0 0.90 -1.0                    <- Zone: min_pts, max_days, R2_threshold, slope_threshold
80.0 160.0                          <- Ballistic number bounds [BN_min, BN_max] (G2 floor may extend zone 1 downward)
20 200 40 16 0.8 0.01 0.123         <- GA: pop, gen, bits_e, bits_BN, Pc, Pm, seed (pop=20 required — see v1.15)
2 0 0                               <- Force model: geo_deg, sun_deg, moon_deg
0 7.2921150d-5 3.35281066d-3 1.0    <- Drag: IDRAG(0=off,1=on), WE, EPS_f, FR
0 0.0 0.0 0                         <- SRP: IPSR(0=off,1=on), CR, AM, ISHAD
```

### Step 3: Run

```bash
fpm run orem --compiler ifx --flag "-heap-arrays" -- input/orem_42928.cfg
# or, after a manual (non-fpm) build per §4:
./orem.exe input/orem_42928.cfg
```

### Step 4: Read output

Console output shows per-zone results; a full report is written to `output/OREM_<NORAD>_<DATE>.txt`:
```
Zone  Epoch (JD)      e_opt     BN(kg/m2)  Re-entry (JD)   Re-entry (UTC)   RPE(%)    status
   1   2458152.5397   0.27666     76.78   2458735.7260   2019  9  9      48.41   ok
   ...
PRIMARY estimate (latest zone, Z 4): JD  2458526.0047  ( 2019  2 11 )
  latest-zone RPE:  -15.65 %
Ensemble ( 4 of  4 zones with a predicted re-entry):
  mean re-entry JD  2458602.1499  ( 2019  4 28 )  +-    94.05 days
```

- **e_opt / BN** — optimal eccentricity and ballistic number (kg/m²) fitted by the GA for this zone
- **PRIMARY estimate** — the latest zone's prediction: the shortest extrapolation and the freshest attitude/altitude regime, and the most accurate single estimator on the validation set (median |RPE| 2.4% at 8 zones)
- **Ensemble mean ± std** — agreement/spread indicator across all predicting zones
- **RPE(%)** — relative prediction error vs observed (if provided); **status** — per-zone diagnostic (`ok`/`boundary`/`nobound`/...)

### Notes

- Set `IDRAG=1` (line 11, first number) to enable atmospheric drag — required for re-entry prediction
- Set `IDRAG=0` for fast testing without drag (orbit won't decay)
- The `input/ATM.DAT` file must be present for drag computation
- Known re-entry date (line 3) is optional — set to `0 0 0 0 0 0.0` for operational prediction mode

### Example config files

| File | Object | Description |
|---|---|---|
| `input/orem_42928.cfg` | PSLV-C39 R/B | IDRAG=0, fast test (no re-entry) |
| `input/orem_42928_drag.cfg` | PSLV-C39 R/B | IDRAG=1, drag enabled, BN=[80,160] |

To run on a different object: copy the config, change lines 1-3 (TLE file, NORAD, re-entry date), and line 6 (BN bounds).

**Ballistic number (BN):** BN = m/(Cd×A) in kg/m². The GA optimizes BN directly, as in the original NPOE research. With the corrected J71 atmosphere (v1.17) and drag phase (v1.18), fitted values on the 7-object validation set fall in ~20–100 kg/m² per zone; the default [80,160] initial range works because the G2 physics floor automatically extends zone 1's search downward when the object's own decay rate warrants it, and later zones inherit trust-gated re-centered ranges.

---

## 6. Running Tests

`fpm test --compiler ifx --flag "-heap-arrays"` runs all suites below in one command. Or,
after a manual (non-fpm) build per §4, each binary individually:

```bash
./test_propagate_ks.exe        # Propagator tests (10 checks)
./test_tle_evolution.exe       # TLE evolution tests (56 checks)
./test_zone_select.exe         # Zone selection tests (68 checks)
./test_ga.exe                  # GA optimizer tests (71 checks)
./test_rsm.exe                 # RSM integration tests (39 checks)
./test_orem.exe                # Driver + diagnostics + G2 + report tests (29 checks)
./test_reentry.exe             # 7-object re-entry validation (35 checks)
./test_e2e.exe                 # End-to-end integration test, IDRAG=1 (20 checks)
./test_gmat.exe                # GMAT cross-validation + exact-model drag reference (14 checks)
./test_sw.exe                  # Space weather + 2-D atmosphere tests (12 checks)
./test_tle_filter.exe           # TLE quality filtering: outliers, maneuvers, gaps (14 checks)
```

**368 checks total**, all passing as of v1.24.

### test_propagate_ks
Two-body energy conservation, orbit closure, multi-revolution propagation, re-entry detection, input preservation.

### test_tle_evolution (56 tests)
- Basic: 47944 SSO (element ranges, epoch ordering, ha>hp, Sun azimuth)
- 42928 PSLV-C39 Zone 0 (i≈19°, e≈0.32, decaying apogee, RAAN regression)
- Error handling: bad file (ierr=1), no NORAD match (ierr=2)
- Geometric: ha+hp+2Re=2a identity, ha=a(1+e)-Re, hp=a(1-e)-Re
- Finite output (NaN check), perigee radius>0, ra>rp
- Sun azimuth physics: varies over time, seasonal shift
- Spot-check: epoch years, inclination values
- Large catalog: 94597-entry file, maxpts cap, field ranges, Vanguard-1 filter
- Repeatability, boundary conditions (maxpts=1)
- Deduplication: no consecutive epochs within 86 sec, duplicate removal count

### test_zone_select (68 tests)
- linfit unit: perfect linear, negative slope, constant, 2-point, noisy, 1-point
- Synthetic: linear decay, flat, oscillating, rising, empty, single, nzones cap
- Real HEO: 42928 PSLV-C39, 35497 Ariane 5, 37151 Long March 3B, 39615 Proton-M
- Zone validity: indices, non-overlapping, slopes<0, R²>0.90, min points
- Parameter sensitivity: R² threshold, min_zone_pts, max_zone_days, slope threshold
- Advanced: two-segment decay, noisy linear, step function, steep vs gradual
- Deep validation: duration limits, epoch sorting, ha>0 (all 4 objects)
- Ha monotonicity within zones (no local spikes)
- max_zone_days enforcement on all objects
- Zone count reasonableness [1,10]
- Boundary: npts=min_zone_pts exact, npts=min_zone_pts-1
- Sparse data: 30-day gaps, two clusters with gap
- Independent R² verification (manual SS_res/SS_tot)
- Degenerate: identical epochs, 2 points, very steep decay
- Repeatability, robustness (nzones_max=0, large nzones_max)

### test_ga (74 tests)
- TWOINT bilinear interpolation: constant, linear, corners, center, edges, quadratic, boundary
- Chromosome decode: all-zeros, all-ones, single-bit, asymmetric bits (60+20), non-zero lower bound
- RNG: range [0,1), different seeds, reproducibility
- Synthetic optimization: known-optimum (e=0.32, A=120), different seeds, wide bounds, 1000 generations
- ga_tinterp: t=0, t=1, midpoint, t>1 clamp, single-point surface
- Convergence: deterministic (same seed), more gen improves, fitness>0, NaN check, 5gen/1gen bounds
- Surface edge cases: corner optima, steep e-gradient, narrow bounds, few observations
- Robustness: flat surface, bounds checking, fewer generations, non-negative RMS
- High eccentricity: e~0.68 (39615), e~0.63 (35497), e~0.56 (37151), extreme asymmetric sensitivity, GEO-scale apogee, finite/bounds checks
- GA internals: ga_stats (sum/avg/max/min/ibest, equal fitness), ga_mutate (bit flip, pm=0/1), ga_iflip (prob=0), ga_irnd (range, variation), ga_fitness (perfect match)

### test_rsm (39 tests)
- jd2cal: J2000, 2017-09-22, leap year 2000-02-29, 2019-03-03 (re-entry), fractional hours, midnight
- Grid construction: e_mid/a_mid, BN=m/(Cd*A), perigee-preserving SMA adjustment
- Error handling: nzone<2 → ierr=1
- 42928 Zone 0 integration: 9 propagation runs (two-body+J2), surface physicality, center nearest obs, tobs/apobs extraction
- Surface quality: higher e → higher ha, all finite (NaN check), physical range [5k-20k km], center nearest, repeatability
- RSM→GA integration: feed real RSM surfaces into ga_optimize, e_opt/a_opt in bounds, rms valid, e_opt near TLE ecc, fitness>0.5

### test_orem (29 tests) — includes issues #12/#13
- compute_rpe: perfect RPE=0, 10-day late RPE~1.9%, mean/std (Mode 2), zero predictions
- Error handling: bad TLE file, wrong NORAD ID
- 42928 integration: full pipeline (TLE→zone→RSM→GA→propagation), e_opt physical, bn_opt physical (positive/finite — G2 floor + the corrected J71 table put fit-consistent BN below the caller's 80), rms valid, zone epochs valid
- Failure recovery/diagnostics (#12): D15 propagator-divergence skip (BN=0 forces a division-by-zero in the drag term → NaN altitude → `zone_status=1`), D16 GA boundary detection with a [20,30] window pinned *below* the real BN (`zone_status=2`; a window pinned above gets un-pinned by the G2 floor), D17 all-zones-fail doesn't crash the driver loop (`nzones_valid=0`)
- G2 physics-based BN floor (#12): 37151's zone 1 with the default `bn_min_init=80` — floor estimate extends `bn_lo` well below 80, letting `bn_opt(1)` land there (structurally impossible before this change)
- Prediction report (#13): R1–R4 real-run report (header/zone table/ensemble/legend), R5–R7 synthetic-array report exercising the with-re-entry path (PRIMARY = latest zone, latest-zone RPE line)

### test_e2e (20 tests) — Issue #16 (closed v1.20)
Full pipeline with IDRAG=1, **full force model** (geo=20, sun=2, moon=3, SRP on: Cr=1.2, A/m=0.01 m²/kg, conical shadow); GA minimizes trajectory RMS:
- E1–E5: 42928 PSLV-C39 R/B (re-entry 2019-03-03): pipeline, zones, e_opt/bn_opt physical, re-entry in all 4 zones
- E6–E10: 42928 zone-0 (14 TLEs, e≈0.32, epoch 2017-09-22): bn_opt physical, re-entry detected
- E11–E15: 39615 Proton-M Briz-M (re-entry 2017-09-15): pipeline, zones, e_opt, per-zone BN in [50,500]; E15 informational — the dedicated zone-1 file (e=0.68) is a weak-signal window whose honest fit can predict beyond the 5-year cap
- E16–E20: 35497 Ariane 5 ESC-A (re-entry 2016-10-31): pipeline, zones, e_opt, per-zone BN; E20 informational (same rationale)
- BN search range is fixed at `[bn_min_init, bn_max_init]` for every zone (no per-zone narrowing/widening) — a v1.21 trust-gated narrow/widen scheme was tried and, separately, a zone-to-zone propagated-trajectory IC-seeding scheme was also tried; both were measured against the 7- and 30-object campaigns and the propagated-IC version measurably worsened mean \|RPE\| (22.9%→27.4% on 7 objects, 22.0%→32.9% on 30, two objects blowing up over 100 points) by letting one zone's own fit uncertainty compound into the next zone's starting geometry. Reverted to always-wide search + independent per-zone SGP4-osculating seeding (`tle_find_osc`, unchanged since #31) for every zone — see issue tracking for the full investigation
- RPE printed as diagnostic; the enforced accuracy evidence lives in the 7-object campaigns (`scratch_rpe/`): latest-zone RPE median 2.4% / mean 4.1% / max 10.4% at 8 zones

### test_gmat (14 tests) — Issue #11
Cross-validates propagate_ks against GMAT R2026a reference runs (`scratch_gmat/gmat_xval_42928z0.script`, run via `GmatConsole.exe`) on 42928 PSLV-C39 R/B, Zone 0 (2017-09-24), plus a first-principles drag-magnitude reference:
- N1-N3: BN monotonicity — higher BN → less apogee decay (each of 3 e-rows)
- N4-N6: e monotonicity — higher e → higher initial apogee (each of 3 BN columns)
- N7-N9: BN sensitivity ratio decay(BN=80)/decay(BN=160) > 1.5 (propagate_ks ~2.0; GMAT ref ~1.45-1.54)
- N10: No divergence across all 9 RSM grid runs
- N11: IDRAG=0 gives < 0.5 km drop in 7 days
- N12: All drops negative for IDRAG=1
- N13-N14: Decay magnitude within ±10% of an **exact RK4 integration of propagate_ks's own drag model** at matched duration/atmosphere (`scratch_gmat/drag_ref.py`; propagate_ks agrees to ~1%). GMAT's 7-day magnitudes are printed as context only — they span 64.1 revolutions vs the 35 tested (the historical mismatch that manufactured issue #25's apparent "2× deficit") and carry J2-aliased apogee sampling and diurnal-bulge geometry a static-atmosphere model cannot reproduce
- **Key finding**: propagate_ks's drag physics is validated at the revolution level; the historical RPE bias was the ATM.DAT profile (fixed v1.17) and an arc-level drag-phase defect (fixed v1.18), not BN physics

### test_reentry (35 tests)
7 objects × 5 checks each: pipeline completion, zone detection, e_opt physical, BN physical (positive/finite — no longer bounds-checked against the caller's [80,160] input range as of v1.12, since the G2 physics-based BN floor can legitimately push `bn_opt` below 80), rms valid
- 42928 PSLV-C39 (i=19.2°, e=0.33, re-entry 2019-03-03)
- 35497 Ariane 5 ESC-A (i=5.7°, e=0.63, re-entry 2016-10-31)
- 37151 Long March 3B (i=24.9°, e=0.56, re-entry 2015-12-03)
- 39615 Proton-M Briz-M (i=48.5°, e=0.68, re-entry 2017-09-15)
- 27526 Ariane 5 R/B (i=17.7°, e=0.59, re-entry 2012-05-09)
- 32007 GSLV R/B (i=25.9°, e=0.29, re-entry 2010-06-06)
- 37819 Proton-M R/B (i=63.4°, e=0.47, re-entry 2013-09-12)

---

## 7. KSROP Dependency

`Subrouts.F`, `Legendre.F`, `TLEread.F` are no longer copied into this repo — they come from
[hari251086/KSROP](https://github.com/hari251086/KSROP) as a real `fpm` dependency, pinned to a
release tag in `fpm.toml`:

```toml
[dependencies]
ksrop = { git = "https://github.com/hari251086/KSROP", tag = "v2.1.0" }
```

To pick up a newer KSROP release, bump the `tag` here (and re-run the full test suite —
`fpm test --compiler ifx --flag "-heap-arrays"` — to confirm nothing regressed before committing).

`propagate_ks.F` (`src/propagate_ks.F`) is **not** part of this dependency — it's OREM's own fork
of KSROP's `driver_KS.F` core with re-entry-specific logic added, and evolves independently.

---

## 8. Version History

**0.1 — 2026-06-23**
- Initial repo: `propagate_ks` refactored from KSROP `driver_KS.F`

**0.2 — 2026-06-23**
- Batch TLE processing (`tle_evolution.F`)
- Epoch deduplication
- 56 tests

**0.3 — 2026-06-23**
- Zone selection (`zone_select.F`, `linfit`)
- 4 HEO TLE histories added
- `max_zone_days` bug fix
- 68 tests

**0.4 — 2026-06-24**
- GA optimizer (`ga.F`), refactored from GENESIS
- Handles high-e orbits + piecewise internals
- 71 tests

**0.5 — 2026-06-24**
- RSM surface generation (`rsm.F`), 9× `propagate_ks` calls per zone
- ATM.DAT reader fix
- RSM→GA integration verified
- 39 tests

**0.5.1 — 2026-06-24**
- Fix `propagate_ks` drag crash (KSROP #16): `ALT_atm` range guard, `H_dg÷0` safety, exp overflow clamp
- 234 total tests

**0.6 — 2026-06-24**
- OREM driver (`orem.F`) + `compute_rpe` (#6, #7)
- Full pipeline run on 42928 (4 zones)
- 7 test objects added from research data
- 14 new tests, 248 total

**0.7 — 2026-06-24**
- 7-object re-entry validation (#8), covering all orbit regimes (i=5.7°–63.4°, e=0.29–0.68)
- 35 new tests, 283 total

**0.8 — 2026-06-27**
- Fix RSM mean anomaly + time coupling: MA read from TLE (was hardcoded 0), surfaces interpolated at observation JDs
- Drag-enabled pipeline
- First re-entry detection on 42928
- 283 tests

**0.9 — 2026-06-27**
- Revert to original BN-based estimation (mass as variable, Cd=1, A=1)
- Config uses BN bounds [80,160] directly
- RSM zone-length propagation only
- 283 tests

**1.0 — 2026-07-04**
- E2E integration test with IDRAG=1 (#16): TLE→zone→RSM→GA→re-entry→RPE proven end-to-end on 42928
- Fix `test_propagate_ks` T2/T6 (per-rev dump)
- Skip re-entry propagation when IDRAG=0
- 298 total tests

**1.1 — 2026-07-04**
- NPOE cross-validation (#11): 14 tests confirm `propagate_ks` correctly models BN sensitivity (ratio ~2.0 vs NPOE's 2.02) and apogee decay direction
- Magnitude is ~50% of Jacchia-70 (ATM.DAT vs Jacchia model)
- RPE inaccuracy diagnosed as a short-zone/noise issue, not a propagator bug
- 312 total tests

**1.2 — 2026-07-04**
- Fix NaN in RSM propagation:
  - `car2oe` clamps all `dacos()` arguments to [-1,1] — floating-point overflow at orbital perigee caused NaN true-anomaly → NaN drag → NaN state in ie=2,3 RSM surfaces
  - `rsm_generate` hardcodes IDRAG=1 — without drag, all 9 RSM surfaces were identical and the GA had no BN signal
- 312 tests still pass

**1.3 — 2026-07-04**
- Add zone-0 E2E run (E6–E10) in `test_e2e.F` using `example_42928_zone0.tle.txt` (14 TLEs, e≈0.32, epoch 2017-09-22)
- Zone-0 RPE = −55.5% vs −87 to −96% for late zones — confirms improved accuracy when propagating from an early orbit
- 317 total tests

**1.4 — 2026-07-04**
- Fix GA array-dimension mismatch bug: `ga_optimize`/`ga_fitness` declared `surfaces` with leading dimension `nsurf_pts` (≈nobs≈26) but callers allocated `surfaces(max_surf=5000,...)` — all surface reads were reading wrong memory, so the GA always returned the lower bound regardless of the fitness landscape
- Fix: add `ld_surf` parameter to `ga_optimize`/`ga_fitness`; callers pass `max_surf`
- Add E11–E20 tests for 39615 and 35497 with zone-specific TLE files
- Zone-0 RPE improves from −55.5% to −16.1%
- 327 total tests

**1.5 — 2026-07-04**
- Multi-zone campaign fitting (#12): replace per-zone independent GA with a single campaign GA that finds one BN consistent across all zones simultaneously
- `ga.F`: add `ga_campaign` (1-D BN-only GA) and `ga_camp_fitness` (mean RMS across valid zones)
- `orem.F`: three-phase Step 5 — Phase 1 generates RSM for all zones on the full BN range [bn_min, bn_max] (no iterative narrowing), Phase 5b runs the campaign GA, Phase 5c propagates re-entry from each zone using the shared campaign BN
- Eccentricity fixed at the zone TLE midpoint
- 327/327 tests pass

**1.6 — 2026-07-04**
- Slope-based BN estimation (#12): replace instantaneous-apogee GA with dha/dt slope fitting
- Lunisolar oscillations (~2-day period) cancel in the linear-regression slope over 10-day zones; drag contribution is secular and accumulates (BN=80 → −2.2 km/day vs BN=160 → −1.1 km/day, 2× signal)
- `ga.F`: add `ga_slope_optimize` (2-D GA matching apogee-rate slope) and `ga_slope_fitness` (bilinear-interpolates `surf_slopes` at (e,BN), compares to obs slope)
- `orem.F`: Step 5 now computes `surf_slopes(3,3)` via linear regression of each RSM surface column, `obs_slope` from TLE apogee history, then calls `ga_slope_optimize` per zone with BN narrowing
- Campaign fitting (v1.5) removed — per-zone BN is physically correct since BN encodes attitude (BN = M/CdA, and A depends on attitude, which varies per zone)
- 327/327 tests pass

**1.7 — 2026-07-05**
- Revert to the original trajectory-matching fitness (`genpoen1.f` algorithm): Step 5 calls `ga_optimize` (not `ga_slope_optimize`), matching RMS of the propagated apogee trajectory against all TLE observations in the zone — identical to the published NPOE research fitness function
- Slope-fitting (v1.6) caused the GA to saturate at the BN lower bound; trajectory matching restores correct BN identification (Z1: BN≈151 vs 80 in v1.6)
- Zone-0 RPE −16% confirmed
- `ga_slope_optimize` kept in `ga.F` for reference
- 327/327 tests pass

**1.8 — 2026-07-12**
- Sync `ksrop/` with the KSROP GMAT validation campaign fixes (KSROP #18–#23):
  - `aLegP` buffer-overflow rewrite (`Legendre.F`) — old version ignored its degree argument and wrote ~50× out of bounds on every call
  - `aleg`/`sleg`/`oleg` off-by-one — force/time-element formulas need degree n+1, so `aLegP(n+1,...)` at both `propagate_ks.F` call sites (live in the pipeline at ngeo_deg=50)
  - `Tau_geo` sign/scale fix (missing `amue`; thesis eq. 2.56) — epoch labeling only
  - `third_body_aux`'s `deg` was implicitly declared `double` instead of `integer` while callers pass integers, so its power-series loop ran zero times and **the third-body force was silently exactly zero**
  - `qsun`/`qmoon` rewritten in the correct KS-elements EOM convention `shape·u + r·Lᵀ(∇shape)`, verified vs `KSJLSDNP.F` to machine precision
  - `solarnpv`/`lunarpv` upgraded to Montenbruck & Gill analytic series — Sun 0.6%→0.097%, Moon 3.6%→0.109% vs DE405
- KSROP-side GMAT validation: Sun-only GTO 1.2 m/rev, Moon-only 0.46 km/rev, full conservative 1.9 km/2 revs
- OREM pipeline currently runs `nsun_deg=nmoon_deg=0` (const_new.DAT), so the third-body/ephemeris fixes are dormant until lunisolar is enabled — but the geopotential fixes are active
- 327/327 tests pass

**1.9 — 2026-07-11**
- Replace NPOE with GMAT as the issue #11 cross-validation reference (`test_npoe.F` → `test_gmat.F`) — NPOE's own atmosphere model confounded the comparison (~50% magnitude gap unrelated to force-model correctness) and only proved consistency with NPOE-era heritage research, not physical correctness
- New reference trajectories generated by `scratch_gmat/gmat_xval_42928z0.script` (GMAT R2026a via `GmatConsole.exe`, no GUI), using the same force-model settings the KSROP↔GMAT validation campaign already trusted (EGM2008 zonal Degree=4/Order=0, Sun+Luna point mass, JacchiaRoberts F107=72/Kp=1.0 matching ATM.DAT)
- N13/N14 magnitude tolerance tightened from factor-3 (NPOE) to ±50% (GMAT), based on the observed 69–112% residual across all 9 grid points
- 327/327 tests pass

**1.10 — 2026-07-11**
- Implement issue #12's original scope: RSM/GA failure recovery and diagnostics
- `orem_run` gains two new outputs, `zone_status(nzones_max)` (0=ok, 1=skip_propfail, 2=boundary, 3=nobound, 4=skip_toofewpts) and `nzones_valid` — purely additive, no existing output array or `ierr` semantics changed
- Detects: propagator divergence (already-existing RSM skip, now tagged + warned); GA optimum within 15% of a search bound (warns and **widens** — not narrows — the next zone's BN range, since narrowing further would entrench a likely-wrong search window); RSM envelope not bounding an observation (warn only, still runs GA); zone TLE count below a fixed RSM-reliability floor of 3 (defensive, unreachable until TLE filtering #10 lands)
- Real production runs surface a genuine finding: several zones across the 7-object validation set have the RSM envelope failing to bound most of their observations (e.g. 42928 Z4: 21/22; 37819: ~30/31 in every zone) — direct evidence for the RPE-inaccuracy problem #12's BN-identifiability work has been chasing
- The boundary-widen behavior also causes BN to escalate without a ceiling for objects that keep hitting it (37151: 151→198→268→373 across 4 zones, vs. the old narrow-only 151→167→174→178) — matches the issue's literal spec, flagged as a follow-up consideration rather than capped here
- New tests D15–D17 in `test_orem.F` (BN=0 forces a division-by-zero in the drag term → deterministic propagator divergence, tests both the single-zone skip and the all-zones-fail path)
- 333/333 tests pass

**1.11 — 2026-07-11**
- Widen `test_e2e.F` (issue #16) to the full force model: geopotential degree 4→20, SRP enabled (previously fully off — Cr=1.2, A/m=0.01 m²/kg, conical shadow, PSR=4.56e-6 N/m² at 1 AU)
- Finding: RSM/GA-fitted BN and e (from the short ~7–10 day zone-fitting window) are unchanged to the last decimal across all 4 objects — drag dominates apogee decay at that timescale, and SRP/J5–J20 zonal terms don't move the fit
- The long-duration re-entry propagation (up to 5 years) *is* measurably affected — zone-0's re-entry trajectory shortened from 1551 to 1504 revolutions, predicted re-entry date ~4.7 days earlier (RPE −69.80%→−70.68%)
- The main 4-zone run's own re-entry propagations are short enough (72–314 revolutions) that the perturbation doesn't shift which revolution crosses the 80 km threshold, so those stayed byte-for-byte identical
- Confirms BN fitting is a drag-only problem at this timescale, independent of what else is enabled
- 20/20 tests still pass, 333 total unchanged

**1.12 — 2026-07-11**
- Implement algorithm-review finding G2 (issue #12): physics-based BN floor
- New `estimate_bn_floor` in `orem.F` fits zone 1's own TLE decay rate (`linfit(epz, smaz, ...)`, reusing `zone_select.F`'s existing routine), then numerically calibrates it against one short (~0.5–3 day) trial `propagate_ks` run at a reference BN=100 kg/m² — decay rate ∝ 1/BN, so `bn_floor_est = 100 × |trial slope| / |observed slope|`
- Numerical calibration (not a closed-form formula) chosen specifically so the estimate can't disagree with `propagate_ks`'s own internal drag/unit conventions
- Applied as a **floor-only safety net**: only ever extends `bn_lo` downward from the caller's `bn_min_init` for zone 1 (`bn_lo = min(bn_lo, 0.5×estimate)`), never touches `bn_hi`, never raises the floor — zones 2+ keep the existing v1.10 narrow/widen carryover untouched
- Validated against real data: fires correctly for 4/7 `test_reentry` objects (e.g. 37151: floor extends to 23.10, `bn_opt(1)` reaches 48.79 — previously impossible below the hardcoded 80 floor); for 35497, the object-level "zone 1" and the issue-referenced "zone 2" TLE file give different physics estimates (266.7 vs 46.5, only the latter fires) — a real, expected consequence of BN varying with attitude/altitude regime across an object's life, not a calibration bug
- Caught and fixed a bug during implementation: the floor logic initially reused the narrow/widen step's `if (bn_lo < 10) bn_lo = 10` safety clamp, which unconditionally raised `bn_lo` and broke the v1.10 BN=0 divergence tests (D15/D17) — removed, since a floor-only guarantee must never raise what the caller passed
- Updated `test_reentry.F`'s stale `bn_opt ∈ [80,160]` assertions (no longer valid once the floor can legitimately go below 80) to a physical sanity check
- New tests: `test_orem.F` G2 section (2 checks)
- 335/335 tests pass

**1.13 — 2026-07-12**
- Ground-truth validation harness for issue #12 (analysis only, no pipeline code changed)
- New `scratch_legacy_validation/compare_bn.F` runs `orem_run` (BN search widened to [10,300], single-zone mode) against the original 2017–2021 GA fitting run's own raw TLE slice, for all 31 zones across the 7-object validation set (sourced from `E:\Research\1. R&D\Re-entry\KSROP-DataPrint`)
- Compared against the legacy `GA/gene9.txt` search bounds, `GA/genesis.dat` best-fit checkpoints, and the independent non-GA `Non-Opt KS/Non-opt.txt` baseline
- Finding: only 6 distinct BN values appear across the 29 zones that returned a result — `205.6167` alone in 20 of 29, tracking `zone_status=nobound` (and, in 2 cases, even `ok`) almost exactly
- With a fixed GA seed and search range, a flat/uninformative RSM fitness landscape reproduces the same decoded chromosome regardless of input; only 2 of 29 fits landed inside the legacy GA's own search box
- Generalizes the existing "RSM doesn't bound observations" note (previously 2 objects) to the dominant outcome across the full 7-object set
- Recommendation: before choosing between G3 (BN floor/range tuning) and G4 (zone distribution), investigate why `rsm_generate`'s fitness landscape goes flat this often — neither fixes a search with no signal
- Caveat: this test forces single-zone mode with a much wider range than OREM's normal chained [80,160]-narrowing operation, so the harness itself as a contributing factor wasn't yet ruled out
- No test count change (no pipeline code touched)

**1.14 — 2026-07-13**
- Root-caused and fixed v1.13's "flat RSM landscape" (issue #12)
- `rsm_generate`'s local propagator scratch buffers (`traj_jd`, `traj_x`, `traj_xd`, `traj_ha`, `traj_t` in `rsm.F`) are declared `SAVE` and were never cleared before each `propagate_ks` call — a grid point whose trajectory is shorter than a previous one (an earlier grid point in the same 3×3 loop, or a previous zone/object entirely) inherits that previous run's leftover trajectory tail
- The end-of-data sentinel only excluded stale points *later* than the current zone by more than `zone_dur+0.5` days — never stale points *earlier* than it — so old data silently bled into the new zone's apogee envelope
- Proved with a new scratch diagnostic (`scratch_legacy_validation/diag_rsm.F`, not part of the build) that calls `rsm_generate` directly and shows order-dependence: run after object 37151, object 42928's envelope showed `smax` pinned at exactly 14833 km (37151's apogee, not 42928's own ~6500 km) at 8 of 10 observations; run first, the contamination vanished
- Fix: zero `traj_jd(1:max_traj)` immediately before every `propagate_ks` call inside `rsm_generate`'s grid loop, so the sentinel can never see a prior call's data regardless of its epoch
- Re-ran the v1.13 31-zone harness after the fix: `zone_status=nobound` dropped from 20/29 zones (69%) to **0/29** — the RSM envelope now bounds every observation in the validation set — and RMS on the recomputed 27526 zones improved by 1–2 orders of magnitude (e.g. Z3: 24.6→0.49, Z4: 40.6→0.48)
- The remaining spread (8/29 landing on the same GA-decoded BN=205.6167, 20/29 flagged `boundary`, mostly near the widened 300 kg/m² ceiling) now looks like real residual signal rather than corruption — informs the G3-vs-G4 decision directly
- 335/335 tests pass (no test-visible behavior change on the existing suite — the bug required a specific multi-call/short-trajectory sequence the unit tests didn't happen to trigger)

**1.15 — 2026-07-13**
- GA population raised from 4 to 20 at all pipeline call sites (issue #12) — experiments on the 31-zone ground-truth harness proved pop=4's output is a *range-invariant seed artifact*, not an optimum:
  - Parameterized `scratch_legacy_validation/compare_bn.F` with optional args (`bn_hi`, `popsize`, output CSV); defaults byte-for-byte reproduce the committed v1.14 baseline
  - Smoking gun: raising the BN ceiling from 300 to 600 at pop=4 moved the pinned value 267.7896 → 534.4685 and the decoy 205.6167 → 407.9788 — both are the *identical binary chromosome* decoded over the wider range (fraction 0.888930 / 0.674540 of either interval). A 4-individual population searching a 56-bit chromosome converges on a seed-determined decode with zero influence from the data; boundary flags got worse (20/29 → 26/29)
  - At pop=20, [10,300]: all 29 zones land on distinct data-driven BN values, `ok` zones 9→20, median GA RMS 0.164→0.058 (better in 24/29 zones), and 42928 Z0 fits BN=151.04 — matching the v1.7 chained result (151.11) and the heritage research
  - Ceiling stays at [10,300]/[80,160]: at pop=20 the wider [10,600] range *degrades* the search (42928 Z0 drifts to 397 with worse RMS) — drag fitness goes flat at high BN, so extra range is noise for the same generation budget
- Changes: `test_e2e.F` (4 sites), `test_orem.F` (7), `test_reentry.F` (1, + widened `F5.2`→`F7.2` BN print format that overflowed at BN≥100), `test_rsm.F` (1), both `input/*.cfg` files, `ga.F` doc comment. `test_ga.F` deliberately stays at pop=4 — its tests exercise GA mechanics, not the production config
- 335/335 tests pass with **zero assertion changes** — including the two sensitive ones: D16 (boundary detection, [200,205] pinned window still traps the optimum at the edge) and G2 (37151 floor: pop-20 GA lands bn_opt(1)=56.79, still below 80, near the physics estimate of 46)
- 7-object fits are now coherent per object: 37151 = 48–91, 27526 = 73–113, 32007 = 77–126, 39615 = 126–200 kg/m²; 42928/35497/37819 start pinned at the top of [80,160] and escalate zone-to-zone (up to 300–438) via the v1.10 boundary-widen carryover — their true zone-1 BN is at or above the caller's ceiling
- **Key negative result: RPE is essentially unchanged** (42928 zone-0: −70.68% → −73.78%; 4-zone best: −72.35%; 39615: −97.68%). The broken optimizer was *masking* the real remaining problem, not causing it — with within-zone fits now excellent (RMS ~0.06) yet long-horizon predictions still 70–97% early, the RPE error must come from downstream of the fit: a BN fitted on a ~10-day zone under-predicting the months-long decay (attitude/regime drift), the static J70 density vs. the real solar cycle, or the re-entry propagation config itself — that's the reframed #12 investigation

**1.16 — 2026-07-14**
- GMAT re-entry cross-check (issue #12, analysis only — no pipeline code changed): new `scratch_gmat/gmat_reentry_42928z0.script` propagates 42928 zone-0's v1.15-fitted state (e=0.3216, BN=117.60 kg/m², same epoch/elements as the issue #11 grid script) all the way to a 102 km altitude threshold in GMAT, twice: JacchiaRoberts with constant F10.7=72/Kp=1.0 (the exact static weather ATM.DAT was generated with), and JacchiaRoberts with the actual 2017–2019 flux history (GMAT's shipped CSSI `SpaceWeather-All-v1.2.txt`)
- Result (`gmat_reentry_42928z0_results.txt`): OREM predicts re-entry at epoch +135.3 d (−74% vs observed); GMAT const-flux at +736.4 d (+40%); GMAT real-flux at +723.9 d (+38%); observed 2019-03-03 = +524.1 d
- **Verdict 1 — the ATM.DAT density profile is the dominant RPE error.** Same state, same BN, same static weather: lifetime differs 5.4× purely on the density model. Within the ~160–170 km zone-fitting window the models roughly agree (v1.9 measured 69–112% over 7 days), so the fit can't see the problem — but as perigee descends, the J70 table's density rises far faster than JacchiaRoberts', and the arc-integrated decay runs ~5× too fast. Also explains the v1.8 regression (zone-0 RPE −16% → −70% when ATM.DAT switched to J70): the old table was too *thin* low down, partially cancelling the profile-shape error
- **Verdict 2 — solar-activity history is irrelevant for this object** (736.4 vs 723.9 d, a 1.7% shift): 2017–2019 was solar minimum with actual F10.7 ≈ 68–70, right at the static 72 assumption. #14's "dynamic solar-activity scaling" would matter for objects decaying across a solar maximum, but it is not the current bottleneck
- **Verdict 3 — the zone-fitted-BN-extrapolation hypothesis is refuted as the main driver**: if the ~10-day-window BN were the problem, GMAT would also have predicted early; instead it lands +38–40% *late* (JacchiaRoberts' own bias for this object; the observed date sits between the two models)
- Two GMAT gotchas documented in the script: JacchiaRoberts throws below 100 km altitude (stop threshold set to 102 km — terminal decay, days from 80 km), and the CSSI file's daily flux steps violate RKV89 Accuracy=1e-12 (real-flux propagator relaxed to 1e-10, MinStep 1e-6)
- Next step for #12: fix the ATM.DAT low-altitude profile — regenerate with proper J70 lower-boundary handling (90–125 km fixed-temperature region) or replace the tabulated-exponential lookup with an analytic Jacchia implementation in `propagate_ks`, then re-fit and re-measure RPE
- No test-count change (335)

**1.17 — 2026-07-14**
- **Regenerated `input/ATM.DAT` with the real Jacchia-71 profile** (issues #12/#14): new `KSROP/gen_atm_jr71.F` replaces `gen_atm_j70.F`, whose hand-rolled single-exponential temperature profile (forced 12 K/km gradient at 90 km) ran ~127 K too warm through the 90–125 km region (458 K at 125 km vs J71's Tx=328 K) and inflated the hydrostatic column above it — quantified by a new GMAT density probe (`scratch_gmat/gmat_density_probe.script`): the old table was **3.3–3.5× denser than JacchiaRoberts across the 140–200 km perigee band** (and 0.4× at 300 km), while agreeing at ~100 km
- New generator implements the J71 structure with the Roberts-1971 polynomial anchors (quartic T(z) 90–125 km with zero gradient at the 90 km minimum; exponential-asymptotic above; δij species anchors at 125 km; ζ total-density anchor at 100 km; Aa mean-molecular-mass barometric 90–100 km), cross-checked coefficient-by-coefficient against the SatelliteToolboxAtmosphericModels.jl reference implementation and validated in a Python prototype to 0.6–6% against GMAT JR over 102–300 km. The generated table tracks GMAT JR at 0.80–0.95 across 102–300 km (residual = diurnal factor; table is nighttime-minimum static, T∞=626.3 K); SCH column is now the local density scale height −dz/d ln ρ (what the King-Hele model actually needs), not R·T/(M̄g)
- **Fits transformed**: 42928 now fits BN=67–75 across all four zones (was 153→299 escalating), 35497 74–98 (was 160→402), 37819 86–98 (was 160→438) — the v1.10 boundary-widen escalation is gone because the fit-consistent BN now sits inside the searchable range. RPE improves but less than the v1.16 arc cross-check predicted: zone-0 −73.78% → **−53.40%**, 4-zone best −72.35% → −72.58%, 39615 −92.52%
- **New finding, filed as issue #25**: `propagate_ks`'s King-Hele drag produces ~2× less 7-day decay than an exact RK4 integration of its own stated model at matched ρ_p/H (`scratch_gmat/drag_ref.py`: −16.2 vs −37.0 km on the N13 case), yet over the full re-entry arc behaves ~3× stronger than that factor implies (242-day zone-0 prediction where the matched-drag GMAT scaling predicts ~700–800) — a non-constant model deficit that does *not* cancel through the BN fit and is now the dominant remaining RPE error. v1.9's N13/N14 "agreement" with GMAT was this deficit cancelling against the old too-dense table
- Assertion updates forced by the corrected table (all documented in-code): D12/E4/E9 BN checks → physical-sanity (G2 floor + thinner table put fit-consistent BN below the caller's 80); D16 boundary-detection window moved from [200,205] (pinned above the old table-consistent BN) to [20,30] (pinned below — the G2 floor un-pins windows from above); N13/N14 tolerance → honest [0.15,0.60]×GMAT band encoding the known factors, with a do-not-rewiden pointer to #25
- 335/335 tests pass. KSROP-side: `gen_atm_jr71.F` + regenerated `input/ATM.DAT` committed to KSROP `HS-dev` separately

**1.18 — 2026-07-14**
- **Fixed the #25 drag-phase defect; RPE collapses from −72..−97% to bracketing zero.** 42928 4-zone RPE now +52.2/+14.6/**+3.2**/−13.8% (ensemble mean re-entry within **+11%** of observed on a 527-day horizon, ±100 d spread); zone-0 +37.1% (was −53.4%)
- Root cause (two-part, and *not* what #25 originally claimed):
  - The issue's "2× deficit at the 7-day window" was a **test artifact**: N13 compared a 35-revolution `propagate_ks` drop against GMAT references spanning 7 days = 64.1 revolutions of that orbit. Duration-matched, `propagate_ks` agrees with an exact RK4 integration of its own drag model to **~1%** (BN=80: −16.24 vs −16.45 km; BN=160: −8.12 vs −8.21; `scratch_gmat/drag_ref.py`, now with the F co-rotation factor and oblate-perigee density matched)
  - The *arc-level* distortion was real: the old analytic eccentric-anomaly sweep (`DE_dg = (VIPP·π − EA₀)/istep`) advances the drag-density phase at **half rate whenever a revolution starts past EA=π** (the `VIPP=4` branch targets 4π over one rev's steps) — intermittently dephasing the density peak from the true perigee passage along every long decay arc as revolution boundaries drift. A phase error that comes and goes with orbit geometry cannot be absorbed by the fitted BN, which kept RPE pinned deep-negative through four generations of upstream fixes
- Fix: the drag density now reads the **true eccentric anomaly from the state** (`pek(7)`, refreshed by the per-stage `car2oe`) instead of the analytic sweep — stage-accurate, covers both wings of the perigee density bump, no steps-per-rev assumption. Bit-identical on the N13 window (perigee-anchored revs never trigger `VIPP=4`), transformative on arcs
- Fits: per-object BN in coherent physical bands across all 7 objects (42928: 45–75, declining with zone; 39615: 58–65; 37819: 64–73; 32007: 38–66). 39615's dedicated zone-1 file (e=0.68, weak drag signal) fits BN≈139 and honestly predicts >5-year lifetime vs the actual 2.2 — the true #12 identifiability limit, no longer masked; E15 made informational with the E20/35497 rationale
- N13/N14 rebased onto the first-principles reference at matched duration (±10% bands around −16.45/−8.21 km); GMAT 7-day magnitudes demoted to context output (they also carry J2-aliased osculating-apogee sampling and diurnal-bulge geometry a static-atmosphere model cannot reproduce)
- 335/335 tests pass. KSROP's `driver_KS.F` carries the same heritage sweep (KSJLSDNP2 lineage) — porting this fix there is flagged KSROP-side

**1.19 — 2026-07-14**
- **Structured prediction report (issue #13)**: new `report.F` / `orem_report` writes `output/OREM_<NORAD>_<DATE>.txt` from `orem_run`'s outputs — config echo, per-zone fit/prediction table (epoch, e_opt, BN, re-entry JD + UTC date, RPE, `zone_status` label), and the headline **ensemble block**: mean re-entry ± std, relative spread (% of the zone-1→mean horizon), and ensemble RPE when an observed date is provided. Wired into `main_orem.F`; new R1–R4 tests in `test_orem.F`. **339 total tests**
- Fixed a latent `main_orem.F` bug found during wiring: its `orem_run` call was never updated for v1.10's `zone_status`/`nzones_valid` arguments (positional mismatch — `rpe` landed in `zone_status`'s slot). Latent only because `orem.exe` had not been rebuilt since; now threaded and rebuilt
- **First full 7-object drag-enabled RPE campaign post-v1.18** (`scratch_rpe/rpe_campaign.F`, full force model, results in `rpe_campaign.csv`): ensemble RPE per object — 42928 **+15.3%**, 35497 +238.7%, 37151 −7.4%, 39615 **+8.4%**, 27526 +20.4%, 32007 **+0.7%**, 37819 −17.7%. **Six of seven objects within ±21%** (median |ensemble RPE| 15.3%); the outlier 35497 is the known i=5.7° solar-apsidal-resonance case that issue #9 (3-variable optimization with inclination) was written for — its zone-4 alone predicts −1.1%, while early zones run +170..+520%
- Design signal for #12/G4: **the latest zone is consistently the sharpest single predictor** (35497 Z4 −1.1%, 37151 Z4 the only zone predicting at all, 42928 Z3/Z4 best) — drag signal concentrates as perigee decays, motivating recency-weighted ensembles or late-zone selection over the current uniform mean

**1.20 — 2026-07-14**
- **Latest-zone primary estimate (closes #16's <10% accuracy target)**: offline evaluation of five ensemble schemes against the 7-object campaign (`scratch_rpe/ensemble_eval.py`) — uniform mean, latest-zone, index-weighted, inverse-remaining-lifetime-weighted, median. **Latest-zone wins decisively: median |RPE| 8.2%, mean 7.6%, worst object 14.4%** (vs uniform mean's 45.3% mean / 238.7% max). Every object within ±15%; even the 35497 resonance outlier lands at −0.2%. Physical basis: the latest zone has the shortest extrapolation and the freshest attitude/altitude regime
- `orem_report` ensemble block now leads with **"PRIMARY estimate (latest zone, Z n)"** + its RPE, followed by the uniform mean ± std as the spread/agreement indicator
- New tests R5–R7 (synthetic zone arrays, exercising the with-re-entry report path without propagation). **342 total tests**
- #16 closed: the E2E chain is proven and the <10% target is met by the primary estimator (5/7 objects <10%, mean 7.6%). Remaining accuracy work continues under #12 (weak-signal zones: 37151 −14.4%, 27526 +10.8%) and #9 (35497's inclination resonance)

**1.21 — 2026-07-14**
- **Trust-gated BN-range carryover + 8-zone operation: latest-zone RPE now median 2.4% / mean 4.1% / max 10.4% across all 7 objects**
- Measurement first: re-running the campaign with `nzones_max=8` (zone_select returns the top-R² candidates, so a higher cap admits *later* zones) sharpened the latest-zone estimator wherever signal exists (42928 −4.4→0.0%, 39615 8.9→2.7%, 32007 6.4→0.8%) but regressed 37151 (−14.4→−38.0%): its Z1–Z7 all predict no re-entry (weak-signal fits), yet each re-centered the v1.10 BN-range carryover, marching the search from [12.5,160] down to [17.2,24.9] and imprisoning Z8 — the only zone with real signal
- Fix in `orem.F`: the carryover now chains **only from zones whose fit carries real signal** — with drag on, a zone that actually predicted a re-entry; with drag off, an unflagged (`zone_status=0`) zone. Untrusted zones leave the range unchanged. Objects whose zones all predict (42928, 35497, 37819) chain exactly as before — bit-identical e2e results
- Gated 8-zone campaign (`rpe_campaign_8zone_gated.csv`; 4-zone and ungated-8-zone runs preserved alongside): 42928 0.0%, 35497 0.6%, 39615 2.0%, 32007 2.4%, 37819 −5.3%, 37151 **−8.1%** (recovered), 27526 10.4% — **all seven at or under ~10%**
- Shipped configs raised to `nzones_max=8`. 342/342 tests pass (test suites unchanged — they run 4-zone IDRAG=0 paths whose chains are gated on `zone_status` and unaffected in practice)

**1.22 — 2026-07-17**
- **New-object case study: 33587 (1989-039EF, i=65.2° Molniya-class fragment; analysis only, no pipeline code changed).** 306 TLEs spanning 2009 → 2022-12-03, observed decay 2025-04-22 — the prediction must extrapolate 2.4 years past the last TLE, across the 2023–2025 solar maximum. Artifacts: `input/example_33587.tle.txt`, `input/orem_33587.cfg`, `scratch_rpe/*33587*`
- Pipeline result (honest): five zones found (Jul–Oct 2022), **no zone predicts re-entry within the 5-year cap** — under both the shipped F10.7=72 table and a regenerated F10.7=150 variant (`scratch_rpe/ATM_F150.DAT`, T∞=879 K)
- Diagnosis, part 1 — **no drag signal in any fit window**: zone perigees sit at 440–640 km, where a ≤10-day window carries essentially no BN information (fits are table-invariant GA noise; only Z5, hp≈444 km, responded to the 1.5× table change). This object's decay was driven by **lunisolar perigee cycling** (hp 616→341 km in the last five TLE months — far too fast for drag), which hands the orbit to the atmosphere only after the TLE record ends
- Diagnosis, part 2 — **static-atmosphere lifetime error across a solar maximum**: direct long propagation from the last TLE state (`scratch_rpe/prop_33587.F`, 60k-rev cap, full force, F150 table) does re-enter — the modeled lunisolar cycling works — but 4.6–5.7× too slowly: 4044–4929 days vs the observed 870 (RPE +365% to +467% across BN 40–120). The 2023–2025 arc averaged F10.7 ≈ 160–180 with major geomagnetic storms; a static quiet-condition table cannot represent it
- Conclusion: 33587 is **out of scope for the static-atmosphere OREM** — it is the concrete motivating case for #14 (dynamic space-weather along the arc) and exercises the object class where zone selection needs a drag-signal criterion (hp-aware zone quality). The 7-object validation set's accuracy (median 2.4%) is unaffected: those objects' windows are drag-dominated and their arcs mostly avoid solar-max crossings

**1.23 — 2026-07-17**
- **Epoch-resolved space weather implemented (issue #26)**: `input/SW-All.csv` (CelesTrak daily history 1957→present + monthly predicts to 2041; refresh via curl) + `input/ATM2D.DAT` (new `KSROP/gen_atm2d_jr71.F`: the J71 profile over a 550–1500 K T∞ grid, 291×39; profile functions shared via new `KSROP/jr71_profile.F`, 1-D generator verified bit-identical after the split)
- Runtime: `sw_tinf` (JD → T∞, binary-searched — predicted era is monthly) and `atm2d_interp` (bilinear ρ/H in legacy scaled units) live **inside `ksrop/propagate_ks.F`**, hooked into the per-revolution drag reference; loaders (`sw_load`/`atm2d_load`, new `swx.F`) are linked only by opt-in executables, so all legacy builds and results are bit-unchanged. `orem.exe` auto-detects both files and states ENABLED/DISABLED loudly
- New `test_sw.F` (12 checks incl. the hand-verified 2024-05-11 G5-storm T∞=1216.56 K and the W12 smoke test: 55.4 vs 31.1 km 7-day decay at storm vs minimum epochs). **354 tests total**, 342 legacy checks unchanged
- **33587 verdict overturned on attribution**: weather-enabled arcs are nearly identical to the F150-static arcs (T∞=879 K was already a fair proxy for the 2023–2025 average) — the 5× lifetime error is **not density**. The in-record diagnostic (`scratch_rpe/prop_33587_hp.F`) proves it: over the record's last 131 days the observed perigee descends 616→341 km (lunisolar cycle) while the modeled perigee stays flat (±20 km, wrong direction). Secular third-body eccentricity evolution is missing for this critical-inclination orbit — per-rev GMAT validations could never see it. Filed as **#27** (P1) with the ±30% acceptance transferred; GMAT hp(t) comparison specified as the decisive next experiment
- Follow-up before #26 closes: 7-object weather-mode regression campaign

**1.24 — 2026-07-18**
- **Three array-bound bugs fixed**, exposed by validation objects with much longer TLE histories than the original 7-object set (11550: 46-year record, 10,143 deduped TLEs): `orem.F` `maxpts` 10000→15000 (silently truncated the TLE fill loop for long-history objects, no error); `zone_select.F` `max_cand` 100→2000 (Pass-1's candidate buffer silently *stopped recording* — not just deprioritizing — once 100 windows were found; for 11550 the entire 2010–2025 terminal-decay region, including the zone where mean perigee collapses to 96 km at the last TLE, was invisible to the algorithm); `main_orem.F` `mxz` 10→50 (workaround, not a structural fix — recency-biased zone selection is a real open design gap for long-history objects). `report.F` gained an Epoch (UTC) column. 354/354 tests pass.
- **Three new real-decay validation cases** (11550/59347/40943, all real ground-truth decay dates except 11550's inferred one): 59347 and 40943 landed within 1–3 days of the true date (RPE 0.17%/0.23%); 11550 is the **first clear counterexample** to the latest-zone-as-primary heuristic (v1.20) — it underperformed the ensemble mean on both a full-history run and a genuine out-of-sample run withholding all of 2025.
- **Issue #26 (7-object weather regression campaign) run — mixed result, kept open.** 4/7 objects improve under epoch-resolved weather (32007 17.3%→5.6%), 2/7 regress notably (42928 0.02%→−18.3%, 37151 10.0%→31.9%), aggregate mean/median |RPE| both tick up slightly. Not a clean "no regression" pass. Also surfaced: the *static* baseline itself no longer reproduces the 1.21-era published numbers (median 2.4%→10.0%) against the current tree — most visibly 37151, whose zone-8 selection changed once the `max_cand` fix stopped hiding candidate windows. That drift is real and unexplained; flagged, not resolved, here.
- **Issue #27 BN-refit experiment run — deepens rather than resolves the question, kept open.** Fit `(e, BN)` directly against 33587's real TLE observations across the exact 131-day in-record window the GMAT decisive experiment used, wide BN search `[1,150]`. Result: `bn_opt=126` (unremarkable, not boundary-pinned) with RMS ≈ 468 km — 10–100× worse than the validation set's good fits — i.e. **no physically plausible BN explains the observed collapse via King-Hele drag either**, on top of the prior GMAT finding already ruling out third-body truncation. New lead: the record's last TLE has a mean-motion derivative ~143× the first TLE's — possibly a TLE-quality/fitting-artifact question (→ #10), not a physics gap at all.

**1.25 — 2026-07-19**
- **CI pipeline added (issue #22, closes it)**: `.github/workflows/ci.yml` + `test_all.sh` run all 11 suites on every push/PR to `main`/`HS-dev`. Installs Intel oneAPI `ifx` on the Linux runner rather than `gfortran` — the first real CI run found `ksrop/propagate_ks.F` doesn't compile under gfortran at all (implicit-real array dimensions, a function/array name collision on `R`; filed as **#28**, P4, since ifx is the only toolchain this project is actually validated against). A second real CI-only failure: `input/const_new.DAT` vs. the source's `'input/const_new.dat'` — silently fine on case-insensitive Windows, hard failure on case-sensitive Linux (the identical bug KSROP already hit and fixed in its own copy). `test_tle_evolution.F`'s T42-T50 (large-catalog checks) now skip gracefully via an `INQUIRE` guard when the 13 MB gitignored catalog file isn't present, rather than hard-failing on a fresh checkout.
- **Issue #11 closed**: verified `test_gmat.exe` already runs with lunisolar enabled (`geo=4, sun=2, moon=2`) and all 14/14 checks pass, including the tight ±10% exact-model magnitude check (N13/N14) — the "pending lunisolar re-check" this issue was left open for is satisfied.
- **TLE quality filtering implemented (issue #10)**: new `tle_filter.F` — trailing-window outlier rejection (apogee altitude vs local mean) and maneuver detection (eccentricity vs local linear trend, `linfit`), with gap-aware window clipping so points right after a real epoch gap aren't judged against stale pre-gap context. Real data tuning was the actual work: the issue's own suggested 3-sigma/5-10-point defaults flagged ~20% of well-tracked validation objects (42928, 35497) as false positives — real orbital eccentricity has more local curvature than a short strict linear fit tolerates. Tuned to `sigma_ha=4, sigma_e=6, nwin=20` → ~1.5% false-positive rate across 3 real objects. 14 new tests (`test_tle_filter.F`). **Not wired into `orem_run`'s default path** — that's a separate decision needing its own 7-object re-validation.
- **Issue #14 reconciled, re-triaged P2→P3**: its 4 scope items are substantially covered by work done since it was filed (J71-vs-GMAT-JacchiaRoberts validation, v1.17; epoch-resolved weather, v1.23/#26; the 10-object RPE campaign) via a different path than originally specified (GMAT cross-check rather than NRLMSISE-00). The literal NRLMSISE-00 comparison and a dedicated density-sensitivity study remain undone but are no longer a P2 accuracy blocker.
- **368/368 tests pass** (354 + 14 new `test_tle_filter.F` checks), verified both locally (`ifx` on Windows) and on CI (`ifx` on Linux).

**1.26 — 2026-07-19**
- **Issue #23 (Production Roadmap) closed, reconciled rather than executed.** A 2026-06-23 planning doc, substantially stale by the time it was revisited: 17/22 referenced issues were already done, but via a different architecture than planned (the ops layer became a separate repo, `OREM-Watchlist`, not embedded scripts) and with major unplanned accuracy work (#25/#26/#27) the original plan never anticipated. The project reached operational status — `OREM-Watchlist` produced a real `IMMINENT` prediction for a real TIP-listed satellite the same day. README's Version History (this document) is the authoritative record going forward, not a forward-looking roadmap.
- **`tle_filter` wired into `orem_run`'s default pipeline (issue #10, closes it).** The deferred decision from #10's original implementation. Re-ran the 7-object campaign to check for regression: result is a genuine, substantial improvement, not just "no regression" — **mean latest-zone \|RPE\| across the 7-object set drops from 17.6% to 8.3%** (more than halved). 27526 goes from a barely-useful 59.4% to a good 3.5%; 4/7 objects improve, 2/7 unaffected, 1/7 (42928) regresses mildly on an already-near-perfect object (0.02%→9.15%). Every executable linking `orem.F` needed `tle_filter.F` added to its build command — CI caught the missing link immediately. `scratch_rpe/rpe_campaign_prefilter_backup.csv` preserves the before state.
- **Cross-repo: `OREM-Watchlist` issue #9 fixed** using this repo's own tooling (`scratch_rpe/zone_windows_37398.F`) — a real TIP-listed candidate (37398) was getting `NO_PREDICTION` not because of a stale `nzones_max`, but because its TLE tracking cadence dropped from ~92/60days (2015-2016) to ~6/60days recently, so `zone_select`'s default `max_zone_days=10` couldn't form any qualifying window in the last 4 years at any R2 threshold — a real-world illustration of a design gap the 8-zone/10-day defaults have for sparsely-tracked objects outside this repo's densely-tracked 7-object validation set.
- **Issue #14, remaining items dispositioned.** Perigee-altitude density-sensitivity study (`scratch_rpe/density_sensitivity_42928.F`): 42928's own real first-tracked state (2017-08-31, hp=166 km) at ATM.DAT scaled 0.5x/1.0x/2.0x shows real, monotonic sensitivity (~30-60 days shift in time-to-a-fixed-decay-state per density-doubling), though none of the three fully re-entered within the tested window. The literal NRLMSISE-00 comparison remains genuinely infeasible with available tooling (NASA CCMC's Instant Run tool needs interactive form submission, not a GET-able endpoint) — recommended formally accepting the existing GMAT JacchiaRoberts comparison (v1.17) as the practical substitute instead.

**1.27 — 2026-07-22**
- **RPE campaign extended 20→30 objects.** The original ROCKET-BODY-only satcat filter returned zero new candidates against the same catalog snapshot; broadened to also include DEBRIS (precedent: object 20, 48259) to find the next 10. Result deepens issue #29 rather than resolving it: only 6/10 new objects form any valid zone at all, and mean \|RPE\| across all 30 roughly doubles vs. the 20-object baseline (20.0%→43.0%).
- **Issue #31 (GA/RSM mean-vs-osculating fitness mismatch) found, fixed, and empirically characterized.** `tle_evolution.F`'s `ha_out` fed `rsm.F`'s `apobs` straight from raw TLE mean elements, while `propagate_ks`'s `surfaces` are genuinely osculating — `ga.F` was RMS-comparing two different physical bases every zone, every object, since the pipeline's inception. Fixed via a new `tle_find_osc` helper (SGP4 `tle2sv` + `car2oe`) that converts each zone-point to a true osculating state, applied *after* `zone_select`'s own R² linearity check (which deliberately stays on the smooth mean-element series — using the noisier osculating series there dropped several already-marginal test zones below threshold). Also fixes `rsm_generate`'s propagation initial condition, previously seeded from a raw mean anomaly mixed with other mean elements.
- **Investigated the fix's real-world impact before shipping it, not after.** A direct git-stash A/B on 37151 (all 8 zones, same window, only the fitness basis changed) shows per-zone fit RMS gets uniformly worse — 2× to 32× — not better. Root cause: TLE mean elements are themselves a least-squares smoothing of several days of tracking; converting each TLE *independently* to osculating recovers a correct-in-isolation short-period correction, but consecutive TLEs are independently-fit snapshots at essentially random, uncorrelated orbital phases. The resulting `apobs` series carries point-to-point noise with a different statistical character than `propagate_ks`'s own smoothly-evolving continuous trajectory — comparing them is *more* correct in kind but noisier in practice. 7-object mean \|RPE\| moves 8.3%→14.2%. Kept anyway, as a deliberate call: the fix corrects a real, previously-undocumented basis error that had been silently present in every fitted BN in this project's history, and 14.2% is treated as the honest current number rather than reverting to an implementation known to be comparing the wrong things. A noise-matched implementation (e.g. a short local propagation across the zone instead of independent per-TLE conversion) is the natural follow-up, not yet done.
- **358/358 tests pass** (full suite re-verified clean after the fix, including previously-latent failures the fix exposed and fixed along the way: `orem_run`'s `ierr=2` early return wasn't initializing `nzones_valid`/`zone_status`, invisible until zone-selection behavior changed enough to actually hit that path in existing tests).
- **30-object campaign regenerated against the fix (same day).** Result is genuinely mixed, not a uniform win or loss: full-30 mean \|RPE\| ticks slightly better (43.0%→39.4%), median ticks slightly worse; the curated original 7 gets worse (matches the direct 37151 finding above) while the less-tuned 20/30-object set is closer to neutral-to-improved (8 objects improved, 7 regressed, 4 flat). Object counts unchanged (19/30 predicting), confirming `zone_select` stability at scale. Pre-fix baseline preserved at `scratch_rpe/rpe_campaign_30obj_prefix31_backup.csv`.
- **Issue #30 (`ga_optimize` hang/crash on a NaN/Inf-contaminated RSM surface) fixed and closed.** A degenerate trial propagation (e.g. from the same saturated terminal-TLE `ndot` pattern found on 33587/issue #27) can leave NaN/Inf in `rsm.F`'s `surfaces` array; `ga_fitness` propagated that into a NaN fitness, which corrupted `ga_preselect`'s stochastic-remainder selection two distinct ways — `int(NaN)` is undefined behavior (reproduced directly as an out-of-bounds array write), and even after flooring the population average away from literal NaN, an all-identical-fitness population still leaves every stochastic remainder exactly 0, spinning the mating-pool loop forever (reproduced directly as the originally-reported indefinite hang on object 42985, same zone, same symptom). Three layered fixes in `ga.F`: `ga_fitness` treats non-finite fit as the worst-possible chromosome instead of passing NaN downstream; `ga_stats` floors the population average away from exactly 0; `ga_preselect` gets a bounded attempt cap with a round-robin fallback, mirroring the same hard-cap pattern `Kesolve` already uses elsewhere in this codebase. 3 new regression tests (G72-G74). **371/371 tests pass.** Verified against a direct reproduction of the original 42985 hang.
- Progress sentinel added to all three `scratch_rpe/*campaign*.F` multi-object loop scripts: `[PROGRESS] NN% (io of nobj objects done)` after each object, letting a log-tailing Monitor report percent-complete without polling.

**1.28 — 2026-07-23**
- **Zone-to-zone IC chaining + always-wide BN search tried, measured, one reverted.** Motivated by the global RPE investigation's (#32) literature finding that ballistic coefficient shouldn't be treated as independent-per-zone. Two changes to `orem_run`'s per-zone loop: (a) propagate a trusted zone's fitted (e, BN) forward via `propagate_ks` to the next zone's first TLE epoch and use that as the next zone's RSM initial-condition seed, replacing a fresh `tle_find_osc` lookup; (b) remove the v1.21 trust-gated BN-range narrow/widen-on-boundary logic, so every zone always searches the full original `[bn_min_init, bn_max_init]` range. **Both regressed RPE**: (a)+(b) combined moved mean \|RPE\| 22.9%→27.4% (7-object) and 22.0%→32.9% (30-object), with two objects (60328, 61734) blowing up over 100 points. Root cause: chaining lets one zone's own fit uncertainty compound into the next zone's starting geometry, instead of every zone independently re-anchoring to fresh SGP4 data. **(a) reverted**, back to unconditional `tle_find_osc` for every zone. **(b) shipped anyway** per explicit user direction to isolate and measure it standalone — still a regression on its own (22.9%→26.2% / 22.0%→23.2%, notably smaller than combined and without the catastrophic per-object blowups) but not yet root-caused as thoroughly as (a); not blocking. `zone_chain.F`/`test_zone_chain.F` (the chaining implementation, verified correct via a 9-test unit suite before the campaign regression was found) were deleted along with (a). Also found and documented (not fixed, out of scope): `propagate_ks`'s own `car2ks`/`ks2car` round-trip at initialization introduces a small, stable, one-time inclination/RAAN/AOP offset from the analytical input elements — harmless for this codebase's magnitude-based apogee fitting, previously undetected since nothing else ever checked orientation fidelity out of a `propagate_ks` call. Full writeup: issue #33. **354 tests passed during development (345 after `zone_chain.F`'s removal), zero regressions, zero new pipeline failures on either campaign.**

**1.29 — 2026-07-24**
- **RPE investigation Phase 2: cheap, correlational error-budget decomposition** (`scratch_rpe/phase2_error_budget.py`, issue #32). Reused the existing 30-object campaign — `rpe_campaign.F` now also logs `zepoch`/`rms_fit` per zone (both already computed by `orem_run`, so this is a zero-new-propagation, additive change; the rerun reproduced every `e_opt`/`bn_opt`/`rpe_pct` bit-for-bit, confirming no regression). Three findings:
  1. **Fitted BN and each TLE's own published BSTAR are strongly anti-correlated within an object** (median r=-0.70 across 20 objects with ≥3 zones; survives a zone-to-zone first-difference check guarding against a shared-trend artifact, median Δr=-0.58). Physically expected — BSTAR≈Cd·A/m while BN≈m/(Cd·A), approximate reciprocals — and confirms the GA's independent per-zone refit is tracking real, independently-derived drag-term variation (SGP4's own BSTAR fit vs. OREM's own apogee-trajectory fit), not just noise. Supports #33's literature premise (BC genuinely isn't constant) while reinforcing that *how* to exploit it matters — trajectory-chaining already failed; using each zone's nearby BSTAR as an informative prior on the BN search range (not yet built) is a more promising next candidate.
  2. **RPE is not a horizon-length artifact.** \|RPE\| vs. prediction horizon shows ~zero correlation (r=-0.06; r=-0.10 vs. 1/horizon) across 101 zone-level predictions, and mean \|RPE\| is statistically indistinguishable for short (<100 d, 67.9%) vs. long (≥100 d, 74.2%) horizons. The large scatter in RPE reflects genuine per-zone fit/extrapolation problems, not a metric artifact — no reason to replace RPE as the project's core accuracy metric.
  3. **Fit RMS predicts boundary-pinned degeneracy but not re-entry-date accuracy.** Boundary-pinned zones (`zstat=2`) have ~9x worse median RMS than trusted zones (0.282 vs. 0.032) — RMS is a good structural-degeneracy signal. But `r(rms_fit, |RPE|) = +0.003` across all 101 predicting zones — essentially zero. A zone can fit its own short observed window very tightly and still be wildly wrong about the extrapolated re-entry date: the window is fundamentally underdetermined for long extrapolation. This directly implies the long-standing "G4 signal-weighted zone selection" idea (#12/#16 history — weight zones by fit quality) would not help RPE, since fit RMS carries no information about extrapolation accuracy.

**1.30 — 2026-07-24**
- **RPE investigation Phase 2 closed out: the last two items (TLE noise, density-model error).** Both cheap — no new propagation.
  1. **TLE noise has a real floor, not fully averaging-reducible.** Extended `scratch_rpe/phase2_error_budget.py` to join each zone's fit RMS against its TLE point count (already printed per-zone, just not previously extracted): median RMS falls 0.070→0.035→0.023 km as points-per-zone go <2000→2000-4000→4000-6000, then **flattens at ~0.02 km beyond ~4000-6000 points** — more tracking data stops helping. Matches #31's own diagnosis that independent per-TLE SGP4 conversion carries a systematic, not purely random, noise component. Consistent with the literature (Frueh & Schildknecht 2012): intrinsic TLE self-consistency noise for HEO (0.8-1.4 km/24h) is ~30-40x smaller than the ~35 km near-epoch error vs. independent optical truth — most real-world TLE error is SGP4/SDP4 model bias, not OD noise, and no amount of extra data removes it.
  2. **Density-model error is already well-characterized and is likely the largest remaining physics gap — but fixing it is a project, not a diagnostic.** Three independent literature sources (Sharma 1997a: `propagate_ks.F` confirmed missing the diurnal density-bulge term, validated to <1% error if added; Swinerd & Boulton 1983: real-tracking-fit scale height runs ~11% high vs. J77 with a 7% fit-residual cost from dropping oblateness, though this is a near-circular validation that may not transfer proportionally to OREM's e=0.3-0.85 regime; ISO/CD 27852: industry-standard statement that atmospheric error dominates propagator fidelity, plus a direct critique of the exponential density form `propagate_ks.F:345` uses) all point the same direction. Posted to #14 rather than implemented here — out of scope for Phase 2's cheap/correlational mandate.
- **Phase 2 is now fully closed (all 5 items).** Next candidates: Phase 3 (#27's Wang & Gurfil resonance test) or the BSTAR-as-prior BN-search idea from v1.29 finding 1.

**1.31 — 2026-07-24**
- **RPE investigation Phase 3: Wang & Gurfil solar-apsidal-resonance test against 33587 (issue #27).** Both previously-tested hypotheses for 33587's unexplained 131-day hp collapse (616→341 km) were already ruled out (GMAT third-body predicts flat hp like OREM's own model; a direct BN refit across the whole window finds no physically plausible constant BN). New `scratch_rpe/phase3_resonance_33587.py` computes `a(t)/e(t)/i(t)` directly from all 306 real TLEs (Kepler's third law on the TLE's own mean motion, no propagation) and evaluates Wang & Gurfil's (2017, *Adv. Space Res.* 59, Eq. 12) resonance criterion — the GTO's combined RAAN+perigee J2 drift rate becoming commensurate with the Sun's apparent motion, at which point the normal ~180-day periodic perigee oscillation turns monotonic. **Result: a genuine crossing on 2022-10-25**, 11 days before the previously-documented multi-point acceleration onset (2022-11-05) and 39 days before the record's last TLE — a real, well-timed candidate mechanism, neither of the two already-ruled-out explanations. **Caveats, stated plainly**: 33587's ~65° inclination sits well outside the paper's own validated low-inclination (6°) numerical examples and near the classical critical inclination, where the criterion's underlying term changes sign (handled via the correct real odd-root branch, `abs(x)**(2/7)`, but still an extrapolation); and the match found is on crossing *timing*, not a magnitude reproduction (would need Wang & Gurfil's SAOD model actually implemented and integrated to check that). An exploratory sweep across all 30 campaign objects finds 8/30 show a crossing at all, several landing within days of the object's own decay (44187: 0-1 days; 60328: 0-10 days; 59347: 28-30 days) while several other near-65° objects show no crossing at all (so it isn't firing on inclination alone) — a real, testable pattern, not yet independently confirmed against a naive-baseline check. **Not proposing to close #27** on this alone; recommend keeping it open with this as a documented, credible lead. Full writeup: #27 and #32 comments.

**1.32 — 2026-07-24**
- **Issue #26 baseline drift resolved; weather regression re-measured against a fresh current-tree baseline.** The "static numbers don't reproduce v1.21" question flagged (not chased) in the 1.28 era wasn't actually mysterious — `#10`/`#12`/array-bound fixes/`#31`/`#33` all landed since and were each validated individually but never re-consolidated into one current number. Rebuilt and reran the static (`rpe_campaign_7obj_issue31`) and weather (`rpe_campaign_weather`) campaigns together fresh against `HEAD`: current static baseline is mean\|latest-zone RPE\|=32.6%, mean\|ensemble RPE\|=26.2%. Weather vs. that baseline: latest-zone (the primary estimator) **improves** 32.6%→30.8%; ensemble mildly regresses 26.2%→28.0%, isolated entirely to one zone (37151 Z7) newly predicting a badly-off re-entry it previously lacked entirely under static — not a systemic density-model problem. `#26`'s "no regression" DoD is arguably met on the primary metric; left open per explicit direction pending a deliberate close decision. Full writeup: #26 comment.

**1.33 — 2026-07-24**
- **Issue #32 (G3): BSTAR-informed BN search-range prior, first accuracy-positive result from this investigation since #10's TLE filter.** Phase 2 (v1.29) found fitted BN anti-correlates with each TLE's own published BSTAR (median r=-0.70) but flagged trajectory-chaining as already-tried-and-failed (#33) — this builds the not-yet-tried alternative: use each zone's own median BSTAR as an informative prior narrowing the GA's BN search range, instead of trying to exploit BC variability through IC continuity. Calibrated (not first-principles) via a pooled log10(BN)=slope·log10(BSTAR)+intercept regression over the 30-object campaign's zone-level data (`scratch_rpe/phase2b_bstar_bn_calibration.py`, n=147, R²=0.49, rmse=0.31 dex); new `estimate_bn_bstar_prior`/`tle_find_bstar` in `orem.F` narrow `[bn_lo,bn_hi]` by intersection only (±2·rmse window, ~4x each side) — never widens beyond the caller's own range, same floor-only spirit as the existing G2 mechanism, falls back cleanly when a zone's BSTAR is unavailable or the prior window doesn't overlap the caller's range at all. Binds for a genuine subset of zones (10/56 in the 7-object set), not a no-op and not universal. **Measured**: mean\|latest-zone RPE\| 32.6%→29.6% (real improvement, driven mainly by 37819 Z1-8 converging on a tighter, more accurate BN), mean\|ensemble RPE\| 26.2%→26.6% (flat). 371/371 tests pass. Full writeup: #32 comment.

**1.34 — 2026-07-25**
- **4-core-capped parallel partitioning added to the 30-object campaign** (`scratch_rpe/rpe_campaign.F`, per `GitHub\CLAUDE.md` §1's repo-wide max-4-cores rule). `orem.F`/`rsm.F`/`propagate_ks.F` use `SAVE`d arrays and common blocks throughout, so this is not safe to thread within one process — instead, optional `<io_start> <io_end>` CLI args partition the object loop across up to 4 independent OS processes, each writing its own `rpe_campaign_part_<start>_<end>.csv`. No-args behavior (every prior invocation) is unchanged.
- **G3's improvement confirmed to generalize to the full 30-object set — directly answers #29's "doesn't generalize beyond the curated 7" concern in G3's favor.** Ran via 4 concurrent processes (objects 1-8/9-16/17-23/24-30), merged: mean\|latest-zone RPE\| 48.94%→43.29%, mean\|ensemble RPE\| 23.22%→22.99%, object counts unchanged (22/30 predicting, no new failures). Driven mainly by two large single-object wins — 40943 and 37819 (matching the 7-object finding) — with most other objects flat and a few small mixed movements; not a uniform win but a real net positive. Pre-G3 baseline preserved at `scratch_rpe/rpe_campaign_30obj_preg3_backup.{csv,log}`. **Correction (v1.35): the originally-posted 39.04%/34.33% figures had a methodology bug — a zone with no predicted re-entry (`reentry_jd=0`) was wrongly counted as a perfect 0%-error prediction whenever it was an object's most recent zone; only the 30-object latest-zone metric was affected (ensemble RPE and the 7-object set were not), same direction, wrong magnitude.** Full writeup: #32 and #29 comments.

**1.35 — 2026-07-25**
- **Issue #33 revisit: v1.21's trust-gated BN-range narrow/widen logic restored alongside G3.** #33 shipped an always-wide BN search only because the alternative (narrow/widen the *global* range toward the previous zone's own fitted BN, gated on that zone actually carrying real drag signal) measurably regressed RPE when tested *before* G3 existed. G3 changes the calculus: it narrows each zone's range from that zone's own BSTAR, independent of zone-to-zone history, so the two mechanisms compose rather than conflict — restored logic adjusts the persistent `bn_lo`/`bn_hi` between zones, G3 then intersects its own per-zone prior against whatever that currently is.
- **Real, substantial further improvement on top of G3 alone**, corrected methodology throughout:

  | Metric | Pre-G3 | G3-only | G3+trust-gate |
  |---|---|---|---|
  | 7-obj mean\|latest-zone RPE\| | 32.6% | 29.6% | **14.0%** |
  | 7-obj mean\|ensemble RPE\| | 26.2% | 26.6% | **22.9%** |
  | 30-obj mean\|latest-zone RPE\| | 48.9% | 43.3% | **39.4%** |
  | 30-obj mean\|ensemble RPE\| | 23.2% | 23.0% | **22.1%** |

  Every metric on both sets moves the right direction, none regress. 371/371 tests pass. Full writeup: #33 comment.

**1.36 — 2026-07-25**
- **Issue #31 revisited: the noise-matched-apobs follow-up, retested alongside G3+trust-gate.** The earlier attempt (smooth linfit-trend correction to the mean/osculating bias, replacing the fully independent per-point SGP4 conversion) measured flat against the pre-G3 baseline and was reverted. Retested unchanged against the new G3+trust-gate baseline, on the theory that a narrower, more history-informed BN search might respond differently to the observation-side noise character. Result is genuinely mixed, not a repeat of either extreme:

  | Metric | G3+trust-gate | +noise-matched apobs |
  |---|---|---|
  | 7-obj latest-zone | 14.0% | **11.6%** |
  | 7-obj ensemble | 22.9% | 23.4% |
  | 30-obj latest-zone | 39.4% | **38.7%** |
  | 30-obj ensemble | 22.1% | 22.2% |

  Small, consistent improvement on the primary (latest-zone) metric on both sets; small, consistent regression on the secondary (ensemble) metric on both sets. Nowhere near #33's magnitude either direction. **Kept per explicit user direction** — net positive on the metric this project treats as primary. 371/371 tests pass. Full writeup: #31 comment.

**1.37 — 2026-07-25**
- **Issue #26 revisited against the current shipped state (G3 + trust-gate + noise-matched apobs).** The earlier weather-mode regression (37151 Z7 newly predicting a badly-off re-entry, dragging that object's ensemble RPE from 4.3%→14.9%) no longer appears at meaningful magnitude — the BN-search changes since then evidently changed which zones converge where. Fresh static-vs-weather comparison: mean\|latest-zone RPE\| 11.59%→11.08%, mean\|ensemble RPE\| 23.39%→23.33%, both flat-to-slightly-better, no per-object outlier this time (largest ensemble delta is 35497, the already-known unrelated resonance case, at +6.2 points). **"No regression" DoD now cleanly met on both metrics**, not just the primary one as before. Full writeup: #26 comment.

**1.38 — 2026-07-25**
- **Issue #32 Phase 4: re-ran Phase 2's error-budget decomposition (`scratch_rpe/phase2_error_budget.py`, unchanged) against the current shipped state (G3 + trust-gate + noise-matched apobs).** Pure re-analysis of the already-current 30-object campaign CSV, no code change. Two findings shifted, two held up:
  - **Changed: BSTAR↔BN correlation weakened substantially** (median r=-0.70 pre-G3 → **-0.43** now) — expected mechanically, since G3 explicitly pulls fitted BN toward the BSTAR-predicted value where it binds, so the *residual* correlation left over is naturally smaller. Not a red flag; flags that G3 has partially "used up" the very signal it was built to exploit.
  - **Changed: RPE now shows a real horizon-dependence Phase 2 didn't find.** Mean\|RPE\| for horizon<100d is 33.2% vs 69.6% for horizon≥100d, a large gap — Phase 2 (pre-G3) found these "statistically indistinguishable" (67.9% vs 74.2%). Raw linear correlations stay weak (r=+0.01, r=-0.14 vs 1/horizon), so this reads as a threshold/bucket effect rather than a clean linear one, but suggests the improved fitting disproportionately helps short-horizon (later-zone, more-TLE-history) predictions while long-horizon extrapolation stays fundamentally hard — reinforces the latest-zone-as-primary design logic rather than undermining it.
  - **Unchanged: fit RMS still doesn't predict extrapolation accuracy** (r(rms_fit,\|RPE\|)=+0.01, boundary-pinned zones still ~9x worse median RMS) and **TLE noise floor still flattens around 4000-6000 points/zone** — both Phase 2 conclusions hold up against the new pipeline.
  - Full output: `scratch_rpe/phase4_error_budget_current.log`. Full writeup: #32 comment.
- **Correction to finding 1 above**: the "-0.70→-0.43 weakened" claim compared two different statistics — Phase 2's per-object within-zone correlation vs. G3's actual calibration basis, a pooled cross-object regression, never recomputed against current data before writing the original comment. Redone correctly (`scratch_rpe/phase2b_bstar_bn_calibration.py`, unchanged, rerun against current data): **R² 0.491→0.509, rmse 0.3096→0.2689 dex — the relationship held up, if anything slightly tighter, not weaker.** Also verified directly: splitting the current data by whether G3 actually bound for that zone, both subsets remain strongly correlated (bound r=-0.64 n=42, unbound r=-0.70 n=105) — no sign of degradation.
- **G3's window tightened accordingly** (`RMSE_DEX` 0.3096→0.2689 in `estimate_bn_bstar_prior`, center/`SLOPE`/`INTCPT` unchanged — single-variable test of the narrower window alone). Measured, another clean improvement:

  | Metric | Wide window | Tight window |
  |---|---|---|
  | 7-obj latest-zone | 11.59% | **9.89%** |
  | 7-obj ensemble | 23.39% | 23.53% (flat) |
  | 30-obj latest-zone | 38.71% | **35.49%** |
  | 30-obj ensemble | 22.20% | **21.53%** |

  3 of 4 metrics improve, the 4th is flat (noise-level). 371/371 tests pass. Shipped.

**1.39 — 2026-07-25**
- **Issue #32 Phase 5: first GMAT independent cross-check of today's fit.** Extracted 42928's zone-8 (latest zone, primary estimator) exact SGP4-osculating IC via a temporary debug dump (reverted, not shipped) and built the perigee-preserving GA-fitted Keplerian state (SMA=7525.205157 km, ECC=0.133748, TA resolved from the zone's MA via Kepler's equation using the fitted eccentricity — exactly what `oe2car` does inside `rsm_generate`), then independently propagated it through GMAT's own JacchiaRoberts drag model (`scratch_gmat/gmat_xval_42928z8_2026.script`), same force model as the campaign (EGM2008 zonal deg=20, Sun+Luna point mass, SRP on, static F107=72/Kp=1 matching the static campaign this fit came from).
  - **Result: a real, substantial disagreement.** OREM predicts re-entry 86.4 days after zone start (RPE -4.0% vs the observed 90.0-day horizon — an excellent fit by OREM's own metric). GMAT's independent propagation of the *identical* fitted state reaches the same 102 km threshold in only 48.3 days — 38 days (42%) sooner.
  - **Correct interpretation, not "OREM is wrong": BN is a `propagate_ks`-specific calibrated parameter, not a transferable physical ballistic coefficient.** OREM's GA fits BN specifically so `propagate_ks`'s own density profile reproduces the observed TLE decay — the same numeric BN plugged into a *different* atmosphere model (JacchiaRoberts) doesn't preserve equivalent drag acceleration if the two models' density-vs-altitude profiles differ. This cross-check quantifies, for the first time with a real number, the density-model gap already flagged in Phase 0/2 literature (Swinerd & Boulton 1983: real-tracking-fit scale height ~11% high vs J77; ISO/CD 27852's direct critique of `propagate_ks.F`'s exponential density form) — `propagate_ks`'s implied density is apparently understated enough that JacchiaRoberts decays the same nominal state 42% faster over this horizon.
  - Not a regression or a bug to fix here — a credibility/context finding for the whole RPE campaign: OREM's internal self-consistency (fitting BN to match TLEs, then extrapolating with the same model) is what the RPE metric measures, and this cross-check is evidence that self-consistency is doing real work masking a real, sizeable density-model discrepancy from any *external* physical ballistic-coefficient interpretation. Full writeup: #32 comment.
- **Extended to a second object (39615, Proton-M Briz-M) — the gap is consistent, not object-specific.** Same methodology (`scratch_gmat/gmat_xval_39615z8_2026.script`), a very different case: inclination 48.4° vs 42928's 19.1°, eccentricity 0.28 vs 0.13, BN 53.1 vs 93.5, six-month-different epoch. OREM predicts 150.2 days to re-entry (RPE +9.4% vs the 137.3-day observed horizon); GMAT's independent propagation of the identical fitted state reaches 102 km in 80.3 days — **46.5% sooner**, remarkably close to 42928's 44.1% gap despite the objects sharing almost nothing else. Strong evidence this is a systematic `propagate_ks`-vs-JacchiaRoberts density-model offset, not an object- or altitude-specific artifact (n=2, both perigee ~130-140 km at fit time — a third check at a different perigee regime would further test that boundary). Full writeup: #32 comment.
- **Extended to the full curated 7-object set — a more nuanced picture than n=2 suggested, one real reversal found.** Same methodology on the remaining 5 (35497, 37151, 27526, 32007, 37819):

  | norad | perigee alt (km) | inc (deg) | OREM (days) | GMAT (days) | GMAT vs OREM |
  |---|---|---|---|---|---|
  | 27526 | 118.5 | 17.6 | 11.6 | 3.7 | -67.9% |
  | 37151 | 130.5 | 24.8 | 137.5 | 47.6 | -65.4% |
  | 39615 | 130.7 | 48.4 | 150.2 | 80.3 | -46.5% |
  | 42928 | 140.6 | 19.1 | 86.4 | 48.3 | -44.1% |
  | 35497 | 109.0 | 5.7 | 34.2 | 20.6 | -39.6% |
  | 32007 | 169.4 | 25.8 | 133.2 | 85.0 | -36.2% |
  | 37819 | 139.7 | **63.2** | 96.3 | 101.3 | **+5.2%** |

  6 of 7 objects still show GMAT decaying faster (direction consistent), but the *magnitude* ranges widely (-36% to -68%), not the tight ~45% band the first two suggested — the earlier framing was accurate as far as it went but didn't yet have enough data to see the real spread. **37819 reverses the direction entirely** — notably the only object with inclination above 50° (63.2° vs everyone else's <50°), a plausible but unconfirmed candidate explanation (n=1 for that regime, not established). No clean single-variable trend against perigee altitude either (169km gives the smallest gap among the 6 consistent cases, but 109km gives a smaller gap than 118-130km, not a monotonic relationship). **Honest summary: real, substantial, mostly-one-directional density-model disagreement between `propagate_ks` and JacchiaRoberts, but not a single fixed correction factor — inclination (and possibly other factors) modulates it, not yet characterized precisely.** Full writeup: #32 comment.
- **Issue #32 Phase 5, second candidate: ISO/CD 27852 standard comparison — read the primary source in full (19 pages, `E:\Research\References\40 March 2019\Space systems — Determining orbit lifetime, 2007.pdf`), not just the earlier lighter literature-survey pass.** Real correction found, worth stating plainly:
  - **Scope mismatch, stated explicitly**: this standard targets 25-year *post-mission LEO-crossing disposal* lifetime for collision-risk compliance (IADC's 25-year guideline) — a generic, long-duration, forward-forecast problem. OREM predicts a specific near-term reentry *date* for an already-decaying real HEO/GTO object using already-observed TLE history. Different problem class; not everything in the standard transfers directly.
  - **Correction to the earlier framing**: the previous literature note said the standard "critiques `propagate_ks.F`'s exponential density form." Re-reading closely, the standard's own atmosphere-model comparison (§5.1, Fig. 7-9) treats "Exponential" as one specific, crude, standalone named model (runtime 1.4s/1.8M evals, visibly diverges from every other model at low/high altitude in their own Fig. 9) — a different thing from `propagate_ks`'s King-Hele/Sharma-lineage *locally*-exponential drag form driven by a real empirical reference-density/scale-height table (J71-based `ATM2D`, now with #26's epoch-resolved F10.7/Ap variation). The standard's own words on the Jacchia-1971-class tier `propagate_ks`'s lineage actually belongs to: "not ideal, [but] can work well for long-duration orbit lifetime studies." **OREM's atmosphere model sits in the standard's own acceptable-if-not-best-fidelity tier, not its explicitly-named "avoid" category.**
  - **Where OREM already complies**: the standard's core atmospheric-fidelity concern is "models that do not accommodate solar activity variations" — #26's epoch-resolved weather (real historical F10.7/Ap, not a static table) already addresses this directly.
  - **Genuine remaining gaps**: OREM doesn't offer MSIS2000/Jacchia-Bowman (the standard's explicitly-preferred best-fidelity models) as an option; OREM's gravity model (degree-20 zonal) is simpler than the standard's Method-3 ceiling (up to 300×300 spherical harmonics) — though zonal-only is a common, defensible simplification for these eccentric-orbit regimes, not obviously wrong for OREM's own purpose.
  - **Ballistic-coefficient treatment**: the standard's own fast method (Method 2, ~45s/30yr case) uses one constant average BC; only its slowest, highest-fidelity method (Method 3, ~1700s/30yr case) supports a detailed attitude/angle-of-attack-dependent BC. OREM's per-zone independently-refit BN sits between these two tiers — more adaptive than a single global constant, short of full attitude modeling — and qualitatively aligns with the broader literature consensus (Walsh 2012, Russell et al. 2012) that constant BC underperforms.
  - **Not applicable to OREM at all**: the standard's §5.2 (bulk of its remaining technical content) is about *forecasting* solar/geomagnetic activity uncertainty over an unknown 20-30-year future window via random historical-triad draws — OREM's short-horizon predictions use already-observed real history, never a multi-decade forecast, so this entire section doesn't transfer.
  - **Net assessment**: no new actionable gap found beyond what #26/#14 already track — the main value of this exercise was correcting an overstated prior claim (OREM's atmosphere model isn't in the standard's "avoid" category) and confirming #26's weather work already satisfies the standard's central concern. This closes out Phase 5's second candidate. Full writeup: #32 comment.

**1.40 — 2026-07-25**
- **RPE campaign extended 30 → 50 objects (issue #29).** 20 new objects filtered from a fresh `satcat.csv` pull (APOGEE>8000, PERIGEE<3000, ROCKET BODY/DEBRIS, decay after 2015, excluding all 30 already used) — 41 raw candidates, screened down using the exact same criteria every prior expansion round established (already-known saturated-ndot/too-few-TLE/multi-decade-gap exclusions, plus 3 new anomalous-apogee CZ-5-class exclusions likely to be trans-lunar-injection stages, not GTO). A large SL-12 debris cluster shared near-identical 2025-05-07/08 decay dates (likely one breakup/catalog-processing batch, not independent events) — kept only a representative subset rather than the whole cluster. **Reproducibility verified**: the original 30 objects reproduce bit-identical `bn_opt` across all 148 predicting zones — the expansion adds data without perturbing anything existing.
  - **Result is genuinely mixed, not a clean win or loss**:

    | Set | mean\|latest-zone RPE\| | mean\|ensemble RPE\| | predict rate |
    |---|---|---|---|
    | 30-object (unchanged) | 35.5% | 21.5% | 63% (19/30) |
    | 20 new objects only | 47.4% (worse) | 15.6% (**better**) | 45% (9/20, lower) |
    | 50-object combined | 39.3% | 19.0% | 56% (28/50) |

  - The new 20 predict *less often* (45% vs 63%) — consistent with #29's core generalization concern — but the ones that *do* predict do notably better on the ensemble metric than the curated 30, while doing worse on latest-zone. Not a uniform "broader sets are worse" story on every axis. Pre-expansion 30-object baseline preserved at `scratch_rpe/rpe_campaign_30obj_pre50_backup.{csv,log}`. Full writeup: #29 and #32 comments.

**1.41 — 2026-07-25**
- **`zone_select` Pass 2 gains a recency-guaranteed slot (issues #29/#32).** Root cause: Pass 2 selected the top-`nzones_max` candidate windows purely by R², with no recency preference — for a long TLE history, numerous older high-R² windows could crowd out the physically critical terminal-decay zone entirely before the cap was reached (documented for object 11550; previously worked around only by raising `mxz` to 50, never fixed at the selection-criterion level). Fix: since Pass 1 already stores candidates in chronological order (each acceptance jumps the scan past the whole window before the next start is tried), the single most recent candidate is now always guaranteed one of the `nzones_max` slots before the R²-greedy pass fills the rest. No new tunable parameter — same "floor/prior, never touches the well-behaved case" style as the existing G2 BN floor and G3 BSTAR prior. 3 new unit tests (Z66-Z68), synthetic case reproducing the exact crowding pattern; **374/374 tests pass**.
  - **Confirmed on the motivating case**: at `nzones_max=8` (this campaign's own setting), 11550's terminal zone (ending 32 days before its last tracked TLE, perigee 96 km) was previously excluded entirely, crowded out by eight 2006-2009 windows — now included, displacing one redundant near-duplicate old window.
  - **50-object campaign regenerated** (pre-fix baseline preserved at `scratch_rpe/rpe_campaign_prezoneselectfix_backup.csv`): net improvement on every aggregate metric tracked — curated-7 mean\|ens RPE\| 23.53%→21.32% (median 11.47%→9.69%, max 106.41%→86.99%, driven mainly by 35497's 106.4%→87.0%); all-50 mean 19.02%→18.58% (median 9.28%→8.13%).
  - **Not a uniform win**: 17/50 objects changed at all. 42928 regresses +7.5pp (12.0%→19.6%) — root-caused to a real second-order effect, not the terminal zone itself: displacing any zone from the sequence to make room perturbs the *trust-gated BN-range carryover* (which chains sequentially from each zone's fit into the next zone's search bounds, `orem.F:646-672`), shifting the range fed into whichever zones immediately follow the displaced slot. Two of 42928's zones show notably worse fits purely from this carryover disruption, on otherwise-identical TLE windows; the chain re-converges a couple of zones later. Accepted despite this cost, since the aggregate effect across every tracked metric is positive and no object catastrophically fails. Full writeup: #29 and #32 comments.

**1.42 — 2026-07-25**
- **Diurnal atmospheric density bulge added to `propagate_ks`, auto-enabled with epoch-resolved weather mode (issue #32 Phase 5 follow-up).** `propagate_ks` previously had no diurnal bulge at all — density depended only on altitude (and, since #26, exospheric temperature), never on where the satellite sits relative to the Sun, a limitation `ALGORITHM.md` §9 already flagged as the likely largest remaining physics gap after Phase 5's GMAT cross-check found `propagate_ks` reads 36-68% "thin" vs GMAT JacchiaRoberts at low-perigee terminal altitudes on 6/7 curated objects (JacchiaRoberts has a bulge; OREM didn't). Implemented Sharma 1997a's exact missing term (`rho = rho_0*(1+F*cos(phi))*exp{...}`, read in full this session) with every constant literature-cited, none guessed: `lambda=37deg` bulge lag and the `T_max`/`T_min` extrema formula (`0.15`/`0.24`/`K_max=0.5458`/`K_min=-0.5722`) both come from Swinerd & Boulton 1983's real-tracking fit against 52 orbits of satellite 1963-27A (also read in full), whose own fit found the amplitude within 2% of this exact formula. Two new subroutines in `propagate_ks.F` — `jacchia_bulge_extrema` (temperature extrema from `T_mean`+Sun declination) and `bulge_cosphi` (geometric factor from the satellite's and Sun's position, both already computed for existing physics, zero new ephemeris calls) — reuse the existing `atm2d_interp` table to map `T_max`/`T_min` to a density-ratio amplitude, so the feature only activates when weather mode (#26) is already loaded; when it isn't, `f_bulge=0` is an exact no-op. No signature change to `propagate_ks` and no call-site updates anywhere (avoids the class of interface-change bug that bit `main_orem.F` at v1.10). 6 new unit tests (`test_sw.F` W13-W18: exact closed-form extrema at the equinox, exact geometric cases for `bulge_cosphi`); **380/380 tests pass**, zero regressions with weather mode off.
  - **Isolated-physics validation** (42928 Z8, same osculating IC, bypassing the BN-carryover chain via `scratch_gmat/bulge_xval_42928z8.F`): bulge makes the same zone's re-entry propagation decay ~4% faster (28.49→27.36 days to 102 km) — the *correct direction* (GMAT independently predicts faster decay than OREM), so the physics isn't backwards. Both weather-mode figures are already far below the original static-flux comparison (86.4d OREM / 48.3d GMAT) for reasons predating bulge entirely (real historical F10.7/Kp vs constant-flux, issue #26) — this test can't attribute closing "the" 36-68% gap to bulge alone, only confirm it nudges the right way, by a small amount, at this one zone.
  - **Full curated-7 weather campaign regenerated** (pre-bulge baseline: `scratch_rpe/rpe_campaign_weather_prebulge_backup.csv`): **net regression** — mean|ensemble RPE| 20.00%→23.03%, median 10.70%→14.49%, 6/7 objects worse (35497 88.7%→98.2% the largest; only 37819 improves, 10.1%→9.3%). Root cause: the same trust-gated BN-range carryover sequencing from v1.41 — a small BN shift from adding real physics in one zone changes the search range inherited by the next, amplifying non-uniformly rather than reflecting the bulge term's own (small, correctly-directed) effect.
  - **Accepted and kept auto-enabled per explicit user decision**, same reasoning already applied to v1.41's 42928 regression: the underlying mechanism is literature-grounded and directionally validated in isolation; the campaign-level cost is a characterized, understood interaction with existing BN-carryover logic, not evidence the physics itself is wrong. Full writeup: #32 comment.
  - **Ported to KSROP** (`KSROP/app/driver_KS.F` HS-dev `0e0110e`) same session — identical drag-model lineage, same `DENS` formula, already had `sw_tinf`/`atm2d_interp`/`solarnpv` from the #26 port. 6 new tests there too (`test_sw.F` W12-W17); 425/425 KSROP tests pass.

**1.43 — 2026-07-25**
- **`max_zone_days` shipped default raised 10→20 days (issue #29).** Root-caused the 50-object campaign's 12 zero-valid-zone objects into two distinct failure modes (`scratch_rpe/diag_zero_zone.F`/`diag_zero_zone2.F`, applying `tle_filter` first to match the real pipeline exactly): **7/12 were a pure density problem** — a real, clean decay signal exists, but never inside any 10-day window (rescued at 15-30 days: 41553, 44187, 44591, 44892, 46429, 48259, 28572); **5/12 have no clean decay signal in their tracked history at any window width up to 60 days** (23647, 41693, 43670, 44563, 57804) — a genuine data-availability gap (matches the already-documented "last-TLE-to-decay gap" pattern), not fixable by this parameter at all.
  - **Tested `max_zone_days` ∈ {20, 30} against the curated-7 campaign first**: both roughly halve mean\|ens RPE\| (23.03%→13.75%/13.83%) and cut the max from 98.15% to 42.27%/30.08% — driven mainly by 35497 (this project's long-standing worst curated-7 case) finally fitting well. But 30 days introduces two severe individual regressions 20 days doesn't have (39615: +1.4%→-27.9%, 37151: 14.5%→30.1%) — **20 days chosen** as the safer, more balanced value; no object regresses by more than ~1.4pp at 20 days.
  - **Confirmed on the actual generalization target — the full 50-object campaign** (pre-change baseline preserved at `scratch_rpe/rpe_campaign_maxzonedays10_backup.csv` / `..._weather_maxzonedays10_backup.csv`; `max_zone_days=30` alternative preserved at `scratch_rpe/rpe_campaign_weather_mzd30.csv` for reference): genuine predict rate **58%→72%** (29→36 of 50 objects, using `reentry_jd>0` as the real gate — not the `ens_rpe_pct=0` sentinel, a methodology bug this project has been bitten by before). Curated-7 median 9.69%→7.30%, mean 21.32%→12.03%, max 86.99%→39.33% — a clean, substantial win. Among all 28 objects predicting under both settings, median improves 14.82%→8.82%; mean is flat (22.24%→22.28%) only because of one outlier (61734, +156pp — see below).
  - **Not a clean win**: two real costs, both already-familiar failure classes this session. (1) **61734** (already predicting before, not one of the rescued objects) regresses 25.4%→181.5% — the same zone-composition-perturbs-the-BN-carryover-chain sensitivity as v1.41's 42928 regression and v1.42's campaign regression, a third instance of the identical mechanism. (2) **44187** (newly rescued) produces a technically-passing but wildly untrustworthy fit (565% RPE) from a marginal R²≈0.90 window — widening doesn't just add coverage, it can also admit noisy near-threshold windows that pass the gate but shouldn't be trusted. Neither is fixed by this change; both are documented as known costs, same acceptance reasoning as v1.41/v1.42 (net positive across every tracked aggregate metric, individual costs characterized rather than hidden).
  - Shipped: all 7 `input/orem_*.cfg` files and the 3 tracked campaign harnesses (`rpe_campaign.F`, `rpe_campaign_weather.F`, `rpe_campaign_7obj_issue31.F`) updated to `max_zone_days=20.0`. No source-code change — `zone_select.F`'s own logic is unchanged and already handles arbitrary `max_zone_days` correctly (parametrized tests Z28/Z29/Z31/Z32 already cover this), so this is purely a validated configuration-value change. 380/380 tests pass (unaffected, as expected).
  - **Trust-gated BN-carryover chain has now caused three separate, unrelated changes this session to regress on at least one object** (v1.41's zone_select fix, v1.42's diurnal bulge, this change's 61734) while each change was independently correct/beneficial in aggregate. If a fourth mechanism shows the same pattern, the carryover chain itself — not whatever new mechanism triggers it — is probably what needs fixing next.

**1.44 — 2026-07-26**
- **`orem.F`'s BN-carryover chain: widen-on-boundary now anchors to the existing search-range center instead of the boundary-pinned zone's own value.** Investigated the chain flagged three times in v1.41-v1.43 by tracing 61734's exact BN escalation under `max_zone_days=20` (BN 16.3(boundary)→58.6→89.4(boundary)→136.0(boundary), RPE 25.4%→181.5%). Root cause: a boundary-pinned zone's own fitted BN is definitionally unreliable — it sits at whatever edge the *current* search window had, not necessarily the true value — but the old widen-on-boundary logic still re-centered the *next* zone's window on that unreliable value before widening, letting a run of boundary hits walk/escalate the window in one direction each time instead of ever escaping it. Fix: widen around the window's existing center, not the boundary-pinned zone's own value — preserves the original "give the next zone more room" intent without letting the center itself drift on an untrustworthy point. 380/380 tests pass (unaffected — no existing test asserts exact BN values sensitive to this).
  - **Mechanism confirmed, but doesn't fix the case that motivated it.** Re-tracing 61734 under the fix: zone 3 stops being boundary-pinned and the BN progression genuinely differs (58.6→75.0→100.98→124.0 vs the old escalating 58.6→89.4→136.0) — the fix demonstrably changes chain propagation as designed. But 61734's own ensemble RPE is unchanged (181.5%→181.7%), because **zone 1 itself — the actual defect — is the first zone in the chain, with nothing upstream for a carryover fix to correct.** Zone 1's own GA search independently lands on a boundary-pinned, wildly wrong BN=16.3 (RPE 759.7%), a different problem than carryover propagation.
  - **Aggregate effect across the 50-object campaign: a wash**, not a win. Curated-7 mean\|ens RPE\| 12.03%→12.01% (unchanged), all-50 mean 37.54%→37.52% (unchanged), max unchanged at 565.16% (44187, a separate marginal-fit problem, untouched by this). Only 6 objects move by more than 1pp, split roughly evenly: 42928 improves (8.94%→5.08%, the original v1.41 regression finally reversing) but 32007 gets notably worse (7.30%→11.55%).
  - **Kept per explicit user decision** despite the wash: the fix is real, mechanistically correct for the specific pathology it targets (boundary-value escalation), harmless to every aggregate metric, and does help at least one previously-regressed object. Not shipped as "the" fix for the carryover chain's sensitivity — that remains open. 61734's own defect (why zone 1's fit is so bad under `max_zone_days=20`) is a distinct, uninvestigated question.
  - Pre-fix baselines preserved at `scratch_rpe/rpe_campaign_precarryoverfix_backup.csv` / `..._weather_precarryoverfix_backup.csv`.

- **v1.44 — 61734's zone-1 defect root-caused (not a G3 bug); ensemble secondary metric switched from arithmetic mean to median (`orem.F` `compute_rpe`).** Investigated the case flagged unresolved by v1.41–v1.43 (`scratch_rpe/diag_61734_zone1.F`/`diag_61734_zone1_data.F`): zone 1 sits at ~35,000 km apogee (e=0.72, near-original GTO), a genuinely clean 19-day linear fit (R²=0.994) only 138 days before the object's true re-entry. **Disproved the G3-miscalibration hypothesis directly**: a test build with G3's BSTAR-floor bypassed still lands zone 1 at a non-boundary BN=26.3 with RPE unchanged at 763.9% — the search range was never the problem. Root cause: a single constant-BN KS propagation from this early, high-apogee window cannot linearly extrapolate the ~1000+ days such a fit implies, because the real object's decay accelerates far faster than that over its remaining ~138 days. This is a long-horizon-extrapolation limitation intrinsic to zone 1's own position in the decay history, not a fixable search/calibration bug.
  - **Fix targets the ensemble's vulnerability to this failure mode, not the zone-1 fit itself**: `compute_rpe`'s secondary "ensemble" indicator now reports the **median**, not the arithmetic mean, of valid per-zone predictions (`t_mean` parameter name kept for call-site compatibility; `report.F` labels updated to "median re-entry JD"). One catastrophically-wrong early zone can drag a mean far off; a median is far more robust to exactly one such outlier. Re-verified the underlying premise first (`scratch_rpe/ensemble_eval_50obj.py`, re-running the v1.20 `ensemble_eval.py` comparison against the *current* 50-object campaign rather than assuming the old figures still hold): median beats uniform mean on every aggregate stat today, though both remain behind the already-shipped latest-zone **primary** estimate (v1.20) — this change only affects the secondary spread/agreement indicator.
  - **Result: one of the cleanest wins this investigation has produced.** Curated-7 static mean\|ensemble RPE\| 12.01%→**4.05%**, median 7.27%→**2.78%**, max 38.81%→**12.23%**. Curated-7 weather+bulge mean 14.27%→**4.50%**, median 11.22%→**3.40%**, max 41.82%→**10.00%** (6/7 objects improve; only 37819 regresses, 4.48%→9.40%, still small). All-50 (45 comparable objects) mean 30.02%→**24.26%**, median 7.27%→**3.70%**; max unchanged at 565.16% (44187 has only one valid zone — median of one value is that value, so this single-zone marginal-fit case is untouched by construction, exactly as expected). 12 objects improve by >0.5pp, 4 regress by >0.5pp (all small), 29 unchanged. **61734 itself: 181.68%→-6.77%**, matching this investigation's own hand-derived prediction exactly. Predict rate unaffected (36/45, as expected — a pure reporting-statistic change touches no fitting/propagation code). 382/382 tests pass (2 new: D6 encodes 61734's exact scenario as a synthetic regression test).
  - **Does not touch the fragile BN-carryover chain at all** — no RSM/GA/zone_select changes, so none of the three regression patterns that have hit every other change this session apply here.
  - Pre-fix baselines preserved at `scratch_rpe/rpe_campaign_premedian_backup.csv` / `..._weather_premedian_backup.csv`.

---

## 9. References

- Sellamuthu, H. (2019) Regularized Astrodynamics Using Kustaanheimo-Stiefel Space, Ph.D. Thesis, Karunya Institute of Technology and Sciences
- Sellamuthu, H., Sharma, R.K. & Arumugam, S. Optimal re-entry time prediction of RSO from HEO, Advances in Space Research (submitted)
- Stiefel, E.L. & Scheifele, G. (1971) Linear and Regular Celestial Mechanics, Springer-Verlag
