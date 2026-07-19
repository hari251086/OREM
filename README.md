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
├── ksrop/                          Propagator engine (from KSROP repo)
│   ├── propagate_ks.F              KS propagator as callable subroutine
│   ├── Subrouts.F                  Coordinate transforms, I/O, utilities
│   ├── Legendre.F                  Zonal Legendre polynomial evaluation
│   └── TLEread.F                   TLE reader + SGP4/SDP4 conversion
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
├── tle_evolution.F                 Batch TLE → orbital evolution
├── zone_select.F                   Zone selection — linear apogee decay, top-R² candidates
├── ga.F                            Binary-coded GA optimizer (pop=20; trajectory-matching fitness)
├── rsm.F                           RSM surface generation (9 surfaces per zone)
├── orem.F                          OREM driver: pipeline + G2 BN floor + zone diagnostics
│                                   + trust-gated BN carryover + compute_rpe
├── report.F                        Prediction report writer (latest-zone primary + ensemble)
├── main_orem.F                     Standalone runner (reads orem.cfg, writes the report)
│
├── test_propagate_ks.F             Propagator tests (10)
├── test_tle_evolution.F            TLE evolution tests (56)
├── test_zone_select.F              Zone selection tests (68)
├── test_ga.F                       GA optimizer tests (71)
├── test_ga_sensitivity.F           GA parameter sensitivity study (not in test suite)
├── test_rsm.F                      RSM integration tests (39)
├── test_orem.F                     Driver + diagnostics + G2 + report tests (29)
├── test_reentry.F                  7-object re-entry validation (35)
├── test_e2e.F                      End-to-end integration, IDRAG=1 + full force (20)
├── test_gmat.F                     GMAT cross-validation + exact-model drag reference (14)
└── README.md                       (342 tests total)
```

---

## 3. Propagator Interface

The core propagator is `propagate_ks` — a callable subroutine refactored from KSROP's `driver_KS.F`:

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

Requires **Intel oneAPI Fortran** (`ifx`) or **GNU Fortran** (`gfortran`).

### Windows (Intel oneAPI ifx 2025.0)

`/heap-arrays /F:16777216` (16 MB stack) is required for every executable that links `rsm.F`/`ga.F` — the `surfaces(5000,3,3)` arrays overflow the default stack without it. `test_orem.exe` and `orem.exe` additionally link `report.F` (since v1.19).

```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx /heap-arrays /F:16777216 test_propagate_ks.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F /exe:test_propagate_ks.exe
ifx /heap-arrays /F:16777216 test_tle_evolution.F tle_evolution.F ksrop/TLEread.F ksrop/Subrouts.F ksrop/Legendre.F /exe:test_tle_evolution.exe
ifx /heap-arrays /F:16777216 test_zone_select.F zone_select.F tle_evolution.F ksrop/TLEread.F ksrop/Subrouts.F ksrop/Legendre.F /exe:test_zone_select.exe
ifx /heap-arrays /F:16777216 test_ga.F ga.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F /exe:test_ga.exe
ifx /heap-arrays /F:16777216 test_rsm.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F ga.F /exe:test_rsm.exe
ifx /heap-arrays /F:16777216 test_orem.F orem.F report.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F /exe:test_orem.exe
ifx /heap-arrays /F:16777216 test_reentry.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F /exe:test_reentry.exe
ifx /heap-arrays /F:16777216 test_e2e.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F /exe:test_e2e.exe
ifx /heap-arrays /F:16777216 test_gmat.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F ga.F /exe:test_gmat.exe

ifx /heap-arrays /F:16777216 test_sw.F swx.F orem.F report.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_sw.exe

REM Standalone runner -- swx.F required since v1.23 (main_orem.F calls sw_load/atm2d_load)
ifx /heap-arrays /F:16777216 main_orem.F orem.F report.F swx.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:orem.exe
```

### Unix / gfortran

Same source lists as above with `gfortran ... -o <name>.exe` (no `/heap-arrays` equivalent needed if the default stack suffices; otherwise `ulimit -s unlimited`). Example for the runner:

```bash
gfortran main_orem.F orem.F report.F swx.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o orem.exe
```

---

## 5. How to Run (Quick Start)

### Step 1: Compile

```bat
REM Windows — Intel oneAPI
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx /heap-arrays /F:16777216 main_orem.F orem.F report.F swx.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:orem.exe
```

```bash
# Unix — gfortran
gfortran main_orem.F orem.F report.F swx.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o orem.exe
```

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
```

**354 checks total**, all passing as of v1.23.

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

### test_ga (71 tests)
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
- BN range carryover between zones is trust-gated (v1.21): only zones that actually predicted a re-entry re-center the search range
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

## 7. KSROP Source Files

Files in `ksrop/` are copied from [hari251086/KSROP](https://github.com/hari251086/KSROP). To update after KSROP changes:

```bash
cp ../KSROP/Subrouts.F ksrop/
cp ../KSROP/Legendre.F ksrop/
# propagate_ks.F is a refactored version of driver_KS.F — manual sync
```

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

---

## 9. References

- Sellamuthu, H. (2019) Regularized Astrodynamics Using Kustaanheimo-Stiefel Space, Ph.D. Thesis, Karunya Institute of Technology and Sciences
- Sellamuthu, H., Sharma, R.K. & Arumugam, S. Optimal re-entry time prediction of RSO from HEO, Advances in Space Research (submitted)
- Stiefel, E.L. & Scheifele, G. (1971) Linear and Regular Celestial Mechanics, Springer-Verlag
