# OREM — Optimal Regularized re-Entry estimation Method

Optimal re-entry time prediction for resident space objects from highly elliptical orbits, using Response Surface Methodology (RSM) and Genetic Algorithm (GA) optimization with the KSROP regularized orbit propagator.

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com

---

## 1. Overview

OREM predicts re-entry times of HEO debris (GTO, Molniya, SSTO upper stages) by:

1. Processing TLE history for a target NORAD ID
2. Selecting optimal TLE zones based on solar apsidal resonance
3. Generating mean apogee surfaces via RSM (varying eccentricity and ballistic coefficient)
4. Optimizing initial conditions with GA to minimize RMS error between propagated and observed TLE apogee trajectory
5. Propagating with KSROP until re-entry (altitude < 80 km)

Target accuracy: **< 5% relative prediction error** (RPE) validated against real re-entries.

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
│   └── orem_42928.cfg             Example config file for PSLV-C39
│
├── test_propagate_ks.F             Tests for refactored propagator
│
├── tle_evolution.F                 Batch TLE → orbital evolution (56 tests)
├── zone_select.F                   Zone selection — linear apogee decay (68 tests)
├── test_tle_evolution.F            TLE evolution tests
├── test_zone_select.F              Zone selection tests
├── ga.F                            Binary-coded GA optimizer (71 tests); ld_surf fix v1.4; trajectory-matching fitness
├── test_ga.F                       GA optimizer tests
├── test_ga_sensitivity.F           GA parameter sensitivity study (pop/gen/Pm sweep, not in test suite)
├── rsm.F                           RSM surface generation (39 tests)
├── test_rsm.F                      RSM integration tests
├── main_orem.F                     Standalone runner (reads orem.cfg)
├── orem.F                          OREM driver + compute_rpe (14 tests)
├── test_orem.F                     OREM driver tests
├── test_reentry.F                  7-object re-entry validation (35 tests)
├── test_e2e.F                      End-to-end integration test, IDRAG=1 (20 tests: E1–E10 42928 + E11–E20 39615/35497)
├── test_gmat.F                      GMAT cross-validation: BN sensitivity (14 tests)
└── README.md
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

```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx test_propagate_ks.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F /exe:test_propagate_ks.exe
ifx test_tle_evolution.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_tle_evolution.exe
ifx test_zone_select.F zone_select.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_zone_select.exe
ifx test_ga.F ga.F /exe:test_ga.exe
ifx /heap-arrays /F:16777216 test_rsm.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_rsm.exe
ifx /heap-arrays /F:16777216 test_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_orem.exe
ifx /heap-arrays /F:16777216 test_reentry.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:test_reentry.exe

REM Standalone runner
ifx /heap-arrays /F:16777216 main_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:orem.exe
```

### Unix / gfortran

```bash
gfortran test_propagate_ks.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F -o test_propagate_ks.exe
gfortran test_tle_evolution.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_tle_evolution.exe
gfortran test_zone_select.F zone_select.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_zone_select.exe
gfortran test_ga.F ga.F -o test_ga.exe
gfortran test_rsm.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_rsm.exe
gfortran test_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_orem.exe
gfortran test_reentry.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_reentry.exe

# Standalone runner
gfortran main_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o orem.exe
```

---

## 5. How to Run (Quick Start)

### Step 1: Compile

```bat
REM Windows — Intel oneAPI
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx /heap-arrays /F:16777216 main_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F /exe:orem.exe
```

```bash
# Unix — gfortran
gfortran main_orem.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o orem.exe
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
4                                    <- Max number of zones
8 10.0 0.90 -1.0                    <- Zone: min_pts, max_days, R2_threshold, slope_threshold
80.0 160.0                          <- Ballistic number bounds [BN_min, BN_max]
4 200 40 16 0.8 0.01 0.123          <- GA: pop, gen, bits_e, bits_BN, Pc, Pm, seed
2 0 0                               <- Force model: geo_deg, sun_deg, moon_deg
0 7.2921150d-5 3.35281066d-3 1.0    <- Drag: IDRAG(0=off,1=on), WE, EPS_f, FR
0 0.0 0.0 0                         <- SRP: IPSR(0=off,1=on), CR, AM, ISHAD
```

### Step 3: Run

```bash
./orem.exe input/orem_42928.cfg
```

### Step 4: Read output

The output shows per-zone results:
```
Zone    Epoch (JD)     e_opt  A (m2)  Re-entry (JD)   Re-entry (UTC)    RPE(%)
   1    2458152.54  0.273541   0.500  2458543.20       2019-03-02        -0.30
   2    2458235.23  0.240538   2.180  2458545.80       2019-03-05         1.20
```

- **e_opt** — optimal eccentricity found by GA for this zone
- **A (m2)** — optimal cross-sectional area (m²)
- **Re-entry (JD/UTC)** — predicted re-entry date
- **RPE(%)** — relative prediction error vs observed (if provided)

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

**Ballistic number (BN):** BN = m/(Cd×A) in kg/m². Typical range for GTO/HEO debris: 30–200. The GA optimizes BN directly, as in the original NPOE research.

---

## 6. Running Tests

```bash
./test_propagate_ks.exe        # Propagator tests
./test_tle_evolution.exe       # TLE evolution tests (56 checks)
./test_zone_select.exe         # Zone selection tests (68 checks)
./test_ga.exe                  # GA optimizer tests (71 checks)
./test_rsm.exe                 # RSM integration tests (39 checks)
./test_orem.exe                # OREM driver tests (14 checks)
./test_reentry.exe             # 7-object re-entry validation (35 checks)
./test_e2e.exe                 # End-to-end integration test, IDRAG=1 (20 checks)
./test_gmat.exe                # GMAT cross-validation: propagator BN sensitivity (14 checks)
```

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

### test_orem (22 tests) — includes issue #12
- compute_rpe: perfect RPE=0, 10-day late RPE~1.9%, mean/std (Mode 2), zero predictions
- Error handling: bad TLE file, wrong NORAD ID
- 42928 integration: full pipeline (TLE→zone→RSM→GA→propagation), 4 zones, e_opt/a_opt/rms valid, zone epochs valid
- Failure recovery/diagnostics (#12): D15 propagator-divergence skip (BN=0 forces a division-by-zero in the drag term → NaN altitude → `zone_status=1`), D16 GA boundary detection (`zone_status=2`), D17 all-zones-fail doesn't crash the driver loop (`nzones_valid=0`)
- G2 physics-based BN floor (#12): 37151's zone 1 with the default `bn_min_init=80` — floor estimate extends `bn_lo` well below 80, letting `bn_opt(1)` reach 48.79 (structurally impossible before this change)

### test_e2e (20 tests) — Issue #16
Full pipeline with IDRAG=1, **full force model** (geo=20, sun=2, moon=3, SRP on: Cr=1.2, A/m=0.01 m²/kg, conical shadow — widened from geo=4/SRP-off in v1.11); GA minimizes trajectory RMS (v1.7):
- E1–E5: 42928 PSLV-C39 R/B (re-entry 2019-03-03): pipeline, zones, e_opt, BN in [80,160] per zone, re-entry in all 4 zones
- E6–E10: 42928 zone-0 (14 TLEs, e≈0.32, epoch 2017-09-22): per-zone BN, re-entry detected, zone-0 RPE ≈ −16%
- E11–E15: 39615 Proton-M Briz-M (re-entry 2017-09-15): pipeline, zones, e_opt, per-zone BN in [50,500], re-entry
- E16–E20: 35497 Ariane 5 ESC-A (re-entry 2016-10-31): pipeline, zones, e_opt, per-zone BN; no re-entry predicted (BN>50, within-zone fit; informational)
- BN narrows across zones (each zone reduces the search range by 50%); e_opt per zone reflects zone TLE eccentricity
- RPE printed as diagnostic (not enforced — open as issue #12)
- The full force model leaves short-window (≈7–10 day) RSM/GA-fitted BN and e values unchanged to the last decimal (drag dominates at that timescale) but measurably shifts long-duration (multi-year) re-entry propagations — e.g. zone-0's re-entry trajectory shortened by 47 revolutions and its predicted date moved ~4.7 days earlier

### test_gmat (14 tests) — Issue #11
Cross-validates propagate_ks against GMAT R2026a reference runs (`scratch_gmat/gmat_xval_42928z0.script`, run via `GmatConsole.exe`) on 42928 PSLV-C39 R/B, Zone 0 (2017-09-24). Replaces the earlier NPOE-based comparison — GMAT is the trusted ground truth established by the KSROP↔GMAT validation campaign, and matching the NPOE-era heritage research is no longer a goal for this test:
- N1-N3: BN monotonicity — higher BN → less apogee decay (each of 3 e-rows)
- N4-N6: e monotonicity — higher e → higher initial apogee (each of 3 BN columns)
- N7-N9: BN sensitivity ratio decay(BN=80)/decay(BN=160) > 1.5 (propagate_ks ~2.0; GMAT ref ~1.45-1.54)
- N10: No divergence across all 9 RSM grid runs
- N11: IDRAG=0 gives < 0.5 km drop in 7 days
- N12: All drops negative for IDRAG=1
- N13-N14: Magnitude within 50% of GMAT for BN=80 and BN=160 (propagate_ks is 69-112% of GMAT magnitude across all 9 grid points; tightened from the old factor-3 NPOE tolerance)
- **Key finding**: propagate_ks correctly models BN physics; RPE inaccuracy is due to short zone windows and TLE noise, not a propagator bug

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

---

## 9. References

- Sellamuthu, H. (2019) Regularized Astrodynamics Using Kustaanheimo-Stiefel Space, Ph.D. Thesis, Karunya Institute of Technology and Sciences
- Sellamuthu, H., Sharma, R.K. & Arumugam, S. Optimal re-entry time prediction of RSO from HEO, Advances in Space Research (submitted)
- Stiefel, E.L. & Scheifele, G. (1971) Linear and Regular Celestial Mechanics, Springer-Verlag
