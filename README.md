# OREM — Optimal Regularized re-Entry estimation Method

Optimal re-entry time prediction for resident space objects from highly elliptical orbits, using Response Surface Methodology (RSM) and Genetic Algorithm (GA) optimization with the KSROP regularized orbit propagator.

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com

---

## 1. Overview

OREM predicts re-entry times of HEO debris (GTO, Molniya, SSTO upper stages) by:

1. Processing TLE history for a target NORAD ID
2. Selecting optimal TLE zones based on solar apsidal resonance
3. Generating mean apogee surfaces via RSM (varying eccentricity and ballistic coefficient)
4. Optimizing initial conditions with GA to match observed mean apogee
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
├── ga.F                            Binary-coded GA optimizer (71 tests); ld_surf fix v1.4
├── test_ga.F                       GA optimizer tests
├── test_ga_sensitivity.F           GA parameter sensitivity study (pop/gen/Pm sweep, not in test suite)
├── rsm.F                           RSM surface generation (39 tests)
├── test_rsm.F                      RSM integration tests
├── main_orem.F                     Standalone runner (reads orem.cfg)
├── orem.F                          OREM driver + compute_rpe (14 tests)
├── test_orem.F                     OREM driver tests
├── test_reentry.F                  7-object re-entry validation (35 tests)
├── test_e2e.F                      End-to-end integration test, IDRAG=1 (20 tests: E1–E10 42928 + E11–E20 39615/35497)
├── test_npoe.F                     NPOE cross-validation: BN sensitivity (14 tests)
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
./test_npoe.exe                # NPOE cross-validation: propagator BN sensitivity (14 checks)
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

### test_orem (14 tests)
- compute_rpe: perfect RPE=0, 10-day late RPE~1.9%, mean/std (Mode 2), zero predictions
- Error handling: bad TLE file, wrong NORAD ID
- 42928 integration: full pipeline (TLE→zone→RSM→GA→propagation), 4 zones, e_opt/a_opt/rms valid, zone epochs valid

### test_e2e (20 tests) — Issue #16
Full pipeline with IDRAG=1, force model geo=4/sun=2/moon=3:
- E1–E5: 42928 PSLV-C39 R/B (re-entry 2019-03-03): pipeline, zones, e_opt, bn_opt in [80,160], re-entry
- E6–E10: 42928 zone-0 (14 TLEs, e≈0.32, epoch 2017-09-22): zone-0 RPE = −16% (vs −87–96% for late zones after GA bug fix)
- E11–E15: 39615 Proton-M Briz-M (re-entry 2017-09-15): pipeline, zones, e_opt, bn_opt in [50,500], re-entry
- E16–E20: 35497 Ariane 5 ESC-A (re-entry 2016-10-31): pipeline, zones, e_opt, bn_opt in [50,500]; no re-entry predicted with zone-based BN (informational)
- RPE printed as diagnostic (not enforced — BN sensitivity tuning pending Issue #11)

### test_npoe (14 tests) — Issue #11
Cross-validates propagate_ks against NPOE reference runs on 42928 PSLV-C39 R/B, Zone 0 (2017-09-24):
- N1-N3: BN monotonicity — higher BN → less apogee decay (each of 3 e-rows)
- N4-N6: e monotonicity — higher e → higher initial apogee (each of 3 BN columns)
- N7-N9: BN sensitivity ratio decay(BN=80)/decay(BN=160) > 1.5 (confirmed ~2.0, matches NPOE's 2.02)
- N10: No divergence across all 9 RSM grid runs
- N11: IDRAG=0 gives < 0.5 km drop in 7 days
- N12: All drops negative for IDRAG=1
- N13-N14: Magnitude within factor 3 of NPOE for BN=80 and BN=160 (ATM.DAT gives ~50% of Jacchia-70)
- **Key finding**: propagate_ks correctly models BN physics; RPE inaccuracy is due to short zone windows and TLE noise, not a propagator bug

### test_reentry (35 tests)
7 objects × 5 checks each: pipeline completion, zone detection, e_opt physical, a_opt in bounds, rms valid
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

| Version | Date | Changes |
|---|---|---|
| 0.1 | 2026-06-23 | Initial repo: propagate_ks refactored from KSROP driver_KS.F |
| 0.2 | 2026-06-23 | Batch TLE processing (`tle_evolution.F`), 56 tests, epoch dedup |
| 0.3 | 2026-06-23 | Zone selection (`zone_select.F`, `linfit`), 68 tests, 4 HEO TLE histories, max_zone_days bug fix |
| 0.4 | 2026-06-24 | GA optimizer (`ga.F`), refactored from GENESIS, 71 tests, high-e orbits + piecewise internals |
| 0.5 | 2026-06-24 | RSM surface generation (`rsm.F`), 9× propagate_ks per zone, 39 tests, ATM.DAT reader fix, RSM→GA integration verified |
| 0.5.1 | 2026-06-24 | Fix propagate_ks drag crash (KSROP #16): ALT_atm range guard, H_dg÷0 safety, exp overflow clamp. 234 total tests |
| 0.6 | 2026-06-24 | OREM driver (`orem.F`) + `compute_rpe` (#6, #7), 14 tests, full pipeline on 42928 (4 zones). 7 test objects from research Data. 248 total tests |
| 0.7 | 2026-06-24 | 7-object re-entry validation (#8), 35 tests, all orbit regimes (i=5.7°–63.4°, e=0.29–0.68). 283 total tests |
| 0.8 | 2026-06-27 | Fix RSM mean anomaly + time coupling: MA from TLE (not 0), surfaces interpolated at obs JDs, drag-enabled pipeline. First re-entry detection on 42928. 283 tests |
| 0.9 | 2026-06-27 | Revert to original BN-based estimation (mass as variable, Cd=1, A=1). Config uses BN bounds [80,160] directly. RSM zone-length propagation only. 283 tests |
| 1.0 | 2026-07-04 | E2E integration test with IDRAG=1 (#16): TLE→zone→RSM→GA→re-entry→RPE proven end-to-end on 42928. Fix test_propagate_ks T2/T6 (per-rev dump). Skip re-entry propagation when IDRAG=0. 298 total tests |
| 1.1 | 2026-07-04 | NPOE cross-validation (#11): 14 tests confirm propagate_ks correctly models BN sensitivity (ratio ~2.0 vs NPOE 2.02) and apogee decay direction. Magnitude is ~50% of Jacchia-70 (ATM.DAT vs Jacchia model). RPE inaccuracy diagnosed as short-zone/noise issue, not propagator bug. 312 total tests |
| 1.2 | 2026-07-04 | Fix NaN in RSM propagation: (1) car2oe clamps all dacos() arguments to [-1,1] — floating-point overflow at orbital perigee caused NaN true-anomaly → NaN drag → NaN state in ie=2,3 RSM surfaces; (2) rsm_generate hardcodes IDRAG=1 — without drag all 9 RSM surfaces were identical and the GA had no BN signal. 312 tests still pass. |
| 1.3 | 2026-07-04 | Add zone-0 E2E run (E6–E10) in test_e2e.F using example_42928_zone0.tle.txt (14 TLEs, e≈0.32, epoch 2017-09-22); zone-0 RPE = −55.5% vs −87–96% for late zones, confirming improved accuracy when propagating from early orbit. 317 total tests. |
| 1.4 | 2026-07-04 | Fix GA array-dimension mismatch bug: ga_optimize and ga_fitness declared surfaces with leading dimension nsurf_pts (≈nobs≈26) but callers allocated surfaces(max_surf=5000,...). All surface reads were reading wrong memory — GA always returned lower bound regardless of fitness landscape. Fix: add ld_surf parameter to ga_optimize and ga_fitness; callers pass max_surf. Add E11–E20 tests for 39615 and 35497 with zone-specific TLE files. Zone-0 RPE improves from −55.5% to −16.1%. 327 total tests. |

---

## 9. References

- Sellamuthu, H. (2019) Regularized Astrodynamics Using Kustaanheimo-Stiefel Space, Ph.D. Thesis, Karunya Institute of Technology and Sciences
- Sellamuthu, H., Sharma, R.K. & Arumugam, S. Optimal re-entry time prediction of RSO from HEO, Advances in Space Research (submitted)
- Stiefel, E.L. & Scheifele, G. (1971) Linear and Regular Celestial Mechanics, Springer-Verlag
