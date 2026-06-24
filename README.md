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
│   └── ATM.DAT                     Atmosphere density table (60-630 km)
│
├── output/                         Runtime-generated files
│
├── test_propagate_ks.F             Tests for refactored propagator
│
├── tle_evolution.F                 Batch TLE → orbital evolution (56 tests)
├── zone_select.F                   Zone selection — linear apogee decay (28 tests)
├── test_tle_evolution.F            TLE evolution tests
├── test_zone_select.F              Zone selection tests
├── ga.F                            Binary-coded GA optimizer (71 tests)
├── test_ga.F                       GA optimizer tests
├── rsm.F                           RSM surface generation (39 tests)
├── test_rsm.F                      RSM integration tests
├── orem.F                          (planned) Main OREM driver
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
```

### Unix / gfortran

```bash
gfortran test_propagate_ks.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F -o test_propagate_ks.exe
gfortran test_tle_evolution.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_tle_evolution.exe
gfortran test_zone_select.F zone_select.F tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_zone_select.exe
gfortran test_ga.F ga.F -o test_ga.exe
gfortran test_rsm.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_rsm.exe
```

---

## 5. Running Tests

```bash
./test_propagate_ks.exe        # Propagator tests
./test_tle_evolution.exe       # TLE evolution tests (56 checks)
./test_zone_select.exe         # Zone selection tests (68 checks)
./test_ga.exe                  # GA optimizer tests (71 checks)
./test_rsm.exe                 # RSM integration tests (39 checks)
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

---

## 6. KSROP Source Files

Files in `ksrop/` are copied from [hari251086/KSROP](https://github.com/hari251086/KSROP). To update after KSROP changes:

```bash
cp ../KSROP/Subrouts.F ksrop/
cp ../KSROP/Legendre.F ksrop/
# propagate_ks.F is a refactored version of driver_KS.F — manual sync
```

---

## 7. Version History

| Version | Date | Changes |
|---|---|---|
| 0.1 | 2026-06-23 | Initial repo: propagate_ks refactored from KSROP driver_KS.F |
| 0.2 | 2026-06-23 | Batch TLE processing (`tle_evolution.F`), 56 tests, epoch dedup |
| 0.3 | 2026-06-23 | Zone selection (`zone_select.F`, `linfit`), 68 tests, 4 HEO TLE histories, max_zone_days bug fix |
| 0.4 | 2026-06-24 | GA optimizer (`ga.F`), refactored from GENESIS, 58 tests, high-e orbits validated |
| 0.5 | 2026-06-24 | RSM surface generation (`rsm.F`), 9× propagate_ks per zone, 39 tests, RSM→GA integration verified |

---

## 8. References

- Sellamuthu, H. (2019) Regularized Astrodynamics Using Kustaanheimo-Stiefel Space, Ph.D. Thesis, Karunya Institute of Technology and Sciences
- Sellamuthu, H., Sharma, R.K. & Arumugam, S. Optimal re-entry time prediction of RSO from HEO, Advances in Space Research (submitted)
- Stiefel, E.L. & Scheifele, G. (1971) Linear and Regular Celestial Mechanics, Springer-Verlag
