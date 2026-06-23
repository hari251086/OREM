# OREM — Application & Physical Architecture

## 1. System Overview

**OREM** (Optimal Regularized re-Entry estimation Method) predicts the atmospheric re-entry time of resident space objects (RSO) decaying from highly elliptical orbits (HEO). It combines a regularized orbit propagator (KSROP) with response surface methodology (RSM) and genetic algorithm (GA) optimization to compensate for the low accuracy of Two-Line Element (TLE) catalog data.

### Problem Statement

RSO in HEO (GTO, Molniya, SSTO) experience complex orbital evolution under luni-solar gravity, oblateness, and atmospheric drag. Their re-entry times are highly sensitive to:
- **Initial conditions** — TLE accuracy is limited (~km-level for HEO)
- **Ballistic coefficient** — unknown tumbling state, cross-sectional area uncertainty
- **Eccentricity** — small errors amplify through luni-solar resonance dynamics

OREM treats re-entry prediction as an optimization problem: find the (eccentricity, cross-sectional area) pair that best fits the observed TLE orbital evolution, then propagate forward to re-entry.

### Target Accuracy

< 5% relative prediction error (RPE), validated against 4 known re-entry cases spanning GTO, Molniya-like, and medium-eccentricity orbits.

---

## 2. Application Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         OREM DRIVER (orem.F)                     │
│                                                                   │
│  Input: TLE file, NORAD ID, config file                         │
│  Output: re-entry predictions per zone, optimal e/A, RPE         │
│                                                                   │
│  ┌──────────────┐  ┌───────────────┐                             │
│  │ TLE Evolution │→ │ Zone Selection │                             │
│  │ tle_evolve()  │  │ zone_select() │                             │
│  └──────────────┘  └───────┬───────┘                             │
│                            │                                       │
│              ┌─────────────▼──────────────────────┐               │
│              │  FOR EACH ZONE (iterative loop):    │               │
│              │                                      │               │
│              │  ┌──────────────────────────────┐   │               │
│              │  │ RSM Surface Generation        │   │               │
│              │  │ 9× propagate_ks              │   │               │
│              │  │ 3 ecc × 3 area values        │   │               │
│              │  │ → 9 mean-apogee surfaces      │   │               │
│              │  └──────────────┬───────────────┘   │               │
│              │                 │                     │               │
│              │  ┌──────────────▼───────────────┐   │               │
│              │  │ GA Surface Search             │   │               │
│              │  │ (NO propagation)              │   │               │
│              │  │ TWOINT bilinear interpolation │   │               │
│              │  │ → optimal e_opt, A_opt        │   │               │
│              │  └──────────────┬───────────────┘   │               │
│              │                 │                     │               │
│              │  ┌──────────────▼───────────────┐   │               │
│              │  │ Re-entry Propagation (×1)     │   │               │
│              │  │ propagate_ks(e_opt, A_opt)    │   │               │
│              │  │ → reentry_jd(iz)              │   │               │
│              │  └──────────────┬───────────────┘   │               │
│              │                 │                     │               │
│              │  Narrow A bounds for next zone       │               │
│              │  A_range = A_range * 0.5             │               │
│              └─────────────────┬────────────────────┘               │
│                                │                                     │
│              ┌─────────────────▼────────────────────┐               │
│              │ RPE Computation                       │               │
│              │ Mode 1: vs observed re-entry          │               │
│              │ Mode 2: vs ensemble mean (no obs)     │               │
│              └───────────────────────────────────────┘               │
└─────────────────────────────────────────────────────────────────────┘
```

### Propagation Call Budget

`propagate_ks` is called at two distinct stages with different purposes:

1. **RSM surface generation**: 9 × N_zones short propagations (within each zone's time span only). These produce the mean-apogee surfaces that the GA searches.
2. **Per-zone re-entry prediction**: 1 × N_zones long propagation with optimal (e, A) from the GA, running from the zone epoch forward until altitude < 80 km.

The GA **never calls the propagator** — it only evaluates bilinear interpolation of the pre-computed surfaces. This is what makes OREM computationally feasible.

### Module Descriptions

| Module | File | Purpose | Calls propagate_ks? |
|--------|------|---------|---------------------|
| **TLE Evolution** | `tle_evolution.F` | Process TLE history → orbital evolution with epoch dedup | No |
| **Zone Selection** | `zone_select.F` | Identify linear apogee decay intervals (R² test + linfit) | No |
| **RSM Surfaces** | `rsm.F` | Generate 9 mean-apogee surfaces per zone | **Yes — 9× per zone** |
| **GA Optimizer** | `ga.F` | Search pre-computed surfaces for optimal (e, A) | **No — TWOINT interpolation only** |
| **KSROP Propagator** | `ksrop/propagate_ks.F` | KS regular elements orbit propagation | (called by RSM and re-entry prediction) |
| **RPE Metric** | `rpe.F` | Relative prediction error (two modes) | No |
| **OREM Driver** | `orem.F` | Orchestrates iterative zone loop | **Yes — 1× per zone (re-entry)** |

---

## 3. Physical Architecture

### 3.1 File Structure

```
OREM/
├── orem.F                        Main driver program
├── tle_evolution.F               Batch TLE → orbital evolution (56 tests)
├── zone_select.F                 Zone selection + linfit (68 tests)
├── rsm.F                         Response surface methodology
├── ga.F                          Genetic algorithm optimizer
├── rpe.F                         Relative prediction error
│
├── ksrop/                        KSROP propagator engine (from KSROP repo)
│   ├── propagate_ks.F            Refactored KS propagator subroutine
│   ├── Subrouts.F                Coordinate transforms, I/O, utilities
│   ├── Legendre.F                Zonal Legendre polynomials
│   └── TLEread.F                 TLE reader + SGP4/SDP4 + frame transforms
│
├── input/
│   ├── const_new.dat             Physical constants (mu, R_Earth, AU, etc.)
│   ├── ATM.DAT                   Atmosphere density table
│   ├── example_35497.tle.txt     Ariane 5 ESC-A (i=5.7°, GTO)
│   ├── example_37151.tle.txt     Long March 3B (i=24.9°, GTO)
│   ├── example_39615.tle.txt     Proton-M Briz-M (i=48.5°, HEO)
│   ├── example_42928.tle.txt     PSLV-C39 (i=19.2°, HEO)
│   ├── example_42928_zone0.tle.txt  42928 Zone 0 subset (14 TLEs)
│   ├── example_47944.tle.txt     SSO LEO test case
│   └── example_multi.tle.txt     Multi-satellite catalog (94597 entries)
│
├── test_propagate_ks.F           Propagator subroutine tests
├── test_tle_evolution.F          TLE evolution tests (56)
├── test_zone_select.F            Zone selection tests (68)
└── README.md
```

### 3.2 Data Flow

```
TLE File (input)
    │
    ▼
┌──────────────────┐
│ tle_evolve()     │  ← tle_evolution.F
│                  │  Reads TLE, filters by NORAD ID,
│                  │  deduplicates epochs (<86 sec),
│                  │  computes SMA, ha, hp, Λ_S
└────────┬─────────┘
         │  epochs, a, e, inc, raan, aop, ha, hp, Λ_S  (npts points)
         ▼
┌──────────────────┐
│ zone_select()    │  ← zone_select.F
│                  │  Sliding-window R² test for linear
│                  │  apogee decay. Parameters: min_zone_pts,
│                  │  max_zone_days, r2_threshold, slope_threshold
└────────┬─────────┘
         │  zone_start(1:nzones), zone_end(1:nzones)
         ▼
┌──────────────────────────────────────────────────────────┐
│ FOR EACH ZONE iz = 1, nzones:                             │
│                                                            │
│   Build 3×3 grid:                                          │
│     e-axis: e_mid ± δe  (δe from TLE scatter in zone)     │
│     A-axis: [A_min, A_mid, A_max]  (from driver bounds)   │
│                                                            │
│   STAGE 1: RSM Surface Generation                          │
│   ┌────────────────────────────────────────────────┐      │
│   │ rsm_generate()                                  │      │
│   │                                                  │      │
│   │  FOR i=1,3 (eccentricity):                      │      │
│   │    FOR j=1,3 (cross-sectional area):            │      │
│   │      propagate_ks(e_i, A_j, fixed_mass)         │      │
│   │      → mean apogee curve eXbY(t)                │      │
│   │                                                  │      │
│   │  Output: 9 surfaces + tleobs (observed apogee)   │      │
│   └────────────────────┬───────────────────────────┘      │
│                        │                                    │
│   STAGE 2: GA Surface Search (NO propagation)              │
│   ┌────────────────────▼───────────────────────────┐      │
│   │ ga_optimize()                                   │      │
│   │                                                  │      │
│   │  Binary-coded GA (40 bits e, 16-28 bits A)      │      │
│   │  Pop=4, Gen=500, Pc=0.8, Pm=0.01               │      │
│   │  Fitness: TWOINT bilinear interpolation of      │      │
│   │    9 surfaces vs observed TLE mean apogee       │      │
│   │  RMS = sqrt(Σ((interp - obs)/100)² / N)        │      │
│   │                                                  │      │
│   │  Output: e_opt(iz), A_opt(iz)                   │      │
│   └────────────────────┬───────────────────────────┘      │
│                        │                                    │
│   STAGE 3: Re-entry Propagation (×1 per zone)             │
│   ┌────────────────────▼───────────────────────────┐      │
│   │ propagate_ks(e_opt, A_opt, fixed_mass)          │      │
│   │ Long propagation → until altitude < 80 km       │      │
│   │ → reentry_jd(iz)                                │      │
│   └────────────────────┬───────────────────────────┘      │
│                        │                                    │
│   Narrow A bounds for next zone:                           │
│     A_range = A_range * 0.5                                │
│     A_min = A_opt - A_range/2                              │
│     A_max = A_opt + A_range/2                              │
│                                                            │
└──────────────────────────┬─────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│ RPE Computation                                           │
│                                                            │
│ Mode 1 (validation — observed re-entry known):             │
│   RPE(iz) = (reentry_jd(iz) - t_obs) /                   │
│             (t_obs - t_zone(iz)) × 100%                    │
│                                                            │
│ Mode 2 (operational — no observed re-entry):               │
│   t_mean = mean(reentry_jd(1:nzones))                     │
│   t_std  = std(reentry_jd(1:nzones))                      │
│   RPE(iz) = (reentry_jd(iz) - t_mean) /                  │
│             (t_mean - t_zone(iz)) × 100%                   │
└──────────────────────────────────────────────────────────┘
         │
         ▼
    Results: reentry_jd per zone, e_opt/A_opt per zone, RPE, t_mean, t_std
```

### 3.3 Subroutine Interfaces

#### tle_evolve (implemented, 56 tests)
```fortran
subroutine tle_evolve(tle_file, target_norad, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
c  Input:
c    tle_file       : path to TLE file
c    target_norad   : NORAD ID (-1 for all)
c    maxpts         : dimension of output arrays
c  Output:
c    epochs(maxpts)  : Julian dates (deduplicated, <86 sec)
c    a_out..hp_out   : SMA, ecc, inc, RAAN, AOP, ha, hp (from TLE mean elements)
c    lambda_s_out    : Sun azimuth angle (degrees)
c    npts            : actual count
c    ierr            : 0=ok, 1=file error, 2=no match
```

#### zone_select (implemented, 68 tests)
```fortran
subroutine zone_select(epochs, ha, npts,
     &                 min_zone_pts, max_zone_days,
     &                 r2_threshold, slope_threshold,
     &                 nzones_max,
     &                 zone_start, zone_end, nzones)
c  Input:
c    epochs(npts), ha(npts) : from tle_evolve
c    min_zone_pts           : minimum TLE observations per zone (e.g. 8)
c    max_zone_days          : maximum zone span in days (e.g. 10)
c    r2_threshold           : R² linearity threshold (e.g. 0.90)
c    slope_threshold        : negative slope cutoff km/day (e.g. -1.0)
c    nzones_max             : max zones to return
c  Output:
c    zone_start/end(nzones_max) : index pairs into epochs/ha
c    nzones                     : actual count (sorted, non-overlapping)
```

#### linfit (implemented, tested within zone_select)
```fortran
subroutine linfit(x, y, nn, slope, intercept, r2)
c  Least-squares linear fit with R² (coefficient of determination)
```

#### rsm_generate (planned — Issue #5)
```fortran
subroutine rsm_generate(
     &   epochs_zone, ha_zone, e_zone, nzone,
     &   cal0, x0, xd0,
     &   e_min, e_max, a_min, a_max,
     &   fixed_mass, cd_drag,
     &   propagator_params...,
     &   surfaces, tobs, apobs, nsurf_pts,
     &   ierr)
c  Generates 9 mean-apogee surfaces for one zone.
c  e-axis: [e_min, e_mid, e_max] — centered on zone TLE eccentricity
c  A-axis: [A_min, A_mid, A_max] — cross-sectional area (m²)
c  Calls propagate_ks 9× with fixed mass, varied (e, A).
c  Output: 9 apogee time-series + observed TLE apogee for GA.
```

#### ga_optimize (planned — Issue #4)
```fortran
subroutine ga_optimize(
     &   surfaces, tobs, apobs, nsurf_pts,
     &   e_grid, a_grid,
     &   e_bounds, a_bounds,
     &   pop_size, ngen, nbits_e, nbits_a, pc, pm,
     &   e_opt, a_opt, rms_opt)
c  Binary-coded GA searching pre-computed RSM surfaces.
c  Fitness: TWOINT bilinear interpolation vs observed apogee.
c  No propagation — pure surface evaluation.
c  Output: optimal (e, A) and RMS fitness.
```

#### compute_rpe (planned — Issue #7)
```fortran
subroutine compute_rpe(
     &   reentry_jd, zone_epochs, nzones,
     &   t_obs,
     &   rpe, t_mean, t_std,
     &   ierr)
c  Mode 1 (t_obs > 0): RPE against known re-entry time
c  Mode 2 (t_obs = 0): mean/std from ensemble, RPE vs mean
```

#### propagate_ks (implemented)
```fortran
subroutine propagate_ks(
     &   x0, xd0, cal0,
     &   nrev, istep, tole,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   BN, IDRAG, WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)
c  exit_code: 0=normal, 1=reentry (alt<80km), 2=divergence (NaN)
```

---

## 4. Physical Models (inherited from KSROP)

| Perturbation | Model | Source |
|---|---|---|
| Earth gravity | EGM2008 zonal harmonics, configurable degree (0–2190) | Subrouts.F |
| Luni-solar | Third-body Legendre expansion (degree 2–3) | Subrouts.F |
| Atmospheric drag | Oblate co-rotating exponential, tabulated density (ATM.DAT) | propagate_ks.F |
| Solar radiation pressure | Cannonball + cylindrical/conical shadow | propagate_ks.F |
| TLE conversion | SGP4 (near-Earth) + SDP4 (deep-space) → J2000 | TLEread.F |
| Sun geometry | Analytical ephemeris + orbital-plane azimuth angle | Subrouts.F |

---

## 5. Optimization Architecture

### 5.1 Design Variables

| Variable | Symbol | Range | Source of uncertainty |
|---|---|---|---|
| Eccentricity | e | e_mid ± δe (TLE scatter) | SGP4/SDP4 reconstruction error |
| Cross-sectional area | A (m²) | Iteratively narrowed | Unknown tumbling state, shape |

**Note:** Original NPOE/KS research varied mass (Cd=1, A=1 m²). OREM fixes mass (known from launch records) and varies cross-sectional area A to achieve the same ballistic coefficient range. B = Cd × A / m.

### 5.2 RSM 3×3 Grid Construction

For each zone, 9 propagation runs map the (e, A) → mean_apogee response:

```
         A₁       A₂       A₃
    e₁  [e1b1]   [e1b2]   [e1b3]     ← propagate_ks called 9 times
    e₂  [e2b1]   [e2b2]   [e2b3]     ← short propagation (zone duration)
    e₃  [e3b1]   [e3b2]   [e3b3]     ← each produces a mean-apogee curve
```

**Eccentricity axis:**
- e₂ = zone's first TLE eccentricity (center)
- δe = half the TLE eccentricity scatter within the zone
- e₁ = e₂ - δe, e₃ = e₂ + δe
- Typical δe ≈ 0.0003–0.0005 (corresponds to ~3.5–5 km in apogee altitude)

**Cross-sectional area axis:**
- A₁ = A_min, A₂ = (A_min + A_max)/2, A₃ = A_max
- Fixed mass m (from launch records), fixed Cd
- Ballistic coefficient B = Cd × A / m

### 5.3 Iterative A-bound Narrowing Across Zones

The area (ballistic coefficient) search range shrinks as successive zones refine the estimate:

```
Zone 0:  A=[A_min, A_mid, A_max]  wide range   (no prior knowledge)
Zone 1:  narrowed by ~50% around Zone 0 GA result
Zone 2:  narrowed further
  ...
Zone N:  converged — A range typically 25% of initial
```

Verified from 42928 research data (mass-equivalent):
- Zone 0: range=80, Zone 2: range=40, Zone 12: range=40, Zone 13: range=20

RSM accepts A_min/A_max as **inputs** from the driver. GA output from zone N feeds zone N+1.

### 5.4 RSM-GA Coupling

**Phase A — RSM Surface Generation (expensive, done once per zone):**
- 9 calls to `propagate_ks`, each with a different (e, A) from the 3×3 grid
- Each call propagates from zone start epoch through zone duration
- Output: 9 mean-apogee radius time-series (eXbY format)
- Also extracts observed TLE mean apogee at same time epochs (tleobs)

**Phase B — GA Surface Search (cheap, iterative):**
- GA evaluates candidate (e, A) pairs via TWOINT bilinear interpolation across the 9 surfaces
- Fitness = RMS of (interpolated - observed) mean apogee
- **No propagation calls** — pure arithmetic
- 500 generations × 4 population = 2000 fitness evaluations

**Phase C — Re-entry Propagation (one per zone):**
- Long propagation with optimal (e_opt, A_opt) from zone start until altitude < 80 km
- Produces reentry_jd(iz) — independent re-entry prediction from each zone

### 5.5 Objective Function (Fitness)

```
RMS = sqrt( (1/N) × Σ ((ha_interp(e, A, t_i) - ha_obs(t_i)) / 100)² )
```

`ha_interp` uses TWOINT 2D bilinear interpolation across the 9 pre-computed surfaces. The GA minimizes RMS.

### 5.6 GA Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Variables | 2 (e, A) | Eccentricity + cross-sectional area |
| Population size | 4 | Small — only 2 variables |
| Generations | 500 | Sufficient for convergence |
| Bit encoding | 40 bits (e) + 16-28 bits (A) | ~10⁻¹² resolution for e |
| Crossover probability | 0.8 | Standard |
| Mutation probability | 0.01 | Low — maintain diversity |

### 5.7 Zone Selection Strategy

Implemented in `zone_select.F` using sliding-window R² test:

| Parameter | Default | Purpose |
|---|---|---|
| min_zone_pts | 8 | Minimum TLE observations per zone |
| max_zone_days | 10 | Maximum zone span (days) |
| r2_threshold | 0.90 | Minimum R² for linearity |
| slope_threshold | -1.0 | Minimum negative slope (km/day) |

---

## 6. RPE Metric

### Mode 1: Validation (observed re-entry known)
```
RPE(iz) = (t_pred(iz) - t_obs) / (t_obs - t_zone(iz)) × 100%
```
Target: |RPE| < 5% for all zones.

### Mode 2: Operational (no observed re-entry)
```
t_mean = mean(t_pred(1:nzones))
t_std  = std(t_pred(1:nzones))
RPE(iz) = (t_pred(iz) - t_mean) / (t_mean - t_zone(iz)) × 100%
```
t_mean is the best estimate. t_std quantifies ensemble uncertainty.

---

## 7. Validation Cases

| Object | NORAD | Type | i (°) | e₀ | TLE entries | Known re-entry |
|---|---|---|---|---|---|---|
| Ariane 5 ESC-A | 35497 | GTO, low-i | 5.7 | 0.63 | 3093 | 2016-10-31 |
| Long March 3B | 37151 | GTO, mid-i | 24.9 | 0.56 | 3948 | 2015-12-04 |
| Proton-M Briz-M | 39615 | HEO, high-i | 48.5 | 0.68 | 2853 | 2017-09-16 |
| PSLV-C39 | 42928 | HEO, low-e | 19.2 | 0.32 | 2693 | 2019-03-03 |

---

## 8. Development Status

| Issue | Component | Status | Tests |
|---|---|---|---|
| #1 | Batch TLE processing (`tle_evolution.F`) | **Done** | 56 |
| #2 | Mean orbital elements | **Closed** — TLE mean elements used directly | — |
| #3 | Zone selection (`zone_select.F` + `linfit`) | **Done** | 68 |
| #4 | Genetic Algorithm (`ga.F`) | Open — unblocked | — |
| #5 | Response Surface Methodology (`rsm.F`) | Open — unblocked | — |
| #6 | OREM driver (`orem.F`) | Open — blocked by #4, #5 | — |
| #7 | RPE metric (`rpe.F`) | Open — blocked by #6 | — |
| #8 | Test suite (4 re-entry cases) | Open — blocked by #6 | — |

---

## 9. Configuration File Format (`orem.cfg`)

```
input/example_35497.tle.txt       ! TLE file path
35497                             ! Target NORAD ID
2016 10 31 00 00 0.0              ! Observed re-entry epoch (0 if unknown)
4                                 ! Max number of zones
8 10.0 0.90 -1.0                  ! Zone params: min_pts, max_days, R2, slope
0.5 5.0                           ! Area bounds (m²) — A_min, A_max
1000.0                            ! Fixed spacecraft mass (kg)
2.2                               ! Drag coefficient Cd
4 500 40 16 0.8 0.01              ! GA: pop, gen, bits_e, bits_A, Pc, Pm
50 2 2                            ! Geo degree, Sun degree, Moon degree
1 7.2921150d-5 3.35281066d-3 1.0  ! IDRAG, WE, EPS, FR
1.2 0.01 0 1                      ! SRP: CR, AM, IPSR, ISHAD
```

---

## 10. Build & Run

### Compile (Intel oneAPI ifx 2025.0)
```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\Fortran\compiler\2025.0\env\vars.bat"

ifx orem.F tle_evolution.F zone_select.F rsm.F ga.F rpe.F ^
    ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F ^
    /exe:orem.exe
```

### Run
```bash
./orem.exe input/orem.cfg
```

### Test
```bash
./test_propagate_ks.exe        # Propagator tests
./test_tle_evolution.exe       # TLE evolution tests (56)
./test_zone_select.exe         # Zone selection tests (68)
```

---

## 11. KSROP Linkage

The `ksrop/` directory contains files copied from the KSROP repo. When KSROP is updated:

1. Copy updated files: `cp $KSROP/{Subrouts.F,Legendre.F,TLEread.F} ksrop/`
2. If `driver_KS.F` changes, re-apply the refactoring to `propagate_ks.F`
3. Run `test_propagate_ks.exe` to verify compatibility

The common block `/xy/` (pi, d2r, r2d, amue, AU, R_Earth) is the interface contract between KSROP files and OREM modules. `init_constants()` must be called before any KSROP subroutine.
