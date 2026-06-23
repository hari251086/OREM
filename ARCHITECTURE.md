# OREM — Application & Physical Architecture

## 1. System Overview

**OREM** (Optimal Regularized re-Entry estimation Method) predicts the atmospheric re-entry time of resident space objects (RSO) decaying from highly elliptical orbits (HEO). It combines a regularized orbit propagator (KSROP) with response surface methodology (RSM) and genetic algorithm (GA) optimization to compensate for the low accuracy of Two-Line Element (TLE) catalog data.

### Problem Statement

RSO in HEO (GTO, Molniya, SSTO) experience complex orbital evolution under luni-solar gravity, oblateness, and atmospheric drag. Their re-entry times are highly sensitive to:
- **Initial conditions** — TLE accuracy is limited (~km-level for HEO)
- **Ballistic coefficient** — unknown tumbling state, surface area uncertainty
- **Eccentricity** — small errors amplify through luni-solar resonance dynamics

OREM treats re-entry prediction as an optimization problem: find the (eccentricity, ballistic coefficient) pair that best fits the observed TLE orbital evolution, then propagate forward to re-entry.

### Target Accuracy

< 5% relative prediction error (RPE), validated against 4 known re-entry cases spanning GTO, Molniya-like, and medium-eccentricity orbits.

---

## 2. Application Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         OREM DRIVER (orem.F)                     │
│                                                                   │
│  Input: TLE file, NORAD ID, config file                         │
│  Output: predicted re-entry epoch, optimal e/Bcoeff, RPE         │
│                                                                   │
│  ┌──────────────┐  ┌───────────────┐  ┌────────────────────┐    │
│  │ TLE Evolution │→ │ Zone Selection │→ │ RSM Surface Gen.   │    │
│  │ (tle_evol.F) │  │ (zone_sel.F)  │  │ (rsm.F)            │    │
│  └──────────────┘  └───────────────┘  └────────┬───────────┘    │
│                                                  │                │
│                                        ┌─────────▼─────────┐    │
│                                        │ GA Optimization    │    │
│                                        │ (ga.F)             │    │
│                                        │ Design vars: e, Bc │    │
│                                        └─────────┬─────────┘    │
│                                                  │                │
│                                        ┌─────────▼─────────┐    │
│                                        │ KSROP Propagation  │    │
│                                        │ (propagate_ks.F)   │    │
│                                        │ → until re-entry   │    │
│                                        └─────────┬─────────┘    │
│                                                  │                │
│                                        ┌─────────▼─────────┐    │
│                                        │ RPE Computation    │    │
│                                        │ (rpe.F)            │    │
│                                        └───────────────────┘    │
└─────────────────────────────────────────────────────────────────┘
```

### Module Descriptions

| Module | File | Purpose | I/O |
|--------|------|---------|-----|
| **TLE Evolution** | `tle_evolution.F` | Process TLE history for single NORAD ID → osculating orbital evolution time-series | In: TLE file, NORAD ID. Out: epoch[], a[], e[], I[], Ω[], ω[], ha[], hp[], Λ_S[] |
| **Mean Elements** | `mean_elements.F` | Sliding-window average of osculating elements → mean apogee/perigee | In: osculating time-series. Out: mean ha[], mean hp[] |
| **Zone Selection** | `zone_select.F` | Identify TLE epoch intervals with quasi-linear mean apogee decay | In: mean ha[], Λ_S[]. Out: zone boundaries (start/end indices) |
| **RSM Surfaces** | `rsm.F` | Generate 9 mean-apogee surfaces per zone from (e, B_coeff) grid | In: zone, TLE elements, e-bounds, B-bounds. Out: surfaces, polynomial fit |
| **GA Optimizer** | `ga.F` | Binary-coded genetic algorithm to find optimal (e, B_coeff) | In: fitness function, bounds, GA params. Out: optimal e, B_coeff |
| **KSROP Propagator** | `ksrop/propagate_ks.F` | KS regular elements orbit propagation until re-entry | In: initial state, force config, atm tables. Out: trajectory, exit_code |
| **RPE Metric** | `rpe.F` | Compute relative prediction error against observed re-entry | In: predicted epoch, observed epoch, zone epoch. Out: RPE value |
| **OREM Driver** | `orem.F` | Main program orchestrating the full pipeline | In: config file. Out: results to stdout + output file |

---

## 3. Physical Architecture

### 3.1 File Structure

```
OREM/
├── orem.F                        Main driver program
├── tle_evolution.F               Batch TLE → orbital evolution
├── mean_elements.F               Sliding-window mean elements
├── zone_select.F                 Zone selection algorithm
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
│   ├── EGM2008_to2190_TideFree   Geopotential coefficients (optional)
│   └── orem.cfg                  OREM configuration file
│
├── test_propagate_ks.F           Propagator subroutine tests
├── test_orem.F                   OREM integration tests (4 re-entry cases)
├── lint_check.sh                 Code lint checks
├── test_all.sh                   Unified test runner
├── Makefile                      Build system
└── README.md                     Documentation
```

### 3.2 Data Flow

```
TLE File (input)
    │
    ▼
┌──────────────────┐
│ read_tle()       │  ← TLEread.F
│ tle2sv()         │  ← SGP4/SDP4 conversion
│ sun_azimuth()    │  ← Subrouts.F
└────────┬─────────┘
         │  Osculating elements: epoch, a, e, I, Ω, ω, ha, hp, Λ_S
         ▼
┌──────────────────┐
│ mean_elements()  │  Sliding window over 1 orbital period
└────────┬─────────┘
         │  Mean: ha_mean[], hp_mean[]
         ▼
┌──────────────────┐
│ zone_select()    │  Detect linear apogee decay regions
└────────┬─────────┘
         │  Zones: [(start_idx, end_idx, epoch_start, epoch_end)] × N_zones
         ▼
┌──────────────────────────────────────────────────┐
│ FOR EACH ZONE:                                    │
│                                                    │
│   ┌────────────────┐                              │
│   │ RSM Surface Gen │  9 propagation runs per zone │
│   │                 │  3 values of e × 3 of Bcoeff │
│   └────────┬───────┘                              │
│            │  9 mean-apogee curves                  │
│            ▼                                        │
│   ┌────────────────┐                              │
│   │ GA Optimization │  Fitness: |propagated - TLE| │
│   │ pop=4, gen=500  │  Design vars: e, Bcoeff     │
│   └────────┬───────┘                              │
│            │  Optimal: e_opt, Bcoeff_opt            │
│            ▼                                        │
│   ┌────────────────────────────────┐              │
│   │ propagate_ks(optimal ICs)      │              │
│   │ → propagate until re-entry     │              │
│   │   (altitude < 80 km)           │              │
│   └────────┬───────────────────────┘              │
│            │  Predicted re-entry epoch              │
│            ▼                                        │
│   ┌────────────────┐                              │
│   │ RPE Computation │                              │
│   └────────────────┘                              │
└──────────────────────────────────────────────────┘
         │
         ▼
    Results: predicted epoch, optimal e/Bcoeff, RPE per zone
```

### 3.3 Subroutine Interfaces

#### tle_evolution
```fortran
subroutine tle_evolution(tle_file, norad_id, maxpts,
     &                   epochs, a_osc, e_osc, inc_osc,
     &                   raan_osc, aop_osc,
     &                   ha_osc, hp_osc, lambda_s,
     &                   npts, ierr)
c  Input:
c    tle_file    : character*120 — path to TLE file
c    norad_id    : integer — NORAD catalog number to filter
c    maxpts      : integer — dimension of output arrays
c  Output:
c    epochs(maxpts)    : Julian dates
c    a_osc..lambda_s   : osculating elements + Sun azimuth
c    npts              : actual number of points
c    ierr              : 0=ok
```

#### mean_elements
```fortran
subroutine mean_elements(epochs, ha_osc, hp_osc, npts,
     &                   ha_mean, hp_mean)
c  Sliding-window average over one orbital period.
c  Input:  epochs(npts), ha_osc(npts), hp_osc(npts)
c  Output: ha_mean(npts), hp_mean(npts)
```

#### zone_select
```fortran
subroutine zone_select(epochs, ha_mean, lambda_s, npts,
     &                 nzones_max,
     &                 zone_start, zone_end, nzones)
c  Input:  time-series + max zones desired
c  Output: zone_start(nzones), zone_end(nzones) — index pairs
c          nzones — actual number found
```

#### rsm_generate
```fortran
subroutine rsm_generate(
     &   e_bounds, bc_bounds,
     &   zone_start_idx, zone_end_idx,
     &   epochs, e_osc, a_osc, inc_osc, raan_osc, aop_osc,
     &   ha_osc, npts,
     &   force_config, atm_tables,
     &   surfaces, nsurfaces)
c  Generates 9 mean-apogee surfaces per zone.
c  Calls propagate_ks internally for each (e, Bcoeff) combination.
```

#### ga_optimize
```fortran
subroutine ga_optimize(fitness_func, nvars,
     &                 lower_bounds, upper_bounds,
     &                 pop_size, ngen, nbits, pc, pm,
     &                 optimal_vars, optimal_fitness)
c  Binary-coded GA with generic fitness function interface.
c  Input:  bounds, GA parameters
c  Output: optimal variable values, fitness value
```

#### propagate_ks
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
c  Already implemented. exit_code: 0=normal, 1=reentry, 2=divergence
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
| Eccentricity | e | ±0.003 around TLE value | SGP4/SDP4 reconstruction error |
| Ballistic coefficient | B_coeff | 30–160 kg/m² | Unknown tumbling state, shape |

### 5.2 Objective Function (Fitness)

For a given zone with N_obs TLE observations:

```
fitness(e, Bcoeff) = 1 / (1 + RMSE)

where RMSE = sqrt( (1/N_obs) * Σ (ha_propagated_i - ha_observed_i)² )
```

The GA maximizes fitness → minimizes RMSE between propagated and observed mean apogee altitude within the zone.

### 5.3 Response Surface Model

For each zone, 9 propagation runs map the (e, B_coeff) → mean_apogee_altitude response:

```
         B₁    B₂    B₃
    e₁  [run1] [run2] [run3]
    e₂  [run4] [run5] [run6]
    e₃  [run7] [run8] [run9]
```

Each run produces a mean-apogee curve over the zone duration. The observed TLE mean apogee must be bounded within the 9 surfaces. First-order polynomial fit for linear zones.

### 5.4 GA Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Population size | 4 | Small — only 2 variables |
| Generations | 500 | Sufficient for convergence |
| Bit encoding | 80 bits | 40 bits per variable → ~10⁻¹² resolution |
| Crossover probability | 0.8 | Standard |
| Mutation probability | 0.01 | Low — maintain diversity |

### 5.5 Zone Selection Strategy

Zone selection depends on the orbital dynamics regime:

| Regime | Indicator | Zone Placement |
|---|---|---|
| Low-inclination GTO (i < 15°) | Solar apsidal resonance (Λ_S u-turn within 180°) | After the u-turn, when apogee decays linearly |
| Medium-inclination HEO (15° < i < 50°) | Λ_S completes one revolution | After first revolution, during linear apogee decay |
| High-inclination (i > 50°) | No resonance pattern | During any sustained linear apogee decay period |
| Low-eccentricity (e < 0.4) | Continuous decay from launch | Anywhere in the evolution |

---

## 6. Development Plan

### Phase 1: Foundation (Issues #1, #2) — Unblocked

| Task | Issue | Effort | Description |
|---|---|---|---|
| Batch TLE processing | #1 | 3 days | `tle_evolution.F`: read TLE, filter by NORAD, convert via SGP4/SDP4, compute ha/hp/Λ_S |
| Mean elements | #2 | 1 day | `mean_elements.F`: sliding-window average subroutine |

**Deliverable:** Given a TLE file + NORAD ID → complete orbital evolution time-series with mean elements.

### Phase 2: Optimization Core (Issues #3, #4) — Unblocked (#4), #3 needs Phase 1

| Task | Issue | Effort | Description |
|---|---|---|---|
| Genetic Algorithm | #4 | 3 days | `ga.F`: binary-coded GA with generic fitness interface |
| Zone selection | #3 | 2 days | `zone_select.F`: detect linear apogee decay regions from mean elements + Λ_S |

**Deliverable:** GA optimizer + zone boundaries from orbital evolution data.

### Phase 3: Surface Generation (Issue #5) — Needs Phase 1+2

| Task | Issue | Effort | Description |
|---|---|---|---|
| RSM surfaces | #5 | 3 days | `rsm.F`: generate 9 mean-apogee surfaces per zone, polynomial fit, bound checking |

**Deliverable:** Response surfaces bounding observed TLE data.

### Phase 4: Integration (Issue #6) — Needs Phase 1–3

| Task | Issue | Effort | Description |
|---|---|---|---|
| OREM driver | #6 | 3 days | `orem.F`: main program orchestrating full pipeline (TLE → zones → RSM → GA → propagation → re-entry) |

**Deliverable:** Complete working OREM tool.

### Phase 5: Validation (Issues #7, #8) — Needs Phase 4

| Task | Issue | Effort | Description |
|---|---|---|---|
| RPE metric | #7 | 0.5 day | `rpe.F`: compute relative prediction error |
| Test suite | #8 | 2 days | `test_orem.F`: 4 validated re-entry cases (NORAD 35497, 37151, 39615, 42928), target RPE < 5% |

**Deliverable:** Validated OREM with published accuracy.

### Timeline Summary

```
Phase 1 ████████░░░░░░░░░░░░  (4 days)  — Foundation
Phase 2 ░░░░████████░░░░░░░░  (5 days)  — Optimization core
Phase 3 ░░░░░░░░░░░█████░░░░  (3 days)  — Surface generation
Phase 4 ░░░░░░░░░░░░░░░█████  (3 days)  — Integration
Phase 5 ░░░░░░░░░░░░░░░░░░██  (2.5 days)— Validation
                                Total: ~18 days
```

### Dependency Graph

```
#1 (TLE evolution) ──┐
                     ├──→ #3 (Zone selection) ──┐
#2 (Mean elements) ──┘                          ├──→ #5 (RSM) ──┐
                                                │                ├──→ #6 (OREM) ──→ #7 (RPE)
#4 (GA optimizer) ──────────────────────────────┘                │                  #8 (Tests)
```

---

## 7. Configuration File Format (`orem.cfg`)

```
input/example_35497.tle.txt       ! TLE file path
35497                             ! Target NORAD ID
2016 10 31 00 00 0.0              ! Observed re-entry epoch (for validation)
4                                 ! Number of zones
7.0                               ! Mean apogee altitude bound (km) for RSM
30.0 160.0                        ! Ballistic coefficient bounds (kg/m²)
4 500 80 0.8 0.01                 ! GA params: pop, gen, bits, Pc, Pm
50 2 2                            ! Geo degree, Sun degree, Moon degree
50.0 1 7.2921150d-5 3.35281066d-3 1.0  ! Drag: BN, IDRAG, WE, EPS, FR
1.2 0.01 0 1                      ! SRP: CR, AM, IPSR, ISHAD
```

---

## 8. Build & Run

### Compile
```bash
ifx orem.F tle_evolution.F mean_elements.F zone_select.F rsm.F ga.F rpe.F \
    ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F \
    /exe:orem.exe
```

### Run
```bash
./orem.exe input/orem.cfg
```

### Test
```bash
bash test_all.sh
```

---

## 9. KSROP Linkage

The `ksrop/` directory contains files copied from the KSROP repo. When KSROP is updated:

1. Copy updated files: `cp $KSROP/{Subrouts.F,Legendre.F} ksrop/`
2. If `driver_KS.F` changes, re-apply the refactoring to `propagate_ks.F`
3. Run `test_propagate_ks.exe` to verify compatibility
4. TLEread.F updates (e.g., SDP4) copy directly

The common block `/xy/` (pi, d2r, r2d, amue, AU, R_Earth) is the interface contract between KSROP files and OREM modules. `init_constants()` must be called before any KSROP subroutine.
