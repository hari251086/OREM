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
│  ┌──────────────┐  ┌───────────────┐                             │
│  │ TLE Evolution │→ │ Zone Selection │                             │
│  │ (tle_evol.F) │  │ (zone_sel.F)  │                             │
│  └──────────────┘  └───────┬───────┘                             │
│                            │                                       │
│                   ┌────────▼────────────────────────────┐         │
│                   │ RSM Surface Generation (rsm.F)      │         │
│                   │                                      │         │
│                   │  Calls propagate_ks.F 9× per zone   │         │
│                   │  (3 values of e × 3 of Bcoeff)      │         │
│                   │                                      │         │
│                   │  ┌──────────────────────────┐       │         │
│                   │  │ propagate_ks (×9/zone)   │       │         │
│                   │  │ short propagation within  │       │         │
│                   │  │ zone duration only        │       │         │
│                   │  └──────────────────────────┘       │         │
│                   │                                      │         │
│                   │  Output: 9 mean-apogee curves        │         │
│                   │        + polynomial surface fit       │         │
│                   └────────────────┬────────────────────┘         │
│                                    │                               │
│                   ┌────────────────▼────────────────────┐         │
│                   │ GA Optimization (ga.F)               │         │
│                   │                                      │         │
│                   │  Searches pre-computed surfaces      │         │
│                   │  (NO propagation — surface interp.)  │         │
│                   │  Fitness: |surface_ha - TLE_ha|      │         │
│                   │                                      │         │
│                   │  Output: optimal e, Bcoeff           │         │
│                   └────────────────┬────────────────────┘         │
│                                    │                               │
│                   ┌────────────────▼────────────────────┐         │
│                   │ Final Propagation (propagate_ks ×1) │         │
│                   │                                      │         │
│                   │  Uses optimal ICs from GA            │         │
│                   │  Long propagation → until re-entry   │         │
│                   │  (altitude < 80 km)                  │         │
│                   └────────────────┬────────────────────┘         │
│                                    │                               │
│                   ┌────────────────▼────────────────────┐         │
│                   │ RPE Computation (rpe.F)              │         │
│                   └─────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────┘
```

### Propagation Call Budget

`propagate_ks` is called at two distinct stages with different purposes:

1. **RSM surface generation**: 9 × N_zones short propagations (within each zone's time span only). These produce the mean-apogee surfaces that the GA searches.
2. **Final re-entry prediction**: 1 long propagation with optimal (e, B_coeff) from the GA, running from the zone epoch forward until altitude < 80 km.

The GA **never calls the propagator** — it only evaluates the polynomial fit of the pre-computed surfaces. This is what makes OREM computationally feasible: expensive propagation is done upfront in RSM, and GA is a cheap surface search.

### Module Descriptions

| Module | File | Purpose | Calls propagate_ks? |
|--------|------|---------|---------------------|
| **TLE Evolution** | `tle_evolution.F` | Process TLE history → osculating orbital evolution time-series | No |
| **Mean Elements** | `mean_elements.F` | Sliding-window average → mean apogee/perigee | No |
| **Zone Selection** | `zone_select.F` | Identify TLE intervals with quasi-linear mean apogee decay | No |
| **RSM Surfaces** | `rsm.F` | Generate 9 mean-apogee surfaces per zone | **Yes — 9× per zone** |
| **GA Optimizer** | `ga.F` | Search pre-computed surfaces for optimal (e, B_coeff) | **No — surface interpolation only** |
| **KSROP Propagator** | `ksrop/propagate_ks.F` | KS regular elements orbit propagation | (called by RSM and final prediction) |
| **RPE Metric** | `rpe.F` | Compute relative prediction error | No |
| **OREM Driver** | `orem.F` | Orchestrates full pipeline; calls final propagate_ks ×1 | **Yes — 1× final prediction** |

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
│   STAGE 1: RSM Surface Generation                  │
│   ┌────────────────────────────────────────┐      │
│   │ rsm_generate()                         │      │
│   │                                        │      │
│   │  FOR i=1,3 (eccentricity values):      │      │
│   │    FOR j=1,3 (Bcoeff values):          │      │
│   │      ┌─────────────────────────┐       │      │
│   │      │ propagate_ks(e_i, Bc_j) │       │      │
│   │      │ short prop within zone  │       │      │
│   │      │ → mean apogee curve_ij  │       │      │
│   │      └─────────────────────────┘       │      │
│   │                                        │      │
│   │  Output: 9 mean-apogee curves          │      │
│   │        + polynomial surface fit         │      │
│   └───────────────────┬────────────────────┘      │
│                       │                            │
│   STAGE 2: GA Surface Search (NO propagation)      │
│   ┌───────────────────▼────────────────────┐      │
│   │ ga_optimize()                          │      │
│   │                                        │      │
│   │  Fitness: evaluate polynomial surface  │      │
│   │           vs observed TLE mean apogee  │      │
│   │  (interpolation only — no propagation) │      │
│   │                                        │      │
│   │  Output: optimal e_opt, Bcoeff_opt     │      │
│   └───────────────────┬────────────────────┘      │
│                       │                            │
│   STAGE 3: Final Re-entry Propagation (×1)         │
│   ┌───────────────────▼────────────────────┐      │
│   │ propagate_ks(e_opt, Bcoeff_opt)        │      │
│   │                                        │      │
│   │ Long propagation from zone epoch       │      │
│   │ forward until altitude < 80 km         │      │
│   └───────────────────┬────────────────────┘      │
│                       │                            │
│   ┌───────────────────▼────────────────────┐      │
│   │ RPE Computation                        │      │
│   │ ε = (predicted - observed) /           │      │
│   │     (observed - last_zone_epoch)       │      │
│   └────────────────────────────────────────┘      │
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

### 5.2 RSM-GA Coupling (Computational Architecture)

The optimization has two distinct computational phases:

**Phase A — RSM Surface Generation (expensive, done once per zone):**
- 9 calls to `propagate_ks`, each propagating from the zone start epoch through the zone duration
- Each call uses a different (e, B_coeff) combination from the 3×3 grid
- Output: 9 mean-apogee altitude curves + a polynomial surface fit
- This is the **only** stage where orbit propagation occurs during optimization

**Phase B — GA Surface Search (cheap, iterative):**
- GA evaluates candidate (e, B_coeff) pairs by interpolating the pre-computed polynomial surface
- Fitness = how well the interpolated mean apogee matches the observed TLE mean apogee
- **No propagation calls** — pure arithmetic on the polynomial coefficients
- 500 generations × 4 population = 2000 fitness evaluations, all via surface interpolation

This separation makes OREM computationally feasible: 9 propagations per zone (seconds each) rather than thousands.

### 5.3 Objective Function (Fitness)

For a given zone with N_obs TLE observations, the GA evaluates:

```
fitness(e, Bcoeff) = 1 / (1 + RMSE)

where RMSE = sqrt( (1/N_obs) * Σ (ha_surface(e,Bc,t_i) - ha_observed_i)² )
```

`ha_surface(e, Bc, t_i)` is the polynomial interpolation of the 9 pre-computed RSM surfaces at the candidate (e, Bc) point and observation time t_i. No propagation is involved.

The GA maximizes fitness → minimizes RMSE between the RSM surface prediction and observed TLE mean apogee altitude within the zone.

### 5.4 Response Surface Model

For each zone, 9 propagation runs map the (e, B_coeff) → mean_apogee_altitude response:

```
         B₁    B₂    B₃
    e₁  [run1] [run2] [run3]     ← propagate_ks called 9 times
    e₂  [run4] [run5] [run6]     ← short propagation (zone duration)
    e₃  [run7] [run8] [run9]     ← each produces a mean-apogee curve
```

Each run produces a mean-apogee curve over the zone duration. The observed TLE mean apogee must be bounded within the 9 surfaces. A first-order polynomial is fitted for linear zones (higher-order for curved evolution).

The polynomial surface `ha(e, Bc, t)` is then used by the GA for cheap fitness evaluation without further propagation.

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
