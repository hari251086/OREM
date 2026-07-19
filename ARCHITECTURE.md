# OREM — Application & Physical Architecture

*(current as of v1.23, 2026-07-17; the Version History in README.md is the authoritative changelog)*

## 1. System Overview

**OREM** (Optimal Regularized re-Entry estimation Method) predicts the atmospheric re-entry time of resident space objects (RSO) decaying from highly elliptical orbits (HEO). It combines a regularized orbit propagator (KSROP) with response surface methodology (RSM) and genetic algorithm (GA) optimization to compensate for the low accuracy of Two-Line Element (TLE) catalog data.

### Problem Statement

RSO in HEO (GTO, Molniya, SSTO) experience complex orbital evolution under luni-solar gravity, oblateness, and atmospheric drag. Their re-entry times are highly sensitive to:
- **Initial conditions** — TLE accuracy is limited (~km-level for HEO)
- **Ballistic number** BN = m/(Cd·A) — unknown tumbling state, cross-sectional area uncertainty, and attitude regime that changes over the object's life
- **Eccentricity** — small errors amplify through luni-solar resonance dynamics

OREM treats re-entry prediction as an optimization problem: per decay zone, find the (eccentricity, ballistic number) pair that best fits the observed TLE apogee evolution, then propagate each zone's fit forward to re-entry. The **latest zone's prediction is the primary estimate** — it carries the shortest extrapolation and the freshest attitude/altitude regime.

### Achieved Accuracy (v1.21) and the v1.22–v1.23 Frontier

Latest-zone RPE **median 2.4%, mean 4.1%, worst object 10.4%** across the 7-object validation campaign (full force model, 8 zones, `scratch_rpe/rpe_campaign_8zone_gated.csv`). This required, in sequence: a working GA (population 20, v1.15), a correct atmosphere table (Jacchia-71, v1.17), a correct drag phase (v1.18), the latest-zone estimator (v1.20), and a trust-gated BN-range carryover (v1.21).

That 7-object set is drag-dominated and moderate-inclination. Testing a new object outside that regime (33587, i=65° critical inclination, hp descending 616→341 km over 131 days, v1.22) exposed the next accuracy frontier: **epoch-resolved space weather (v1.23, §4.1)** replaced the static-atmosphere assumption, but a decisive GMAT cross-check (2026-07-17, issue #27) showed the object's *observed* TLE mean elements lose ~1900 km of semi-major axis with simultaneous circularization over the same window — the signature of substantial drag, not a missing lunisolar term — while both OREM and an independent GMAT full-force propagation predict almost none at the object's currently-assumed BN=55 kg/m². The open question is therefore **BN identifiability for this object/window** (§7), not gravity-model fidelity; the 7-object validation numbers above are unaffected.

---

## 2. Application Architecture

```
┌───────────────────────────────────────────────────────────────────┐
│                   OREM DRIVER (orem.F: orem_run)                    │
│                                                                     │
│  Input : TLE file, NORAD ID, config                                 │
│  Output: per-zone (e_opt, bn_opt, reentry_jd, rpe, zone_status),    │
│          ensemble t_mean/t_std, nzones_valid                        │
│                                                                     │
│  ┌──────────────┐  ┌───────────────┐  ┌──────────────────────────┐ │
│  │ TLE Evolution │→ │ Zone Selection │→ │ G2 BN floor (zone 1 only)│ │
│  │ tle_evolve()  │  │ zone_select() │  │ estimate_bn_floor()      │ │
│  └──────────────┘  └───────────────┘  └──────────────────────────┘ │
│                                                                     │
│  FOR EACH ZONE (up to nzones_max, 8 recommended):                   │
│    1. rsm_generate(): 9× propagate_ks over the zone span            │
│       (3 ecc × 3 BN grid → 9 mean-apogee surfaces, pre-interpolated │
│        at the TLE observation times; scratch buffers zeroed per     │
│        call — v1.14)                                                │
│    2. Diagnostics: RSM envelope bounds observations? (zone_status)  │
│    3. ga_optimize(): pop=20 GA on the surfaces (TWOINT bilinear     │
│       interpolation, NO propagation) → e_opt, bn_opt                │
│       + GA-boundary detection (15% of range → zone_status=2)        │
│    4. propagate_ks(e_opt, bn_opt): long propagation → reentry_jd    │
│    5. Trust-gated BN carryover (v1.21): re-center the next zone's   │
│       BN range on this zone's fit ONLY if this zone actually        │
│       predicted a re-entry (drag on) / is unflagged (drag off);     │
│       widen ×1.5 if the fit sat at a search boundary, else ×0.5     │
│                                                                     │
│  compute_rpe(): per-zone RPE + ensemble t_mean/t_std                │
└───────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────┐
│  REPORT (report.F: orem_report — called by main_orem.F)             │
│  output/OREM_<NORAD>_<DATE>.txt:                                    │
│    config echo · per-zone table (epoch, e, BN, re-entry, RPE,       │
│    status) · PRIMARY estimate (latest zone) + its RPE ·             │
│    ensemble mean ± std + relative spread · status legend            │
└───────────────────────────────────────────────────────────────────┘
```

### Propagation Call Budget

`propagate_ks` is called at three stages:

1. **G2 BN floor** (zone 1 only): one short trial run at BN=100 to calibrate the physics-based floor numerically against the propagator's own conventions.
2. **RSM surface generation**: 9 × N_zones short propagations (each spanning only its zone).
3. **Per-zone re-entry prediction**: 1 × N_zones long propagation with the fitted (e, BN), until altitude < 80 km or a 5-year cap.

The GA **never calls the propagator** — it evaluates bilinear interpolation of the pre-computed surfaces. This is what makes OREM computationally feasible.

### Module Descriptions

| Module | File | Purpose | Calls propagate_ks? |
|--------|------|---------|---------------------|
| **TLE Evolution** | `tle_evolution.F` | TLE history → orbital evolution with epoch dedup | No |
| **Zone Selection** | `zone_select.F` | Linear apogee-decay windows, top-R² candidates | No |
| **RSM Surfaces** | `rsm.F` | 9 mean-apogee surfaces per zone, interpolated at obs times | **Yes — 9× per zone** |
| **GA Optimizer** | `ga.F` | Search surfaces for optimal (e, BN); pop=20 | **No — TWOINT only** |
| **KSROP Propagator** | `ksrop/propagate_ks.F` | KS regular-elements propagation | (called by G2/RSM/re-entry) |
| **OREM Driver** | `orem.F` | Zone loop, G2 floor, diagnostics, trust-gated carryover, `compute_rpe` | **Yes — G2 trial + 1× per zone** |
| **Report** | `report.F` | Prediction report with latest-zone primary estimate | No |
| **Runner** | `main_orem.F` | Reads `orem.cfg`, runs pipeline, writes report | (via orem_run) |

---

## 3. Data Flow

```
TLE File → tle_evolve() → epochs, a, e, i, Ω, ω, ha, hp, Λs  (deduplicated)
         → zone_select() → up to nzones_max windows with clean linear
                           apogee decay (R² ≥ 0.90 over ≤ 10 days),
                           ranked by R², sorted by epoch
         → estimate_bn_floor() → may extend zone 1's bn_lo downward
                           (floor-only; never raises; G2, v1.12)

FOR EACH ZONE:
  3×3 grid:  e-axis = e_mid ± δe (TLE scatter in the zone)
             BN-axis = [bn_lo, mid, bn_hi] (kg/m²)
  rsm_generate() → surfaces(i, ie, ibn) = mean apogee at obs time i
  zone_status: 1=propagator failure (>6 of 9 grid points diverged)
               3=RSM envelope fails to bound an observation
  ga_optimize() → (e_opt, bn_opt); zone_status=2 if within 15% of a bound
  propagate_ks(e_opt, bn_opt) → reentry_jd or 0
  trust gate → BN range for the next zone

compute_rpe():
  Mode 1 (validation): RPE(iz) = (t_pred(iz) − t_obs)/(t_obs − t_zone(iz)) × 100%
  Mode 2 (operational): t_mean/t_std over predicting zones; RPE vs t_mean
```

### Key Subroutine Interfaces (as implemented)

```fortran
call orem_run(
     &   tle_file, norad_id, t_obs_cal,
     &   nzones_max, min_zone_pts, max_zone_days,
     &   r2_thresh, slope_thresh,
     &   bn_min_init, bn_max_init, idrag_flag,
     &   ipopsize, maxgen, nbits_e, nbits_a, pcross, pmute, ga_seed,
     &   ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD, PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   reentry_jd, e_opt_out, bn_opt_out, rms_out,
     &   zone_epoch, nzones_used,
     &   zone_status, nzones_valid,          ! v1.10 — keep in every call site
     &   rpe_out, t_mean, t_std, ierr)

call orem_report(
     &   rep_file, norad_id, t_obs_cal,
     &   bn_min_init, bn_max_init, idrag_flag,
     &   reentry_jd, e_opt, bn_opt, rms_o,
     &   zone_epoch, nzones_used, zone_status, nzones_valid,
     &   rpe, t_mean, t_std, nzmax, ierr_rep)

call propagate_ks(x0, xd0, cal0, nrev, istep, tole,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   BN, IDRAG, WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD, PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump, traj_jd, traj_x, traj_xd, exit_code)
c  exit_code: 0=normal, 1=reentry (alt<80 km), 2=divergence (NaN)
```

`zone_status` codes: 0=ok, 1=skip_propfail, 2=GA_boundary, 3=RSM_nobound, 4=skip_toofewpts.

---

## 4. Physical Models

| Perturbation | Model | Notes |
|---|---|---|
| Earth gravity | EGM2008 zonal harmonics, configurable degree | `geo_coeff` reads J2..Jn from `EGM2008_to2190_TideFree` |
| Luni-solar | Third-body Legendre expansion (degree 2–3) | M&G analytic ephemerides (KSROP sync v1.8) |
| Atmospheric drag | Per-revolution King-Hele: ρ_p at the oblate perigee altitude, exp(−βae(1−cosE)) along the rev, co-rotation factor F | Density phase keyed to the **true eccentric anomaly of the state** (v1.18 — the analytic sweep it replaced dephased along decay arcs) |
| Atmosphere reference | **Static**: Jacchia-71 profile (Roberts-1971 anchors) from `input/ATM.DAT`, F10.7=72, Kp=1.0, nighttime-min T∞=626.3 K. **Epoch-resolved** (opt-in, v1.23): 2-D ρ(h,T∞)/H(h,T∞) table from `input/ATM2D.DAT`, looked up at the exospheric temperature implied by real F10.7/Kp history for the revolution's own epoch | `KSROP/gen_atm_jr71.F` (static, validated 0.80–0.95× GMAT JacchiaRoberts over 102–300 km, v1.17) / `KSROP/gen_atm2d_jr71.F` + `swx.F` (epoch-resolved, §4.1). SCH column = local −dz/dlnρ in both |
| Solar radiation pressure | Cannonball + cylindrical/conical shadow | |
| TLE conversion | SGP4/SDP4 → J2000 | `TLEread.F` |

Validation lineage: two-body/zonal/third-body/SRP validated against GMAT R2026a in the KSROP campaign; OREM-side drag magnitude validated to ~1% against an exact RK4 integration of the same drag model (`scratch_gmat/drag_ref.py`); the re-entry arc cross-checked against GMAT JacchiaRoberts (`scratch_gmat/gmat_reentry_42928z0.script`); the epoch-resolved weather mechanism hand-verified against a G5-storm exospheric temperature and a storm-vs-quiet decay smoke test (`test_sw.F`); the 33587 in-record dynamics cross-checked against an independent GMAT full-force + gravity-only propagation (`scratch_gmat/gmat_hp_33587_issue27_{full,grav}.script`, issue #27).

### 4.1 Epoch-Resolved Space Weather (v1.23, issue #26)

```
input/SW-All.csv (CelesTrak, 1957→present, ~3×/day updates,        ┐
  observed F10.7/Kp daily + PRM/PRD predicted rows to ~2041)         │  sw_load()
input/ATM2D.DAT (KSROP gen_atm2d_jr71.F: ρ(h,T∞)/H(h,T∞) over        │  atm2d_load()
  h=90–1500 km × T∞=550–1500 K grid, 291×39, shares jr71_profile.F   │  (swx.F)
  with the static generator — bit-identical 1-D profile after split) ┘
                    │                              │
                    ▼                              ▼
         common /swdat/                   common /atm2dc/
    (JD, F10.7_ADJ, F10.7_ADJ_81, Kp)      (ρ, H tables, legacy-scaled
    binary-searched by JD; predicted        internal units — drag math
    rows are monthly, Kp defaults to 2.0    downstream is untouched)
    when absent)
                    │                              │
                    └──────────────┬───────────────┘
                                   ▼
              propagate_ks (ksrop/propagate_ks.F, per-rev hook, ~line 288):
                sw_tinf(JD) → T∞ for this revolution's epoch
                atm2d_interp(h_perigee, T∞) → ρ_p, H   [istat=0]
                falls back to the legacy static ALT/DEN/SCH table
                lookup whenever nothing is loaded, T∞≤0, or the point
                falls outside the 2-D table [istat≠0] — computed ONCE
                per revolution, same slot the static path always filled
```

**Opt-in, zero-cost when unused**: `sw_tinf`/`atm2d_interp` live inside `propagate_ks.F` itself (not `swx.F`) so that executables which never call `sw_load`/`atm2d_load` link without `swx.F` and run the legacy single-table path bit-unchanged — every pre-v1.23 build, result, and test is unaffected. `main_orem.F` auto-detects both files at startup and prints ENABLED/DISABLED loudly; only opt-in executables (`orem.exe`, `test_sw.exe`) link `swx.F`.

**Status**: mechanism-validated (`test_sw.F`, 12 checks — G5-storm T∞=1216.56 K hand-verified against the Jacchia formula; W12 smoke test shows storm-epoch decay 55.4 km vs quiet-epoch 31.1 km over 7 days, correct direction and rough magnitude). Not yet regression-tested across the 7-object validation campaign (issue #26 stays open pending that run — expected neutral-to-better since those objects' fit windows sit in solar-quiet periods close to the static table's own F10.7=72 assumption).

---

## 5. Optimization Architecture

### 5.1 Design Variables

| Variable | Symbol | Range | Source of uncertainty |
|---|---|---|---|
| Eccentricity | e | e_mid ± δe (TLE scatter in the zone) | SGP4/SDP4 reconstruction error |
| Ballistic number | BN (kg/m²) | zone 1: [bn_min, bn_max] with G2 floor; later zones: trust-gated carryover | Tumbling state, attitude regime drift |

BN is optimized directly (mass-as-variable convention: Cd=1, A=1 m², DryMass=BN), as in the original NPOE-era research. Fitted values on the validation set fall in ~20–100 kg/m² per zone.

### 5.2 RSM 3×3 Grid

9 short propagations per zone map (e, BN) → mean-apogee history; each surface is pre-interpolated at the zone's TLE observation times, so the GA compares like-for-like. Propagator scratch buffers are zeroed before every grid run (v1.14 — SAVE'd buffers otherwise leak a previous zone's trajectory tail into the envelope).

### 5.3 BN Search Range Across Zones (v1.21)

```
Zone 1 : [bn_min_init, bn_max_init], with bn_lo extended downward when the
         G2 physics floor (zone 1's own TLE decay rate, calibrated against
         one propagate_ks trial) estimates the true BN below the caller floor.
Zone k : IF zone k−1 is TRUSTED (drag on: it predicted a re-entry;
              drag off: zone_status=0):
             re-center on bn_opt(k−1); ×0.5 width if the fit was interior,
             ×1.5 if it sat at a search boundary (true value likely outside)
         ELSE: range unchanged (weak-signal fits must not steer the search —
              37151's seven no-prediction zones once marched the range from
              [12.5,160] down to [17.2,24.9], imprisoning the only real zone)
```

### 5.4 GA Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Variables | 2 (e, BN) | |
| Population size | **20** (= maxpop) | Population 4 converged to a range-invariant seed artifact — the decoded chromosome was independent of the data (v1.15). Do not lower. |
| Generations | 200–500 | |
| Bit encoding | 40 bits (e) + 16 bits (BN) | |
| Crossover / mutation | 0.8 / 0.01 | |
| Fitness | 1/(1+RMS) of TWOINT-interpolated vs observed mean apogee over all zone observations | Trajectory matching (genpoen1.f heritage); slope-fitting variant kept unused in ga.F |

### 5.5 Estimator (v1.20)

**Primary = the latest zone's prediction.** Offline evaluation of five schemes on the 7-object campaign (`scratch_rpe/ensemble_eval.py`): latest-zone median |RPE| 8.2% / mean 7.6% / max 14.4% at 4 zones (2.4%/4.1%/10.4% at 8 zones, gated) vs uniform mean 45% mean error. The uniform ensemble mean ± std is retained as the agreement/spread indicator.

---

## 6. Validation Cases & Results (v1.21, 8-zone gated campaign)

| Object | NORAD | i (°) | e₀ | Known re-entry | Latest-zone RPE |
|---|---|---|---|---|---|
| PSLV-C39 R/B | 42928 | 19.2 | 0.33 | 2019-03-03 | **0.0%** |
| Ariane 5 ESC-A | 35497 | 5.7 | 0.63 | 2016-10-31 | **0.6%** |
| Proton-M Briz-M | 39615 | 48.5 | 0.68 | 2017-09-15 | **2.0%** |
| GSLV R/B | 32007 | 25.9 | 0.29 | 2010-06-06 | **2.4%** |
| Proton-M R/B | 37819 | 63.4 | 0.47 | 2013-09-12 | **−5.3%** |
| Long March 3B | 37151 | 24.9 | 0.56 | 2015-12-03 | **−8.1%** |
| Ariane 5 R/B | 27526 | 17.7 | 0.59 | 2012-05-09 | **10.4%** |

35497's early zones run +170..+520% (solar apsidal resonance at i=5.7° — the motivation for issue #9's 3-variable optimization); its final zone lands at −1.1% per-zone.

---

## 7. Development Status

Core algorithm **complete** (all closed): #1–#8 pipeline, #12 diagnostics/identifiability, #13 report, #16 E2E + accuracy target, #25 drag audit, KSROP #24. Open, priority order: **#27** (P1/critical — see below), **#26** (P1 — space weather implemented, §4.1, open pending 7-object regression), #14 (P2 — dynamic solar activity, largely superseded by #26's mechanism), #10/#11/#22 (P3 — TLE quality filtering, GMAT re-check pending lunisolar, CI), #9 (P4 — inclination as third variable), and the P4 operational backlog (#15, #17–#21, #23–#24). **354 tests across 10 suites** (added `test_sw.F`, 12 checks, v1.23).

### Open problem: critical-inclination decay (#27) — BN identifiability, not gravity fidelity

33587 (i=65°, hp 616→341 km over 131 days in-record) motivated #26, but a decisive GMAT cross-check (2026-07-17) found:

- **GMAT's full-force run** (exact point-mass Sun/Moon — not degree-truncated like OREM's qsun/qmoon — EGM2008 J2–J20, real 2022 CSSI flux, SRP) reproduces OREM's own flat hp(t) (616→~620 km), **not** the observed descent. A gravity-only companion run (drag/SRP off) gives nearly the same curve — third-body/oblateness truncation is not the gap.
- **Independently re-deriving hp(t) from the raw TLE mean elements** (bypassing OREM/GMAT entirely) shows semi-major axis collapsing 11128→9211 km with simultaneous circularization (e: 0.37→0.27) over the same window — the signature of atmospheric drag removing orbital energy, not a conservative lunisolar term (which cannot secularly change SMA).
- This contradicts the issue's original premise that drag is negligible there (that check only compared weather-on vs weather-off at a fixed, likely-wrong BN=55 kg/m² — a value inherited from an earlier weak-signal zone fit, not a validated physical estimate for this object).

**Working hypothesis**: the true ballistic coefficient during this arc is far lower (draggier) than 55 kg/m², and re-fitting BN against the in-record TLE window — rather than any gravity-model change — is the next step (not yet run; flagged in the issue, `#27` comment 2026-07-17). See `scratch_gmat/gmat_hp_33587_issue27_{full,grav}.script`.

---

## 8. Configuration File (`orem.cfg`)

```
input/example_42928.tle.txt          TLE file path
42928                                NORAD ID
2019 3 3 0 0 0.0                     Observed re-entry (0s = operational mode)
8                                    Max zones (8 recommended, v1.21)
8 10.0 0.90 -1.0                     Zone: min_pts, max_days, R², slope
80.0 160.0                           BN bounds kg/m² (G2 floor may extend down)
20 200 40 16 0.8 0.01 0.123          GA: pop, gen, bits_e, bits_BN, Pc, Pm, seed
2 0 0                                Force: geo_deg, sun_deg, moon_deg
0 7.2921150d-5 3.35281066d-3 1.0     Drag: IDRAG, WE, EPS_f, FR
0 0.0 0.0 0                          SRP: IPSR, CR, AM, ISHAD
```

Build/run/test commands: see README §4–§6 (all executables need `/heap-arrays /F:16777216`; `orem.exe` and `test_orem.exe` link `report.F`).

---

## 9. KSROP Linkage

The `ksrop/` directory contains files copied from the KSROP repo. When KSROP is updated:

1. Copy updated files: `cp $KSROP/{Subrouts.F,Legendre.F,TLEread.F} ksrop/`
2. If `driver_KS.F` changes, re-apply the refactoring to `propagate_ks.F` (the two carry the same physics; fixes flow in both directions — e.g. the v1.18 drag-phase fix was ported back as KSROP #24)
3. `input/ATM.DAT` (static) is generated by `KSROP/gen_atm_jr71.F`; `input/ATM2D.DAT` (epoch-resolved, v1.23) by `KSROP/gen_atm2d_jr71.F` — both share the 1-D profile code in `KSROP/jr71_profile.F`, so the two stay physically consistent by construction. Regenerate in KSROP and copy.
4. Run the full OREM suite to verify (10 executables — see README §4/§6 for the `test_sw.exe` build line, added v1.23)

The common block `/xy/` (pi, d2r, r2d, amue, AU, R_Earth) is the interface contract between KSROP files and OREM modules. `init_constants()` must be called before any KSROP subroutine.
