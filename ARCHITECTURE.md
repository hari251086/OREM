# OREM — Application & Physical Architecture

*(current as of v1.45, 2026-07-27; the Version History in README.md is the authoritative changelog)*

## 1. System Overview

**OREM** (Optimal Regularized re-Entry estimation Method) predicts the atmospheric re-entry time of resident space objects (RSO) decaying from highly elliptical orbits (HEO). It combines a regularized orbit propagator (KSROP) with response surface methodology (RSM) and genetic algorithm (GA) optimization to compensate for the low accuracy of Two-Line Element (TLE) catalog data.

### Problem Statement

RSO in HEO (GTO, Molniya, SSTO) experience complex orbital evolution under luni-solar gravity, oblateness, and atmospheric drag. Their re-entry times are highly sensitive to:
- **Initial conditions** — TLE accuracy is limited (~km-level for HEO)
- **Ballistic number** BN = m/(Cd·A) — unknown tumbling state, cross-sectional area uncertainty, and attitude regime that changes over the object's life
- **Eccentricity** — small errors amplify through luni-solar resonance dynamics

OREM treats re-entry prediction as an optimization problem: per decay zone, find the (eccentricity, ballistic number) pair that best fits the observed TLE apogee evolution, then propagate each zone's fit forward to re-entry. The **latest zone's prediction is the primary estimate** — it carries the shortest extrapolation and the freshest attitude/altitude regime.

### Achieved Accuracy (current, v1.45) — and an honest note on where the numbers stand

Across the curated-7 validation set, **static atmosphere**: mean\|latest-zone RPE\| **32.1%**, mean\|ensemble RPE\| **4.1%** (median 2.8%) — see the real per-object table in §6. **With epoch-resolved weather + the diurnal bulge** (§4.1): mean\|latest-zone RPE\| **26.1%**, mean\|ensemble RPE\| **4.5%** (median 3.4%). Across the full 50-object generalization set (§7): mean\|ensemble RPE\| **16.1%**, predict rate **35/50 (70%)**.

These are materially different from — and in the latest-zone metric's case, considerably worse than — the "2.4% median / 4.1% mean" figure this document quoted through v1.26. That older number came from an 8-zone-gated variant of the campaign that hasn't been the shipped default in a long time; the honest current picture is in §6, not a historical best case. **A genuinely useful, current finding**: under today's median-based ensemble metric (§5.5, changed from a mean in v1.44), the ensemble/spread indicator is now *consistently more accurate than the "primary" latest-zone estimator* on this curated-7 set — the opposite of the v1.20 finding that made latest-zone primary in the first place. See §5.5 for why, and treat "primary" as the officially-reported estimate, not necessarily the most accurate one for every object.

The accuracy story since v1.26 has been one of closing generalization gaps rather than one clean number going up: a `zone_select` recency fix (v1.41), a literature-grounded diurnal density bulge (v1.42, §4.1), a widened `max_zone_days` default (10→20 days, v1.43, raising the 50-object predict rate from 58%→72% before the Falcon-9 swap below changed the roster again), a BSTAR-informed per-zone BN prior (**G3**, §5.1/5.3), and the mean→median ensemble-metric switch (v1.44) all measurably helped *something* while each also interacted unpredictably with the BN-range carryover mechanism across zones (§5.3) — four independent, confirmed instances of the same fragility, investigated but not yet structurally resolved (issue #35; see §5.3). Most recently (v1.45), three Falcon-9 R/B objects were removed from the 50-object set and replaced with non-maneuverable substitutes, since `propagate_ks` has no thrust/maneuver modeling and SpaceX is documented to perform active post-separation deorbit burns on Falcon 9 second stages.

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
│    1. G3 BN prior (issue #32): narrow this zone's [bn_lo,bn_hi]     │
│       using its own TLE-published BSTAR (pooled log-log regression │
│       vs. fitted BN, ~4x window) — stateless, per-zone, intersect  │
│       only; falls back cleanly if BSTAR is unavailable             │
│    2. rsm_generate(): 9× propagate_ks over the zone span            │
│       (3 ecc × 3 BN grid → 9 mean-apogee surfaces, pre-interpolated │
│        at the TLE observation times; scratch buffers zeroed per     │
│        call — v1.14)                                                │
│    3. Diagnostics: RSM envelope bounds observations? (zone_status)  │
│    4. ga_optimize(): pop=20 GA on the surfaces (TWOINT bilinear     │
│       interpolation, NO propagation) → e_opt, bn_opt                │
│       + GA-boundary detection (15% of range → zone_status=2)        │
│    5. propagate_ks(e_opt, bn_opt): long propagation → reentry_jd    │
│    6. Trust-gated BN carryover (v1.21, anchor fix v1.44): re-center │
│       the next zone's BN range on this zone's fit ONLY if this      │
│       zone actually predicted a re-entry (drag on) / is unflagged   │
│       (drag off); widen ×1.5 around the window's EXISTING center    │
│       if the fit sat at a search boundary, else ×0.5 re-centered    │
│       on the fit itself. Recursive and path-dependent by            │
│       construction — confirmed (issue #35) to amplify small,        │
│       individually-correct upstream changes unpredictably 4         │
│       separate times; a non-recursive pooled-prior replacement was  │
│       designed, implemented, and campaign-tested but gave genuinely │
│       mixed results and was NOT shipped (see §5.3) — this step is   │
│       still the recursive v1.21/v1.44 logic as of this document     │
│                                                                     │
│  compute_rpe(): per-zone RPE + ensemble t_mean/t_std (median, v1.44)│
└───────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────┐
│  REPORT (report.F: orem_report — called by main_orem.F)             │
│  output/OREM_<NORAD>_<DATE>.txt:                                    │
│    config echo · per-zone table (epoch, e, BN, re-entry, RPE,       │
│    status) · PRIMARY estimate (latest zone) + its RPE ·             │
│    ensemble MEDIAN ± std + relative spread (median, not mean,       │
│    since v1.44 — issue #29) · status legend                         │
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
| **KSROP-lineage Propagator** | `src/propagate_ks.F` | KS regular-elements propagation (refactored from KSROP's `driver_KS.F`; carries the identical drag/gravity physics as a callable subroutine) | (called by G2/RSM/re-entry) |
| **OREM Driver** | `orem.F` | Zone loop, G2 floor, G3 BSTAR prior, diagnostics, trust-gated carryover, `compute_rpe` (median, v1.44) | **Yes — G2 trial + 1× per zone** |
| **Report** | `report.F` | Prediction report with latest-zone primary estimate | No |
| **Runner** | `main_orem.F` | Reads `orem.cfg`, runs pipeline, writes report | (via orem_run) |

---

## 3. Data Flow

```
TLE File → tle_evolve() → epochs, a, e, i, Ω, ω, ha, hp, Λs  (deduplicated)
         → zone_select() → up to nzones_max windows with clean linear
                           apogee decay (R² ≥ 0.90 over ≤ 20 days,
                           max_zone_days 10→20 since v1.43; a Pass-1.5
                           recency guarantee, v1.41, ensures the most
                           recent candidate window always gets a slot),
                           ranked by R², sorted by epoch
         → estimate_bn_floor() → may extend zone 1's bn_lo downward
                           (floor-only; never raises; G2, v1.12)

FOR EACH ZONE:
  estimate_bn_bstar_prior() → narrows [bn_lo,bn_hi] to [bn_lo_zone,
                           bn_hi_zone] using this zone's own TLE BSTAR
                           (G3, issue #32; stateless, intersect-only)
  3×3 grid:  e-axis = e_mid ± δe (TLE scatter in the zone)
             BN-axis = [bn_lo_zone, mid, bn_hi_zone] (kg/m²)
  rsm_generate() → surfaces(i, ie, ibn) = mean apogee at obs time i
  zone_status: 1=propagator failure (>6 of 9 grid points diverged)
               3=RSM envelope fails to bound an observation
  ga_optimize() → (e_opt, bn_opt); zone_status=2 if within 15% of a bound
  propagate_ks(e_opt, bn_opt) → reentry_jd or 0
  trust gate → BN range for the next zone (recursive, path-dependent;
               confirmed fragile 4x, issue #35 — see §5.3)

compute_rpe():
  Mode 1 (validation): RPE(iz) = (t_pred(iz) − t_obs)/(t_obs − t_zone(iz)) × 100%
  Mode 2 (operational): t_mean/t_std over predicting zones; RPE vs t_mean
  (t_mean is the MEDIAN of valid per-zone predictions since v1.44, issue
   #29 — parameter name kept for call-site compatibility; a single
   catastrophically-wrong zone can drag an arithmetic mean far off, a
   median is far more robust to exactly that failure mode)
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
| Diurnal density bulge | ρ = ρ₀(1+F·cos φ)·exp{...} (Sharma 1997a / Swinerd & Boulton 1983), φ from satellite/Sun geocentric angle, bulge lag λ=37° | v1.42, issue #32 Phase 5 follow-up. Auto-enables only when epoch-resolved weather is already loaded (extends §4.1, zero-cost otherwise); literature-cited constants, none guessed |
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
              propagate_ks (src/propagate_ks.F, per-rev hook, ~line 288):
                sw_tinf(JD) → T∞ for this revolution's epoch
                atm2d_interp(h_perigee, T∞) → ρ_p, H   [istat=0]
                falls back to the legacy static ALT/DEN/SCH table
                lookup whenever nothing is loaded, T∞≤0, or the point
                falls outside the 2-D table [istat≠0] — computed ONCE
                per revolution, same slot the static path always filled
```

**Opt-in, zero-cost when unused**: `sw_tinf`/`atm2d_interp` live inside `propagate_ks.F` itself (not `swx.F`) so that executables which never call `sw_load`/`atm2d_load` link without `swx.F` and run the legacy single-table path bit-unchanged — every pre-v1.23 build, result, and test is unaffected. `main_orem.F` auto-detects both files at startup and prints ENABLED/DISABLED loudly; only opt-in executables (`orem.exe`, `test_sw.exe`) link `swx.F`.

**Status**: shipped and closed (issue #26, `test_sw.F` now 18 checks after the v1.42 diurnal-bulge addition — G5-storm T∞=1216.56 K hand-verified against the Jacchia formula; W12 smoke test shows storm-epoch decay 55.4 km vs quiet-epoch 31.1 km over 7 days, correct direction and rough magnitude; W13–W18 hand-verify the bulge extrema formula and geometry). Full curated-7 regression against the static baseline found weather+bulge a real but non-uniform effect (§1) — kept, since the physics is directionally validated and literature-grounded even where it interacts with the BN-carryover chain's own fragility (§5.3).

---

## 5. Optimization Architecture

### 5.1 Design Variables

| Variable | Symbol | Range | Source of uncertainty |
|---|---|---|---|
| Eccentricity | e | e_mid ± δe (TLE scatter in the zone) | SGP4/SDP4 reconstruction error |
| Ballistic number | BN (kg/m²) | zone 1: [bn_min, bn_max] with G2 floor; every zone: narrowed by G3's BSTAR prior; later zones' range also carried over from zone k−1 (recursive, §5.3) | Tumbling state, attitude regime drift |

BN is optimized directly (mass-as-variable convention: Cd=1, A=1 m², DryMass=BN), as in the original NPOE-era research. Fitted values on the validation set fall in ~20–100 kg/m² per zone.

### 5.2 RSM 3×3 Grid

9 short propagations per zone map (e, BN) → mean-apogee history; each surface is pre-interpolated at the zone's TLE observation times, so the GA compares like-for-like. Propagator scratch buffers are zeroed before every grid run (v1.14 — SAVE'd buffers otherwise leak a previous zone's trajectory tail into the envelope).

### 5.3 BN Search Range Across Zones (v1.21, G3 added v1.32-era, anchor fix v1.44)

Two distinct mechanisms narrow each zone's BN search range, one stateless and one recursive:

**G3 (stateless, per-zone)** — `estimate_bn_bstar_prior`: narrows whatever range is in effect for this zone using that zone's own TLE-published BSTAR value, via a pooled log10(BN) = slope·log10(BSTAR) + intercept regression (median r≈−0.7 across the validation set). Intersection-only, never widens, falls back cleanly if BSTAR is unavailable. Carries nothing between zones — every call is independent.

**Trust-gated carryover (recursive, across zones)**:

```
Zone 1 : [bn_min_init, bn_max_init], with bn_lo extended downward when the
         G2 physics floor (zone 1's own TLE decay rate, calibrated against
         one propagate_ks trial) estimates the true BN below the caller floor.
Zone k : IF zone k−1 is TRUSTED (drag on: it predicted a re-entry;
              drag off: zone_status=0):
             re-center on bn_opt(k−1); ×0.5 width if the fit was interior;
             if it sat at a search boundary (true value likely outside),
             widen ×1.5 around the window's EXISTING center (v1.44 — the
             original v1.21 logic re-centered on the boundary-pinned value
             itself, which let a run of boundary hits escalate the window
             in one direction instead of escaping it; see the 61734 trace
             in README v1.43/v1.44)
         ELSE: range unchanged (weak-signal fits must not steer the search —
              37151's seven no-prediction zones once marched the range from
              [12.5,160] down to [17.2,24.9], imprisoning the only real zone)
```

**This recursive mechanism is a confirmed, recurring source of fragility, not a settled design.** Because each zone's window depends only on the immediately preceding zone's single-point result, a small shift in what BN any one zone converges to — from new physics, different zone boundaries, or different input data — gets amplified and redirected unpredictably down the chain. Four independent, individually-correct changes this project has shipped have each triggered a real regression through this exact mechanism: the `zone_select` recency fix (v1.41), the diurnal density bulge (v1.42), `max_zone_days` 10→20 (v1.43, 61734's BN escalated 16.3→58.6→89.4→136.0), and a re-test of the already-shipped noise-matched-apobs fix (issue #34) that found it now regresses under the current stack when it was originally kept for improving it.

**Issue #35 investigated a structural fix directly** — full history in project records, not repeated here in detail:
- *Tested removing the recursive carryover entirely*, on the hypothesis that G3 plus the mean→median ensemble-metric switch (§5.5) might already substitute for it, since removal was last tried (issue #33) before either existed. **Result: clear loss** on the primary metric across all three validation campaigns. Not shipped.
- *Designed and fully implemented a non-recursive replacement* (`estimate_bn_pool_prior`: median ± 3·MAD pooled from all prior trusted, non-boundary-pinned zones, recomputed fresh each zone rather than mutated in place — eliminates the path-dependency by construction). Campaign-tested at three pool sizes (unbounded, last-4, last-1 zones). **Result: genuinely mixed** — improves the primary metric on 2 of 3 campaigns and fixes 37151's chronic regression, but consistently regresses the ensemble metric on one already-known-difficult object (35497) regardless of pool size. Not shipped; production code reverted to the recursive logic above pending an explicit decision. The full implementation survives as a ready-to-reapply patch (`scratch_rpe/issue35_phase2_pool_prior.patch`).

**What this means in practice**: the mechanism shown in the code block above is what actually runs today. Any future change anywhere in the pipeline that could shift a zone's fitted BN even slightly should be re-validated against the full campaign set, not assumed safe from this interaction just because it's individually correct.

### 5.4 GA Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Variables | 2 (e, BN) | |
| Population size | **20** (= maxpop) | Population 4 converged to a range-invariant seed artifact — the decoded chromosome was independent of the data (v1.15). Do not lower. |
| Generations | 200–500 | |
| Bit encoding | 40 bits (e) + 16 bits (BN) | |
| Crossover / mutation | 0.8 / 0.01 | |
| Fitness | 1/(1+RMS) of TWOINT-interpolated vs observed mean apogee over all zone observations | Trajectory matching (genpoen1.f heritage); slope-fitting variant kept unused in ga.F |

### 5.5 Estimator (primary v1.20; ensemble statistic changed mean→median v1.44)

**Primary (officially reported) = the latest zone's prediction.** Originally chosen because an offline evaluation of five schemes on the 7-object campaign (`scratch_rpe/ensemble_eval.py`) found latest-zone median |RPE| 8.2% / mean 7.6% / max 14.4% at 4 zones, decisively ahead of a uniform-mean ensemble's 45% mean error — the latest zone has the shortest extrapolation and the freshest attitude/altitude regime. `report.F` still leads with this as "PRIMARY estimate."

**The ensemble/spread indicator switched from an arithmetic mean to a median in v1.44** (issue #29), specifically because one catastrophically-wrong zone (typically an early zone whose window sits too far from the true re-entry to extrapolate linearly, e.g. 61734's zone 1 at ~35,000 km apogee, 138 days before the real event) can drag a mean far off in a way a median resists.

**Honest current caveat**: on today's curated-7 static campaign, the median-based ensemble metric is now *more* accurate than the "primary" latest-zone estimator for most objects (§1, §6) — the opposite of the relationship that made latest-zone primary in v1.20. This hasn't been re-evaluated as a possible change to which estimate is reported as PRIMARY; it's flagged here as a real, current data point rather than acted on.

---

## 6. Validation Cases & Results (current, v1.45, `scratch_rpe/rpe_campaign_7obj_issue31.csv`, static atmosphere, 8-zone)

| Object | NORAD | i (°) | e₀ | Known re-entry | Latest-zone RPE | Ensemble (median) RPE |
|---|---|---|---|---|---|---|
| PSLV-C39 R/B | 42928 | 19.2 | 0.33 | 2019-03-03 | −41.5% | **3.3%** |
| Ariane 5 ESC-A | 35497 | 5.7 | 0.63 | 2016-10-31 | −76.2% | **2.8%** |
| Proton-M Briz-M | 39615 | 48.5 | 0.68 | 2017-09-15 | 20.0% | **0.2%** |
| GSLV R/B | 32007 | 25.9 | 0.29 | 2010-06-06 | 5.6% | 12.2% |
| Proton-M R/B | 37819 | 63.4 | 0.47 | 2013-09-12 | 35.7% | **−8.1%** |
| Long March 3B | 37151 | 24.9 | 0.56 | 2015-12-03 | −42.3% | **−1.8%** |
| Ariane 5 R/B | 27526 | 17.7 | 0.59 | 2012-05-09 | 3.3% | **0.1%** |
| **mean\|·\|** | | | | | **32.1%** | **4.1%** |
| **median\|·\|** | | | | | **35.7%** | **2.8%** |

The ensemble (median) column beats the "primary" latest-zone column for 6 of 7 objects here (§5.5) — a genuine, current reversal from the relationship that made latest-zone primary, worth knowing when reading any single object's report rather than trusting the label alone. With epoch-resolved weather + the diurnal bulge (`scratch_rpe/rpe_campaign_weather.csv`): mean\|latest\|=26.1%, mean\|ensemble\|=4.5%, median\|ensemble\|=3.4% — a similar pattern, non-uniform per object (§1).

35497 is a persistent outlier on both metrics (i=5.7° solar-apsidal-resonance case, the motivation for issue #9's still-open 3-variable optimization) and is also the one object whose ensemble metric regressed under every pool-size variant tested in the issue #35 carryover-chain investigation (§5.3) — worth investigating on its own terms if that thread is picked up again, since pool size alone doesn't move it.

---

## 7. Development Status

Core algorithm **complete** (all closed): #1–#8 pipeline, #12 diagnostics/identifiability, #13 report, #16 E2E + accuracy target, #25 drag audit, #26 epoch-resolved space weather (§4.1), #27 critical-inclination decay (root-caused as a BN-identifiability issue, not a gravity-model gap — see git/issue history for the GMAT cross-check that resolved it), KSROP #24.

**Open, priority order**: **#32** (P1, critical — the umbrella tracking issue for the whole ongoing RPE-accuracy investigation; most items below are its follow-ons); **#29** (P2 — the shipped `nzones_max=8`/50-object generalization gap; substantially improved via `max_zone_days` 10→20, a `zone_select` recency fix, the median ensemble metric, and the Falcon-9 roster swap, but the underlying BN-carryover fragility it kept surfacing is still open as #35); **#35** (P2 — BN-carryover chain: a non-recursive replacement was designed, implemented, and tested but gave mixed results and was not shipped, §5.3); #14 (P3 — atmospheric density model validation); #34 (P3 — one already-shipped fix (issue #31) was found to regress under the current stack, flagged, not reverted); #9 (P4 — inclination as a third design variable); #28 (P4 — `propagate_ks.F` doesn't compile under gfortran); and the P4 operational backlog (#15, #17, #24). **382 tests across 12 suites** as of the last shipped change (v1.44's median-metric fix); the reverted issue #35 Phase 2 work added 10 more that are not currently part of the shipped suite (patch preserved, §5.3).

---

## 8. Configuration File (`orem.cfg`)

```
input/example_42928.tle.txt          TLE file path
42928                                NORAD ID
2019 3 3 0 0 0.0                     Observed re-entry (0s = operational mode)
8                                    Max zones (8 recommended, v1.21)
8 20.0 0.90 -1.0                     Zone: min_pts, max_days, R², slope
                                     (max_days 10->20 since v1.43, issue #29)
80.0 160.0                           BN bounds kg/m² (G2 floor may extend down)
20 200 40 16 0.8 0.01 0.123          GA: pop, gen, bits_e, bits_BN, Pc, Pm, seed
2 0 0                                Force: geo_deg, sun_deg, moon_deg
0 7.2921150d-5 3.35281066d-3 1.0     Drag: IDRAG, WE, EPS_f, FR
0 0.0 0.0 0                          SRP: IPSR, CR, AM, ISHAD
```

Build/run/test commands: see README §4–§6 (all executables need `/heap-arrays /F:16777216`; `orem.exe` and `test_orem.exe` link `report.F`).

---

## 9. KSROP Linkage

Both OREM and KSROP were restructured into `fpm` (Fortran Package Manager) packages (`src/`/`app/`/`test/` layout). **OREM now consumes KSROP as a real dependency via `fpm.toml` (git + tag), not hand-copied files** — `fpm build`/`fpm test` fetches the pinned KSROP tag into `build/dependencies/ksrop/` automatically; there is no `ksrop/` directory in OREM's own source tree to keep manually in sync. When KSROP changes:

1. KSROP ships a new tagged release on its own `main`.
2. Bump the tag in OREM's `fpm.toml` KSROP dependency entry.
3. `fpm build`/`fpm test --compiler ifx --flag "/heap-arrays /F:16777216"` re-fetches and rebuilds against the new tag automatically.
4. `input/ATM.DAT` (static) and `input/ATM2D.DAT` (epoch-resolved) are still generated in KSROP (`gen_atm_jr71.F`/`gen_atm2d_jr71.F`, sharing `jr71_profile.F` so the two stay physically consistent by construction) and copied into OREM's `input/` — this one step remains manual, since these are data files, not source.
5. Re-run the full OREM suite to verify (`fpm test`, 12 test executables — see README for the exact build/run commands).

Bugs whose root cause is in KSROP-lineage code (e.g. `Subrouts.F`, `TLEread.F`, `Legendre.F`, or physics shared with `driver_KS.F`) get filed and fixed in KSROP first, regardless of which repo surfaced them — OREM's own issue only tracks "pull in KSROP vX.Y.Z." `src/propagate_ks.F` itself is OREM's own file (refactored from KSROP's `driver_KS.F` into a callable subroutine, no file I/O) and is edited directly in OREM; fixes to its shared physics are ported back to `driver_KS.F` by hand (e.g. the v1.18 drag-phase fix, the v1.42 diurnal bulge), since the two are no longer literally the same file.

The common block `/xy/` (pi, d2r, r2d, amue, AU, R_Earth) is the interface contract between KSROP files and OREM modules. `init_constants()` must be called before any KSROP subroutine.
