# OREM Algorithm Review — GA & RSM Gap Analysis

**Version:** OREM v1.7 (commit 96e143d, branch HS-dev)  
**Date:** 2026-07-05  
**Scope:** Gaps in OREM's RSM+GA pipeline compared to the original research (genpoen1.f / NPOE baseline)

---

## ⚑ Resolution status (added 2026-07-14, v1.21) — review CLOSED

All gaps below are resolved or dispositioned. OREM v1.21 **matches the original research's accuracy**: latest-zone RPE median 2.4% / mean 4.1% / max 10.4% on the 7-object campaign (vs the baseline's "within ±3%, none over 8%" on 4 objects).

| Gap | Resolution |
|---|---|
| G1 atmosphere | v1.8 J70 table → **superseded by v1.17**: the hand-rolled J70 generator's temperature profile was itself wrong (~3.3–3.5× too dense at 140–200 km); replaced by the real Jacchia-71 profile (`gen_atm_jr71.F`), validated 0.80–0.95× GMAT JacchiaRoberts |
| G2 BN range | v1.12 physics-based BN floor (`estimate_bn_floor`, zone 1, floor-only) |
| G3 narrowing | v1.10 widen-on-boundary + **v1.21 trust gate**: only zones that actually predicted a re-entry re-center the range |
| G4 zone distribution | v1.20/v1.21: zone_select's top-R² ranking distributes zones naturally once `nzones_max=8`; the **latest-zone prediction is the primary estimate** (drag signal concentrates as perigee decays) |
| G5 e-axis | Unchanged — inherent short-window limitation, shared with the original research |
| G6–G8 | No change needed (as concluded below) |
| §4 operational RPE | v1.19/v1.20 `report.F`: ensemble mean ± std, relative spread, and latest-zone primary reported unconditionally |

**The lesson this review could not see:** the v1.7 −93..−55% RPE was dominated not by algorithm-vs-heritage gaps but by four implementation defects found afterwards — stale SAVE'd RSM buffers leaking trajectories between zones (v1.14), a GA population of 4 whose output was a range-invariant seed artifact (v1.15), the G1 table's wrong temperature profile (v1.17), and a drag-density phase error along decay arcs (v1.18). Each was masked by the others; the heritage comparison only became meaningful after all four were fixed. The section below is preserved as the historical analysis.

---

## 1. Baseline — what the original research achieved

Predictions from DATES.OUT across four objects. RPE is signed: negative = predicted earlier than actual.

| Object | Zone | Zone epoch | Predicted re-entry | Actual re-entry | RPE |
|--------|------|------------|--------------------|-----------------|-----|
| 42928 PSLV-C39 | Z0 | 2017-09-24 | 2019-02-16 | 2019-03-03 | −2.9% |
| 42928 PSLV-C39 | Z2 | 2017-10-03 | 2019-02-22 | 2019-03-03 | −1.7% |
| 42928 PSLV-C39 | Z12 | 2018-01-20 | 2019-03-11 | 2019-03-03 | +1.7% |
| 42928 PSLV-C39 | Z13 | 2018-02-11 | 2019-02-03 | 2019-03-03 | −7.5% |
| 35497 Ariane 5 | Z0 | 2013-10-16 | 2016-11-13 | 2016-10-31 | +1.1% |
| 35497 Ariane 5 | Z2 | 2015-06-06 | 2016-11-19 | 2016-10-31 | +3.7% |
| 35497 Ariane 5 | Z3 | 2015-10-29 | 2016-11-01 | 2016-10-31 | ~0% |
| 35497 Ariane 5 | Z4 | 2016-05-19 | 2016-10-29 | 2016-10-31 | −1.8% |
| 37151 Long March 3B | Z2 | 2014-01-10 | 2015-11-21 | 2015-12-04 | −2.0% |
| 37151 Long March 3B | Z4 | 2015-03-04 | 2015-12-08 | 2015-12-04 | +1.1% |
| 37151 Long March 3B | Z5 | 2015-07-10 | 2015-12-12 | 2015-12-04 | +4.8% |
| 39615 Proton-M Briz-M | Z1 | 2015-07-16 | 2017-09-04 | 2017-09-15 | −1.5% |
| 39615 Proton-M Briz-M | Z2 | 2016-06-17 | 2017-09-13 | 2017-09-15 | −0.4% |
| 39615 Proton-M Briz-M | Z4 | 2017-03-02 | 2017-09-28 | 2017-09-15 | +6.6% |

Most predictions are within ±3%; no zone exceeds 8% RPE. This is the performance target for OREM.

OREM v1.7 on 42928 returns **−93% to −55% RPE** — a factor of 10–40× worse. The gaps below explain why.

---

## 2. Integration step constraint — confirmed

The KS integration step is fixed at **1° per KS regularized angle** (`istep=360` steps per revolution). This is hardcoded in both RSM generation (`orem.F:246`) and re-entry propagation (`orem.F:310`). It is a proven value for KS stability and must not be changed to adaptive stepping. It matches the original research DINP.DAT parameter of `istep=360`.

---

## 3. Gap analysis

### G1 — Atmospheric density model [RESOLVED — commit cf37576] (→ Issue #14)

**What the original used:** Jacchia-70 atmosphere with live F10.7 solar flux input.  
**What OREM previously used:** ATM.DAT — a static tabulated exponential model with ~50% of the correct thermospheric density.

**Previous divergence (`test_npoe` N13–N14):** OREM's `propagate_ks` produced roughly **50% of NPOE's apogee decay** over the same zone.

| BN (kg/m²) | OREM Δha (7 days) | NPOE ref Δha | Ratio (old) | Ratio (new J70) |
|-----------|-------------------|--------------|-------------|-----------------|
| 80 | −29.5 km → −59.1 km | −61.9 km | 48% | **95%** |
| 120 | −19.7 km → −39.4 km | −41.1 km | 48% | **96%** |
| 160 | −14.8 km → −29.6 km | −30.6 km | 48% | **97%** |

**Fix applied:** `input/ATM.DAT` regenerated using Jacchia-70 multi-species diffusive equilibrium (KSROP `gen_atm_j70.F90`, committed 81f470b). Conditions: F10.7=72, F10.7B=72, Kp=1.0, T_inf=640 K — matching the NPOE Zone-0 reference.

Key values at perigee altitude (170 km):
- Old: ρ = 8.378×10⁻¹⁰ kg/m³, H = 24.5 km
- New: ρ = 1.754×10⁻⁹ kg/m³, H = 18.4 km (correct J70 composition at low solar activity)

All 327 OREM tests pass with the new ATM.DAT (OREM commit cf37576).

**Remaining note:** ATM.DAT is static — it represents low solar activity (F10.7=72). For operational use at higher solar activity, density will be underestimated. A future enhancement would dynamically scale ATM.DAT based on the epoch F10.7 value (→ Issue #14).

---

### G2 — BN initial search range is hardcoded [HIGH] (→ Issue #12, #15)

OREM uses `bn_min_init=80, bn_max_init=160` from config for all objects. The original research shows this range is object-specific:

| Object | Zone | Original BN range | OREM hardcoded | Status |
|--------|------|--------------------|----------------|--------|
| 42928 PSLV-C39 | Z0 | [80, 160] | [80, 160] | matches |
| 35497 Ariane 5 upper | Z0 | [80, 120] | [80, 160] | range too wide |
| 39615 Proton-M Briz-M | Z1 | [30, 110] | [80, 160] | lower half missed entirely |

For 39615, the true BN is ≈30–50 kg/m² and the GA can never find it with floor=80. The original BN ranges were hand-chosen from knowledge of each object's mass/area ratio. A physics-based initial estimate from TLE mean-motion decay rate (observed Δn̄/Δt) could replace the hardcoded floor.

---

### G3 — BN narrowing is strictly sequential — original is not [MEDIUM] (→ Issue #12)

OREM narrows BN range by 50% per zone, centred on the previous zone's BN. The original research shows a more complex pattern for 42928:

| Zone | Original BN range | Pattern |
|------|------------------|---------|
| Z0 | [80, 160] range=80 | Initial |
| Z2 | [100, 140] range=40 | Narrowed from Z0 |
| Z12 | [90, 130] range=40 | Further narrowed |
| Z13 | [80, 160] range=80 | **Widened back** — validation run |

Key observation: Zone-13 reverts to the full [80, 160] range — confirming the algorithm is not a strict sequential pipeline. The researcher re-opened the range for cross-validation. Also, OREM's 50% narrowing compounds rapidly: by Zone 4 the range may be only 10 km/m² wide, which excludes the true BN if any early zone gave a biased estimate. A minimum range floor (e.g., 20 kg/m²) should be added to `orem.F:333–336`.

---

### G4 — Zone selection strategy: consecutive vs. distributed [MEDIUM] (→ Issue #12)

OREM selects 4 consecutive zones from recent TLE history. The original research selected zones at well-separated epochs spanning the full decay history:

- **42928:** Z0 (day 1), Z2 (day 9), Z12 (day 118), Z13 (day 140) — span 142 of 524 days
- **35497:** Z0 (day 1), Z2 (day 599), Z3 (day 743), Z4 (day 946) — span 1027 of 1111 days

Four consecutive zones drawn from one narrow epoch window provide correlated BN estimates — the ensemble mean is less robust than estimates from distributed epochs. Distributing zone selection across the observable TLE history (e.g., at 0%, 25%, 50%, 75% of predicted remaining lifetime) would more closely replicate the original intent and improve ensemble independence.

---

### G5 — Eccentricity axis in RSM: range too narrow for independent information [LOW]

The e-grid spans `e_mid ± de`, where `de` is the eccentricity variation within the zone window, typically ±0.0003–0.0007. For a=9,000 km this is Δha ≈ 9 km — within TLE noise (~5 km). The three e-grid RSM surfaces are nearly indistinguishable. The GA's e-optimisation effectively fits TLE noise. The same narrow de appears in the original research — this is an inherent short-window RSM limitation. The BN axis is what carries the physical information.

---

### G6 — imax ≠ itmax in original: RSM points vs. fitness TLE count [LOW]

In the original research for 42928 Z12: `imax=12` (RSM time-steps) but `itmax=9` (TLE observations in GA fitness). OREM uses `nsurf_pts = nobs = nzone` (all zone TLEs for both RSM and fitness). Using all TLEs is the correct general approach. The original's split arose from manually excluding outlier TLEs. TLE quality filtering (Issue #10) achieves the same effect automatically. No structural change needed.

---

### G7 — TWOINT extrapolation at GA boundary candidates [LOW]

When the GA proposes (e, BN) outside the 3×3 grid, `ga_twoint` extrapolates using the boundary interval. This happens frequently in early zones when BN converges near the search boundary. The extrapolated fitness value is a linear extension of the surface gradient, which may be unreliable >5% beyond the grid edge. This is the same behaviour as the original TWOINT — not a correctness bug, but understood as a limitation.

---

### G8 — INTPOL time-axis in original is uninitialized [NOTE: OREM is correct]

In genpoen1.f, the `tnpoe` array (time axis for INTPOL when imax≠itmax) is declared but never populated. OREM avoids this entirely by storing RSM surfaces directly at TLE observation times via `rsm_tinterp`. OREM's approach is architecturally correct. No change needed.

---

## 4. Operational RPE — independent of known re-entry date

The current RPE is a validation metric requiring the observed re-entry date `t_obs`. For future predictions, `t_obs` is undefined. Two modes already exist in `compute_rpe`:

| Mode | Trigger | Reference | Use case |
|------|---------|-----------|----------|
| Mode 1 | `t_obs > 0` | Known re-entry JD | Historical validation only |
| Mode 2 | `t_obs = 0` | Ensemble mean `t_mean` | Operational — no known re-entry |

The missing piece is that the ensemble **spread** is not reported as the primary uncertainty output. The operational metric should be:

```
relative_spread = t_std / (t_mean - t_now) × 100%
```

This is dimensionless and independent of how far away the predicted re-entry is. It quantifies whether the zone ensemble has converged:

- **spread < 5%:** ensemble converged; BN well-constrained across zones
- **spread 5–15%:** moderate uncertainty; more zones or wider observable history needed  
- **spread > 15%:** BN poorly constrained; prediction is indicative only

No structural change to `compute_rpe` is required. The reporting in `orem.F` should additionally print `t_mean`, `t_std`, and `relative_spread` unconditionally, regardless of whether `t_obs` is supplied (→ Issue #13).

---

## 5. Recommended actions — all complete (status added 2026-07-14)

| Priority | Action | Issues | Outcome |
|----------|--------|--------|---------|
| **Critical** | ~~Validate and correct ATM.DAT~~ | #14 | v1.8 J70 table, then **superseded by v1.17's real Jacchia-71 profile** after a GMAT density probe showed the J70 generator was 3.3–3.5× too dense in the perigee band |
| **High** | ~~Physics-based BN floor from TLE decay rate~~ | #12 | **DONE v1.12** (`estimate_bn_floor`, numerically calibrated against propagate_ks itself, floor-only) |
| **Medium** | ~~Fix narrowing over-convergence~~ | #12 | **DONE v1.10 + v1.21**: widen-on-boundary, then the trust gate (only re-entry-predicting zones re-center the range) — a minimum-range floor proved unnecessary once untrusted zones stopped steering |
| **Medium** | ~~Ensemble spread reporting~~ | #13 | **DONE v1.19/v1.20** (`report.F`: mean ± std, relative spread, latest-zone PRIMARY estimate) |
| **Later** | ~~Distributed zone selection~~ | #12 | **Dispositioned v1.21**: `nzones_max=8` + top-R² ranking distributes zones across the decay; the latest-zone estimator makes explicit lifetime-percentage placement unnecessary |
| **Later** | ~~Recalibrate E-series after ATM.DAT fix~~ | — | **DONE v1.17/v1.18** (assertions moved to physical-sanity + first-principles drag reference) |
