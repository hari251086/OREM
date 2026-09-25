# Zone-selection methodology in re-entry prediction: literature review vs. OREM

**Scope note (ambiguous inputs resolved per task instructions):**
- `<PATH_TO_OREM_METHOD_DOC_OR_CODE>` → resolved to `ALGORITHM.md` §3-4
  (OREM's own methodology document, git-tracked per repo convention) plus
  `src/zone_select.F` (the actual implementation) for exact numeric
  thresholds not stated in prose.
- `<OUTPUT_PATH>` → resolved to `scratch_rpe/` (this repo's established
  location for investigation/campaign write-ups, e.g.
  `scratch_rpe/resonance_plots/`, `scratch_rpe/rpe_campaign.F` outputs).

---

## 0. Source document

- **Title**: *Space Debris and Space Situational Awareness Research
  Studies in ISRO*
- **Authors**: Ram Krishan Sharma, A. K. Anil Kumar, Adimurthy Vipparthi
- **Year**: 2024 ("Monograph-2024")
- **Venue/publisher**: Indian Space Research Organisation (ISRO)
- **File**: `E:\Research\References\00_Author_Collections\1714143725793.pdf`
  (10.06 MB, text-extractable, not scanned — `pdftotext -layout` succeeded
  cleanly, 8430 lines of text)
- **Relevant chapter**: Chapter 8, "Re-entry Prediction Studies," pp.
  141–192, 27 subsections (8.1–8.27)
- **Note**: this is the same monograph previously deep-read in an earlier
  OREM session (via a differently-named duplicate PDF,
  `1714143800970.pdf` — same title/content, confirmed by identical
  chapter/page structure). That earlier pass extracted accuracy figures
  and the algorithm's citation lineage; this review is a fresh,
  zone-selection-methodology-specific pass with page-level citations,
  not a re-summary of the earlier one.

## 1. OREM's own zone-selection methodology (reference description)

From `ALGORITHM.md` §3–4 and `src/zone_select.F`:

- **Definition**: a "zone" is a contiguous window of TLE-derived mean
  apogee-altitude points showing a statistically significant *linear
  decay trend* — not a fixed time span or fixed TLE count by design.
- **Selection rule** (`zone_select.F`, Pass 1+2): scan the (quality-
  filtered) apogee-altitude time series; a candidate zone must have
  **R² ≥ `r2_thresh` (default 0.90)**, slope steeper than `slope_thresh`
  (default −1.0 km/day), **at least `min_zone_pts` points (default 8)**,
  spanning **at most `max_zone_days` (default 20, raised from 10 in
  v1.43 — issue #29)**. Up to `nzones_max` candidates are kept, preferring
  later/more-recent-to-decay zones when the cap binds; zones are
  non-overlapping by construction (Pass 2's candidate-removal step) and
  returned sorted by epoch.
- **Pre-filtering**: `tle_filter` drops points >4σ off a local 20-TLE
  windowed trend (30-day max gap) before zone selection ever runs — an
  explicit, quantified outlier-handling step (tuned against 3 real
  objects, ~1.5% false-positive rate).
- **Event dependence**: none, by default. Zones are found purely from
  the TLE apogee-decay trend itself; they are not anchored to a distinct
  physical event (e.g. a resonance epoch) unless a caller explicitly
  restricts the input TLE file to a sub-range (as done in this repo's
  own `scratch_rpe/rpe_campaign_resonance16.F` / the validation-RPE
  campaign this session, both one-off investigative uses, not the
  production default).
- **Propagator**: KSROP's Kustaanheimo-Stiefel (KS) regularized
  numerical propagator (`propagate_ks`, RKG4 integrator), fitted via a
  3×3 response-surface grid in (eccentricity, ballistic number) + a
  genetic algorithm (`rsm.F`/`ga.F`).
- **Re-entry threshold**: perigee altitude < 80 km.
- **Accuracy metric**: RPE = (predicted − actual) / (actual − fit-window
  epoch) × 100%, both per-zone and as an ensemble across zones.

## 2. Part 1 — Deep read: zone/window methodology per study (Chapter 8)

Every re-entry-prediction study Chapter 8 discusses in enough detail to
extract a fitting-window methodology, in the order presented. Page
numbers are the monograph's own printed page numbers (footer). Items
with no described window/zone method (pure background, campaign
narrative, or a different problem class) are listed as "not applicable."

| § | Study (page) | Window/zone definition | Length / count / TLE range | Event dependence | TLE filtering/outliers | Propagator | Accuracy | Test objects |
|---|---|---|---|---|---|---|---|---|
| 8.3 (142–144) | Anil Kumar et al. 2003c/2007, Kalman filter | Not zone-based — sequential recursive filter over the whole TLE stream | Not reported (whole available history, constant gains) | None | Not reported | Simple drag-only LEO propagator, USSA-1976 | Qualitative "close to actual" (3 objects) | 25947, SROSS-C2, Cosmos 1043 |
| 8.3.4 (147–148) | Bandyopadhyay et al. 2001, GSLV-D1 | Not zone-based — full continuous numerical propagation | Whole tracked history, monthly-averaged F10.7 | None | Not reported | NPOE (numerical, J6,6 + luni-solar + MSIS90) | Predicted ~600d vs actual 639d | GSLV-D1 upper stage |
| 8.3.5 (148) | **Mutyalarao & Sharma 2010** | Single fixed-length fitting window before re-entry | **165 days before re-entry epoch** (single window, not multi-zone) | None stated | Not reported | RSM+GA (propagator not specified in text) | Accurate 7 days before actual re-entry | GSLV-F01/CS |
| 8.3.5 (149–151) | **Mutyalarao & Sharma 2011** | **Explicit "zones"**, apogee-altitude-based | **8 zones (A–H)** over **200 days** of orbit evolution, "rough linear variation of mean apogee" | None (fixed calendar span) | Not reported | RSM+GA | <7% abs. error (6 decayed objects, Tables 8.1–8.2) | GSLV-F04/CS (#32051) + 6 more |
| 8.3.5 (152–153) | **Fletcher & Sharma 2014** | **4 intervals (A–D)**, apogee-altitude-based | Near-linear apogee segments, count/width not machine-extractable (table) | None | Not reported | RSM+GA | 144–148d predicted vs 154d actual | GSLV-D5 U/S (#39499) |
| 8.3.5 (153) | Mutyalarao & Raj 2015 | "Selected different time zones", **semi-major-axis**-based (not apogee) | Not reported (count) | None | Not reported | RSM+GA | Good agreement (qualitative) | GSLV-D5/CUS, Phobos-Grunt, ROSAT |
| 8.5 (154–155) | Sharma & Anil Kumar 2005, Kalman filter | Not zone-based — one continuous 11-day measurement stream | **54 TLEs, 17–28 Jan 2005** | None | Not reported | KS-element propagation | Close to actual, low uncertainty bands | Cosmos 2332 |
| 8.6 (156–157) | Anil Kumar & Subba Rao 2002 (BALEST) | Sequential online epochs, not discrete zones | "Last 15 epochs" gave good online predictions | None | Not reported | Simple mean-atmosphere propagator | Qualitative | SROSS-C2, Soyuz-1 stage |
| 8.7 (157–159) | **Sharma et al. 2004b (KSGEN)** | BC estimated across "different epochs", not explicit zones | Not reported (count); F10.7 = 81-day average prior to epoch | None | Not reported explicitly | KSNUM (KS elements, J2–J6, oblate diurnal atmosphere, RKG4, Jacchia-1977) | Close to actual across 7 IADC campaigns | Cosmos 1043, Cosmos 389, 25947, SROSS-C2 |
| 8.9 (159–161) | Chaotic-motion study | Discrete calendar epochs, not linear-decay zones | 6 epochs (Jan 2001–May 2002); separately last-3-days/15 epochs | None | Not reported | NPOE (36×36 gravity + drag + luni-solar, MSIS90) | Qualitative match | Unnamed rocket body |
| 8.10 (161–163) | **Sharma et al. 2009b (KSGEN)** | Fixed trailing-day windows, not linear-trend zones | **SROSS-C2: 14 TLEs / last 5 days. SL-12 R/B: 13 TLEs / last 3 days** | None (terminal-window only) | Not reported | KSGEN (KS + GA) | Max 7.1% (SROSS-C2), 12.4%/4.3min (SL-12) | SROSS-C2, SL-12 R/B |
| 8.11 (163–165) | Mutyalarao & Raj 2012 | "Selected time intervals" | Not reported (count) | None | Not reported | RSM+GA | Qualitative | Phobos-Grunt, ROSAT |
| 8.11 (165) | Mutyalarao & Raj 2014a | Trailing-day window | **Last 5 days before re-entry**, 5 epochs (Table 8.10) | None | Not reported | RSM+GA fit, HPOP(STK) propagation | Not reported numerically here | GSLV-D5/CUS |
| 8.12 (164–166) | **Anil Kumar et al. 2014** | "Selected TLE time intervals" | Not reported (count; tables are graphics) | None | Not reported | RSM+GA fit, HPOP(STK) | <6.9% abs. error | Phobos-Grunt, ROSAT (LEO, e<0.2) |
| 8.13 (166–167) | **STK OPTIM (Anil Kumar et al. 2014)** | **Sparse 3-point sampling**, NOT a decay-trend window | Exactly **3 TLEs: one at start, one at middle, one at end of a short recent span** | None | Not reported | HPOP (STK) | — | (method description) |
| 8.14 (167–168) | Re-entry Time Prediction of Phobos-Grunt/ROSAT | Trailing-day window, multiple predictions | **4 predictions from last 4 days**; 15 TLEs same-day for first case | None (terminal-only) | Not reported | HPOP (STK), random search over BC 340–400 | <6% (Phobos-Grunt), <3.4% (ROSAT) | Phobos-Grunt, ROSAT |
| 8.15 (168–169) | STK LTOptim | "A set of TLEs" propagated across a BC range | Not reported (count) | None | Not reported | HPOP (STK) | Table 8.15 (GOCE) | GOCE |
| 8.17 (169–171) | Anil Kumar et al. 2017 | 3-method comparison (RSM-GA / STKOptim / STKLTOptim) | Not reported (count) | None | Not reported | HPOP (STK) + RSM-GA | RSM-GA lowest error at last prediction | CZ-2C R/B (LEO, e<0.2) |
| 8.18 (172–173) | Mutyalarao & Anil Kumar 2018 | Linear-apogee-based intervals | Not reported (count) | None | Not reported | RSM+GA | Accurate 1 day before re-entry | GSLV-F09/CUS |
| **8.19 (173–175)** | **Sellamuthu & Sharma 2018 — KSROP** | **Explicit "zones"**, apogee-altitude-based, **directly analogous to OREM's own zone_select** | **5 zones** (Fig. 8.13); per-zone TLE count in Table 8.20 (not machine-extractable, table is a graphic) | None stated | Not reported | **KSROP** (same KS-regularized lineage OREM's `propagate_ks` descends from), oblate atmosphere (density scale height vs. altitude), J6 zonals | Max 5.78% (zone 5, Table 8.21); overall <6% | Several HEO RSOs |
| 8.20 (174–176) | Lawrence & Sharma 2019 | Not zone-based — trajectory/breakup/survival-rate tool, different problem class | Not applicable | Not applicable | Not reported | DRAMA 2.0 + MATLAB + NPOE | Re-entry time within ~22 min; landing-site offset ~250km lat/long | OSIRIS 3U |
| 8.21 (175–177) | Sellamuthu et al. 2019 | "Selected different time zones" | Not reported (count; tables are graphics) | None | Not reported | NPOE (not KSROP here) | <5% (6 objects) | 6 Molniya satellites |
| 8.22 (176–178) | Dutt et al. 2020 | RSM-GA minimizes error "in the specified time interval" (single window per instance) | Not reported (count) | None | Not reported | HPOP (STK) | Percentage error comparison, 3 methods | Electron R/B (#43166), PSLV-C39/IRNSS-1H (#42928) |
| 8.24 (178–179) | Dutt et al. 2021c | "Different time zones/time instances" | Table 8.24 gives per-time-zone accuracy (count not extractable) | None | Not reported | HPOP (STK), 4 methods incl. ABPro | ≤2.0% avg. w/ last TLE (3 methods) | Starlink-26, CZ-5B + 6 more (3 GTO, 3 LEO) |
| **8.25 (179–180)** | **Dutt et al. 2021b** | STKOptim/EBP fitting, manual-intervention pain points identified | Not reported (count) | None | **Only explicit "outlier removal" mention in the whole chapter** — listed as a manual step the paper automates away, no algorithm given | HPOP (STK) | Success rate / function evals / exec time (Table 8.25) | Various |
| 8.26 (180–182) | **Dutt et al. 2023 (peer-reviewed, ASR)** | Same as 8.24, extended to long-horizon comparison | Same test objects; window/horizon varies by method | **Method choice IS horizon-dependent**: ABPro best near-decay, RSM-GA best ~6 months out, STKOptim best at intermediate horizons (GTO); different split for LEO | Not reported | HPOP (STK), 4 methods | Table 8.26 | Starlink-26, CZ-5B, +6 more |

**Not applicable / no fitting-window methodology to extract**: §8.1
(general intro), §8.2 (IADC campaign background), §8.3.1–8.3.3 (physics
discussion, no fitted zones), §8.4 (material degradation — thermal/
experimental, not orbit-fitting), §8.8 (Bandyopadhyay et al. 2006 —
general RSM lifetime-vs-inclination study, no explicit zone count
given), §8.16 (STK Lifetime tool, launch-time sensitivity, no zone
count), §8.23 (Sellamuthu et al. 2020 — analytical SRP theory
validation, not a re-entry-time-from-TLE method), §8.27 (results
summary table of 19 IADC campaigns, not a methodology description).

## 3. Part 2 — Comparative review: OREM vs. literature

### 3.1 Comparison table

| Method/study | Zone definition | Start/end rule | TLE count/zone | Spread/coverage | Event dependence | Propagator | Accuracy metric | Reference (page) |
|---|---|---|---|---|---|---|---|---|
| **OREM (this repo)** | Statistically-significant linear apogee-decay window | R²≥0.90, slope<−1 km/day, min 8 pts, max 20 days; auto-detected | 8–~29 typical (min_zone_pts=8 floor, informally capped by 20-day window) | Up to `nzones_max`, spread across full history, later zones preferred | None by default (event-anchored only in one-off investigations) | KSROP (KS-regularized, `propagate_ks`) | RPE (per-zone + ensemble), 80 km re-entry threshold | `ALGORITHM.md` §3-4, `src/zone_select.F` |
| Mutyalarao & Sharma (2011) | Explicit "zones," rough linear apogee variation | Manual/visual segmentation over a fixed calendar span | Not reported | 8 zones over 200 days | None | Not specified | <7% (6 objects) | p.149–151 |
| Fletcher & Sharma (2014) | 4 near-linear apogee intervals | Not reported | Not reported | 4 intervals, span not reported | None | Not specified | 144–148d vs 154d actual | p.152–153 |
| **Sellamuthu & Sharma (2018)** | Explicit "zones," apogee-altitude-based | Not reported (visual/manual per fig. 8.13) | Per-zone count in Table 8.20 (not extractable) | 5 zones | None | **KSROP** (same lineage as OREM) | Max 5.78% (zone 5) | p.173–175 |
| Sellamuthu et al. (2019) | "Time zones," not detailed | Not reported | Not reported | Not reported | None | NPOE | <5% (6 objects) | p.175–177 |
| Sharma et al. (2004b/2009b), KSGEN | Fixed trailing-day windows, not decay-trend zones | Terminal N-day/N-TLE window before re-entry | 13–14 TLEs / 3–5 days | Single terminal window per prediction instant | None (terminal-proximity only) | KSNUM (KS + oblate diurnal atmosphere) | 7.1–12.4% max | p.157–163 |
| STK OPTIM (Anil Kumar et al. 2014) | **Not a window at all** — 3 discrete points (start/mid/end) | Manual pick within a short recent span | 3 TLEs | Single short span, sparse sampling | None | HPOP (STK) | Table 8.16–8.19 series | p.166–167 |
| Dutt et al. (2020/2021c/2023) multi-method | Comparative — each sibling method (RSM-GA/STKOptim/STKLTOptim/ABPro) uses its own window definition | Varies by method | Not reported | Multiple time instances/horizons per object | **Horizon-dependent method choice** (ABPro near-decay, RSM-GA ~6mo out, STKOptim intermediate) | HPOP (STK), 4 methods | % error vs. actual, method-vs-horizon table | p.176–182 |

### 3.2 Written analysis (~750 words)

**Common practices.** Every zone/window-based re-entry study surveyed
here — spanning three decades of ISRO work (Sharma et al. 2004b through
Dutt et al. 2023) — converges on the same two-stage architecture OREM
also implements: (1) segment the TLE-derived apogee (or semi-major-axis)
history into one or more "clean" fitting windows, and (2) recover a
small number of unobservable physical parameters (almost universally
eccentricity and an effective ballistic coefficient/number) by matching
a propagated trajectory against that window's observations, usually via
response-surface-methodology-plus-genetic-algorithm (RSM+GA) or an
equivalent nonlinear optimizer. This is not a coincidence: as OREM's own
`ALGORITHM.md` §1 already documents, this exact zones+RSM+GA lineage
traces directly back to Sharma, Bandyopadhyay & Adimurthy (2006) and
Mutyalarao & Sharma (2010/2011) — two of the very studies extracted
above — so the "common practice" found here is in a real sense OREM's
own ancestry, not an independent convergence.

A second point of near-universal agreement: **nobody in this literature
treats zone/window selection as a fully specified, reproducible
algorithm with numeric thresholds** the way OREM's `zone_select.F` does.
Mutyalarao & Sharma (2011) describe 8 zones from "rough linear
variation" of apogee; Fletcher & Sharma (2014) describe 4 "near linear"
intervals; Sellamuthu & Sharma (2018) show 5 zones in a figure without
a stated selection rule. In every one of these cases the zones read as
*manually or visually* segmented, with the RSM+GA machinery applied
per-zone afterward — the segmentation step itself is not the object of
methodological scrutiny the paper is written to defend. OREM, by
contrast, makes zone selection itself a first-class, parameterized,
auto-detecting algorithm (R²/slope thresholds, min-points/max-days
bounds) — a genuine methodological advance over its own literature
lineage, not just an implementation detail.

**Where OREM differs, and likely strengths/limitations of each.**
Three concrete numeric divergences stand out. First, the terminal
re-entry altitude threshold: OREM uses 80 km; Sharma et al. (2004b/
2009b, KSGEN) use 90 km; Anil Kumar et al. (2014/2017) use a *mean
perigee altitude* of 10 km. None of the three is obviously "more
correct" — they trade off differently against SGP4/propagator validity
near the edge of the atmosphere versus how close to actual breakup the
number is meant to represent — but a reader comparing RPE figures
across these papers and OREM's own results should not assume the
percentages are computed on a like-for-like re-entry definition.
Second, several of the closest OREM ancestors (KSGEN-era studies, STK
OPTIM/LTOptim, the multi-method Dutt et al. comparisons) rely on
**fixed trailing-day or fixed-TLE-count terminal windows** ("last 3
days," "last 5 TLEs," "3 points: start/mid/end") rather than a
decay-trend-quality test. This is simpler and guarantees a window
exists whenever there are recent TLEs at all, but it cannot
distinguish a genuinely quasi-linear decay segment from a noisy or
non-monotonic one — exactly the failure mode OREM's R²/slope gate is
designed to reject. Third, and most significant for OREM's own current
research thread: **none of the studies surveyed here condition zone
placement on a distinct physical event** (a resonance crossing, a
maneuver, a density anomaly) — every "zone" or "window" in this
literature is chosen purely from the TLE data's own statistical
behavior or its calendar proximity to re-entry. OREM's own recent
solar-apsidal-resonance validation work (issue #56, this session's
340-zone campaign) is, as far as this document search found, not
paralleled anywhere in the ISRO monograph's Chapter 8 — this is a real
gap OREM is currently the only party (in this literature set) actively
investigating.

**Gaps OREM does not address.** The multi-method comparisons (Dutt et
al. 2020/2021c/2023) are the single most important finding for framing
OREM's own limitations honestly: ISRO's own conclusion, reached across
three papers culminating in a peer-reviewed 2023 Advances in Space
Research publication, is that **no single method (including their own
RSM-GA) dominates across the full re-entry horizon** — different
methods (ABPro, STKOptim, STKLTOptim, RSM-GA) win at different
time-to-decay regimes, motivating automated method-switching rather
than single-method refinement. OREM implements RSM-GA exclusively; this
is a literature-backed, institution-level precedent that a purely
single-method RPE ceiling should be expected, not a sign that OREM's
own zone-selection or GA tuning is under-optimized. Separately, the
explicit outlier-handling gap: only one study in the whole chapter
(Dutt et al. 2021b, §8.25) names "outlier removal" as a step at all,
and even there it is listed as a *manual* pain point being automated
away, not a described algorithm. OREM's `tle_filter` (4σ/20-TLE-window/
30-day-gap, quantified false-positive rate) is measurably more rigorous
on this specific point than anything found in this survey.

**Points suitable for citing in a methodology/discussion section.**
(1) OREM's zones+RSM+GA architecture's direct citation lineage
(Sharma/Mutyalarao/Sellamuthu, already in `README.md`/`ALGORITHM.md`).
(2) The re-entry-altitude-threshold inconsistency across the field (80
km / 90 km / 10 km perigee) as a caveat when comparing RPE numbers
cross-study. (3) Dutt et al. (2023)'s no-single-method-dominates finding
as the field's own explanation for an RSM-GA-only system's expected
accuracy ceiling. (4) The near-total absence of a *quantified*,
reproducible zone-selection algorithm elsewhere in this literature, as
the specific, citable point of methodological novelty for OREM's own
`zone_select.F`.

### 3.3 Reference list

1. Anil Kumar, A. K. and Subba Rao, P. V. (2002), "Re-entry prediction
   accuracy improvement using Genetic Algorithm," 34th COSPAR
   Scientific Assembly, COSPAR 02-A01372.
2. Anil Kumar, A. K., Raj, M. X. J., Mutyalarao, M., Gupta, S. and
   Kumari, T. R. S. (2014), "Models for Re-entry Time Prediction of
   Risk Objects," National Conference on Space Debris Management and
   Mitigation Techniques, ASI-SDMMT.
3. Anil Kumar, A. K., Raj, M. X. J., Mutyalarao, M., Dutt, P.,
   Bhanumathy, P. B., Kumari, T. R. S. and Selvan, A. (2017),
   "Performance of prediction models on re-entry time prediction of
   CZ-2C rocket body," First International Conference on Recent
   Advances in Aerospace Engineering, ICRAAE 2017.
4. Bandyopadhyay, P., Sharma, R. K. and Adimurthy, V. (2001), "The
   Orbiting Third Stage of GSLV-D1 as Space Debris," Workshop on Space
   Debris Software Tools, VSSC/AERO/TR-001/2001, 86–98.
5. Bandyopadhyay, P., Sharma, R. K., Mutyalarao, M. and Adimurthy, V.
   (2006), "Lifetime Estimation of Upper Stages Re-entering from GTO
   with different inclinations," IAC-06-B6.2.11, 57th International
   Astronautical Congress, Valencia.
6. Dutt, P., Mutyalarao, M., Bhanumathy, P., Kumari, T. R. S., Negi, D.,
   Anil Kumar, A. K., Kumar, A. and Ashok, V. (2020), "Assessment of
   in-house developed re-entry prediction methods," IAA-UT Space
   Traffic Management Conference, IAA-UT-STM-03-04.
7. Dutt, P., Negi, D., Kumar, A. and Ashok, V. (2021b), "Assessment of
   Nonlinear Optimization Algorithms on Weighted Least-Square-Based
   Method for Re-Entry Prediction," Advances in Multidisciplinary
   Analysis and Optimization, Lecture Notes in Mechanical Engineering,
   Springer.
8. Dutt, P., Mutyalarao, M., Bhanumathy, P., Kumari, T. R. S., Negi, D.,
   Anil Kumar, A. K., Kumar, A. and Ashok, V. (2021c), "Assessment of
   in-house re-entry prediction algorithms," Stardust-R Global Virtual
   Workshop II.
9. Dutt, P., Mutyalarao, M., Bhanumathy, P., Kumari, T. R. S., Negi, D.,
   Anil Kumar, A. K., Kumar, A. and Ashok, V. (2023), "Assessment of
   in-house algorithms on re-entry time prediction of uncontrolled
   space objects," Advances in Space Research, 72, 2535–2551.
10. Fletcher, J. and Sharma, R. K. (2014), "Lifetime estimation of the
    upper stage of GSAT-14 in Geostationary transfer orbit,"
    International Scholarly Research Notices 2014, Article ID 864953.
11. Lawrence, C. S. and Sharma, R. K. (2019), "Re-Entry of Space Objects
    from Low Eccentricity Orbits," International Journal of Astronomy
    and Astrophysics, 9, 200–216.
12. Mutyalarao, M. and Sharma, R. K. (2010), "Optimal re-entry time
    estimation of an upper stage from geostationary transfer orbit,"
    Journal of Spacecraft and Rockets, 47, 686–690.
13. Mutyalarao, M. and Sharma, R. K. (2011), "On prediction of re-entry
    time of an upper stage from GTO," Advances in Space Research, 47,
    1877–1884.
14. Mutyalarao, M. and Raj, M. X. J. (2012), "Optimal re-entry
    prediction of space objects from low Earth parking orbits using
    RSM with GA," 39th COSPAR Scientific Assembly.
15. Mutyalarao, M. and Raj, M. X. J. (2014a), "Prediction of re-entry
    time of an upper stage from geostationary transfer orbit using RSM
    with GA," Proceedings of 59th Congress of ISTAM.
16. Mutyalarao, M. and Raj, M. X. J. (2015), "Performance of Response
    surface method using genetic algorithm on re-entry time prediction
    problem," Canadian Journal of Basic and Applied Sciences, 3,
    182–194.
17. Mutyalarao, M. and Anil Kumar, A. K. (2018), "Optimal re-entry time
    prediction of an uncontrolled spent stage of GSLV," Proceedings of
    National Conference on Essence of Mathematics and Engineering
    Applications (EMEA 17), 11–15.
18. Sellamuthu, H. and Sharma, R. K. (2018), "Optimal re-entry time
    prediction of resident space objects from highly elliptical
    orbits," 42nd COSPAR Scientific Assembly, Abstract PEDAS.1-18-18.
19. Sellamuthu, H., Sharma, R. K., Solomon, B. J., Ravichandran, S. and
    Khadri, S. P. M. S. (2019), "Optimized re-entry time prediction of
    Molniya orbit objects," IAC-19-A6.9.10, 70th International
    Astronautical Congress.
20. Sellamuthu, H., Sharma, R. K., Pushparaj, N., Ravichandran, S. and
    Khadri, S. P. M. S. (2020), "Regularized analytical orbit theory
    with solar radiation pressure," IAC-20-C1.VP.x60338, 71st
    International Astronautical Congress.
21. Sharma, R. K. and Anil Kumar, A. K. (2005), "Kalman Filter Approach
    for Re-Entry Predictions of Risk Objects with K-S Element
    Equations," 56th International Astronautical Congress,
    IAC-05-B6.3.05.
22. Sharma, R. K., Anil Kumar, A. K. and Raj, M. X. J. (2004b), "An
    integrated approach for risk object re-entry predictions in terms
    of KS elements and genetic algorithm," 35th COSPAR Scientific
    Assembly, 4146.
23. Sharma, R. K., Anil Kumar, A. K., Raj, M. X. J. and Sabarinath, A.
    (2009b), "On Re-Entry Prediction of Near Earth Objects with
    Genetic Algorithm Using KS Elements," 5th European Conference on
    Space Debris, ESA SP-672.
24. Sharma, R. K., Bandyopadhyay, P. and Adimurthy, V. (2004a),
    "Consideration of lifetime limitations for spent stages," Advances
    in Space Research, 34, 1227–1232.
25. Sharma, R. K., Anil Kumar, A. K. and Adimurthy, V. (2024), *Space
    Debris and Space Situational Awareness Research Studies in ISRO*,
    Monograph-2024, Indian Space Research Organisation. [this document]

---
*Prepared per the credit-saving session rules: PDF extracted once via
`pdftotext -layout` to a scratch text file, then read with targeted
`Read`/`grep` passes — no re-extraction, no broad repo exploration
beyond `ALGORITHM.md`/`src/zone_select.F`.*
