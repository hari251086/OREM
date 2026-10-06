# Combined resonant-azimuth theory for high-inclination HEO re-entry (draft v1, 2026-10-02)

Follow-up to OREM #56 and `E:\Research\References\11_Molniya_Resonances\literature_review.md`.
Goal: a theory in the spirit of Wang & Gurfil's solar-azimuth U-turn, but built on the drivers that
operate **above 46.4°**, where the solar apsidal resonance cannot occur.

Evidence tags: **[L]** literature (see review for tag detail), **[D]** derived here, **[N]** numerical result
from the scripts in this folder (TLE pool of 253 GTO/HEO objects; i ≥ 46.4° subset), **[?]** not established.

---

## 1. The observation that organises everything

For i ≥ 46.4° the apogee radius stays almost the same *relative to* the perigee radius (project owner's
observation). Quantitatively:

- **Mechanism [L+D].** In the quadrupole model averaged over the satellite's mean anomaly and the perturber's
  motion, the semi-major axis is an exact constant (Celletti et al. 2016, full text: "M is an ignorable variable,
  its conjugated action L (or equivalently the semi-major axis a) is a constant"). Hence δr_a = −δr_p = a·δe and
  the relative sensitivities satisfy (δr_p/r_p)/(δr_a/r_a) = −(1+e)/(1−e).
- **Check on TLEs [N]** (`check_apogee_perigee.py`; 91 objects, i ≥ 46.4°, first 70 % of each record, median e = 0.71):
  spread of log a = 0.12 %, of log r_a = 2.1 %, of log r_p = 7.5 %; measured ratio r_p/r_a spread **5.6** vs
  predicted (1+e)/(1−e) = **5.9**; correlation of increments dr_a vs dr_p = **−0.91**.

**Consequence:** the state variable that carries the physics is the perigee radius q = a(1−e), and the dynamics
is one degree of freedom in (e, angle) at fixed a. This is why a single "azimuth" angle can be enough.

## 2. The azimuth per perturber

For a circular perturber P (Sun, Moon) the doubly-averaged quadrupole disturbing function is [L: Lidov 1962,
Kozai 1962; same model as Celletti 2016, Alessi 2021]

    R_P = C_P [ (2 + 3e²)(3cos²ψ_P − 1) + 15 e² sin²ψ_P cos 2ω′_P ],   C_P = μ_P a² / (16 a_P³)

with ψ_P the mutual inclination between the satellite orbit and the perturber's orbital plane and **ω′_P the
argument of perigee measured from their mutual node**. Lagrange's equation gives [D]

    de/dt = (15/8)(μ_P / a_P³ n) · e √(1−e²) · sin²ψ_P · sin 2ω′_P ,    dq/dt = −a · de/dt .

So **ω′_P plays the role of the "solar azimuth"**: while 2ω′_P circulates, the drive oscillates and averages out;
when 2ω′_P turns round / lingers (ω̇′_P → 0) the drive keeps one sign and q changes monotonically — the
U-turn mechanism, now in a form valid at any inclination and for both perturbers. The perturber planes enter
through unit normals: Sun = ecliptic pole (fixed in the equatorial frame); Moon = ecliptic-referenced node Ω_L
regressing at −0.05295°/day (18.6 yr), inclination 5.145°, rotated by the obliquity (`quadrupole_drive.py`).

Wang & Gurfil's U-turn is the **singly-averaged, semi-secular** relative of this: its angle is (Ω+ω) − λ_sun and
needs d(Ω+ω)/dt = n_sun > 0 with d(Ω+ω)/dt ∝ (5cos²i − 2cos i − 1), i.e. i < 46.4° [L: Wang 2017 Eqs. 8–9,
full text; D]. Above 46.4° that drift is retrograde, so it cannot lock to the Sun.

## 3. The family of resonant angles (what ω′_P is made of)

Expanding 2ω′_P in the satellite's own angles gives the quadrupole family [L: Celletti 2016, (2−2p)ω̇ + mΩ̇ + sΩ̇_k = 0]

    Ψ_{k,m,s} = k·ω + m·Ω + s·Ω_L ,   Ψ̇ = λ(a,e)·G(i) + s·Ω̇_L ,   G = k(5cos²i − 1)/2 − m·cos i,
    λ(a,e) = 1.5 J2 n (R/p)² ,   p = a(1−e²)     (J2 secular rates ω̇ = λ(5c²−1)/2, Ω̇ = −λ c)

- **s = 0 (inclination-only) members reproduce the critical inclinations [N]** (`loci_table.py`, matches the
  literature values verified in full text): 46.38° (k=2,m=2), 56.06° (m=1), 63.43° (m=0), 69.01° (m=−1),
  73.15° (m=−2), 90°, and the retrograde pairs 106.85°, 110.99°, 116.57°, 123.94°, 133.62°.
- **s ≠ 0 (lunar-node) members** are stationary when λ·G = −s·Ω̇_L = s·0.0530°/day. With λ ≈ 0.1–0.3°/day in
  the HEO regime (e.g. 0.260°/day at a = 26,562 km, e = 0.70) these are reachable for |G| ~ 0.2–2, i.e. at
  essentially every inclination above 46.4° and a ≈ 20,000–45,000 km [N]. The Sun's semi-secular member
  (2ω + 2Ω − 2λ_sun, rate 2n_sun ≈ 1.97°/day) is reachable only for GTO-like a and e [N]; the Moon's
  semi-secular members (≥ 26°/day) are not.
- **Molniya-like state (a = 26,562 km, e = 0.70, i = 63.4°) [N]:** Ψ̇ for 2ω ≈ +0.0006°/day, and
  2ω+Ω−2Ω_L, 2ω−Ω+2Ω_L within ~0.01°/day: several drivers are simultaneously near stationary. The family is
  dense, so "resonance overlap" is the generic situation at these inclinations, in line with the chaos
  literature [L: Daquin 2016, Gkolias 2016, Rosengren 2019; abstracts].
- Members with k = 0 carry no ω, so ∂R/∂ω = 0 and they do **not** pump e directly at quadrupole order; they act
  through i and Ω and couple to e only via the precession rates.

## 4. Numerical results so far

| # | Test | Result | Reading |
|---|---|---|---|
| R1 | apogee vs perigee spread, 91 objects (§1) | 5.6× vs predicted 5.9×; corr(dr_a, dr_p) = −0.91 | confirms the one-DOF picture |
| R2 | s=j=0 family members → critical inclinations | 46.38, 56.06, 63.43, 69.01, 73.15, 90 deg | family is consistent with the literature |
| R3 | window-mean of the exact Sun+Moon quadrupole drive vs observed dq/dt (180-day windows; 62 objects, 3,845 windows; windows with a-drift > 0.3 % or perigee height < 400 km dropped) | per-object **median Pearson +0.993** (IQR 0.91–0.997), **Spearman +0.977**, **slope +1.00 (IQR 0.99–1.01)**, sign agreement 1.00; lag-shift null rejects chance for **89 %** of objects | the parameter-free drive explains the observed perigee trend |
| R4 | coherence κ = \|mean\|/rms of the drive; slowest driver in each window | κ median 0.997 — 180-day windows are shorter than the angle's own period; the "slowest of ~40 drivers" is always slow (median 0.004–0.007°/day) | **inconclusive**: windows too short, family too dense for argmin labelling |
| R5 | strength E = A/\|ω̇′\| (amplitude/frequency) vs observed q-range in 540-day windows | Spearman +0.54 pooled (log–log corr 0.64, slope 0.45); +0.85 at 46–56°, +0.18 at 56–70°, +0.54 above 70°; per-object median +0.30 | **partial**: E is the asymptotic excursion; windows shorter than a half-cycle saturate at ≈ A·W |

Caveats on R3 (important): (i) the pooled Pearson is dominated by a few large-amplitude objects (top 5 carry
82 % of Σobs²) — the per-object numbers above are the ones to use; (ii) TLE elements are mean elements fitted
with an SDP4-type model that itself contains lunisolar terms, so R3 shows the physics is *consistent with the
tracked orbit evolution*, not that TLEs are independent truth; (iii) all in-sample, no free parameters fitted;
(iv) only 62 of the 188 high-inclination objects had enough clean, non-drag windows.

## 5. What this gives OREM (proposals, not validated)

1. **A parameter-free secular perigee predictor.** R3 says Δq over a window ≈ ∫(Sun+Moon drive) dt using only the
   current mean elements — usable as a cheap cross-check of any OREM perigee trajectory at i > 46.4°, and as the
   dynamical core of a one-degree-of-freedom propagator (a constant, e and ω′ evolving).
2. **A resonance flag in the spirit of the U-turn:** persistence time of sign(drive) — long persistence ⇒
   ω̇′_P ≈ 0 ⇒ monotone q change ahead. Still to be defined and validated (see §6).
3. **Domain gate for the existing detector** (done in #56): i < 46.4° → solar-azimuth U-turn; i ≥ 46.4° → this
   framework.
4. Use as a **flag / dispersion widener**, not as a "predict from here" trigger — the #56 finding (post-resonance
   restriction degraded accuracy) still applies.

## 6. Open items / next steps

- **Closed form of ω̇′_P** — done (§8.3): exact chain rule, verified; but not a usable stand-alone U-turn locus (the azimuth rate swings sign with the nodal phase).
- **Persistence-time detector** — superseded by the forecast detector of §8.4, which predicts the sign change of the drive (the perigee peak) and the fall below 300 km years ahead.
- **[?] Generality check:** the same machinery should also reproduce the drive at i < 46.4° (including the 7
  confirmed solar U-turn objects); not yet run.
- **[?] Beyond quadrupole:** Moon octupole grows with a (Celletti 2016, full text); Sun semi-secular terms near
  GTO; J2-modified Kozai integral.
- **[?] Link to re-entry:** Δq-to-decay relation — when does the drive-driven perigee drop cross the drag
  threshold, and how does the lead time distribute (cf. 0 d – 8 yr spread found for the U-turn)?
- Figures (drive vs observed; ω′_P tracks with stationary episodes) — to be made with the dataviz skill.

## 7. Files

`drivers.py` (family, J2 rates, loci) · `loci_table.py` · `quadrupole_drive.py` (exact drive, ω′_P) ·
`check_apogee_perigee.py` (§1) · `validate_drive.py`, `robust_stats.py` (R3) · `identify_drivers.py` (R4) ·
`validate_strength.py` (R5). Data: `apogee_perigee_stats.csv`, `drive_windows.json`.

---

## 8. The Molniya perigee cycle: what drives the rise before decay (second pass, 2026-10-02)

**Question (from the OREM validation report plots):** for Molniya-type orbits the perigee rises by 1–3 thousand km over
~5–7 years and then falls to decay, while the apogee moves the opposite way. What drives it, and is it a known resonance?

### 8.1 Observation (OREM_Validation_Report.pdf, evolution pages) [N]
MOLNIYA 1-28 (NORAD 7480): perigee altitude 650 → 1,890 km (peak at 6.7 yr) → decay at ~11 yr; apogee 39,700 → 38,450 km
over the first 7 yr — the constant-a signature of §1 (δr_p ≈ −δr_a) on a multi-year scale.

### 8.2 The driver is the lunisolar quadrupole [N]  (`molniya_cycle.py`, `forecast_test.py`, `make_figures.py`)
19 Molniya-type objects (i ≈ 62–66°, e ≈ 0.7, a ≈ 26,000–27,000 km: Molniya satellites, SL-6 R/Bs, Cosmos debris):
- **Reconstruction.** Integrating the exact Sun+Moon quadrupole drive along the observed elements reproduces the perigee
  rise in 18 of 19 objects: correlation ≈ 1.000, median normalised RMSE 0.009; rise amplitudes 1,100–3,300 km with
  modelled/observed ratio 1.00. The Moon supplies ≈ 70 % (65–73 %), the Sun ≈ 30 %. (Exception: COSMOS 1030 DEB, 12907,
  peak at 32 yr, outside the regime.)
- **The peak is a sign change of the Moon's drive.** At the observed perigee maximum, 2ω′_Moon lies within a few degrees of
  180° for every object (175°–180°), i.e. sin 2ω′_Moon changes sign and the perigee turns.
- **Forecast, no free parameters** (secular propagator `secular_propagator.py`: J2 + Sun + Moon quadrupole, a = const, no
  drag, started from TLE mean elements): from the *first* TLE (horizon ≈ 10.7 yr) the perigee-peak epoch has median error
  0.11 yr and the epoch of perigee < 300 km 0.19 yr (bias +0.00; 17 objects); started 4 yr / 2 yr before the peak, 0.11 / 0.08 yr.
  A naive "every object lives the median lifetime" baseline has mean error 1.67 yr vs 0.22 yr (skill 0.87); predicted and
  observed horizons (9–15 yr) correlate at 0.991.
- **Ablation (which driver is required):** J2 alone gives no rise; J2+Sun gets the peak epoch wrong by 7.7 yr; J2+Moon gets the
  peak right (0.23 yr) but the fall wrong by ≈ 5.8 yr; **both perturbers are needed** — the Moon makes the rise and its peak,
  the Sun the descent (Fig. 1, middle panel; Fig. 2, right).
- **Beyond the Molniya clones** (`forecast_pool.py`): 46 pool objects with i_mean ≥ 46.4°, e > 0.3, record reaching 300 km:
  median relative error 3 %, 86 % within 10 %, 90 % within 25 %, mean abs error 0.89 yr vs 4.45 yr for the median-lifetime
  baseline (skill 0.80), log–log corr 0.944. Mix: 30 objects at 62–66°, 3 at 46–56°, 6 at 66–90°. Failures: 4 of 46 give no
  drop within 45 yr (e.g. 10370, inclination varies 36 to >46°) and two are off by a factor 2–3.

### 8.3 Which resonance? — and a correction about the closed form
- **Known-resonance reading [N+L].** These orbits sit at the **critical inclination (2ω̇ ≈ 0, 63.43°)**: the J2 apsidal drift almost
  cancels, so the quadrupole drive ∝ sin 2ω′ keeps one sign for years. At a Molniya state (a = 26,562 km, e = 0.70, i = 63.4°)
  Ψ̇ = +0.0006°/day for 2ω and ≈ ±0.01°/day for 2ω+Ω−2Ω_L, 2ω−Ω+2Ω_L (§3). This is the lunisolar critical-inclination resonance
  of the literature, here in its Kozai-type "apsidal slow-passage" form relative to the Moon/Sun planes.
- **Azimuth rate, exact chain rule (verified).** ω̇′ = ω̇ − (∂ω_N/∂Δ)(Ω̇ − Ω̇_P) − (∂ω_N/∂i)·i̇ − (∂ω_N/∂I)·İ_P with
  tan ω_N = sin I sin Δ /(sin i cos I − cos i sin I cos Δ) (`omega_prime_rate.py`): matches the vector geometry to 3×10⁻¹⁴° and
  finite differences along the propagated trajectory to 4 digits.
- **But it is not a stand-alone U-turn locus.** (i) For the Sun the J2-only stationarity condition f(i, Ω) = (5c²−1)/2 + c·∂ω_N/∂Δ = 0
  has a node-dependent critical inclination, i*(Ω) = 57.7° (Ω = 180°) … 69.6° (Ω = 0°) [N]. (ii) Along a real Molniya orbit the
  azimuth rate is **not** slow and signed: it swings from −18 to +21°/yr within ~6 yr because Ω circulates at ≈ −41°/yr
  (period ≈ 8.8 yr) and ∂ω_N/∂Δ depends on the nodal phase. A straight-line "observed sweep rate" therefore disagrees with the
  J2-only closed form (median ratio 0.36, spread −0.95…+4.6) and with the full-quadrupole chain rule (median 0.34); this is a
  comparison artefact of averaging a strongly oscillating rate, not an error in the formula. **Physical reading:** ω is nearly
  frozen (critical inclination) while the nodal phase modulates the drive with the 8.8-yr period; the drive is positive for
  ≈ 6–7 yr and then negative, which is the multi-year perigee cycle. Hence a single inequality in (a, e, i) is the wrong detector;
  the quantity to follow is the sign and phase of the full drive, i.e. the forecast below.

### 8.4 A U-turn-style time detector for i > 46.4°
Detector = **epoch at which the Moon+Sun quadrupole drive changes sign (perigee maximum), and the epoch the perigee then drops below
300 km**, obtained by integrating the secular model from the latest TLE mean elements. Performance (Molniya-class, §8.2): peak epoch
±0.04–0.11 yr, fall below 300 km ±0.08–0.19 yr, 2–10.7 yr ahead. As a real-time flag: "perigee is rising and the Moon drive is still
positive" ⇒ time to peak and to the 300 km crossing are forecastable, giving years of lead before decay; after the peak the
remaining lifetime is bounded by the same forecast.

### 8.5 Caveats (read before using)
- **Mean elements.** TLEs are SDP4-type mean elements; agreement shows consistency with tracked evolution, not independence from the TLE model.
- **Sample.** The 17 Molniya-class objects are near-clones (similar a, e, i, ω; different RAAN/epoch); the 46-object pool is 65 % Molniya-class.
  Out-of-class generality rests on 9 objects (46–56° and 66–90°).
- **Definitions.** E300 (perigee < 300 km, smoothed) is chosen to stay above the drag regime; no drag is modelled, so the final months are not
  predicted. No parameter was fitted, but event definitions and the 22-year/45-year horizon caps were set before running; results are in-sample on this pool.
- **Model scope.** Quadrupole only (octupole matters for the Moon at larger a); Moon and Sun on circular orbits; a held constant.
- **Status.** Analysis tooling in `scratch_rpe/`; no change to OREM's prediction pipeline.

### 8.6 Next
(1) Run the forecast as an *independent cross-check* beside OREM's Molniya predictions (OREM median RPE 0.98 % overall on the 28-object
campaign); (2) add drag at low perigee to extend the forecast to re-entry itself; (3) test the same forecast on the 7 confirmed solar-U-turn
objects (i < 46.4°); (4) quantify sensitivity to initial-element noise (propagate perturbed initial states) to turn the forecast into a distribution.

**Figures:** `figures/fig1_molniya_7480_forecast.png`, `figures/fig2_forecast_skill.png`. **Scripts:** `molniya_cycle.py`,
`secular_propagator.py`, `forecast_test.py`, `forecast_baseline.py`, `forecast_pool.py`, `omega_prime_rate.py`, `make_figures.py`.

---

## 9. Trends across the 28 plots of OREM_Validation_Report.pdf (third pass, 2026-10-02)

Read all-object evolution plots (pages = summary + 1 of each 4-page object block; 11 viewed directly, all 28 re-extracted from the
same TLE data by `trend_all28.py`; `trend_phase.py`; Fig. 3 `figures/fig3_report_trends.png`). The report's "resonance angle" is the
solar azimuth ψ = (Ω + ω − L_sun) mod 360° (Wang & Gurfil's angle).

| Group | n | Inclination | Solar azimuth ψ | Perigee | Apogee |
|---|---|---|---|---|---|
| High-i | 23 | 63.5–68° | circulates steadily, ψ̇ = −1.08…−1.18°/day (period ≈ 0.9 yr), same value in the last 15 % of the record | hump +1,100…+3,300 km (median ≈ 1,500), peak 4.4–7.6 yr after the first TLE, then fall to decay | mirror image (a const) |
| Low-i | 5 | 5.7–25.5° | ψ̇ ≈ −0.3°/day and shrinking; **sign flips in the last 30 % for 3 of 5** (35497, 27526, 37151) | stays 120–260 km (small ±20–35 km annual swings; FALCON 9 R/B +147 km) | falls steadily 40,000 → < 10,000 km (apogee-led decay) |

- **High-i: the solar U-turn is never engaged [N+D].** ψ̇ = (ω̇ + Ω̇) − n_sun with (ω̇ + Ω̇) ∝ (5cos²i − 2cos i − 1) < 0 above 46.4°, so ψ̇ ≈ −(0.14 + 0.99)°/day
  ≈ −1.1°/day for ever, as measured. The perigee hump therefore has nothing to do with the plotted angle; it is the lunisolar quadrupole of §8
  (Moon ≈ 70 %, Sun ≈ 30 %), reproduced to ~1 % by the model.
- **Low-i: the U-turn is visible in the plots [N].** |ψ̇| falls and crosses zero within ~1–2 yr of the end for 3 of 5 objects, as the solar-resonance picture predicts;
  the perigee, however, stays low (drag regime), so the signature is a late apogee-led collapse, not a perigee hump. (40943 and 59347 do not reverse: ψ̇ −0.17 / −0.85°/day.)
- **Post-peak lifetime is bimodal and tied to the starting nodal phase [N].** For the high-i group the time from the perigee peak to the last TLE is ≈ 3.8–4.7 yr for
  first-TLE Ω ≈ 120–225° and ≈ 7.4–10.2 yr for Ω ≈ 20–115° (cos/sin regression R² = 0.63, n = 22); the time of the peak (R² = 0.19) and the rise amplitude
  (R² = 0.10) depend only weakly on Ω. Consistent with the 8.8-yr nodal modulation of §8.3 and reproduced by the forecast of §8.2 (error ≈ 0.2 yr). Correlation, small n,
  first-TLE epochs differ; 9269 (Ω₀ = 205°, 8.0 yr) and 27902 (Ω₀ = 356°, 12.8 yr) do not follow the grouping.
- **Exceptions in the high-i group:** COSMOS 1109 DEB (32977) starts near its peak (rise 47 km); COSMOS 1030 DEB (12907) and 1261 DEB (27902) show multi-decade cycles with
  perigee 3,500–5,800 km (the model over-predicts 12907); FREGAT R/B (35009) starts at 196 km perigee and rises 900 km.

**Recommendation for the report / Object Detail page:** for i > 46.4° the top panel (solar azimuth) is a featureless sawtooth and hides the driver. Replace or
supplement it with the Moon azimuth 2ω′_Moon (or the Sun+Moon drive dq/dt, Fig. 1c), whose sign change marks the perigee peak.

**Implemented 2026-10-02:** the recommendation above is done for the validation report (OREM-Watchlist `lunisolar_drive.py` + `plot_object_evolution(tle_path=...)`; 23 of 28 objects now show the Sun+Moon drive in the top panel; `OREM_Validation_Report.pdf` replaced, 135 pages). The dashboard Object Detail page is unchanged until its evolution CSV carries i, Ω, ω.

---

## 10. Cross-check against OREM's Molniya predictions (2026-10-03)

Question: do the parameter-free secular forecast (§8) and OREM's validation predictions agree, and where do they diverge?
Script `cross_check_orem.py` (outputs `cross_check_orem_{300,200,150}.csv`, `figures/fig4_cross_check_orem.png`).

**Setup, and why it is not a like-for-like contest [read first].**
- OREM: validation-campaign hindcast, TLE history cut ≈ 90 d before decay, full multi-zone RSM/GA; error = predicted − SATCAT decay (days).
- Forecast: J2+Sun+Moon secular propagator from the **first TLE only** (≈ 12–13 yr ahead), no drag; it predicts the epoch the perigee altitude first
  falls below a threshold after the perigee peak. To express it as a decay date, that epoch is shifted by the **leave-one-out median** of
  (decay − event epoch) over the *other* objects — one empirical constant, never the object's own decay date.
- The two use different information (OREM sees the end of the life; the forecast sees only the start), so they test different things.
- Objects: the Molniya-type members of the 28-object campaign with a usable event (17 / 16 / 11 at 300 / 200 / 150 km; 9506 never reaches 300 km after its peak).
- Event definition: first crossing *after* the smoothed perigee peak, same on observed and model series. A first run used the first crossing anywhere and
  gave identical numbers; a 30-yr propagation window was tried and rejected (the model's global peak lands on a later hump; the earlier tests used 22 yr).
  Three thresholds are reported because 300 km turned out to be a poor decay proxy for three objects (below); none was selected after the fact.

| Event: perigee < | n | OREM median \|err\| (bias) | Forecast median \|err\| (bias) | Forecast mean / max \|err\| | median \|OREM − forecast\| | corr. of errors |
|---|---|---|---|---|---|---|
| 300 km | 17 | 41 d (−41) | 155 d (−39) | 360 / 1,639 d | 146 d | +0.02 |
| 200 km | 16 | 43 d (−43) | 123 d (−38) | 278 / 1,602 d | 122 d | −0.07 |
| 150 km | 11 | 53 d (−53) | 121 d (−37) | 330 / 1,494 d | 138 d | −0.72 (n = 11, unstable) |

As a share of the history span (decay − first TLE): OREM ≈ 1 %, forecast ≈ 3 %.

**Agreement.** At 300 km the forecast decay date is within 300 d of the actual for 14 of 17 objects (within 150 d for 8), from the first TLE, a decade ahead; the two
methods' predicted decay dates differ by a median of 122–146 d. Both are biased early by a similar amount (OREM ≈ −41 d, forecast ≈ −39 d) — an observation only; no cause established.
The error correlation is ≈ 0 at 300 and 200 km (n ≈ 16–17), i.e. the errors look independent; the −0.72 at 150 km rests on 11 objects and is not stable.

**Divergence — one clear failure mode.** For NORAD 8844, 9269 and 11073 the forecast is off by −1,273 / −1,639 / −1,264 d while OREM is within 14–66 d. In the TLE
data the perigee of these objects dips below 300 km on a *marginal first minimum*, then lunisolar forcing raises it again (second hump), and decay only happens 3.5–5 yr later.
A threshold event on a no-drag model flags the first approach as decay; OREM, which sees the real tail, does not. The same sensitivity (a minimum near the
threshold) is the weak point of using this forecast as a decay-date predictor; it is better read as locating the **approach windows** (perigee peak, perigee minima)
than as a decay date unless the minimum depth is resolved.

**Reading.** Complementary, not competing: OREM is the near-term tool (data up to ~90 d before decay, median ≈ 41 d); the forecast gives a decade-ahead window with
median ≈ 4–5 months for these objects and flags where a perigee minimum is marginal. **Untested idea:** use the forecast as a prior/flag for "second-hump risk" when
OREM is run live near a marginal dip — the situation in which an extrapolating method is most likely to call decay too early.

**Next test.** Mid-life cutoff for 8844 / 9269 / 11073: run OREM with the TLE history cut at the first dip below 300 km, and the forecast from the same epoch, and compare both with the real decay.

**Caveats.** Small n (17/16/11), Molniya-type only, forecast event and LOO gap constructed after seeing a first run (disclosed above), TLE mean elements, no drag.

---

## 11. Mid-life cutoff test: OREM vs the forecast at the first perigee dip (2026-10-03)

Follow-up to §10. Idea being tested: at a marginal first perigee minimum an extrapolating method might call decay too early, and the lunisolar forecast —
which sees the perigee recover — could flag that. Script `midlife_test.py` (outputs `midlife_cutoffs.csv`, `midlife_forecast.csv`, `midlife_summary.csv`).

**Design.** Cutoff = first epoch the smoothed perigee falls below 300 km after its peak. Truncated TLE histories (830–1,472 TLEs) were run through **production OREM**
(v1.48 exe, Watchlist cfg conventions, `write_orem_cfg` defaults; nothing written to the Watchlist DB or cache) and the secular forecast was started from the
elements at the same epoch (15-yr propagation, event = first model crossing below a perigee threshold, shifted by the leave-one-out median gap decay − event from the
*other* objects' full lives). Compared with the real SATCAT decay. Objects: the 3 second-hump failures of §10 (8844, 9269, 11073; decay 1,318–1,818 d after the cutoff) and
4 ordinary controls (7480, 8833, 14297, 7903; decay 96–572 d after). **n = 7 — indicative only.**

| NORAD | Group | Decay after cutoff | OREM | Forecast err, 150 km event | Forecast err, 200 km event |
|---|---|---|---|---|---|
| 8844 | second hump | 1,525 d | no re-entry date | +75 d | +107 d |
| 9269 | second hump | 1,818 d | no re-entry date | +121 d | −1,422 d |
| 11073 | second hump | 1,318 d | no re-entry date | +147 d | −887 d |
| 7480 | control | 336 d | no re-entry date | +39 d | +5 d |
| 8833 | control | 183 d | no re-entry date | +75 d | +108 d |
| 14297 | control | 96 d | PRIMARY 1998-05-24 (+13 d) | +212 d | +237 d |
| 7903 | control | 572 d | no re-entry date | −215 d | −241 d |

**Findings.**
- **The hypothesis was not supported.** OREM did *not* call decay too early at the dip. For 6 of 7 objects it reported "no zone predicted a re-entry within the
  propagation cap" (including all three second-hump objects and three of four controls); only 14297, 96 d from decay, got a date (error +13 d). The reason is structural:
  OREM extrapolates the apogee decay trend, which has not started yet when the perigee first dips; its predictions become informative only in the last phase.
- **So the forecast's value here is information, not correction of a wrong OREM call:** it gives a decay date at the dip for all seven objects — median |error| 121 d at the
  150 km event (mean 126 d; second-hump objects +75/+121/+147 d, 3.6–5 yr ahead) — in the situation where OREM is silent. For comparison the first-TLE forecast of §10 was off by
  −1,273/−1,639/−1,264 d for the same three objects; started from the dip the model sees the perigee behaviour and recovers.
- **Threshold sensitivity again:** the 200 km event fails for two of the three second-hump objects (−1,422 and −887 d) because the model perigee dips just below 200 km on its first
  minimum; the 150 km event does not. The 150 km result for these objects was seen after the 200 km one, so it should be treated as a post-hoc observation until confirmed on
  more objects (§10 reports all three thresholds; this test reports two).
- Controls: forecast errors −215…+212 d (median |err| 75 d at 150 km); OREM silent on three of four (it needs the final apogee collapse).

**Caveats.** n = 7; Molniya-type only; the leave-one-out gap borrows the other objects' full lives (as in §10); TLE mean elements; no drag; "OREM silent" is a statement about this
configuration at this cutoff, not about OREM's accuracy near decay (§10: median 41 d at 90 d before decay).

**Next.** Confirm the 150 km (vs 200 km) behaviour on all Molniya-type objects with a mid-life cutoff; add drag so the event is decay itself rather than a threshold; consider surfacing the forecast
as an "approach window" annotation alongside OREM's PRIMARY when OREM reports no re-entry (the situation of this test).

---

## 12. Mid-life cutoff test, extended: all 17 Molniya-type objects + 2 extra second-hump objects (2026-10-03)

Confirmation of §11 (script `midlife_test.py`, `midlife_extra.py`, `scan_second_hump.py`; `midlife_summary.csv`). §11's seven objects were examined before the
150-vs-200 km choice, so the other 10 of the 17 serve as an out-of-sample check; the pool was also scanned for further second-hump objects.

**Results (cutoff = first dip below 300 km after the perigee peak; decay 96–1,818 d later).**

| Set | n | OREM gave a date | Forecast 200 km event | Forecast 150 km event |
|---|---|---|---|---|
| Seen in §11 | 7 | 1 / 7 (14297, +13 d) | median \|err\| 237 d, mean 430, >1 yr off: 2 | median 121 d, mean 126, within 250 d: 7/7 |
| **New, out-of-sample** | 10 | **0 / 10** | median 68 d, mean 66, within 250 d: 10/10 | **median 48 d, mean 46, within 250 d: 10/10** |
| All | 17 | 1 / 17 | median 96 d, mean 215, within 250 d: 15/17 | median 69 d, mean 79, within 250 d: 17/17 |

- **OREM is silent at this cutoff on 16 of 17 objects** ("no zone predicted a re-entry within the propagation cap"), including the ten new ones whose decay came 159–433 d later. Production
  configuration, OREM v1.48; this concerns this mid-life cutoff, not OREM near decay (§10: median 41 d at 90 d before decay).
- **The forecast dates all of them.** On the ten new objects every 150 km estimate is within 85 d of the real decay (−61…+85 d). Errors are mostly positive (+19…+85 d): the leave-one-out gap
  slightly over-shifts; no correction was applied.
- **Second-hump objects.** The pool contains only two further objects with a perigee dip below 300 km, recovery and > 2 yr of remaining record (42908, i = 47.7°, decay 4.03 yr after the
  dip, model max perigee 513 km; 6231, i = 65.9°, 2.69 yr, 415 km). The forecast from the dip gives +171 / +175 d at 200 km and +151 / +160 d at 150 km — both thresholds work on them. With the
  original three (150 km: +75 / +121 / +147 d) that is **five second-hump objects, all within 160 d at 150 km**; at 200 km two of the five (9269, 11073) fail by 2.4–3.9 yr.

**What is and is not confirmed.**
- Confirmed (out of sample, n = 10 + 2): the forecast from a mid-life dip gives a usable decay date where OREM gives none, at 150 km and at 200 km.
- *Not* confirmed: that 150 km is better than 200 km for second-hump objects. The 200 km failures are two in-sample objects (9269, 11073, marginal first minima); the two out-of-sample second-hump
  objects succeed at both thresholds. At 150 km the forecast has never been worse than 215 d on any of 19 objects, but the preference rests on those two. Note 42908 shows the approach also works at i = 47.7°.
- The ten "new" objects are all ordinary (decay 159–433 d after the dip), so they test the ordinary case, not the second-hump case.

**Caveats.** Molniya-type and two further high-inclination objects only; leave-one-out gap borrows other objects' full lives; TLE mean elements; no drag; event = threshold crossing, not decay itself;
OREM-silent statement is configuration- and cutoff-specific.

**Next.** Add drag so the event is decay itself (removes the gap conversion); show the forecast as an "approach window" alongside OREM's PRIMARY when OREM reports no re-entry; extend to objects outside the
Molniya class once more second-hump examples exist.

---

## 13. Adding drag: the event becomes decay itself (2026-10-03)

Follow-up to §10-§12, whose forecasts predicted a perigee-threshold epoch and converted it to a decay date with an empirical leave-one-out gap. Here the secular propagator is extended
with orbit-averaged drag so the model reaches decay on its own (`secular_drag.py`, `secular_drag_prop.py`, `drag_test.py`, `test_secular_drag.py`; `drag_test.csv`; `figures/fig5_forecast_with_drag.png`).

**Model.**
- State (a, e, i, Ω, ω); J2 + Sun + Moon quadrupole rates as before (`secular_propagator.rates`) plus orbit-averaged drag. Drag acceleration a_D = −½ B ρ v² v̂ (B = C_d A/m; atmosphere at rest,
  co-rotation neglected). Instantaneous da/dt = −(a²/μ) B ρ v³ and dh/dt = −½ B ρ v h (h = angular momentum); de/dt from h² = μ a (1−e²); both averaged over the orbit by quadrature in eccentric anomaly.
- Density: OREM's own `input/ATM.DAT` (60-500 km scale height and density — a static, **low-activity** atmosphere, 5×10⁻¹² kg/m³ at 300 km), extended above 500 km with a growing scale height.
  No solar-cycle variation, no diurnal bulge, no epoch-resolved space weather (OREM proper uses all three).
- **B is not fitted**: C_d, A, m from the Watchlist object parameters (DISCOS for 16 of the 17 objects: BN = m/(C_d A) ≈ 50-61 kg/m²; 14297 uses a generic default, BN 23), with ×0.5 and ×2 sensitivity.
- Adaptive step (perigee altitude changes ≤ 10 km per step). **Decay** = perigee altitude < 80 km (OREM's re-entry threshold) or apogee altitude < 150 km (an orbit wholly below this circularizes in dense atmosphere and
  decays within days). Sanity tests: circular-orbit limit da/dt = −Bρ√(μa) within 2 %, density anchors and monotonicity, drag acting at perigee.
- Physical behaviour seen in the model (7480): perigee falls to ~160 km, then the apogee collapses over ~130 d (40,000 → 12,000 km) and the orbit circularizes near 93 km — the same shape as the observed last phase.

**Disclosure of a false start [read before the numbers].** A first run stopped when e → 0 and returned "no decay" for most objects; its summary statistics, computed only over the runs that returned a date, looked good and were
**survivorship-biased**, so they were discarded. The termination rule was corrected on physical grounds (a circularized orbit at low altitude *is* decay) after inspecting the 7480/7641 trajectories (which also showed the model decay
within ~35 d of the real one — so the corrected rule was fixed having seen two objects' outcomes). Statistics below count every object; all 17 produce a date.

**Results — decay-date error in days (17 Molniya-type objects).**

| Start | Method | median \|err\| | mean | bias | within 250 d | > 1 yr off |
|---|---|---|---|---|---|---|
| Mid-life cutoff (first dip < 300 km; decay 96-1,818 d ahead) | no drag + gap, 150 km event (§12) | 69 | 79 | +55 | 17/17 | 0 |
| | **with drag, B × 0.5 / × 1 / × 2** | **26 / 30 / 35** | 42 / 47 / 58 | −26 / −28 / −33 | 17/17 | 0 |
| First TLE (≈ 10-15 yr ahead) | no drag + gap, 300 km event (§10) | 155 | 360 | −39 | 12/17 | 3 |
| | with drag, B × 0.5 / × 1 / × 2 | 200 / 202 / 156 | 277 / 271 / 268 | −40 / −41 / −41 | 12/17 | 3 / 4 / 4 |

- **From a mid-life dip, drag helps and removes the empirical gap**: median 30 d (vs 69 d), 17/17 within 250 d, and the answer is insensitive to B over a factor of 4. The three second-hump objects now come out
  −35 / −79 / −20 d (they were +75 / +121 / +147 d with the 150 km event), 3.6-5 yr ahead. The model decays on average ~28 d early.
- **From the first TLE, drag does not help**: errors of the two methods are similar (median 202 d vs 155 d; 12/17 within 250 d each). Drag fixes some objects (8844, 9269, 11073) and worsens others (7480 +670 d, 7903 +1,220 d);
  at that horizon the error is dominated by the timing of lunisolar minima (whether a marginal perigee minimum is crossed), not by the drag law.
- For comparison, production OREM at this mid-life cutoff gave a date for 1 of 17 objects (§12).

**Caveats.** Static low-activity atmosphere (no solar cycle; real decays at solar maximum come earlier and at minimum later — the −28 d bias is not explained by this and not investigated); atmosphere co-rotation neglected;
B from DISCOS/default values; TLE mean elements (mean vs osculating perigee not distinguished); n = 17, Molniya-type, in-sample for the stopping-rule choice (HA_DECAY = 150 km, fixed after the false start, not tuned on errors);
no independent validation set beyond these objects.

**Next.** Epoch-resolved solar activity (OREM's `SW-All.csv` / ATM2D); then surface the forecast as an "approach window" next to OREM's PRIMARY when OREM reports no re-entry; confirm on non-Molniya high-inclination objects.

## 14. Window coverage: does max(250 d, 0.15 H) hold the true decay date? (2026-10-03)

The Watchlist "approach window" is the forecast decay date +- max(250 d, 0.15 x forecast horizon H). Two coverage tests, forecast started 0.5 / 1 / 2 yr before the real decay, "no decay date within the horizon" counted as a miss (`window_coverage_test.py`, `window_coverage_pool.py`).

- **In-sample (the 17 Molniya-type objects):** 100 % coverage, worst error 247 d. Not independent: the window rule was chosen on these.
- **Out-of-sample (40 decayed, i >= 46.4 deg, e >= 0.3 objects from the 253-object pool, none of the 17; B from DISCOS/SATCAT RCS, 3 default):**

| lead (yr) | n | decay date produced | median abs err (d) | 90th pct | max | coverage |
|---|---|---|---|---|---|---|
| 0.5 | 36 | 31/36 | 30 | 78 | 112 | 0.86 |
| 1.0 | 39 | 34/39 | 33 | 127 | 779 | 0.82 |
| 2.0 | 38 | 34/38 | 39 | 160 | 676 | 0.82 |

**Reading.** The median error is as good as in-sample, but 11-14 % of objects that really decayed within 2 yr get no forecast decay date (the model keeps the perigee above the decay threshold), and a few forecasts that do produce a date miss by ~2 yr. Where a date is produced the window holds the truth in 91-100 % of cases; the shortfall is mainly the no-date objects. Consequence for the dashboard: a forecast without a decay date is shown blank and never as "long-lived"; the caption and Definitions quote both tests. Not investigated: why the no-date objects keep their perigee up (second-hump behaviour as in 8844/9269/11073, an extreme B, or the static atmosphere).

**Status (2026-10-04).** The secular + drag forecast is a testing and verification tool only. A Watchlist "approach window" built on it was implemented, published for one day and then withdrawn at the user's direction: it is not deployed, and no forecast-derived dates are shown to dashboard users. This code stays here in `scratch_rpe/combined_resonance/`.

## 15. KSROP full-force check of the secular + drag forecast, all test cases (2026-10-06)

Question: how much of the forecast's error is the secular/averaged model, as opposed to the starting elements, B and atmosphere? Every test case of sections 13-14 (mid-life cutoff 17, first TLE 17, leads 0.25-4 yr on the 17 Molniya-type objects 83, leads 0.5/1/2 yr on the 40 pool objects 113: **230 cases, 57 objects**) was re-run with KSROP's `driver_KS` (KS-regularised, v2.11.2 sources, fresh ifx build) from the identical start elements (`ksrop_check_cases.py`; the secular dates recomputed there reproduce the stored errors exactly). KSROP set-up: EGM2008 degree 20, Sun (degree 2) and Moon (degree 3), drag with OREM's static `ATM.DAT` (no `SW-All.csv`/`ATM2D.DAT`, so the same atmosphere family as the secular model), co-rotation on, no SRP, B = Cd A/m as before; the TLE mean elements are used as osculating elements at the start epoch (both models do the same). Decay = the driver's own stop at altitude < 80 km (`ksrop_check_run.py`, `ksrop_check_analyze.py`, `ksrop_check_*.csv`).

| set | n | KSROP date | secular date | abs err vs real, KSROP (med / p90 / max, d) | abs err vs real, secular | abs(secular - KSROP), both dated (med / p90 / max, d) |
|---|---|---|---|---|---|---|
| mid-life | 17 | 17 | 17 | 45 / 133 / 239 | 30 / 101 / 207 | 16 / 37 / 42 |
| first TLE | 17 | 17 | 17 | 157 / 768 / 1330 | 202 / 674 / 1220 | 70 / 302 / 478 |
| leads, 17 Molniya | 83 | 83 | 83 | 29 / 166 / 305 | 23 / 90 / 247 | 10 / 64 / 424 |
| leads, 40 pool | 113 | 96 | 99 | 29 / 112 / 776 | 33 / 117 / 779 | 16 / 40 / 240 |

**Reading.**
1. The secular + drag model reproduces the full-force KS trajectory's decay date well: over the 212 cases where both give a date the median difference is 14 d for horizons up to 2 yr (p90 42 d), bias -5 d overall; it degrades with horizon (first-TLE cases, 10-15 yr: median 70 d, up to 478 d). So the secular averaging is not what limits the forecast at the horizons that matter operationally.
2. Full force is **not** closer to the real decay: KSROP is closer than the secular forecast in only 114 of 212 cases, and its errors are as large (med 29-45 d at short leads; 157 d from the first TLE). The error budget is therefore dominated by what both share: start elements (TLE mean elements), B, the static atmosphere, no solar-cycle activity.
3. The failures to give a date are the same objects: of the 17 cases where KSROP reaches no decay (13 no-decay, 4 diverged), 13 also have no secular date (objects 27906, 21935, 9892, 27963, 45349, 26463, ...). For 27906/21935/9892/27963/45349 the real decay was 0.5-2 yr away and neither model gets there, so the 11-14 % "no date" rate of section 14 is not an artefact of the averaging; it points at inputs (B, atmosphere, elements). The 4 diverged runs stop within 6-124 d of the real decay (steep final decay, KSROP's fixed per-revolution stepping), not counted as dates.
4. The largest secular-vs-KSROP gaps are first-TLE cases (10605, 11073, 9269, 7480) and a few short ones (10605 at 0.5 yr: KSROP +304 d, secular +15 d; 7480 at 0.5-1 yr: KSROP +140 d, secular +12 d), where the full-force run decays later than the secular one.

**Not attributed.** KSROP and the secular model differ in several terms at once (full geopotential incl. tesseral vs J2, full lunisolar vs quadrupole, co-rotating drag vs none, orbit averaging); no ablation arm was run, so the 129-424 d outliers are not assigned to a term. SRP and epoch-resolved solar activity were not used on either side.
