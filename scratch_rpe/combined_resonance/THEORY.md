# Combined resonant-azimuth theory for high-inclination HEO re-entry (draft v0, 2026-10-02)

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

- **[?] Closed form of ω̇′_P** from J2 + the perturber's plane precession, and the U-turn locus in (a, e, i) per
  perturber — the direct analogue of Wang & Gurfil's λ_crit(i). The family in §3 gives its Fourier content; the
  single-angle form should give the exact condition.
- **[?] Persistence-time detector** with windows of 2–4 years (R4 showed 180 days is too short), validated on
  decay timing.
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
