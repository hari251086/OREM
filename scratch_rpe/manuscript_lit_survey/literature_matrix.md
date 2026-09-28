# SOTA Literature Matrix — OREM Manuscript Introduction (part of issue #46)

Companion to `references.bib` and `intro_paragraphs.tex` in this folder.
Scope: orbital decay, atmospheric drag dissipation, and lifetime prediction
of RSOs from HEO/GTO — lunisolar third-body perturbations, solar apsidal
resonance, drag/lifetime prediction, and the regulatory/SDA context. This
is the broader background survey; the specific novelty paragraph
distinguishing OREM from Mutyalarao & Sharma / Sellamuthu & Sharma / Sharma
& Raj is already drafted separately in issue #46's own comment thread and
is not repeated here.

**Verification note**: every row below traces to either (a) a primary
source already read in full during this project's earlier literature
passes (`reference_orem_reentry_literature.md`,
`reference_wider_reentry_materials_survey.md`), or (b) a direct
publisher/DOI hit confirmed this session. Two rows are flagged where a
detail could not be independently confirmed — see their Notes column
rather than treating them as settled.

## a) Historical / Analytical Baselines

| Author(s) & Year | Journal / Venue | Primary Methodology / Model | Key Contribution / Focus Area |
|---|---|---|---|
| Shute (1965) | NASA TN X-643-65-402 | Analytical perturbation theory | Earliest dedicated treatment of highly-eccentric-satellite lifetime; identifies luni-solar terms as a major perturbation alongside drag |
| Shute & Chiville (1966) | Planetary and Space Science 14(4) | Analytical, luni-solar + drag | Quantifies the luni-solar effect on orbital lifetime specifically for highly eccentric orbits |
| Kozai (1962) | The Astronomical Journal 67 | Secular perturbation theory | The Kozai(-Lidov) mechanism — eccentricity/inclination coupling under third-body forcing; foundational to all later critical-inclination/resonance work |
| King-Hele (1987) | *Satellite Orbits in an Atmosphere* (Blackie) | Analytical air-drag contraction theory | Definitive monograph-form synthesis of King-Hele's own drag-contraction theory; **substituted for the requested "King-Hele 1982"** — no specific 1982 paper could be verified with confidence, see `references.bib` note |
| Sharma & Raj (1988) | Earth, Moon, and Planets 42 | KS uniformly-regular canonical elements + oblateness | Establishes KS-regularized elements (this work's own propagator lineage) as a long-term orbit-computation framework including Earth oblateness |

## b) Lunisolar Perturbations & Solar Apsidal Resonances

| Author(s) & Year | Journal / Venue | Primary Methodology / Model | Key Contribution / Focus Area |
|---|---|---|---|
| Hough (1981) | Celestial Mechanics 25 | Hamiltonian perturbation theory | Rigorous critical-inclination treatment showing lunisolar forcing controls equilibrium stability for a > 1.4 R⊕ (no drag) |
| Colombo (2015) | AAS 15-395 | Semi-analytical, drag-free | Independent (ESA/Politecnico di Milano) corroboration of the critical-inclination/Kozai-Lidov eccentricity-libration mechanism for HEO |
| Wang & Gurfil (2016) | Acta Astronautica 128 | Full-perturbation numerical (drag+oblateness+lunisolar) | Dynamical modeling/lifetime analysis establishing the GTO long-term-evolution baseline used by the 2017 resonance papers |
| Wang & Gurfil (2017a) | Advances in Space Research 59(8) | Full-perturbation numerical + resonance-condition derivation | Defines and characterizes *solar apsidal resonance* — a computable, per-object scalar criterion converting periodic perigee oscillation into monotonic collapse/rise |
| Wang & Gurfil (2017b) | Journal of Spacecraft and Rockets | Full-perturbation numerical | Deorbit-timing strategy exploiting the same luni-solar resonance for deliberate lifetime reduction |
| Luo & Wang (2019) | Acta Astronautica 165 | Full-perturbation numerical | Generalizes the solar-apsidal-resonance condition to inclined GTOs (Wang & Gurfil's original treatment was near-equatorial) |
| Sellamuthu & Sharma (2021) | Astrodynamics 5 | KS-regularized numerical | Applies KS-regularized (not Cartesian) dynamics to luni-solar gravity on RSOs — same propagator family as this work |

## c) Atmospheric Drag & Lifetime Prediction

| Author(s) & Year | Journal / Venue | Primary Methodology / Model | Key Contribution / Focus Area |
|---|---|---|---|
| Sharma (1997) | Proc. Royal Society of London A 453 | KS-regularized analytical, diurnal bulge | Exact analytical diurnal-density-bulge term for KS-regularized drag theory, oblate atmosphere |
| Swinerd & Boulton (1983) | Proc. Royal Society of London A 389(1796) | Analytical, near-circular oblate diurnal atmosphere | Real-satellite (tracking-data) validation of the diurnal-bulge density model later adopted by the Sharma KS lineage |
| Vallado & Finkleman (2014) | Acta Astronautica 95 | Comparative survey | Authoritative, standards-author assessment of satellite-drag/atmospheric-density-model accuracy limits across the field |
| Sharma, Bandyopadhyay & Adimurthy (2006) | IAC-06-B6.2.11 | Response-Surface Methodology + Genetic Algorithm | Originating RSM+GA cost-function formulation for fitting eccentricity/ballistic-coefficient to observed decay — this work's own methodological ancestor |
| Mutyalarao & Sharma (2010) | Journal of Spacecraft and Rockets 47(4) | RSM+GA, single-zone | Optimal reentry-time estimation for a GTO upper stage; short-term (<3-month) horizon |
| Mutyalarao & Sharma (2011) | Advances in Space Research 47 | RSM+GA, explicit decay "zones" | First explicit "zone"-based RSM+GA decay-fitting formulation — closest direct textual/methodological ancestor of this work's own zone-selection architecture |
| Gupta & Anilkumar (2015) | Journal of Spacecraft and Rockets 52(1) | Integrated multi-objective EBP fit | Independent, structurally similar ballistic-parameter fit from TLE apogee/perigee, validated on terminal-decay-phase real objects |
| Dutt et al. (2023) | Advances in Space Research 72 | Multi-method benchmark (RSM-GA + 3 sibling methods) | Finds no single re-entry-prediction method dominates across the full decay horizon — motivates multi-method/ensemble consideration |
| ISO/CD 27852 (2007) | ISO committee draft | Standards-body comparative review | Authoritative statement that atmospheric-density-model error, not propagator fidelity, dominates orbit-lifetime prediction error industry-wide |

## d) Space Debris Regulations & SDA Infrastructure

| Author(s) & Year | Journal / Venue | Primary Methodology / Model | Key Contribution / Focus Area |
|---|---|---|---|
| IADC (2021) | IADC-02-01, Rev. 3 | International guideline | Current international consensus space-debris-mitigation guidelines (post-mission disposal, collision-avoidance practice) |
| FCC (2022) | FCC 22-74 (37 FCC Rcd 11818) | U.S. federal regulation | Adopts the binding "5-year rule": LEO-crossing satellites must deorbit within 5 years of end-of-mission, replacing the legacy 25-year guideline for FCC-licensed operators |
| Federal Register (2024) | 2024-17093 | U.S. federal regulation (codification) | Formal Federal Register publication of the FCC's orbital-debris mitigation rule, confirming its regulatory force |

**Not included** (found during this search, left out rather than guessed
per the accuracy requirement): a 2019 Celestial Mechanics and Dynamical
Astronomy GTO dynamical-lifetime survey and a 2023 AIAA JSR statistical
long-term-lifetime-prediction paper (DOI `10.2514/1.A35706`) — both
plausible additions to category (c), but their author lists could not be
confirmed with confidence from available sources this session. Verify
directly against AIAA/Springer before adding.
