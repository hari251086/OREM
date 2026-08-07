# DRAMA-based re-entry validation

Uses ESA's DRAMA suite (specifically OSCAR, its orbital lifetime/decay tool)
as an independent cross-check on OREM's own re-entry-date predictions --
the same role `scratch_gmat/` plays for orbit-propagation cross-validation
against GMAT.

## Setup

DRAMA itself and its Java prerequisite are installed outside this repo (they're
large third-party binaries, not something to git-track):

- **Java 17** (Eclipse Temurin, portable): `E:\Java\jdk-17`. Only needed for
  SARA-RISK; OSCAR and ARES in DRAMA 4.1.4 are native executables and don't
  need Java.
- **DRAMA 4.1.4**: `E:\DRAMA` (installed unattended via
  `setup_DRAMA-4.1.4_windows.exe --mode unattended --prefix E:\DRAMA`).

pyDRAMA (ESA's official Python wrapper, package name `drama`) is installed
into a dedicated venv at `drama/.venv` from the copy bundled inside the
DRAMA install:

```
py -3.13 -m venv drama/.venv
drama\.venv\Scripts\pip install E:\DRAMA\TOOLS\drama_python_package
drama\.venv\Scripts\pip install sgp4   # for TLE parsing in scripts/tle_utils.py
```

## Layout

- `scripts/tle_utils.py` -- parses an OREM `input/example_*.tle.txt` file's
  first TLE entry into the classical elements DRAMA's tools expect.
- `scripts/run_oscar_reentry.py` -- runs an OSCAR natural-decay analysis for
  one OREM object and writes the result to `output/`. Working example:
  `drama\.venv\Scripts\python drama\scripts\run_oscar_reentry.py 21670`.
- `scripts/run_ares.py`, `scripts/run_sara.py` -- scaffolded stubs for
  DRAMA's collision-probability (ARES) and re-entry casualty-risk (SARA)
  tools. Installed and importable, not yet exercised against real data --
  they answer different questions than OREM's re-entry-date predictions.
- `scripts/compare_orem_drama.py` -- stub for the future full comparison
  campaign (OREM's `scratch_rpe/rpe_campaign.csv` predictions vs. OSCAR's).
- `input/` -- DRAMA-format input/config files, kept separate from the
  repo's own `input/` (which is OREM/TLE-format).
- `output/` -- OSCAR run results (JSON); only `.gitkeep` is tracked, actual
  run output is gitignored, same convention as the repo's own `output/`.

## Why NORAD 21670 for the working example

It's an object OREM's own RPE campaign already fit
(`scratch_rpe/rpe_campaign.csv`, zone 1): `e_opt=0.739988` matches the TLE's
own eccentricity (0.7399072) almost exactly, and OREM's fitted ballistic
number `BN_opt=10.246 kg/m^2` is what `run_oscar_reentry.py` uses to derive
OSCAR's spacecraft mass (holding drag coefficient and cross-section area at
OSCAR's own defaults) -- so the OSCAR run is traceable to what OREM already
believes about this object rather than an arbitrary guess. OSCAR predicted
a 0.699-year lifetime from the Aug 1991 epoch (re-entry ~ Apr/May 1992),
consistent with the object's real TLE history showing decay through early
1992.

## Status

Scaffolding + one verified working OSCAR run only. A full multi-object
validation campaign (comparable in scope to the KSROP<->GMAT campaign) is a
deliberate follow-on, not done here.

## Campaign results (OREM vs. DRAMA/OSCAR, 53 comparable objects)

Full run of `run_oscar_campaign.py` + `compare_orem_drama.py` against all
92 objects in `scratch_rpe/rpe_campaign.csv`:

| metric | OREM | DRAMA/OSCAR |
|---|---|---|
| median \|RPE%\| | 31.89 | 41.29 |
| mean \|RPE%\| | 49.34 | 236.22 |

Median is the fairer headline comparison -- OSCAR's mean is dragged out by
7 extreme outliers (NORAD [39802, 41679, 41686, 41695, 44482, 45036, 46429], all with
\|orem_rpe - oscar_rpe\| > 500%: OSCAR predicts a multi-decade lifetime
where the object actually decayed within ~1-1.5 years). Two hypotheses were tested and
ruled out: an eccentricity mismatch (OREM's own fitted eccentricity for the same zone is
close to the raw TLE eccentricity used here) and OSCAR's solar/geomagnetic activity
scenario choice (tested directly across all 4 of OSCAR's scenario options for NORAD
39802 -- lifetime only varies 27-37 years across all of them, nowhere near enough to
explain a 30x gap from the true ~1.35-year lifetime). What the outliers actually share:
**6/7** sit within 3&deg; of the Molniya
critical inclination (63.4&deg;, inclination range 25.9-64.9&deg;),
perigee altitude 203-385 km --
and their NORAD IDs cluster tightly (41679/41686/41695), consistent with fragments of the
same debris-generating breakup event. The most likely explanation is that OREM's own
ballistic-number fit for these small-sample, late-zone objects doesn't transfer as a
physically portable quantity into OSCAR's own density/decay model for this specific
regime -- not root-caused further than that here (would need per-fragment attitude/shape
modeling or a direct atmosphere-model diff to close out).

DRAMA/OSCAR's independent prediction landed closer to the true observed
re-entry date than OREM's own primary-zone estimate for
**20/53** objects.

31 objects had no OREM primary-zone estimate to
compare against (DRAMA/OSCAR was still run on these -- see
`oscar_campaign_results.csv`); 17 objects showed
no OSCAR re-entry within the 100-year propagation window;
4 OSCAR runs errored.

This campaign used the "simpler" orbit-state method (raw TLE semi-major
axis/eccentricity; only OREM's fitted ballistic number feeds in, via mass)
rather than perigee-preserving state matching -- see `drama/README.md`'s
setup notes. That choice trades some physical consistency for simplicity
and is a plausible secondary contributor to the outliers above.

Full per-object results: `output/oscar_campaign_results.csv`.
