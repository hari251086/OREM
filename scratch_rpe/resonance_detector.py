"""Hardened Wang & Gurfil (2017) solar-apsidal-resonance detector (issue #56).

Replaces the naive per-TLE regime flip used by resonance_campaign_sweep.py /
prepare_plot_data_fullpool.py, whose manual review (253 objects, 2026-09-25)
found 15 of 27 flags to be false positives. Fixes, one per observed failure
mode:

  1. DOMAIN GATE -- the resonance exists only for i < 46.4 deg (where
     5cos^2i - 2cosi - 1 > 0); higher-inclination orbits are out of domain,
     not "non-resonant". Reported as OUT_OF_DOMAIN, never RESONANT.
  2. MIN_TLES -- records too sparse (interpolated / a few weeks) to show any
     secular behaviour are INSUFFICIENT_DATA.
  3. a(t) OUTLIER FILTER -- single-point TLE glitches (e.g. a: 27000 -> 7000
     -> 27000 km) are dropped against a rolling median before the test.
  4. PERSISTENCE -- a regime change counts only if the new regime holds for
     >= PERSIST_PTS consecutive points spanning >= PERSIST_DAYS, so chatter
     around lambda_crit is not a crossing.
  5. TERMINAL GUARD -- a crossing within the last TERMINAL_PTS TLEs is the
     a/e terminal-collapse degeneracy, not a distinct U-turn; reported as
     TERMINAL_COLLAPSE. (Point-count only: a days criterion wrongly rejected
     the manually confirmed 59347, whose crossing is 28 days / 21 TLEs
     before the last TLE.)
  6. BAD_DATA -- a record with >BAD_PERIGEE_FRAC of points at physically
     impossible negative perigee height (a(1-e) < R_E) is unreliable
     (1970s-era noise); reported as BAD_DATA. Motivated by 5980 only -- one
     example, so treat the threshold as untested elsewhere.

Pure functions on (jd, a_km, e, i_deg) rows; reuses the validated lambda_crit
from resonance_campaign_sweep.py.
"""
import math
import sys

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
from resonance_campaign_sweep import resonance_lambda_crit  # noqa: E402

I_CRIT_DEG = 46.4        # upper inclination bound of the paper's regime
MIN_TLES = 30
OUTLIER_WIN = 9          # rolling-median window (odd)
OUTLIER_REL = 0.05       # |a - median| / median above this => glitch
PERSIST_PTS = 5
PERSIST_DAYS = 10.0
TERMINAL_PTS = 10
BAD_PERIGEE_FRAC = 0.02
R_EARTH = 6378.1363


def _median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def filter_outliers(rows):
    """Drop points whose a deviates >OUTLIER_REL from the centred rolling
    median of a. rows: (jd, a, e, i)."""
    n = len(rows)
    if n < OUTLIER_WIN:
        return list(rows)
    h = OUTLIER_WIN // 2
    keep = []
    for k, r in enumerate(rows):
        lo, hi = max(0, k - h), min(n, k + h + 1)
        med = _median([x[1] for x in rows[lo:hi]])
        if abs(r[1] - med) / med <= OUTLIER_REL:
            keep.append(r)
    return keep


def persistent_crossings(rows, lam_crit):
    """Indices (into rows) where the PRE/POST regime changes AND the new
    regime persists for >= PERSIST_PTS points spanning >= PERSIST_DAYS."""
    regimes = []
    for jd, a, e, i in rows:
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        regimes.append(lam > lam_crit)          # True = PRE
    out = []
    n = len(rows)
    k = 1
    while k < n:
        if regimes[k] != regimes[k - 1]:
            j = k
            while j < n and regimes[j] == regimes[k]:
                j += 1
            run_pts = j - k
            run_days = rows[j - 1][0] - rows[k][0]
            ends_record = (j == n)
            # a trailing run is allowed to be short only in points/days
            # terms if it reaches the end of the record -- the terminal
            # guard handles that case explicitly below.
            if (run_pts >= PERSIST_PTS and run_days >= PERSIST_DAYS) \
                    or ends_record:
                out.append(k)
            k = j
        else:
            k += 1
    return out


def detect(rows):
    """rows: list of (jd, a_km, e, i_deg) sorted by jd.

    Returns dict(status, resonance_jd, crossing_idx, n_used, i_mean).
    status in {OUT_OF_DOMAIN, INSUFFICIENT_DATA, NOT_RESONANT,
               BAD_DATA, TERMINAL_COLLAPSE, RESONANT}.
    """
    if not rows:
        return {'status': 'INSUFFICIENT_DATA', 'n_used': 0}
    i_mean = sum(r[3] for r in rows) / len(rows)
    if i_mean >= I_CRIT_DEG:
        return {'status': 'OUT_OF_DOMAIN', 'i_mean': i_mean,
                'n_used': len(rows)}
    if len(rows) < MIN_TLES:
        return {'status': 'INSUFFICIENT_DATA', 'i_mean': i_mean,
                'n_used': len(rows)}
    n_neg = sum(1 for r in rows if r[1] * (1.0 - r[2]) < R_EARTH)
    if n_neg / len(rows) > BAD_PERIGEE_FRAC:
        return {'status': 'BAD_DATA', 'i_mean': i_mean, 'n_used': len(rows)}
    clean = filter_outliers(rows)
    if len(clean) < MIN_TLES:
        return {'status': 'INSUFFICIENT_DATA', 'i_mean': i_mean,
                'n_used': len(clean)}
    lam_crit = resonance_lambda_crit(math.radians(i_mean))
    idx = persistent_crossings(clean, lam_crit)
    if not idx:
        return {'status': 'NOT_RESONANT', 'i_mean': i_mean,
                'n_used': len(clean)}
    k = idx[-1]
    jd = clean[k][0]
    n_after = len(clean) - 1 - k
    if n_after < TERMINAL_PTS:
        return {'status': 'TERMINAL_COLLAPSE', 'resonance_jd': jd,
                'crossing_idx': k, 'i_mean': i_mean, 'n_used': len(clean)}
    return {'status': 'RESONANT', 'resonance_jd': jd, 'crossing_idx': k,
            'i_mean': i_mean, 'n_used': len(clean)}
