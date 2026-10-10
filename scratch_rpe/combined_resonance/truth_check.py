"""Validity check for a catalogued decay date against the object's own TLE record (issue #56, THEORY section 17).

A decay date is rejected as truth when the tracking contradicts it:
  1. tracked_after_decay - the last TLE is more than MAX_TRACKED_AFTER_D days AFTER the catalogued decay;
  2. final_orbit_too_high - from the last tracked orbit, reaching the 120 km perigee by the catalogued date would need a
     perigee descent faster than MAX_DESCENT_PER_A, normalised by the semi-major axis ((dhp / a) per day; perigee
     change ~ a * de, so this is size-independent).
MAX_DESCENT_PER_A is 2 x the 99th percentile of the fastest 30-day descent observed in each pool object's own history
while its perigee was above 200 km and > 60 d before decay (perigee_descent_rates.py: p99 7.9e-4 / day). The threshold
comes from data, not from the forecast model, so objects the model cannot follow (e.g. high area-to-mass debris) are not
rejected for that reason."""
import numpy as np

RE = 6378.137
MAX_TRACKED_AFTER_D = 30.0
MAX_DESCENT_PER_A = 1.6e-3      # per day; 2 x the pool's p99 (7.92e-4)
HP_REENTRY_KM = 120.0


def check_decay_truth(jd, a, e, decay_jd):
    """(ok, reason, detail) for TLE arrays jd, a (km), e and a catalogued decay date (JD)."""
    after = float(jd[-1] - decay_jd)
    if after > MAX_TRACKED_AFTER_D:
        return False, 'tracked_after_decay', f'last TLE {after:.0f} d after the catalogued decay'
    gap = float(decay_jd - jd[-1])
    if gap <= 0.0:
        return True, '', ''
    k = slice(max(0, len(jd) - 3), len(jd))                 # median of the last 3 TLEs: robust to one bad element set
    a_last = float(np.median(a[k]))
    hp_last = float(np.median(a[k] * (1.0 - e[k]))) - RE
    need = (hp_last - HP_REENTRY_KM) / a_last / gap
    if need > MAX_DESCENT_PER_A:
        return False, 'final_orbit_too_high', (f'last perigee {hp_last:.0f} km, {gap:.0f} d before decay: needs '
                                               f'{need:.2e}/day > {MAX_DESCENT_PER_A:.1e}')
    return True, '', ''
