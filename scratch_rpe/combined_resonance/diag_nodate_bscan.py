"""Follow-up to diag_nodate.py (issue #56, THEORY section 16): for the no-date cases whose B came from the SATCAT RCS fallback
(category area + an assumed 500 kg mass), the B multiplier needed to reach decay at the real date, on a grid extended to x4096,
and the ballistic number BN = m/(Cd A) it implies. Writes diag_nodate_bscan.csv."""
import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
import diag_nodate as dn  # noqa: E402  (module-level run is skipped below via the CASES override)

D = dn.D
diag = pd.read_csv(f'{D}/diag_nodate.csv')
targets = diag[(~diag['dated']) & (diag['B_src'] == 'satcat_rcs')]
rows = []
for _, c in targets.iterrows():
    n, L = int(c['norad']), float(c['lead_yr'])
    o = dn.objects[n]
    rec = dn.load_full_record(f"{dn.BASE}/{o['tle_file']}")
    arr = np.array(dn.clean_rows(rec)); jd, a, e, inc, O, w = arr.T
    decay = dn.cal2jd(*o['decay_ymd'])
    p = dn.get_object_params(n, dn.cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6; BN0 = p.mass_kg / (p.cd * p.area_m2)
    st, jd0 = dn.start_state(jd, a, e, inc, O, w, decay - L * 365.25)
    years = L + 4.0
    first = None; lo = hi = None
    for m in [2.0 ** k for k in range(0, 13)]:                     # x1 ... x4096
        dj, _ = dn.decay_with(st, jd0, B0 * m, years)
        if dj == dj and first is None:
            first = m
        if dj != dj or dj > decay:
            lo = m
        elif hi is None:
            hi = m
            break
    m_best = float('nan')
    if lo is not None and hi is not None:
        for _ in range(8):
            mid = math.sqrt(lo * hi)
            dj, _ = dn.decay_with(st, jd0, B0 * mid, years)
            if dj != dj or dj > decay:
                lo = mid
            else:
                hi = mid
        m_best = math.sqrt(lo * hi)
    rows.append(dict(norad=n, lead_yr=L, BN_assumed=BN0, m_first_dated=first, m_needed=m_best,
                     BN_implied=(BN0 / m_best) if m_best == m_best else float('nan')))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f'{D}/diag_nodate_bscan.csv', index=False)
