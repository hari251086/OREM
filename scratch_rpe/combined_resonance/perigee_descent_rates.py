"""Observed natural perigee-descent rates in the 253-object pool's TLE histories (issue #56, THEORY section 17).

Sets a model-independent threshold for the coverage tests' truth check: how fast can a high-eccentricity object's
perigee fall while it is still above 200 km and more than 60 d before decay? Rates are measured over 30-day spans
(median of the TLEs in a +-3 d window at each end, to suppress single-TLE noise). Writes perigee_descent_rates.csv."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record     # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from validate_drive import clean_rows                                                 # noqa: E402

RE = 6378.137
SPAN, HALF = 30.0, 3.0
rows = []
for o in build_object_list():
    try:
        rec = load_full_record(f"{BASE}/{o['tle_file']}")
    except FileNotFoundError:
        continue
    if len(rec) < 50:
        continue
    jd, a, e, inc, O, w = np.array(clean_rows(rec)).T
    if np.median(e) < 0.3:
        continue
    hp = a * (1 - e) - RE
    decay = cal2jd(*o['decay_ymd'])
    keep = jd < min(decay, jd[-1]) - 60.0
    jd, hp, a = jd[keep], hp[keep], a[keep]
    if len(jd) < 20:
        continue
    t = jd[0]
    best = 0.0; best_e = 0.0
    while t + SPAN <= jd[-1]:
        m0 = np.abs(jd - t) <= HALF
        m1 = np.abs(jd - (t + SPAN)) <= HALF
        if m0.sum() >= 2 and m1.sum() >= 2:
            h0, h1 = np.median(hp[m0]), np.median(hp[m1])
            if h0 > 200.0 and h1 > 200.0:
                best = max(best, (h0 - h1) / SPAN)
                # perigee change ~ a * de, so (dhp / a) per day is the size-independent rate
                best_e = max(best_e, (h0 - h1) / np.median(a[m0 | m1]) / SPAN)
        t += 10.0
    rows.append(dict(norad=o['norad'], name=o.get('name', ''), max_descent_km_per_day=best,
                     max_descent_per_day_over_a=best_e))
df = pd.DataFrame(rows).sort_values('max_descent_per_day_over_a', ascending=False)
df.to_csv('E:/GitHub/OREM/scratch_rpe/combined_resonance/perigee_descent_rates.csv', index=False)
r = df['max_descent_km_per_day']
print(f"objects: {len(df)}; per-object max 30-d perigee descent (km/day): median {r.median():.2f}, "
      f"p90 {r.quantile(.9):.2f}, p99 {r.quantile(.99):.2f}, max {r.max():.2f}")
q = df['max_descent_per_day_over_a']
print(f"per-object max of (dhp/a)/day: median {q.median():.2e}, p90 {q.quantile(.9):.2e}, p99 {q.quantile(.99):.2e}, max {q.max():.2e}")
print(df.head(8).to_string(index=False))
