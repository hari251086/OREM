"""Out-of-sample window-coverage test: the drag forecast started 0.5 / 1 / 2 yr before the real decay on every decayed, high-inclination
(i >= 46.4 deg), strongly eccentric (e >= 0.3) object of the 253-object pool that was NOT one of the 17 Molniya-type validation objects.
B from the Watchlist object parameters when available (else a typical BN of 56 kg/m^2, flagged).  'No decay within the horizon' counts as a miss."""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                                              # noqa: E402
from heowatch.object_info import get_object_params                                   # noqa: E402
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record     # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from validate_drive import clean_rows                                                 # noqa: E402
from secular_drag_prop import propagate_drag                                          # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
cfg = load_config()
seen = set(pd.read_csv(f'{D}/midlife_cutoffs.csv')['norad'])
LEADS = (0.5, 1.0, 2.0)
rows = []
for o in build_object_list():
    n = o['norad']
    if n in seen:
        continue
    try:
        rec = load_full_record(f"{BASE}/{o['tle_file']}")
    except FileNotFoundError:
        continue
    if len(rec) < 300:
        continue
    inc_mean = sum(r[3] for r in rec) / len(rec)
    arr = np.array(clean_rows(rec)); jd, a, e, inc, O, w = arr.T
    if inc_mean < 46.4 or np.median(e[:20]) < 0.3:
        continue
    decay = cal2jd(*o['decay_ymd'])
    if decay < jd[0] or decay - jd[-1] > 400:           # record must reach (nearly) to the decay
        continue
    try:
        p = get_object_params(n, cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6; src = p.source
    except Exception:
        B0 = 2.2 * 5.0 / 900.0 * 1e-6 * 0 + 2.2 / 56.0 * 1e-6; src = 'default BN 56'      # 2.2 / BN
    for L in LEADS:
        t0 = decay - L * 365.25
        i1 = int(np.searchsorted(jd, t0, side='right'))
        if i1 < 12 or jd[i1 - 1] < t0 - 120 or a[i1 - 11:i1].min() <= 0:
            continue
        sl = slice(i1 - 11, i1)
        st = (float(np.median(a[sl])), float(np.median(e[sl])), float(np.median(inc[sl])),
              float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360),
              float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360))
        if st[1] < 0.3 or st[2] < 46.4:
            continue
        dj, _ = propagate_drag(*st, jd[i1 - 1], B0, years=L + 4.0)
        H = dj - jd[i1 - 1] if dj == dj else float('nan')
        rows.append(dict(norad=n, lead_yr=L, i=st[2], e=st[1], hp0=st[0] * (1 - st[1]) - 6378.137, B_src=src,
                         fc_horizon_d=H, err_d=(dj - decay) if dj == dj else float('nan')))
    print(n, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(f'{D}/window_coverage_pool.csv', index=False)
print(f"\nobjects: {df['norad'].nunique()}  (B source: {df['B_src'].value_counts().to_dict()})")
print('lead (yr) |  n | forecast gave a decay date | median |err| d | 90th pct | max |err| | coverage of max(250 d, 0.15 H) (no date = miss)')
for L, g in df.groupby('lead_yr'):
    has = g['err_d'].notna()
    e_ = g.loc[has, 'err_d'].abs()
    cov = float(np.mean(has & (g['err_d'].abs() <= np.maximum(250.0, 0.15 * g['fc_horizon_d']))))
    print(f"{L:8.2f}  | {len(g):2d} | {int(has.sum()):2d}/{len(g)} | {e_.median():10.0f} | {np.percentile(e_, 90):8.0f} | {e_.max():9.0f} | {cov:.2f}")
