"""Operational test for an 'approach window': start the drag forecast at fixed lead times before the real decay and ask how wide the
window must be.  17 Molniya-type objects (DISCOS B).  For lead L in (0.25, 0.5, 1, 2, 4 yr): elements = median of the 11 TLEs ending at decay - L;
forecast decay (drag, a const-free) vs the real decay.  Candidate half-width rules are scored by coverage."""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                      # noqa: E402
from heowatch.object_info import get_object_params           # noqa: E402
from prepare_plot_data_fullpool import load_full_record      # noqa: E402
from validate_drive import clean_rows                         # noqa: E402
from secular_drag_prop import propagate_drag                  # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
cfg = load_config()
cut = pd.read_csv(f'{D}/midlife_cutoffs.csv')
LEADS = (0.25, 0.5, 1.0, 2.0, 4.0)
rows = []
for _, r in cut.iterrows():
    n = int(r['norad']); decay = float(r['decay_jd'])
    p = get_object_params(n, cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6
    rec = np.array(clean_rows(load_full_record(f"{ROOT}/{n}/{n}.tle.txt"))); jd, a, e, inc, O, w = rec.T
    for L in LEADS:
        t0 = decay - L * 365.25
        i1 = int(np.searchsorted(jd, t0, side='right'))
        if i1 < 12:
            continue
        if jd[i1 - 1] < t0 - 120:                # the (truncated) record ends too early for this lead
            continue
        sl = slice(i1 - 11, i1)
        st = (float(np.median(a[sl])), float(np.median(e[sl])), float(np.median(inc[sl])),
              float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360),
              float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360))
        hp0 = st[0] * (1 - st[1]) - 6378.137
        dj, _ = propagate_drag(*st, jd[i1 - 1], B0, years=20.0)
        horizon = decay - jd[i1 - 1]
        rows.append(dict(norad=n, lead_yr=L, horizon_d=horizon, hp0=hp0, err_d=(dj - decay), fc_horizon_d=dj - jd[i1 - 1]))
    print(n, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(f'{D}/window_coverage.csv', index=False)
print('\nlead (yr) |  n | median |err| d | 90th pct | max |err| | coverage of window half-width rules (err within +-w)')
for L, g in df.groupby('lead_yr'):
    e = g['err_d'].abs()
    hz = g['fc_horizon_d']
    cov = {}
    for lab, fn in (('max(120,0.10*H)', lambda H: np.maximum(120, 0.10 * H)), ('max(180,0.12*H)', lambda H: np.maximum(180, 0.12 * H)),
                    ('max(250,0.15*H)', lambda H: np.maximum(250, 0.15 * H)), ('+-250 d', lambda H: 250 + 0 * H)):
        cov[lab] = float(np.mean(e.to_numpy() <= fn(hz.to_numpy())))
    print(f"{L:8.2f}  | {len(g):2d} | {e.median():10.0f} | {np.percentile(e, 90):8.0f} | {e.max():9.0f} | " + '  '.join(f"{k}: {v:.2f}" for k, v in cov.items()))
