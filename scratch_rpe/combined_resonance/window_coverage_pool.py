"""Out-of-sample window-coverage test: the drag forecast started 0.5 / 1 / 2 yr before the real decay on every decayed, high-inclination
(i >= 46.4 deg), strongly eccentric (e >= 0.3) object of the 253-object pool that was NOT one of the 17 Molniya-type validation objects.
B from the Watchlist object parameters when available (else a typical BN of 56 kg/m^2, flagged).  'No decay within the horizon' counts as a miss.

--b-source bstar (section 22): B from the object's own TLEs instead -- BN = 1/(12.741621 B*) (SGP4 reference density
2.461e-8 kg/m^3 x earth radius = 0.15697 kg/m^2/ER, B* = rho0 B/2), median of the positive B* of the TLEs in the
year before the forecast start (causal), falling back to the object parameters when fewer than BSTAR_MIN_N exist.
Writes window_coverage_pool_bstar.csv.
--b-source tlefit (section 22): B fitted to the drag the TLEs show -- the linear da/dt of the mean semi-major axis over
the TLE_FIT_D days before the forecast start (causal), divided by the model's drag da/dt at unit B at that window's
median (a, e) under the chosen atmosphere; BN clamped to [BN_FIT_MIN, BN_FIT_MAX]; object parameters when the window
has fewer than 8 TLEs or a is not decreasing.  Writes window_coverage_pool_tlefit*.csv.
--atmosphere hp (section 22): Harris-Priester mean-activity density instead of the static ATM.DAT (secular_drag);
appends _hp to the output names."""
import argparse
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                                              # noqa: E402
from heowatch.object_info import get_object_params                                   # noqa: E402
from prepare_plot_data_fullpool import BASE, build_object_list, load_bstar, load_full_record  # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from validate_drive import clean_rows                                                 # noqa: E402
from secular_drag import drag_rates, set_atmosphere                                  # noqa: E402
from secular_drag_prop import propagate_drag                                          # noqa: E402
from truth_check import check_decay_truth                                             # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
ap = argparse.ArgumentParser()
ap.add_argument('--b-source', choices=('params', 'bstar', 'tlefit'), default='params')
ap.add_argument('--atmosphere', choices=('atm_dat', 'hp'), default='atm_dat')
args = ap.parse_args()
B_SOURCE = args.b_source
set_atmosphere(args.atmosphere)
OUT = 'window_coverage_pool' + ('' if B_SOURCE == 'params' else '_' + B_SOURCE) + ('_hp' if args.atmosphere == 'hp' else '')
BSTAR_MIN_N = 5
TLE_FIT_D, BN_FIT_MIN, BN_FIT_MAX = 180.0, 0.05, 2000.0
cfg = load_config()
seen = set(pd.read_csv(f'{D}/midlife_cutoffs.csv')['norad'])
LEADS = (0.5, 1.0, 2.0)
rows = []
rejected = []
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
    ok, why, detail = check_decay_truth(jd, a, e, decay)  # section 17: reject truth the tracking contradicts
    if not ok:
        rejected.append(dict(norad=n, name=o.get('name', ''), reason=why, detail=detail))
        print(n, 'REJECTED as truth:', why, detail, flush=True)
        continue
    try:
        p = get_object_params(n, cfg); B_par = p.cd * p.area_m2 / p.mass_kg * 1e-6; src_par = p.source
    except Exception:
        B_par = 2.2 * 5.0 / 900.0 * 1e-6 * 0 + 2.2 / 56.0 * 1e-6; src_par = 'default BN 56'      # 2.2 / BN
    bs = np.array(load_bstar(f"{BASE}/{o['tle_file']}")) if B_SOURCE == 'bstar' else None
    for L in LEADS:
        t0 = decay - L * 365.25
        B0, src = B_par, src_par
        if bs is not None:
            k = (bs[:, 0] > t0 - 365.25) & (bs[:, 0] <= t0) & (bs[:, 1] > 0)
            if k.sum() >= BSTAR_MIN_N:
                B0 = 12.741621 * float(np.median(bs[k, 1])) * 1e-6; src = 'bstar'     # Cd A/m = 1/BN (m^2/kg) -> km^2/kg
            else:
                src = src_par + ' (no bstar)'
        if B_SOURCE == 'tlefit':
            k = (jd > t0 - TLE_FIT_D) & (jd <= t0)
            src = src_par + ' (no fit)'
            if k.sum() >= 8:
                slope = np.polyfit(jd[k] - t0, a[k], 1)[0]                                  # km/day
                unit = drag_rates(float(np.median(a[k])), float(np.median(e[k])), 1e-6)[0] * 86400.0   # at BN 1
                if slope < 0.0 and unit < 0.0:
                    bn = min(max(unit / slope, BN_FIT_MIN), BN_FIT_MAX)
                    B0 = 1e-6 / bn; src = 'tlefit'
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
                         bn=1e-6 / B0,
                         fc_horizon_d=H, err_d=(dj - decay) if dj == dj else float('nan')))
    print(n, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(f'{D}/{OUT}.csv', index=False)
pd.DataFrame(rejected, columns=['norad', 'name', 'reason', 'detail']).to_csv(f'{D}/{OUT}_rejected.csv', index=False)
print(f"rejected as truth: {len(rejected)} object(s)", [r['norad'] for r in rejected])
print(f"\nobjects: {df['norad'].nunique()}  (B source: {df['B_src'].value_counts().to_dict()})")
print('lead (yr) |  n | forecast gave a decay date | median |err| d | 90th pct | max |err| | coverage of max(250 d, 0.15 H) (no date = miss)')
for L, g in df.groupby('lead_yr'):
    has = g['err_d'].notna()
    e_ = g.loc[has, 'err_d'].abs()
    cov = float(np.mean(has & (g['err_d'].abs() <= np.maximum(250.0, 0.15 * g['fc_horizon_d']))))
    print(f"{L:8.2f}  | {len(g):2d} | {int(has.sum()):2d}/{len(g)} | {e_.median():10.0f} | {np.percentile(e_, 90):8.0f} | {e_.max():9.0f} | {cov:.2f}")
