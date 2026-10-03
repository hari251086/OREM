"""Cross-check: parameter-free secular forecast (first TLE only) vs OREM's validation predictions, Molniya-type objects.

OREM (validation campaign, heowatch.validation.compute_validation_rpe): hindcast cutoff ~90 d before decay, full multi-zone RSM/GA.
Secular forecast (secular_propagator.py): J2+Sun+Moon quadrupole, a const, NO drag, started from the first TLE (~10 yr ahead);
it yields the epoch E300 when the perigee altitude falls below 300 km.  To compare in decay-date terms the forecast E300 is
converted with the leave-one-out median gap (decay - E300) of the other objects -- one empirical constant, never the object's own.
"""
import sys, json, math
import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                          # noqa: E402
from heowatch.validation import compute_validation_rpe          # noqa: E402
from heowatch.watchlist_db import init_db                       # noqa: E402
from prepare_plot_data_fullpool import load_full_record         # noqa: E402
from validate_drive import clean_rows                            # noqa: E402
from secular_propagator import propagate, RE                     # noqa: E402
from molniya_cycle import NORADS, ROOT, smooth                   # noqa: E402

THR = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0     # perigee-altitude threshold (km) defining the event
cfg = load_config()
data_root = __import__('pathlib').Path(cfg['data_root'])
conn = init_db(data_root / 'predictions' / 'watchlist.db')
satcat = pd.read_csv(data_root / 'catalog' / 'latest' / 'satcat.csv')
val = compute_validation_rpe(conn, satcat)
print(f'=== event: perigee altitude < {THR:.0f} km after the perigee peak ===')
val = val.set_index('norad_id')


def jd_of(iso):
    import datetime as dt
    return (dt.date.fromisoformat(str(iso)[:10]) - dt.date(1970, 1, 1)).days + 2440587.5


rows = []
for n in NORADS:
    if n not in val.index:
        print(n, 'not in validation set'); continue
    r = val.loc[n]
    rec = np.array(clean_rows(load_full_record(f"{ROOT}/{n}/{n}.tle.txt"))); jd, a, e, inc, O, w = rec.T
    hp = a * (1 - e) - RE
    s = smooth(hp); ipk_obs = int(np.argmax(s))
    idx = np.where((np.arange(len(hp)) > ipk_obs) & (s < THR))[0]       # first crossing AFTER the perigee peak
    if len(idx) == 0:
        print(n, 'record never below 300 km after the peak'); continue
    e300_obs = jd[idx[0]]
    k = 11
    O0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[:k]))))) % 360)
    w0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[:k]))))) % 360)
    e0, i0, aa = float(np.median(e[:k])), float(np.median(inc[:k])), float(np.median(a[:k]))
    traj = propagate(e0, i0, O0, w0, aa, jd[0], years=22.0, dt_days=10.0)
    hpp = aa * (1 - traj[:, 1]) - RE
    sp = smooth(hpp); ipk_p = int(np.argmax(sp[: int(len(sp) * 0.9)]))
    hit = np.where((np.arange(len(hpp)) > ipk_p) & (sp < THR))[0]       # same definition on the model trajectory
    e300_pred = traj[hit[0], 0] if len(hit) else float('nan')
    decay = jd_of(r['actual_utc']); orem = jd_of(r['predicted_utc']); cutoff = jd_of(str(r['run_date'])[:4] + '-' + str(r['run_date'])[4:6] + '-' + str(r['run_date'])[6:8]) if len(str(r['run_date'])) == 8 else float('nan')
    rows.append(dict(norad=n, first_tle=float(jd[0]), decay=decay, orem_pred=orem, cutoff=cutoff, rpe_orem=float(r['rpe_pct']),
                     e300_obs=float(e300_obs), e300_pred=float(e300_pred)))
df = pd.DataFrame(rows).dropna(subset=['e300_pred'])
df['gap'] = df['decay'] - df['e300_obs']
# leave-one-out median gap
gaps = df['gap'].to_numpy()
df['gap_loo'] = [np.median(np.delete(gaps, i)) for i in range(len(gaps))]
df['fc_decay'] = df['e300_pred'] + df['gap_loo']
df['err_orem_d'] = df['orem_pred'] - df['decay']
df['err_fc_d'] = df['fc_decay'] - df['decay']
df['span_yr'] = (df['decay'] - df['first_tle']) / 365.25
df['lead_orem_d'] = df['decay'] - df['cutoff']
df['err_orem_pct_span'] = 100 * df['err_orem_d'] / (df['decay'] - df['first_tle'])
df['err_fc_pct_span'] = 100 * df['err_fc_d'] / (df['decay'] - df['first_tle'])
df.to_csv(f'E:/GitHub/OREM/scratch_rpe/combined_resonance/cross_check_orem_{int(THR)}.csv', index=False)

pd.set_option('display.width', 200)
print(df[['norad', 'span_yr', 'gap', 'err_orem_d', 'err_fc_d', 'rpe_orem', 'err_fc_pct_span']].round(2).to_string(index=False))
print(f"\nobjects: {len(df)}   observed gap decay-E300: median {df['gap'].median():.0f} d, IQR [{df['gap'].quantile(.25):.0f}, {df['gap'].quantile(.75):.0f}]")
for lab, c in (('OREM (cutoff ~%d d before decay)' % df['lead_orem_d'].median(), 'err_orem_d'), ('secular forecast (first TLE, ~%.1f yr ahead)' % df['span_yr'].median(), 'err_fc_d')):
    x = df[c]
    print(f"{lab:<46} median |err| {x.abs().median():6.1f} d   mean |err| {x.abs().mean():6.1f} d   bias {x.median():+6.1f} d   max |err| {x.abs().max():6.1f} d")
for lab, c in (('OREM', 'err_orem_pct_span'), ('forecast', 'err_fc_pct_span')):
    x = df[c]
    print(f"{lab:<9} error as % of (decay - first TLE): median |{x.abs().median():.2f}|%  mean |{x.abs().mean():.2f}|%")
print(f"corr(OREM error, forecast error) = {np.corrcoef(df['err_orem_d'], df['err_fc_d'])[0,1]:+.2f};  corr(|errors|) = {np.corrcoef(df['err_orem_d'].abs(), df['err_fc_d'].abs())[0,1]:+.2f}")
print(f"OREM-vs-forecast predicted decay difference: median |{(df['orem_pred'] - df['fc_decay']).abs().median():.0f}| d")
