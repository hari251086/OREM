"""Compare KSROP full-force decay dates with the secular+drag forecast and the real decay, per case family.  Missing decay dates count as misses."""
import sys
import numpy as np
import pandas as pd
D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
arm = sys.argv[1] if len(sys.argv) > 1 else 'full'
c = pd.read_csv(f'{D}/ksrop_check_cases.csv'); r = pd.read_csv(f'{D}/ksrop_check_results_{arm}.csv')
d = c.merge(r, on='case')
d['ks_decay_jd'] = np.where(d['status'] == 'decay', d['ksrop_end_jd'], np.nan)
d['ks_err_d'] = d['ks_decay_jd'] - d['real_decay_jd']
d['sec_minus_ks_d'] = d['sec_decay_jd'] - d['ks_decay_jd']
print('status counts:', d['status'].value_counts().to_dict(), '| wall s total', round(d['wall_s'].sum()))
def st(x):
    x = x.dropna().abs()
    return f"{len(x):3d} med {x.median():6.0f} p90 {np.percentile(x,90) if len(x) else np.nan:6.0f} max {x.max() if len(x) else np.nan:6.0f}"
groups = {'mid-life (17)': d.family == 'mid', 'first TLE (17)': d.family == 'full', 'lead, 17 Molniya': d.family == 'lead', 'lead, 40 pool': d.family == 'pool'}
print(f"{'set':18s} | n  | KSROP date | secular date | |KS-real| d            | |secular-real| d       | |secular-KS| d (both have date)")
for k, m in groups.items():
    g = d[m]; both = g.dropna(subset=['ks_decay_jd', 'sec_decay_jd'])
    print(f"{k:18s} | {len(g):3d} | {g.ks_decay_jd.notna().sum():3d} | {g.sec_decay_jd.notna().sum():3d} | {st(g.ks_err_d)} | {st(g.sec_err_d)} | {st(both.sec_minus_ks_d)}")
both = d.dropna(subset=['ks_decay_jd', 'sec_decay_jd'])
print('\nall cases with both dates:', len(both), '| secular-KS bias (median, d):', round(both.sec_minus_ks_d.median()), '| |KS-real| <= |sec-real| in', int((both.ks_err_d.abs() <= both.sec_err_d.abs()).sum()), 'of', len(both))
print('\nlargest |secular - KSROP| cases:')
print(both.reindex(both.sec_minus_ks_d.abs().sort_values(ascending=False).index)[['case', 'real_horizon_d', 'ks_err_d', 'sec_err_d', 'sec_minus_ks_d']].head(12).round(0).to_string(index=False))
miss = d[d.ks_decay_jd.isna() | d.sec_decay_jd.isna()][['case', 'status', 'real_horizon_d', 'ks_err_d', 'sec_err_d']]
print('\ncases where either model gave no decay date:'); print(miss.round(0).to_string(index=False))
d.to_csv(f'{D}/ksrop_check_merged_{arm}.csv', index=False)
