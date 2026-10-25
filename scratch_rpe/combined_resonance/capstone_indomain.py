"""Section 21 capstone restricted to the Wang & Gurfil domain (i < 46.4 deg): mean |ensemble RPE| of each object's last
zone, full-history campaign (rpe_campaign.csv) versus post-resonance-only campaign (rpe_campaign_resonance16.csv), over
the objects predicting in both.  Without the restriction this reproduces section 21's 12.0 % -> 15.4 %, 5/8 worse."""
import sys

import numpy as np
import pandas as pd

sys.path[:0] = ['E:/GitHub/OREM/scratch_rpe', 'E:/GitHub/OREM/scratch_rpe/resonance_plots']
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record  # noqa: E402

R = 'E:/GitHub/OREM/scratch_rpe'


def last_zone_abs_ens(path):
    d = pd.read_csv(path, skipinitialspace=True)
    d = d[d.ens_rpe_pct.notna() & (d.ens_rpe_pct.abs() < 1e5)]
    return d.sort_values('zone').groupby('norad').ens_rpe_pct.last().abs()


full, post = last_zone_abs_ens(f'{R}/rpe_campaign.csv'), last_zone_abs_ens(f'{R}/rpe_campaign_resonance16.csv')
objs = {o['norad']: o for o in build_object_list()}
both = sorted(set(full.index) & set(post.index))
inc = {n: float(np.mean([r[3] for r in load_full_record(f"{BASE}/{objs[n]['tle_file']}")])) for n in both}
for lab, ns in (('all', both), ('i < 46.4 deg', [n for n in both if inc[n] < 46.4])):
    print(f'{lab:13s} n={len(ns)}: {full[ns].mean():.1f} % -> {post[ns].mean():.1f} %, '
          f'{int((post[ns] > full[ns]).sum())}/{len(ns)} worse')
    for n in ns:
        print(f'   {n:6d}  i={inc[n]:5.1f}  {full[n]:6.1f} -> {post[n]:6.1f}')
