"""No-date debris cases re-run with the BALLISTIC NUMBER ESTIMATED BY OREM (issue #56, THEORY section 18), instead of the
Watchlist object_info value (BN 758 from the SATCAT-RCS fallback) that sections 14/16 prescribed.

For each case (forecast start t0 = decay - lead): BN from OREM's own fit in rpe_campaign.csv, taking the latest zone whose
fit window (zone start + 20 d, the campaign's max zone length) ended before t0 -- information available at t0 only;
OREM-valid zones (zstat 0) preferred, else the latest zone of any status (flagged). Same start state and secular +
drag model as the coverage test. Writes diag_nodate_orembn.csv."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
import diag_nodate as dn  # noqa: E402

ZONE_DAYS = 20.0
camp = pd.read_csv('E:/GitHub/OREM/scratch_rpe/rpe_campaign.csv', skipinitialspace=True)
pool = pd.read_csv(f'{dn.D}/window_coverage_pool.csv')
cases = pool[pool['err_d'].isna()][['norad', 'lead_yr']]
rows = []
for _, c in cases.iterrows():
    n, L = int(c['norad']), float(c['lead_yr'])
    o = dn.objects[n]
    decay = dn.cal2jd(*o['decay_ymd'])
    t0 = decay - L * 365.25
    z = camp[(camp['norad'] == n) & (camp['zepoch'] + ZONE_DAYS <= t0)]
    if z.empty:
        rows.append(dict(norad=n, lead_yr=L, note='no OREM zone before the forecast start')); print(n, L, 'no OREM zone'); continue
    valid = z[z['zstat'] == 0]
    pick = (valid if not valid.empty else z).sort_values('zepoch').iloc[-1]
    bn = float(pick['bn_opt'])
    prm = dn.get_object_params(n, dn.cfg); bn_info = prm.mass_kg / (prm.cd * prm.area_m2)
    jd, a, e, inc, O, w = np.array(dn.clean_rows(dn.load_full_record(f"{dn.BASE}/{o['tle_file']}"))).T
    st, jd0 = dn.start_state(jd, a, e, inc, O, w, t0)
    dj, traj = dn.decay_with(st, jd0, 1e-6 / bn, L + 4.0)        # B = Cd A / m in km^2/kg = 1e-6 / BN
    rows.append(dict(norad=n, lead_yr=L, zone=int(pick['zone']), zone_valid=bool(pick['zstat'] == 0),
                     zone_age_d=t0 - float(pick['zepoch']), bn_orem=bn, bn_object_info=bn_info, bn_source=prm.source,
                     scope='manoeuvred spacecraft (section 16)' if n == 26463 else '',
                     dated=dj == dj, err_d=(dj - decay) if dj == dj else np.nan,
                     model_hp_min=float(traj[:, 3].min()), t_hp_min_d=float(traj[np.argmin(traj[:, 3]), 0] - jd0),
                     real_decay_after_start_d=decay - jd0))
    print(rows[-1], flush=True)
out = pd.DataFrame(rows)
out.to_csv(f'{dn.D}/diag_nodate_orembn.csv', index=False)
print(out.round(1).to_string(index=False))
