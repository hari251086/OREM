"""Issue #56: what lowers the perigee of the no-date debris objects? (THEORY section 22)

For every object of window_coverage_pool.csv, over the last year of tracking before the catalogued decay:
  - split the observed mean-element perigee change dhp = d[a(1-e)] into a size part (1-e)*da (energy loss: drag) and a
    shape part -a*de (eccentricity forcing at fixed size: lunisolar, SRP);
  - BN implied by the TLEs' own BSTAR (BC = 1/(12.741621 B*) kg/m^2, B* in 1/earth-radii; Vallado), and the area-to-
    mass ratio that implies (Cd = 2.2);
  - the largest perigee change solar radiation pressure could produce in that time at that area-to-mass ratio
    (Gauss, |de/dt| <= 1.5 F sqrt(a(1-e^2)/mu), F = P_sun Cr A/m, Cr = 1.2) -- an upper bound, Sun in the orbit plane.
No-date objects (section 16/20) are compared with the objects that do get a forecast date. Analysis only; writes
diag_nodate_decompose.csv."""
import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, parse_epoch_jd  # noqa: E402
from resonance_campaign_sweep import cal2jd  # noqa: E402
from tle_lines import bstar, tle_pairs  # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
MU, RE = 398600.4418, 6378.137
P_SUN, CR, CD = 4.56e-6, 1.2, 2.2          # N/m^2, -, -
EXCLUDE = {9892, 26463}                     # wrong decay date / manoeuvred (sections 16-17)


def rows_for(path):
    out = []
    for _, l1, l2 in tle_pairs(open(path, encoding='utf-8', errors='replace').read().splitlines()):
        try:
            n = float(l2[52:63]) * 2 * math.pi / 86400.0
            out.append((parse_epoch_jd(l1), (MU / n ** 2) ** (1 / 3), float('0.' + l2[26:33].strip()), bstar(l1)))
        except ValueError:
            continue
    return pd.DataFrame(out, columns=['jd', 'a', 'e', 'bstar']).sort_values('jd').reset_index(drop=True)


pool = pd.read_csv(f'{D}/window_coverage_pool.csv')
nodate = set(pool.loc[pool['err_d'].isna(), 'norad']) - EXCLUDE
objs = {o['norad']: o for o in build_object_list()}
res = []
for norad in sorted(set(pool['norad']) - EXCLUDE):
    o = objs.get(norad)
    if o is None:
        continue
    d = rows_for(f"{BASE}/{o['tle_file']}")
    decay = cal2jd(*o['decay_ymd'])
    w = d[(d.jd >= decay - 365.25) & (d.jd < decay)]
    if len(w) < 20:
        continue
    # robust endpoints: medians of the first / last 10 TLEs of the window, after dropping wild elements
    w = w[(w.e > 0) & (w.e < 1) & (w.a > RE)]
    head, tail = w.iloc[:10], w.iloc[-10:]
    a0, e0, a1, e1 = head.a.median(), head.e.median(), tail.a.median(), tail.e.median()
    T = tail.jd.median() - head.jd.median()
    abar, ebar = 0.5 * (a0 + a1), 0.5 * (e0 + e1)
    dhp = a1 * (1 - e1) - a0 * (1 - e0)
    dhp_size, dhp_shape = (1 - ebar) * (a1 - a0), -abar * (e1 - e0)
    bs = w.bstar[w.bstar > 0].median()
    bn = 1.0 / (12.741621 * bs) if bs == bs and bs > 0 else float('nan')
    am = 1.0 / (CD * bn) if bn == bn else float('nan')
    F = P_SUN * CR * am / 1000.0            # km/s^2
    dedt = 1.5 * F * math.sqrt(abar * (1 - ebar ** 2) / MU) if am == am else float('nan')
    dhp_srp_max = abar * dedt * T * 86400.0
    res.append(dict(norad=norad, group='no date' if norad in nodate else 'dated',
                    window_d=round(T), hp0=round(a0 * (1 - e0) - RE), hp1=round(a1 * (1 - e1) - RE),
                    dhp=round(dhp), dhp_size=round(dhp_size), dhp_shape=round(dhp_shape),
                    da=round(a1 - a0), bstar=bs, bn_bstar=bn, am_bstar=am, dhp_srp_max=dhp_srp_max))

r = pd.DataFrame(res)
r.to_csv(f'{D}/diag_nodate_decompose.csv', index=False)
pd.set_option('display.width', 200)
cols = ['norad', 'window_d', 'hp0', 'hp1', 'dhp', 'dhp_size', 'dhp_shape', 'da', 'bn_bstar', 'am_bstar', 'dhp_srp_max']
print('NO-DATE objects (last year before decay; km, kg/m^2, m^2/kg):')
print(r[r.group == 'no date'][cols].round(3).to_string(index=False))
g = r[r.group == 'dated']
print(f'\nDATED objects (n={len(g)}), medians:')
print(g[cols[1:]].median().round(3).to_string())
for lab, s in (('no date', r[r.group == 'no date']), ('dated', g)):
    frac = (s.dhp_shape.abs() / (s.dhp_size.abs() + s.dhp_shape.abs())).median()
    print(f'{lab:8s}: median share of |dhp| from eccentricity (shape) = {frac:.2f}; median BN from BSTAR = '
          f'{s.bn_bstar.median():.1f} kg/m^2')

# ---- part 2: observed versus modelled rates at each forecast start (window_coverage_pool_bstar.csv cases) ----------
# Over the 120 d after the start, fit da/dt and de/dt to the TLE mean elements and compare them with the secular
# model's own rates at the start state: drag (secular_drag.drag_rates) at the BSTAR B the forecast used, and J2 +
# lunisolar (secular_propagator.rates); drag also under Harris-Priester. SRP changes e but not a (orbit-averaged, shadow neglected), so a shortfall in
# da/dt is missing drag, while a da/dt match with a de/dt shortfall points to a missing eccentricity force.
from prepare_plot_data_fullpool import load_full_record  # noqa: E402
from secular_drag import drag_rates, set_atmosphere  # noqa: E402
from secular_propagator import rates  # noqa: E402
from validate_drive import clean_rows  # noqa: E402

FIT_D = 120.0
pb = pd.read_csv(f'{D}/window_coverage_pool_bstar.csv')
win = np.maximum(250.0, 0.15 * pb['fc_horizon_d'])
pb['case'] = np.where(pb['err_d'].isna(), 'no date', np.where(pb['err_d'].abs() <= win, 'in window', 'outside'))
out2 = []
for norad, g in pb.groupby('norad'):
    o = objs[norad]
    arr = np.array(clean_rows(load_full_record(f"{BASE}/{o['tle_file']}")))
    jd, a, e, inc, O, w = arr.T
    decay = cal2jd(*o['decay_ymd'])
    for _, c in g.iterrows():
        t0 = decay - c['lead_yr'] * 365.25
        k0 = (jd > t0 - 30) & (jd <= t0)
        k = (jd > t0) & (jd <= t0 + FIT_D)
        if k0.sum() < 3 or k.sum() < 8:
            continue
        a0, e0 = np.median(a[k0]), np.median(e[k0])
        i0, O0, w0 = (math.radians(float(np.median(v[k0]))) for v in (inc, O, w))
        da_obs = np.polyfit(jd[k] - t0, a[k], 1)[0]                       # km/day
        de_obs = np.polyfit(jd[k] - t0, e[k], 1)[0]                       # 1/day
        B = 1e-6 / c['bn']
        da_d, de_d = (v * 86400.0 for v in drag_rates(a0, e0, B))
        set_atmosphere('hp')
        da_hp = drag_rates(a0, e0, B)[0] * 86400.0
        set_atmosphere('atm_dat')
        de_ls = rates(np.array([e0, i0, O0, w0]), t0, a0, (True, True, True))[0] * 86400.0
        out2.append(dict(norad=norad, lead_yr=c['lead_yr'], case=c['case'], hp0=round(a0 * (1 - e0) - RE),
                         bn=c['bn'], da_obs=da_obs, da_model=da_d, da_ratio=da_obs / da_d if da_d else float('nan'),
                         da_model_hp=da_hp, da_ratio_hp=da_obs / da_hp if da_hp else float('nan'),
                         de_obs=de_obs, de_model=de_ls + de_d, de_lunisolar=de_ls))
r2 = pd.DataFrame(out2)
r2.to_csv(f'{D}/diag_nodate_rates.csv', index=False)
print(f'\nRATES over the {FIT_D:.0f} d after each forecast start (bstar B): da km/day, de 1e-5/day')
show = r2.copy()
for v in ('de_obs', 'de_model', 'de_lunisolar'):
    show[v] = show[v] * 1e5
print(show[show.case != 'in window'].round(3).to_string(index=False))
for lab, s in r2.groupby('case'):
    print(f'{lab:9s} n={len(s):2d}: median da_obs/da_model = {s.da_ratio.median():6.2f} (ATM.DAT), '
          f'{s.da_ratio_hp.median():6.2f} (Harris-Priester);  median (de_obs - de_model) = '
          f'{(s.de_obs - s.de_model).median() * 1e5:+.2f}e-5/day')

# ---- part 3: solar activity at the forecast start versus over the forecast (tlefit/HP cases) -----------------------
# A fitted B carries the density of its fit window forward; if activity then rises, the drag the forecast needs is
# larger than the fit saw.  F10.7: 81-day centred observed flux from OREM's input/SW-All.csv (as diag_nodate.py).
sw = pd.read_csv('E:/GitHub/OREM/input/SW-All.csv', usecols=['DATE', 'F10.7_OBS_CENTER81'])
sw['jd'] = [cal2jd(int(s[:4]), int(s[5:7]), int(s[8:10])) for s in sw['DATE']]
SWJ, SWF = sw['jd'].to_numpy(), sw['F10.7_OBS_CENTER81'].to_numpy()


def f107(j0, j1):
    k = (SWJ >= j0) & (SWJ <= j1)
    return float(np.nanmean(SWF[k])) if k.any() else float('nan')


pt = pd.read_csv(f'{D}/window_coverage_pool_tlefit_hp.csv')
pt = pt[~pt.norad.isin(EXCLUDE | {26463})]
win = np.maximum(250.0, 0.15 * pt['fc_horizon_d'])
pt['case'] = np.where(pt['err_d'].isna(), 'no date', np.where(pt['err_d'].abs() <= win, 'in window', 'outside'))
rows3 = []
for _, c in pt.iterrows():
    decay = cal2jd(*objs[c['norad']]['decay_ymd'])
    t0 = decay - c['lead_yr'] * 365.25
    rows3.append(dict(norad=c['norad'], lead_yr=c['lead_yr'], case=c['case'], f_fit=f107(t0 - 180.0, t0),
                      f_fc=f107(t0, decay)))
r3 = pd.DataFrame(rows3)
r3['rise'] = r3['f_fc'] / r3['f_fit']
r3.to_csv(f'{D}/diag_nodate_f107.csv', index=False)
print('\nF10.7 (81-d centred) in the 180-d B-fit window vs start-to-decay, tlefit/HP cases:')
print(r3[r3.case != 'in window'].round(2).to_string(index=False))
for lab, s in r3.groupby('case'):
    print(f'{lab:9s} n={len(s):3d}: median F10.7 fit window {s.f_fit.median():5.0f}, forecast {s.f_fc.median():5.0f}, '
          f'median ratio {s.rise.median():.2f}')
