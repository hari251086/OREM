"""Why do some decayed high-inclination objects get no forecast decay date? (issue #56, THEORY section 16)

For every case of window_coverage_pool.csv (forecast started 0.5 / 1 / 2 yr before the real decay), tests the four
candidate causes with the same start state and model as the coverage test:
  B        - multiplier on B needed to reach decay at the real date (m_needed), and a B multiplier FITTED on the year of
             TLE history before the forecast start (m_fit; uses no information after the start), then the forecast with it.
  atmosphere - mean F10.7 (observed, 81-day centred) over [start, real decay] from OREM's input/SW-All.csv.
  elements / lunisolar timing - observed TLE perigee between start and decay vs the model's perigee at the same epochs.
  data   - last TLE epoch vs the catalogued decay date.
Analysis only; writes diag_nodate.csv."""
import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                                              # noqa: E402
from heowatch.object_info import get_object_params                                   # noqa: E402
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record     # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from validate_drive import clean_rows                                                 # noqa: E402
from secular_drag_prop import propagate_drag                                          # noqa: E402
from secular_propagator import RE                                                     # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
cfg = load_config()
pool = pd.read_csv(f'{D}/window_coverage_pool.csv')
nodate_objs = set(pool.loc[pool['err_d'].isna(), 'norad'])
# Every no-date case, plus all cases of the same objects and every lead-1.0 dated case as the comparison group.
cases = pool[(pool['norad'].isin(nodate_objs)) | (pool['lead_yr'] == 1.0)].copy().reset_index(drop=True)
# Optional sharding for up to 4 parallel processes (GitHub/CLAUDE.md section 1): diag_nodate.py <k> <n>
SHARD = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) == 3 else (0, 1)
cases = cases[cases.index % SHARD[1] == SHARD[0]]

sw = pd.read_csv('E:/GitHub/OREM/input/SW-All.csv', usecols=['DATE', 'F10.7_OBS_CENTER81'])
sw['jd'] = [cal2jd(*map(int, d.split('-'))) for d in sw['DATE']]


def f107_mean(jd0, jd1):
    m = (sw['jd'] >= jd0) & (sw['jd'] <= jd1)
    return float(sw.loc[m, 'F10.7_OBS_CENTER81'].mean())


def start_state(jd, a, e, inc, O, w, t0):
    i1 = int(np.searchsorted(jd, t0, side='right'))
    if i1 < 12 or jd[i1 - 1] < t0 - 120:
        return None, None
    sl = slice(i1 - 11, i1)
    st = (float(np.median(a[sl])), float(np.median(e[sl])), float(np.median(inc[sl])),
          float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360),
          float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360))
    return st, float(jd[i1 - 1])


def decay_with(st, jd0, B, years):
    dj, traj = propagate_drag(*st, jd0, B, years=years)
    return dj, traj


def m_needed(st, jd0, B0, decay, years):
    """Smallest-error B multiplier on a log grid, refined by bisection on the sign of (model - real) decay."""
    grid = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64]
    res = []
    for m in grid:
        dj, _ = decay_with(st, jd0, B0 * m, years)
        res.append((m, dj))
    dated = [(m, dj) for m, dj in res if dj == dj]
    if not dated:
        return float('nan'), float('nan')
    m_first = dated[0][0]
    # bracket: late (dj > decay or no date) at lower m, early at higher m
    lo, hi = None, None
    for m, dj in res:
        if dj != dj or dj > decay:
            lo = m
        elif hi is None:
            hi = m
    if lo is None or hi is None:
        return m_first, float('nan')
    for _ in range(8):
        mid = math.sqrt(lo * hi)
        dj, _ = decay_with(st, jd0, B0 * mid, years)
        if dj != dj or dj > decay:
            lo = mid
        else:
            hi = mid
    return m_first, math.sqrt(lo * hi)


def m_fit(jd, a, e, inc, O, w, t0, B0):
    """B multiplier fitted on the year BEFORE t0: model a(t0) from the state at t0 - 365 d vs the observed a(t0)."""
    st_prev, jd_prev = start_state(jd, a, e, inc, O, w, t0 - 365.25)
    st_now, jd_now = start_state(jd, a, e, inc, O, w, t0)
    if st_prev is None or st_now is None:
        return float('nan'), float('nan')
    a_obs = st_now[0]
    da_obs = st_prev[0] - a_obs

    def a_model(m):
        dj, traj = propagate_drag(*st_prev, jd_prev, B0 * m, years=(jd_now - jd_prev) / 365.25 + 1e-3)
        if dj == dj or len(traj) == 0:
            return -1.0                       # decayed inside the fit window: far too much drag
        return float(np.interp(jd_now, traj[:, 0], traj[:, 1]))

    if da_obs < 5.0:                          # less than 5 km of semi-major-axis decay in the year: B not identifiable
        return float('nan'), da_obs
    lo, hi = 0.01, 1000.0
    for _ in range(18):
        mid = math.sqrt(lo * hi)
        if a_model(mid) > a_obs:              # model decays too little: more drag
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi), da_obs


objects = {o['norad']: o for o in build_object_list()}

if __name__ == '__main__':
    rows = []
    for _, c in cases.iterrows():
        n = int(c['norad']); L = float(c['lead_yr']); o = objects[n]
        rec = load_full_record(f"{BASE}/{o['tle_file']}")
        arr = np.array(clean_rows(rec)); jd, a, e, inc, O, w = arr.T
        decay = cal2jd(*o['decay_ymd'])
        try:
            p = get_object_params(n, cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6; src = p.source
            BN = p.mass_kg / (p.cd * p.area_m2)
        except Exception:
            B0 = 2.2 / 56.0 * 1e-6; src = 'default'; BN = 56.0 / 2.2 * 2.2
        t0 = decay - L * 365.25
        st, jd0 = start_state(jd, a, e, inc, O, w, t0)
        if st is None:
            continue
        years = L + 4.0
        dj1, traj = decay_with(st, jd0, B0, years)
        # observed vs model perigee over [start, decay]
        m_obs = (jd > jd0) & (jd <= decay)
        hp_obs = a[m_obs] * (1 - e[m_obs]) - RE
        hp_mod = np.interp(jd[m_obs], traj[:, 0], traj[:, 3]) if len(traj) else np.full(m_obs.sum(), np.nan)
        first, best = m_needed(st, jd0, B0, decay, years)
        mf, da_obs = m_fit(jd, a, e, inc, O, w, t0, B0)
        dj_fit = float('nan')
        if mf == mf:
            dj_fit, _ = decay_with(st, jd0, B0 * mf, years)
        rows.append(dict(
            norad=n, lead_yr=L, B_src=src, BN_kg_m2=BN, inc=st[2], e0=st[1], hp0=st[0] * (1 - st[1]) - RE,
            dated=dj1 == dj1, err_d=(dj1 - decay) if dj1 == dj1 else np.nan,
            m_first_dated=first, m_needed=best, m_fit=mf, da_obs_prev_yr_km=da_obs,
            err_fit_d=(dj_fit - decay) if dj_fit == dj_fit else np.nan, dated_fit=dj_fit == dj_fit,
            f107_mean=f107_mean(jd0, decay),
            hp_obs_min=float(hp_obs.min()) if hp_obs.size else np.nan,
            hp_mod_min=float(np.nanmin(hp_mod)) if hp_mod.size else np.nan,
            hp_bias_med=float(np.nanmedian(hp_mod - hp_obs)) if hp_obs.size else np.nan,
            n_tle_window=int(m_obs.sum()), last_tle_minus_decay_d=float(jd[-1] - decay),
        ))
        print(n, L, 'dated' if dj1 == dj1 else 'NO DATE', f"m_needed={best:.3g} m_fit={mf:.3g}", flush=True)

    out = pd.DataFrame(rows)
    suffix = '' if SHARD[1] == 1 else f'_part{SHARD[0]}'
    out.to_csv(f'{D}/diag_nodate{suffix}.csv', index=False)
    print(out.to_string(max_cols=30))
