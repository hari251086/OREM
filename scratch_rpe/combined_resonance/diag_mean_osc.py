"""TLE mean vs osculating elements at the forecast start (issue #56, THEORY section 19).

The secular forecast starts from TLE MEAN elements, and the perigee it is compared with is the mean-element a(1-e).
For these eccentric orbits SGP4/SDP4's periodic terms can make the perigee the object actually reaches differ from
that. For every case of window_coverage_pool.csv (forecast start t0 = decay - lead, B as in that test):
  dhp      - median over the start TLEs of (minimum SGP4 altitude over one revolution) - (mean-element perigee altitude);
  variant A - start state from the median TLE mean elements (as the coverage test);
  variant B - same, but a and e set from the median SGP4 minimum / maximum radius over one revolution.
Both variants use the same TLEs, read with a loader that keeps every record (prepare_plot_data_fullpool.load_full_record
drops every second record of 3-line files). Writes diag_mean_osc.csv."""
import math
import os
import sys

import numpy as np
import pandas as pd
from sgp4.api import Satrec

sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                                              # noqa: E402
from heowatch.object_info import get_object_params                                   # noqa: E402
from prepare_plot_data_fullpool import BASE, build_object_list                       # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from secular_drag_prop import propagate_drag                                          # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
RE, MU = 6378.137, 398600.4418
cfg = load_config()


def load_tles(path):
    """Every (line1, line2) pair, 2- or 3-line files alike, sorted by epoch; with mean elements."""
    lines = [l.rstrip('\n') for l in open(path)]
    out = []
    for k in range(len(lines) - 1):
        l1, l2 = lines[k], lines[k + 1]
        if l1.startswith('1 ') and l2.startswith('2 ') and l1[2:7] == l2[2:7]:
            sat = Satrec.twoline2rv(l1, l2)
            jd = sat.jdsatepoch + sat.jdsatepochF
            n = float(l2[52:63]) * 2 * math.pi / 86400.0
            a = (MU / n ** 2) ** (1 / 3)
            out.append((jd, a, sat.ecco, math.degrees(sat.inclo), math.degrees(sat.nodeo), math.degrees(sat.argpo), sat))
    out.sort(key=lambda r: r[0])
    return out


def osc_extremes(sat):
    """(min, max) geocentric radius over one revolution from the TLE epoch, by SGP4 (60 s grid, then 1 s refine)."""
    P = 2 * math.pi / sat.no_kozai * 60.0                       # s (no_kozai in rad/min)
    t = np.arange(0.0, P + 60.0, 60.0)
    jd0 = sat.jdsatepoch; fr0 = sat.jdsatepochF + t / 86400.0
    e, r, _ = sat.sgp4_array(np.full_like(fr0, jd0), fr0)
    if np.any(e != 0):
        return float('nan'), float('nan')
    rr = np.linalg.norm(r, axis=1)
    res = []
    for k in (int(np.argmin(rr)), int(np.argmax(rr))):
        tt = np.arange(t[k] - 60.0, t[k] + 60.0, 1.0)
        e2, r2, _ = sat.sgp4_array(np.full_like(tt, jd0), sat.jdsatepochF + tt / 86400.0)
        res.append(np.linalg.norm(r2, axis=1))
    return float(res[0].min()), float(res[1].max())


def circ_mean(deg):
    return float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(deg))))) % 360)


objects = {o['norad']: o for o in build_object_list()}
pool = pd.read_csv(f'{D}/window_coverage_pool.csv')
rows = []
for n, g in pool.groupby('norad'):
    o = objects[n]
    tl = load_tles(f"{BASE}/{o['tle_file']}")
    jd = np.array([r[0] for r in tl])
    decay = cal2jd(*o['decay_ymd'])
    try:
        p = get_object_params(int(n), cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6
    except Exception:
        B0 = 2.2 / 56.0 * 1e-6
    for L in g['lead_yr']:
        t0 = decay - L * 365.25
        i1 = int(np.searchsorted(jd, t0, side='right'))
        if i1 < 12 or jd[i1 - 1] < t0 - 120:
            continue
        win = tl[i1 - 11:i1]
        a_m = np.median([r[1] for r in win]); e_m = np.median([r[2] for r in win])
        i_m = np.median([r[3] for r in win]); O_m = circ_mean([r[4] for r in win]); w_m = circ_mean([r[5] for r in win])
        ext = [osc_extremes(r[6]) for r in win]
        rmin = np.nanmedian([x[0] for x in ext]); rmax = np.nanmedian([x[1] for x in ext])
        hp_mean = a_m * (1 - e_m) - RE
        dhp = float(np.nanmedian([x[0] - r[1] * (1 - r[2]) for x, r in zip(ext, win)]))
        a_o = (rmin + rmax) / 2; e_o = (rmax - rmin) / (rmax + rmin)
        jd0 = float(jd[i1 - 1])
        res = {}
        for tag, st in (('A', (a_m, e_m, i_m, O_m, w_m)), ('B', (a_o, e_o, i_m, O_m, w_m))):
            dj, _ = propagate_drag(*st, jd0, B0, years=L + 4.0)
            res[tag] = (dj - decay) if dj == dj else np.nan
        rows.append(dict(norad=n, lead_yr=L, B_src=g['B_src'].iloc[0], hp_mean=hp_mean, hp_osc=rmin - RE, dhp=dhp,
                         ha_mean=a_m * (1 + e_m) - RE, ha_osc=rmax - RE, err_A=res['A'], err_B=res['B'],
                         dated_A=res['A'] == res['A'], dated_B=res['B'] == res['B'],
                         published_err=float(g.loc[g['lead_yr'] == L, 'err_d'].iloc[0])))
        print(rows[-1], flush=True)
out = pd.DataFrame(rows)
out.to_csv(f'{D}/diag_mean_osc.csv', index=False)
print(out.round(1).to_string(index=False))
