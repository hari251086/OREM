"""TLE osculating vs mean apogee for the 97-object campaign (2026-10-12).

OREM's fit compares the model's osculating apogee (propagated from an SGP4-osculating zone-start state) with the
TLEs' mean apogee. This measures how far apart those two quantities are in the TLE data itself:
  mean apogee       a_mean(1 + e_mean) - R,  a_mean from the TLE mean motion (Kozai n, as tle_evolve does)
  osculating apogee a(1 + e) - R from the SGP4 state at the TLE's own epoch
for every TLE of every campaign object, and separately at each zone's first TLE (the fit's starting point).
Outputs: osc_vs_mean_apogee_objects.csv, osc_vs_mean_apogee_zones.csv, osc_vs_mean_apogee.png."""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sgp4.api import Satrec, WGS72  # noqa: E402

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from tle_lines import tle_pairs  # noqa: E402

MU, RE = 398600.4415, 6378.1363          # OREM input/const_new.dat
SL12 = {41679, 41686, 41688, 41695, 44482, 44555, 45493, 45494, 45495, 45499, 43685, 44892}
BLUE, ORANGE = '#2a78d6', '#eb6834'      # reference palette slots 1-2 (light)
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def tle_rows(norad):
    path = HERE.parent / 'input' / f'example_{norad}.tle.txt'
    if not path.exists():
        return None
    rows = []
    for _, l1, l2 in tle_pairs(path.read_text(encoding='utf-8', errors='replace').splitlines()):
        s = Satrec.twoline2rv(l1, l2, WGS72)
        err, r, v = s.sgp4(s.jdsatepoch, s.jdsatepochF)
        if err:
            continue
        r, v = np.array(r), np.array(v)
        rn, vn = np.linalg.norm(r), np.linalg.norm(v)
        a_osc = 1.0 / (2.0 / rn - vn * vn / MU)
        h = np.cross(r, v)
        e_osc = math.sqrt(max(0.0, 1.0 - np.dot(h, h) / (MU * a_osc)))
        n = s.no_kozai / 60.0                 # rad/s
        a_mean = (MU / n ** 2) ** (1.0 / 3.0)
        rows.append(dict(jd=s.jdsatepoch + s.jdsatepochF, e_mean=s.ecco,
                         ha_mean=a_mean * (1 + s.ecco) - RE, ha_osc=a_osc * (1 + e_osc) - RE,
                         hp_mean=a_mean * (1 - s.ecco) - RE, hp_osc=a_osc * (1 - e_osc) - RE))
    if not rows:
        return None
    d = pd.DataFrame(rows).sort_values('jd').reset_index(drop=True)
    d['d_ha'] = d.ha_osc - d.ha_mean
    d['d_hp'] = d.hp_osc - d.hp_mean
    return d


def main():
    camp = pd.read_csv(HERE / 'rpe_campaign.csv', skipinitialspace=True)
    camp = camp[pd.to_numeric(camp.zone, errors='coerce').notna()].copy()
    camp['zepoch'] = pd.to_numeric(camp.zepoch)
    objs, zones, series = [], [], {}
    for norad in sorted(camp.norad.unique()):
        d = tle_rows(norad)
        if d is None:
            continue
        series[norad] = d
        objs.append(dict(norad=norad, sl12=norad in SL12, n_tle=len(d), e_median=d.e_mean.median(),
                         d_ha_median=d.d_ha.median(), d_ha_p10=d.d_ha.quantile(0.1), d_ha_p90=d.d_ha.quantile(0.9),
                         d_hp_median=d.d_hp.median()))
        for _, z in camp[camp.norad == norad].iterrows():
            k = (d.jd - z.zepoch).abs().idxmin()
            if abs(d.jd[k] - z.zepoch) < 0.01:
                zones.append(dict(norad=norad, zone=int(z.zone), sl12=norad in SL12, e_mean=d.e_mean[k],
                                  hp_mean=d.hp_mean[k], d_ha=d.d_ha[k], d_hp=d.d_hp[k]))
    o, z = pd.DataFrame(objs), pd.DataFrame(zones)
    o.round(3).to_csv(HERE / 'osc_vs_mean_apogee_objects.csv', index=False)
    z.round(3).to_csv(HERE / 'osc_vs_mean_apogee_zones.csv', index=False)

    print(f'{len(o)} objects, {int(o.n_tle.sum())} TLEs, {len(z)} zone-start TLEs')
    for lab, sel in (('all', z), ('SL-12', z[z.sl12]), ('other', z[~z.sl12])):
        a = sel.d_ha
        print(f'  zone starts, {lab:5s} n={len(sel):3d}: osc-mean apogee median {a.median():+6.1f} km, '
              f'|.| median {a.abs().median():5.1f}, p10 {a.quantile(.1):+6.1f}, p90 {a.quantile(.9):+6.1f}, '
              f'|.|>20 km {(a.abs() > 20).mean():.0%}')
    for lo, hi in ((0, 0.3), (0.3, 0.5), (0.5, 0.65), (0.65, 1)):
        s = z[(z.e_mean >= lo) & (z.e_mean < hi)].d_ha
        print(f'  e in [{lo:.2f},{hi:.2f}): n={len(s):3d}  |osc-mean apogee| median {s.abs().median():5.1f} km')

    # ---- figure: one panel per question, no shared/dual axes ----
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2,
                         'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlecolor': INK})
    fig = plt.figure(figsize=(11, 7.2), facecolor='#fcfcfb')
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1], hspace=0.45, wspace=0.28)
    ax0 = fig.add_subplot(gs[0, :])
    other, sl = z[~z.sl12], z[z.sl12]
    ax0.scatter(other.e_mean, other.d_ha, s=22, color=BLUE, edgecolor='#fcfcfb', linewidth=0.8,
                label=f'other campaign objects ({len(other)} zones)', zorder=3)
    ax0.scatter(sl.e_mean, sl.d_ha, s=26, color=ORANGE, edgecolor='#fcfcfb', linewidth=0.8,
                label=f'SL-12 debris ({len(sl)} zones)', zorder=4)
    ax0.axhline(0, color=INK2, linewidth=0.8, zorder=2)
    ax0.set_xlabel('mean eccentricity at zone start')
    ax0.set_ylabel('osculating − mean apogee (km)')
    ax0.set_title('At each fit zone\'s first TLE: osculating apogee (SGP4 state) minus TLE mean apogee',
                  loc='left', fontsize=10)
    ax0.legend(frameon=False, loc='upper left')
    ax0.grid(True, color=GRID, linewidth=0.6)
    ax0.set_facecolor('#fcfcfb')
    for sp in ('top', 'right'):
        ax0.spines[sp].set_visible(False)

    for i, norad in enumerate((42928, 37151, 41688)):
        ax = fig.add_subplot(gs[1, i])
        d = series.get(norad)
        if d is None:
            continue
        t = (d.jd - d.jd.iloc[0]) / 365.25
        ax.plot(t, d.d_ha, color=ORANGE if norad in SL12 else BLUE, linewidth=0.8)
        ax.axhline(0, color=INK2, linewidth=0.8)
        ax.set_title(f'{norad}{" (SL-12 debris)" if norad in SL12 else ""}: every TLE', loc='left', fontsize=9)
        ax.set_xlabel('years since first TLE')
        if i == 0:
            ax.set_ylabel('osculating − mean apogee (km)')
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_facecolor('#fcfcfb')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    fig.savefig(HERE / 'osc_vs_mean_apogee.png', dpi=150, bbox_inches='tight', facecolor='#fcfcfb')
    print('wrote', HERE / 'osc_vs_mean_apogee.png')


if __name__ == '__main__':
    main()
