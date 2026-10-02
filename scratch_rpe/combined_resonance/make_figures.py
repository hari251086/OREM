"""Figures for THEORY.md (dataviz skill: validated categorical slots 1-3 + neutral ink, thin lines, direct labels,
legend for >=2 series, recessive grid, surface #fcfcfb)."""
import sys, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import load_full_record     # noqa: E402
from validate_drive import clean_rows                       # noqa: E402
from secular_propagator import propagate, RE                # noqa: E402

SURFACE, INK, INK2, INK3, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e4e3df'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'      # categorical slots 1-3 (validated)
OUT = 'E:/GitHub/OREM/scratch_rpe/combined_resonance/figures'
ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
import os; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'text.color': INK, 'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE})


def style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.tick_params(length=0)


def fig1(norad=7480):
    rows = np.array(clean_rows(load_full_record(f"{ROOT}/{norad}/{norad}.tle.txt"))); jd, a, e, inc, O, w = rows.T
    t = (jd - jd[0]) / 365.25
    hp, ha = a * (1 - e) - RE, a * (1 + e) - RE
    ser = json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/molniya_series.json'))[str(norad)]
    ipk = ser['ipk']; t_pk = ser['t'][ipk]
    cwrap = lambda x: float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(x[:11]))))) % 360)
    st = (float(np.median(e[:11])), float(np.median(inc[:11])), cwrap(O), cwrap(w), float(np.median(a[:11])))
    variants = [('J2 + Sun + Moon', (True, True, True), C1, '-'), ('J2 + Moon only', (True, False, True), C2, '-'),
                ('J2 + Sun only', (True, True, False), C3, '-'), ('J2 only', (True, False, False), INK3, '-')]
    fig, ax = plt.subplots(3, 1, figsize=(8.2, 9.4), sharex=True, gridspec_kw={'height_ratios': [1.05, 1.15, 0.8]})
    # (a) perigee & apogee, observed vs full model forecast from the first TLE
    traj = propagate(*st[:1], st[1], st[2], st[3], st[4], jd[0], years=12.0, dt_days=10.0)
    tp = (traj[:, 0] - jd[0]) / 365.25
    ax[0].plot(t, hp, color=INK, lw=2, label='observed perigee (TLE)')
    ax[0].plot(tp, st[4] * (1 - traj[:, 1]) - RE, color=C1, lw=2, ls='--', label='forecast from first TLE')
    ax[0].set_ylabel('Perigee altitude (km)'); ax[0].set_ylim(0, 2300)
    ax[0].axvline(t_pk, color=INK3, lw=1, ls=':')
    ax[0].text(t_pk + 0.1, 130, 'observed perigee peak', color=INK2, fontsize=8)
    ax2 = ax[0].twinx() if False else None
    ax[0].legend(loc='upper left', frameon=False, fontsize=8)
    ax[0].set_title(f'NORAD {norad} (MOLNIYA 1-28): forecast from the first TLE only',
                    fontsize=10, loc='left', color=INK)
    # (b) ablation on perigee: which perturber is needed
    ax[1].plot(t, hp, color=INK, lw=2, label='observed')
    for lab, use, col, ls in variants:
        tr = propagate(*st[:1], st[1], st[2], st[3], st[4], jd[0], years=12.0, dt_days=10.0, use=use)
        ax[1].plot((tr[:, 0] - jd[0]) / 365.25, st[4] * (1 - tr[:, 1]) - RE, color=col, lw=1.8, ls=ls, label=lab)
    ax[1].set_ylabel('Perigee altitude (km)'); ax[1].set_ylim(0, 2300)
    ax[1].legend(loc='upper left', frameon=False, fontsize=8, ncol=2)
    ax[1].set_title('Ablation: the Moon makes the rise and its peak; the Sun is needed for the fall', fontsize=9, loc='left', color=INK2)
    # (c) drive: Sun and Moon contributions to dq/dt
    ts = np.array(ser['t'])
    ax[2].axhline(0, color=INK3, lw=0.8)
    ax[2].plot(ts, ser['dm'], color=C2, lw=1.8, label='Moon: dq/dt')
    ax[2].plot(ts, ser['ds'], color=C3, lw=1.8, label='Sun: dq/dt')
    ax[2].set_ylabel('Quadrupole drive dq/dt (km/day)'); ax[2].set_xlabel('Years since first tracked TLE')
    ax[2].axvline(t_pk, color=INK3, lw=1, ls=':')
    ax[2].legend(loc='upper right', frameon=False, fontsize=8)
    ax[2].set_title('The drive changes sign at the perigee peak (azimuth 2w\'_Moon passes 180 deg)', fontsize=9, loc='left', color=INK2)
    for x in ax:
        style(x)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig1_molniya_{norad}_forecast.png', dpi=170); plt.close(fig)


def fig2():
    R = json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/forecast_results.json'))
    P = json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/forecast_pool.json'))
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.4), gridspec_kw={'width_ratios': [1.1, 1]})
    # left: predicted vs observed time to perigee < 300 km
    o = np.array([r['t_obs'] for r in P if math.isfinite(r['t_pred'])]); p = np.array([r['t_pred'] for r in P if math.isfinite(r['t_pred'])])
    inc = np.array([r['i'] for r in P if math.isfinite(r['t_pred'])])
    lim = [2, 45]
    ax[0].plot(lim, lim, color=INK3, lw=1); ax[0].text(3.0, 3.8, '1:1', color=INK2, fontsize=8, rotation=38)
    ax[0].scatter(o, p, s=34, color=C1, edgecolor=SURFACE, linewidth=1.2, zorder=3, label='objects (n=%d)' % len(o))
    ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlim(*lim); ax[0].set_ylim(*lim)
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator
    for axis in (ax[0].xaxis, ax[0].yaxis):
        axis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 40])); axis.set_major_formatter(FormatStrFormatter('%d')); axis.set_minor_locator(NullLocator())
    ax[0].set_xlabel('Observed years from first TLE to perigee < 300 km'); ax[0].set_ylabel('Forecast (J2+Sun+Moon, first TLE only)')
    rel = np.abs(p - o) / o
    ax[0].set_title(f'Median error {100*np.median(rel):.0f}%; {100*np.mean(rel<0.25):.0f}% within 25%', fontsize=9, loc='left', color=INK2)
    ax[0].legend(loc='lower right', frameon=False, fontsize=8)
    # right: ablation, perigee-peak epoch error at first-TLE start
    names = ['J2 + Sun + Moon', 'J2 + Moon only', 'J2 + Sun only', 'J2 only']; keys = ['J2+Sun+Moon', 'J2+Moon', 'J2+Sun', 'J2 only']; cols = [C1, C2, C3, INK3]
    med = []
    for v in keys:
        sel = [abs(r['pk_err_yr']) for r in R if r['start'] == 'first TLE' and r['variant'] == v]
        med.append(float(np.median(sel)))
    y = np.arange(len(names))[::-1]
    ax[1].barh(y, med, color=cols, height=0.55)
    for yi, m in zip(y, med):
        ax[1].text(m + 0.12, yi, f'{m:.2f} yr', va='center', fontsize=8, color=INK2)
    ax[1].set_yticks(y); ax[1].set_yticklabels(names); ax[1].set_xlim(0, 9.5)
    ax[1].set_xlabel('Median |error| of the perigee-peak epoch (years)')
    ax[1].set_title('Which drivers are required (17 Molniya-class objects)', fontsize=9, loc='left', color=INK2)
    for x in ax:
        style(x)
    ax[1].grid(False, axis='y')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig2_forecast_skill.png', dpi=170); plt.close(fig)


if __name__ == '__main__':
    fig1(7480); fig2(); print('ok')
