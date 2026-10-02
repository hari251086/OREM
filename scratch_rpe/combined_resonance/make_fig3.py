"""Fig 3: trends across the 28 objects of OREM_Validation_Report.pdf (dataviz skill: validated slots 1-2, thin lines,
separate panels for different scales, legend for the two-group panel, recessive grid)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE, INK, INK2, INK3, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e4e3df'
C1, C2 = '#2a78d6', '#eb6834'
D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
meta = {r['norad']: r for r in json.load(open(f'{D}/trend_all28.json'))}
ser = json.load(open(f'{D}/trend_all28_series.json'))
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE})
hi = [n for n, r in meta.items() if r['i'] >= 55]
lo = [n for n, r in meta.items() if r['i'] < 46.4]


def style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(length=0)


fig, ax = plt.subplots(3, 1, figsize=(8.4, 10.6))
for n in hi:
    s = ser[str(n)]; te = s['t_end']
    ax[0].plot(np.array(s['t']) - te, s['hp'], color=C1, lw=1.0, alpha=0.55)
ax[0].set_ylim(0, 6000); ax[0].set_xlim(-22, 0.5)
ax[0].set_ylabel('Perigee altitude (km)')
ax[0].set_title(f'High inclination (63-68 deg, n={len(hi)}): perigee rises 1-3 thousand km, then falls to decay', fontsize=9.5, loc='left', color=INK)
for n in lo:
    s = ser[str(n)]; te = s['t_end']
    ax[1].plot(np.array(s['t']) - te, s['hp'], color=C2, lw=1.2, alpha=0.8)
ax[1].set_ylim(80, 450); ax[1].set_xlim(-10, 0.3)
ax[1].set_ylabel('Perigee altitude (km)')
ax[1].set_title(f'Low inclination (6-26 deg, n={len(lo)}): perigee stays 120-260 km; decay is apogee-led', fontsize=9.5, loc='left', color=INK)
ax[2].axhline(0, color=INK3, lw=0.8)
for n in hi:
    s = ser[str(n)]; te = s['t_end']
    ax[2].plot(np.array(s['tt']) - te, s['rr'], color=C1, lw=1.0, alpha=0.55)
for n in lo:
    s = ser[str(n)]; te = s['t_end']
    ax[2].plot(np.array(s['tt']) - te, s['rr'], color=C2, lw=1.2, alpha=0.8)
ax[2].set_xlim(-22, 0.5); ax[2].set_ylim(-1.4, 0.6)
ax[2].plot([], [], color=C1, lw=2, label='high inclination (n=%d)' % len(hi)); ax[2].plot([], [], color=C2, lw=2, label='low inclination (n=%d)' % len(lo))
ax[2].legend(loc='upper left', frameon=False, fontsize=8, ncol=2)
ax[2].set_ylabel('Solar-azimuth rate d(psi)/dt (deg/day)')
ax[2].set_xlabel('Years before the last TLE')
ax[2].set_title("The report's 'resonance angle' psi: steady circulation at high i; slows and reverses at low i (U-turn)", fontsize=9.5, loc='left', color=INK)
for a_ in ax:
    style(a_)
fig.tight_layout(); fig.savefig(f'{D}/figures/fig3_report_trends.png', dpi=170)
print('ok', len(hi), len(lo))
