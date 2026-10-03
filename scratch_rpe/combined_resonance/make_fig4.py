"""Fig 4: per-object |decay-date error| of OREM (hindcast ~90 d before decay) vs the secular forecast (first TLE only, ~13 yr ahead).
dataviz skill: validated slots 1-2, thin marks, log axis, direct labels + legend, recessive grid."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE, INK, INK2, INK3, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e4e3df'
C1, C2 = '#2a78d6', '#eb6834'
D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE})
fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.4), sharey=False)
for ax, thr in zip(axes, (300, 200)):
    df = pd.read_csv(f'{D}/cross_check_orem_{thr}.csv').sort_values('norad').reset_index(drop=True)
    y = np.arange(len(df))[::-1]
    eo, ef = df['err_orem_d'].abs().clip(lower=1), df['err_fc_d'].abs().clip(lower=1)
    for yi, a, b in zip(y, eo, ef):
        ax.plot([a, b], [yi, yi], color=GRID, lw=1.2, zorder=1)
    ax.scatter(eo, y, s=34, color=C1, edgecolor=SURFACE, linewidth=1.2, zorder=3, label='OREM (cutoff ~90 d before decay)')
    ax.scatter(ef, y, s=34, color=C2, edgecolor=SURFACE, linewidth=1.2, zorder=3, label='Secular forecast (first TLE, ~13 yr ahead)')
    ax.set_xscale('log'); ax.set_xlim(1, 3000)
    ax.set_yticks(y); ax.set_yticklabels(df['norad'].astype(int))
    ax.set_xlabel('|Decay-date error| (days, log scale)')
    ax.set_title(f'Forecast event: perigee < {thr} km after the peak (n={len(df)})', fontsize=9.5, loc='left', color=INK)
    ax.grid(True, axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.995))
axes[0].set_ylabel('NORAD ID')
fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(f'{D}/figures/fig4_cross_check_orem.png', dpi=170)
print('ok')
