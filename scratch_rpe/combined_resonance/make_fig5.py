"""Fig 5: decay-date error per object, forecast WITH drag (event = decay) vs earlier no-drag + leave-one-out-gap forecast.
dataviz skill: validated slots 1-2, thin marks, log axis, legend above, recessive grid."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE, INK, INK2, INK3, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e4e3df'
C1, C2 = '#2a78d6', '#eb6834'
D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
df = pd.read_csv(f'{D}/drag_test.csv').sort_values('norad').reset_index(drop=True)
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE})
fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6))
panels = (('Mid-life cutoff (first dip), decay 96-1,818 d ahead', 'mid_prev_150', 'mid_err_x1.0'),
          ('From the first TLE, ~10-15 yr ahead', 'full_prev_300', 'full_err_x1.0'))
for ax, (title, prev, drag) in zip(axes, panels):
    y = np.arange(len(df))[::-1]
    a = df[prev].abs().clip(lower=1); b = df[drag].abs().clip(lower=1)
    for yi, u, v in zip(y, a, b):
        ax.plot([u, v], [yi, yi], color=GRID, lw=1.2, zorder=1)
    ax.scatter(a, y, s=34, color=C1, edgecolor=SURFACE, linewidth=1.2, zorder=3, label='No drag + leave-one-out gap (earlier method)')
    ax.scatter(b, y, s=34, color=C2, edgecolor=SURFACE, linewidth=1.2, zorder=3, label='With drag, event = decay (B from DISCOS)')
    ax.set_xscale('log'); ax.set_xlim(1, 3000)
    ax.set_yticks(y); ax.set_yticklabels(df['norad'].astype(int))
    ax.set_xlabel('|Decay-date error| (days, log scale)')
    ax.set_title(f"{title}\nmedian |err|: gap {df[prev].abs().median():.0f} d, drag {df[drag].abs().median():.0f} d", fontsize=9.5, loc='left', color=INK)
    ax.grid(True, axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
axes[0].set_ylabel('NORAD ID')
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.995))
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(f'{D}/figures/fig5_forecast_with_drag.png', dpi=170)
print('ok')
