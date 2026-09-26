"""Manuscript (issues #45/#47/#48): figures for the 23-object lead-time
hindcast campaign. Reuses this project's own dataviz-skill palette
(imported from OREM-Watchlist's plotting.py, not re-derived) for
visual consistency with the dashboard.

Two figures, per the dataviz skill's "small multiples" rule for a
series count too large for a shared categorical palette (23 objects
can't each get a distinct hue -- they fold into small multiples
instead, one panel per object):

  fig1_predict_rate.png   -- single series: how many of the 23 objects
                             produced ANY prediction, at each lead time.
                             All 23 counted in the denominator (not just
                             the ones with TLE data before the cutoff)
                             so the figure is the honest population-
                             level number, not a survivorship-filtered
                             one -- directly answers the 2019 review's
                             "robustness and general applicability"
                             comment (#48).
  fig2_rpe_by_object.png  -- small multiples, one panel per object
                             (all 23, sorted by median |RPE| among
                             objects that predicted; non-predicting
                             objects sort last and show an explicit
                             "no prediction" panel instead of being
                             dropped) -- ensemble RPE vs lead time,
                             shared y-axis for comparability, a shaded
                             +/-10% reference band, single categorical
                             color (slot 1 blue -- one series per
                             panel, no identity encoding needed within
                             a panel).

Usage: python scratch_rpe/manuscript23/plot_manuscript_figures.py
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("E:/GitHub/OREM-Watchlist/src")))
from heowatch.plotting import _BASELINE, _GRIDLINE, _INK, _INK_MUTED, _INK_SECONDARY, _SURFACE, _BLUE, _RED  # noqa: E402

HERE = Path(__file__).parent
LEADS = [30, 90, 180, 365]


def _style_axes(ax):
    ax.set_facecolor(_SURFACE)
    ax.grid(True, color=_GRIDLINE, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(_BASELINE)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=_INK_MUTED, labelcolor=_INK_SECONDARY, pad=4)
    ax.title.set_color(_INK)


def fig1_predict_rate(df):
    counts = []
    for lead in LEADS:
        sub = df[df["lead_days"] == lead]
        counts.append(int(sub["predicted"].sum()))

    fig, ax = plt.subplots(figsize=(6, 4))
    fig.patch.set_facecolor(_SURFACE)
    x = np.arange(len(LEADS))
    ax.plot(x, counts, color=_BLUE, linewidth=2, marker="o",
            markersize=8, markerfacecolor=_BLUE, markeredgecolor=_SURFACE,
            markeredgewidth=1.5, zorder=3)
    for xi, c in zip(x, counts):
        ax.annotate(f"{c}/23", (xi, c), textcoords="offset points",
                    xytext=(0, 10), ha="center", color=_INK, fontsize=10,
                    fontweight="semibold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}d" for d in LEADS])
    ax.set_ylim(0, 23)
    ax.set_xlabel("Lead time before decay", color=_INK_SECONDARY)
    ax.set_ylabel("Objects with a re-entry prediction (of 23)", color=_INK_SECONDARY)
    ax.set_title("Predict rate vs. lead time -- 23-object resonance set", color=_INK, fontsize=12)
    _style_axes(ax)
    fig.tight_layout()
    out = HERE / "fig1_predict_rate.png"
    fig.savefig(out, dpi=200, facecolor=_SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def fig2_rpe_by_object(df):
    # dropna=False: pivot_table's default silently drops a (norad,name)
    # row when every lead_days value for it is NaN -- exactly the 12
    # never-predicting objects the manuscript specifically must not
    # drop. Confirmed live: default pivot_table returned 11 rows, not
    # 23, before this fix.
    piv = df.pivot_table(index=["norad", "name"], columns="lead_days",
                          values="ens_rpe_pct", aggfunc="first", dropna=False)
    all_pairs = df[["norad", "name"]].drop_duplicates().apply(tuple, axis=1)
    piv = piv.reindex(pd.MultiIndex.from_tuples(sorted(set(all_pairs))))
    assert len(piv) == 23, f"expected 23 objects, got {len(piv)}"
    any_pred = piv.notna().any(axis=1)
    med = piv[any_pred].abs().median(axis=1)
    order_pred = med.sort_values().index.tolist()
    order_nopred = sorted(piv[~any_pred].index.tolist(), key=lambda t: t[0])
    order = order_pred + order_nopred

    n = len(order)
    ncols = 5
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 2.5 * nrows),
                              sharex=True, sharey=True)
    fig.patch.set_facecolor(_SURFACE)
    axes = axes.flatten()

    ymax = max(20.0, np.nanmax(piv.abs().values) * 1.1) if any_pred.any() else 20.0
    x = np.arange(len(LEADS))

    for ax, (norad, name) in zip(axes, order):
        row = piv.loc[(norad, name)]
        ax.axhspan(-10, 10, color=_BASELINE, alpha=0.15, zorder=1)
        ax.axhline(0, color=_BASELINE, linewidth=1, zorder=1)
        vals = [row.get(d, np.nan) for d in LEADS]
        ax.set_ylim(-ymax, ymax)  # shared across every panel, incl. "no prediction"
        if all(np.isnan(v) for v in vals):
            ax.text(0.5, 0.5, "no prediction\nat any lead time",
                    transform=ax.transAxes, ha="center", va="center",
                    color=_INK_MUTED, fontsize=8.5)
        else:
            ax.plot(x, vals, color=_BLUE, linewidth=1.6, marker="o",
                    markersize=5, markerfacecolor=_BLUE,
                    markeredgecolor=_SURFACE, markeredgewidth=0.8, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{d}" for d in LEADS], fontsize=7.5)
        ax.set_title(f"{norad} {name}"[:26], fontsize=8.5, color=_INK, pad=3)
        _style_axes(ax)
        ax.tick_params(labelsize=7.5)

    for ax in axes[n:]:
        ax.axis("off")

    fig.supxlabel("Lead time before decay (days)", color=_INK_SECONDARY, fontsize=10)
    fig.supylabel("Ensemble RPE (%)", color=_INK_SECONDARY, fontsize=10)
    fig.suptitle(
        "Ensemble RPE vs. lead time, all 23 resonance-set objects "
        "(shaded band = |RPE| < 10%)",
        color=_INK, fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0.01, 0.01, 1, 0.97])
    out = HERE / "fig2_rpe_by_object.png"
    fig.savefig(out, dpi=200, facecolor=_SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    df = pd.read_csv(HERE / "manuscript_table_full.csv")
    fig1_predict_rate(df)
    fig2_rpe_by_object(df)


if __name__ == "__main__":
    main()
