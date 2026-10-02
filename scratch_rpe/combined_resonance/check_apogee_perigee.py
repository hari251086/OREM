"""Empirical check of the observation that, for high-inclination HEO objects, the
apogee radius stays nearly constant *relative to* the perigee radius (issue #56 follow-up).

Theory (quadrupole lunisolar, doubly averaged; Celletti et al. 2016, full text): a is an exact
constant of the averaged motion, so d r_a = -d r_p = a d e, and the relative sensitivities are
    (d r_p / r_p) / (d r_a / r_a) = -(1+e)/(1-e).
This script measures, per object, the spread of log r_a, log r_p and log a over the pre-terminal part
of the TLE record (first 70 %, to keep drag-dominated decay out) and compares with the prediction.
"""
import math, sys
import numpy as np

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record  # noqa: E402
from resonance_detector import filter_outliers  # noqa: E402

RE = 6378.1363


def spread(x):
    x = np.asarray(x)
    return float(np.percentile(x, 95) - np.percentile(x, 5))


def main():
    objs = build_object_list()
    out = []
    for o in objs:
        try:
            rows = load_full_record(f"{BASE}/{o['tle_file']}")
        except FileNotFoundError:
            continue
        if len(rows) < 100:
            continue
        i_mean = sum(r[3] for r in rows) / len(rows)
        clean = filter_outliers([(r[0], r[1], r[2], r[3]) for r in rows])
        if len(clean) < 100:
            continue
        cut = int(0.7 * len(clean))
        c = np.array(clean[:cut])
        a, e = c[:, 1], c[:, 2]
        ra, rp = a * (1 + e), a * (1 - e)
        if np.median(rp) < RE + 300 or np.median(e) < 0.2:
            continue                                  # not a strongly eccentric, non-decaying phase
        s_ra, s_rp, s_a = spread(np.log(ra)), spread(np.log(rp)), spread(np.log(a))
        pred = (1 + np.median(e)) / (1 - np.median(e))
        # correlation of detrended increments: theory says d r_a ~ - d r_p (secular, a const)
        dra, drp = np.diff(ra), np.diff(rp)
        corr = float(np.corrcoef(dra, drp)[0, 1]) if len(dra) > 10 else float('nan')
        out.append((o['norad'], i_mean, len(clean), float(np.median(e)), s_ra, s_rp, s_a,
                    s_rp / s_ra if s_ra > 0 else float('nan'), pred, corr))
    out = np.array(out, dtype=float)
    hi = out[out[:, 1] >= 46.4]
    lo = out[out[:, 1] < 46.4]
    for name, d in (('i >= 46.4', hi), ('i < 46.4 ', lo)):
        print(f"{name}: n={len(d)}  median e={np.median(d[:,3]):.3f}  "
              f"spread(log r_a)={np.median(d[:,4]):.4f}  spread(log r_p)={np.median(d[:,5]):.4f}  "
              f"spread(log a)={np.median(d[:,6]):.4f}")
        print(f"          measured spread ratio rp/ra: median {np.nanmedian(d[:,7]):.2f}   "
              f"predicted (1+e)/(1-e): median {np.median(d[:,8]):.2f}   "
              f"corr(d r_a, d r_p): median {np.nanmedian(d[:,9]):.2f}")
    np.savetxt('E:/GitHub/OREM/scratch_rpe/combined_resonance/apogee_perigee_stats.csv', out, delimiter=',',
               header='norad,i_mean,n,e_med,spread_log_ra,spread_log_rp,spread_log_a,ratio_rp_over_ra,pred_ratio,corr_dra_drp',
               fmt='%.5g')


if __name__ == '__main__':
    main()
