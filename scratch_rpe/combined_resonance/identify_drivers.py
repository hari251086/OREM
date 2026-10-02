"""Which slow-angle drivers are stationary when the quadrupole drive is coherent (i.e. when the perigee
is being pumped monotonically)?  Coherence kappa = |window-mean drive| / rms(instantaneous drive):
kappa ~ 0 -> oscillating drive that averages out; kappa -> 1 -> a stationary angle (resonance, U-turn analogue).
For each window we also fit the slope of every family angle Psi = k w + m O + s O_L (unwrapped) from the TLE
series and record the slowest one.
"""
import sys, collections
import numpy as np

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from quadrupole_drive import drive_dqdt                                            # noqa: E402
from drivers import family, theta_L                                                # noqa: E402
from validate_drive import clean_rows, RE, W, STEP                                 # noqa: E402

FAM = [f for f in family() if f[4] == 0 and not (f[2] == 0 and f[3] == 0 and f[1] == 0)]


def main():
    out = []     # (i_mean, kappa, mean_drive, slowest_label, slowest_abs, slowest_pump_label, slowest_pump_abs, norad)
    for o in build_object_list():
        try:
            rows = load_full_record(f"{BASE}/{o['tle_file']}")
        except FileNotFoundError:
            continue
        if len(rows) < 150:
            continue
        i_mean = sum(r[3] for r in rows) / len(rows)
        if i_mean < 46.4:
            continue
        arr = np.array(clean_rows(rows))
        jd, a, e, inc, O, w = arr.T
        if np.median(e) < 0.2:
            continue
        q = a * (1 - e)
        ds, dm = drive_dqdt(jd, a, e, inc, O, w)
        tot = ds + dm
        ok = np.isfinite(tot)
        jd, a, e, q, tot, O, w = jd[ok], a[ok], e[ok], q[ok], tot[ok], O[ok], w[ok]
        if len(jd) < 100:
            continue
        Ou = np.degrees(np.unwrap(np.radians(O))); wu = np.degrees(np.unwrap(np.radians(w)))
        OL = theta_L(jd)
        t0 = jd[0]
        while t0 + W <= jd[-1]:
            sel = (jd >= t0) & (jd < t0 + W)
            t0 += STEP
            if sel.sum() < 8:
                continue
            t = jd[sel]
            if t[-1] - t[0] < 0.6 * W:
                continue
            if (a[sel].max() - a[sel].min()) / np.median(a[sel]) > 0.003:
                continue
            if np.median(q[sel]) - RE < 400.0:
                continue
            span = t[-1] - t[0]
            mean = np.trapezoid(tot[sel], t) / span
            rms = np.sqrt(np.trapezoid(tot[sel] ** 2, t) / span)
            if rms <= 0:
                continue
            kappa = abs(mean) / rms
            best = (1e9, None); best_p = (1e9, None)
            for (lab, k, m, s, j, pump) in FAM:
                psi = k * wu[sel] + m * Ou[sel] + s * OL[sel]
                sl = abs(np.polyfit(t - t[0], psi, 1)[0])
                if sl < best[0]:
                    best = (sl, lab)
                if pump and sl < best_p[0]:
                    best_p = (sl, lab)
            out.append((i_mean, kappa, mean, best[1], best[0], best_p[1], best_p[0], o['norad']))
    print('windows:', len(out))
    kap = np.array([r[1] for r in out])
    lo, hi = np.percentile(kap, 25), np.percentile(kap, 75)
    print(f"kappa quartiles: 25%={lo:.3f}  median={np.median(kap):.3f}  75%={hi:.3f}")
    slow = np.array([r[6] for r in out])
    print(f"slowest e-pumping driver |Psi_dot| (deg/day): coherent(top quartile) median={np.median(slow[kap >= hi]):.4f}  "
          f"incoherent(bottom quartile) median={np.median(slow[kap <= lo]):.4f}")
    from scipy.stats import spearmanr
    print(f"Spearman(kappa, -|Psi_dot|slowest-pump) = {spearmanr(kap, -slow)[0]:+.3f}")
    for name, sel in (('COHERENT (top quartile)', kap >= hi), ('INCOHERENT (bottom quartile)', kap <= lo)):
        cnt = collections.Counter(out[i][5] for i in np.where(sel)[0])
        print(f"\n{name}: slowest e-pumping driver, top 8 of {sel.sum()}")
        for lab, n in cnt.most_common(8):
            print(f"   {lab:>12}  {n:5d}  ({100*n/sel.sum():.0f}%)")
    # by inclination band, coherent windows only
    print("\nCoherent windows by inclination band -> top 3 drivers")
    for lo_i, hi_i in ((46.4, 56.0), (56.0, 63.0), (63.0, 70.0), (70.0, 181.0)):
        idx = [i for i in np.where(kap >= hi)[0] if lo_i <= out[i][0] < hi_i]
        if len(idx) < 15:
            print(f"  {lo_i:.0f}-{hi_i:.0f}: n={len(idx)} (too few)"); continue
        cnt = collections.Counter(out[i][5] for i in idx)
        print(f"  {lo_i:.0f}-{hi_i:.0f} deg: n={len(idx)}  " + "  ".join(f"{l} {100*n/len(idx):.0f}%" for l, n in cnt.most_common(3)))


if __name__ == '__main__':
    main()
