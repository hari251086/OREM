"""Is amplitude/frequency (Talu et al. 2021 ranking idea) the right strength measure for the perigee excursion?
Exact doubly-averaged theory: dq/dt = -A_P(t) sin(2 w'_P).  If w'_P circulates at rate wdot, the perigee
excursion over a half cycle is  E = A / |wdot|  (km).  Resonance (stationary w') => wdot -> 0, E large.
Test: per object, in windows of W2 days, compare the predicted E (Moon+Sun, from window-mean A and the fitted
wdot of the unwrapped angle) with the observed range of q in the window.
"""
import sys
import numpy as np
from scipy.stats import spearmanr
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from quadrupole_drive import angles                                                  # noqa: E402
from validate_drive import clean_rows, RE                                            # noqa: E402

W2, STEP2 = 540.0, 90.0


def main():
    X, Y, band, nid = [], [], [], []
    for o in build_object_list():
        try:
            rows = load_full_record(f"{BASE}/{o['tle_file']}")
        except FileNotFoundError:
            continue
        if len(rows) < 200:
            continue
        i_mean = sum(r[3] for r in rows) / len(rows)
        if i_mean < 46.4:
            continue
        arr = np.array(clean_rows(rows)); jd, a, e, inc, O, w = arr.T
        if np.median(e) < 0.2:
            continue
        q = a * (1 - e)
        (wps, s2s, As), (wpm, s2m, Am) = angles(jd, a, e, inc, O, w)
        ok = np.isfinite(wps) & np.isfinite(wpm)
        jd, a, q, wps, wpm, As, Am = jd[ok], a[ok], q[ok], wps[ok], wpm[ok], As[ok], Am[ok]
        if len(jd) < 150:
            continue
        us = np.degrees(np.unwrap(np.radians(2 * wps))); um = np.degrees(np.unwrap(np.radians(2 * wpm)))  # 2w'
        t0 = jd[0]
        while t0 + W2 <= jd[-1]:
            sel = (jd >= t0) & (jd < t0 + W2); t0 += STEP2
            if sel.sum() < 12:
                continue
            t = jd[sel]
            if t[-1] - t[0] < 0.7 * W2 or (a[sel].max() - a[sel].min()) / np.median(a[sel]) > 0.004:
                continue
            if np.median(q[sel]) - RE < 400.0:
                continue
            rate_m = abs(np.polyfit(t - t[0], um[sel], 1)[0])       # deg/day of 2w'
            rate_s = abs(np.polyfit(t - t[0], us[sel], 1)[0])
            Em = np.mean(np.abs(Am[sel])) / max(np.radians(rate_m), 1e-6)    # km  (A [km/day] / [rad/day])
            Es = np.mean(np.abs(As[sel])) / max(np.radians(rate_s), 1e-6)
            E = Em + Es
            obs = np.ptp(q[sel])
            X.append(E); Y.append(obs); band.append(i_mean); nid.append(o['norad'])
    X, Y, band, nid = map(np.array, (X, Y, band, nid))
    print('windows', len(X), 'objects', len(set(nid)))
    print(f"Spearman(E_pred, observed q-range) = {spearmanr(X, Y)[0]:+.3f}")
    # compare against amplitude-only and rate-only predictors (does the ratio add something?)
    # (recompute quickly from stored arrays is not kept; report only the ratio here)
    lx, ly = np.log10(np.clip(X, 1e-3, None)), np.log10(np.clip(Y, 1e-3, None))
    print(f"log-log slope (obs range vs E) = {np.polyfit(lx, ly, 1)[0]:.2f},  corr(log) = {np.corrcoef(lx, ly)[0,1]:+.3f}")
    for lo, hi in ((46.4, 56.0), (56.0, 70.0), (70.0, 181.0)):
        m = (band >= lo) & (band < hi)
        if m.sum() > 15:
            print(f"  {lo:.0f}-{hi:.0f} deg: n={m.sum():4d}  Spearman={spearmanr(X[m], Y[m])[0]:+.3f}")
    # per-object Spearman
    per = []
    for n_ in set(nid):
        m = nid == n_
        if m.sum() >= 6:
            per.append(spearmanr(X[m], Y[m])[0])
    print(f"per-object Spearman: median {np.nanmedian(per):+.3f}  (n={len(per)} objects with >=6 windows)")


if __name__ == '__main__':
    main()
