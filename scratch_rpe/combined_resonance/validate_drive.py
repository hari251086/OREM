"""Does the window-averaged quadrupole drive predict the observed perigee-radius trend of high-i HEO objects?

For each object (i_mean >= 46.4 deg) and each sliding window (W days), compare
    observed : slope of q = a(1-e) from the TLE series
    predicted: time-average of the instantaneous Sun+Moon quadrupole drive dq/dt (quadrupole_drive.py)
Windows with a-drift > 0.3 % or median perigee height < 400 km are dropped (drag regime, a not constant).
Control: the same windows with the predicted series permuted in time within each object.
"""
import json, math, sys, random
import numpy as np

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from quadrupole_drive import drive_dqdt                                            # noqa: E402

RE = 6378.1363
W = 180.0
STEP = 45.0
random.seed(1); np.random.seed(1)


def clean_rows(rows):
    """Drop single-point a glitches (rolling-median, 5 %)."""
    a = np.array([r[1] for r in rows])
    k = 9
    keep = []
    for idx in range(len(rows)):
        lo, hi = max(0, idx - k // 2), min(len(rows), idx + k // 2 + 1)
        m = np.median(a[lo:hi])
        keep.append(abs(a[idx] - m) / m <= 0.05)
    return [r for r, kk in zip(rows, keep) if kk]


def windows(jd, a, e, q, ds, dm):
    out = []
    t0 = jd[0]
    while t0 + W <= jd[-1]:
        sel = (jd >= t0) & (jd < t0 + W)
        t0 += STEP
        if sel.sum() < 8:
            continue
        t = jd[sel]
        if t[-1] - t[0] < 0.6 * W:
            continue
        aw = a[sel]
        if (aw.max() - aw.min()) / np.median(aw) > 0.003:
            continue
        if np.median(q[sel]) - RE < 400.0:
            continue
        obs = np.polyfit(t - t[0], q[sel], 1)[0]
        span = t[-1] - t[0]
        ps = np.trapezoid(ds[sel], t) / span
        pm = np.trapezoid(dm[sel], t) / span
        out.append((obs, ps, pm, float(np.median(np.degrees(0) + 0 * aw))))
    return out


def main():
    objs = build_object_list()
    recs = []
    for o in objs:
        try:
            rows = load_full_record(f"{BASE}/{o['tle_file']}")
        except FileNotFoundError:
            continue
        if len(rows) < 150:
            continue
        i_mean = sum(r[3] for r in rows) / len(rows)
        if i_mean < 46.4:
            continue
        rows = clean_rows(rows)
        arr = np.array(rows)
        jd, a, e, inc, O, w = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4], arr[:, 5]
        if np.median(e) < 0.2:
            continue
        q = a * (1 - e)
        ds, dm = drive_dqdt(jd, a, e, inc, O, w)
        ok = np.isfinite(ds) & np.isfinite(dm)
        jd, a, e, q, ds, dm = jd[ok], a[ok], e[ok], q[ok], ds[ok], dm[ok]
        if len(jd) < 100:
            continue
        ws = windows(jd, a, e, q, ds, dm)
        if len(ws) < 3:
            continue
        # permutation control: shuffle predicted (sun+moon) across windows of this object
        pred = np.array([w_[1] + w_[2] for w_ in ws])
        recs.append({'norad': o['norad'], 'i': float(i_mean), 'obs': [w_[0] for w_ in ws],
                     'ps': [w_[1] for w_ in ws], 'pm': [w_[2] for w_ in ws]})
    print('objects with usable windows:', len(recs), ' windows:', sum(len(r['obs']) for r in recs))

    def stats(sel, label):
        obs = np.concatenate([np.array(r['obs']) for r in sel]); ps = np.concatenate([np.array(r['ps']) for r in sel])
        pm = np.concatenate([np.array(r['pm']) for r in sel]); pt = ps + pm
        if len(obs) < 10:
            print(f"{label}: too few windows"); return
        c = np.corrcoef(obs, pt)[0, 1]
        beta = float(np.sum(obs * pt) / np.sum(pt * pt))
        sgn = np.mean(np.sign(obs) == np.sign(pt))
        # null: permute predicted within each object
        nulls = []
        for _ in range(300):
            pp = np.concatenate([np.random.permutation(np.array(r['ps']) + np.array(r['pm'])) for r in sel])
            nulls.append(np.corrcoef(obs, pp)[0, 1])
        cm = np.corrcoef(obs, pm)[0, 1]; cs = np.corrcoef(obs, ps)[0, 1]
        print(f"{label:<14} objs={len(sel):3d} win={len(obs):5d}  corr(total)={c:+.3f}  (sun {cs:+.3f}, moon {cm:+.3f})  "
              f"slope={beta:+.2f}  sign-agree={sgn:.2f}  null corr 95%=[{np.percentile(nulls,2.5):+.3f},{np.percentile(nulls,97.5):+.3f}]")

    stats(recs, 'all i>=46.4')
    for lo, hi in ((46.4, 56.0), (56.0, 70.0), (70.0, 90.0), (90.0, 181.0)):
        sel = [r for r in recs if lo <= r['i'] < hi]
        if sel:
            stats(sel, f'{lo:.0f}<=i<{hi:.0f}')
    json.dump(recs, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/drive_windows.json', 'w'))


if __name__ == '__main__':
    main()
