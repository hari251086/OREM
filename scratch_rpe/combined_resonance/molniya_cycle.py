"""Molniya-type objects (i ~ 62-66 deg, e ~ 0.7): what drives the multi-year perigee rise before decay?

For each object: integrate the exact doubly-averaged Sun+Moon quadrupole drive (quadrupole_drive.py) along the
observed mean elements, q_rec(t) = q(t0) + int (dq/dt)_Sun+Moon dt, and compare with the observed perigee radius
over the lunisolar phase (until the perigee altitude first falls below 400 km after the peak, where drag takes over).
Also records the perigee-peak epoch, its lead time before the last TLE, the Sun/Moon shares, and 2w'_P at start/peak.
"""
import sys, json
import numpy as np

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import load_full_record     # noqa: E402
from quadrupole_drive import drive_dqdt, angles              # noqa: E402
from validate_drive import clean_rows, RE                    # noqa: E402

ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
NORADS = [7480, 7641, 7653, 7903, 8187, 8833, 8844, 9269, 9411, 9506, 9574, 10605, 10802, 11073, 11602,
          12907, 13107, 14297, 15030]


def smooth(x, k=15):
    k = min(k, max(3, len(x) // 5) | 1)
    pad = np.pad(x, k // 2, mode='edge')
    return np.convolve(pad, np.ones(k) / k, mode='valid')


def cumtrapz(y, t):
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))
    return out


def main():
    res = []
    series = {}
    for n in NORADS:
        try:
            rows = load_full_record(f"{ROOT}/{n}/{n}.tle.txt")
        except FileNotFoundError:
            print(n, 'no file'); continue
        if len(rows) < 100:
            print(n, 'too few', len(rows)); continue
        arr = np.array(clean_rows(rows)); jd, a, e, inc, O, w = arr.T
        q = a * (1 - e); hp = q - RE; ha = a * (1 + e) - RE
        ds, dm = drive_dqdt(jd, a, e, inc, O, w)
        (wps, _, _), (wpm, _, _) = angles(jd, a, e, inc, O, w)
        ok = np.isfinite(ds) & np.isfinite(dm)
        jd, a, e, inc, q, hp, ha, ds, dm, wps, wpm = [x[ok] for x in (jd, a, e, inc, q, hp, ha, ds, dm, wps, wpm)]
        t_yr = (jd - jd[0]) / 365.25
        qs = smooth(q)
        ipk = int(np.argmax(qs))
        # lunisolar phase ends when perigee altitude first drops below 400 km after the peak (drag takes over)
        after = np.where((np.arange(len(hp)) > ipk) & (smooth(hp) < 400.0))[0]
        iend = int(after[0]) if len(after) else len(hp) - 1
        cum_s, cum_m = cumtrapz(ds, jd), cumtrapz(dm, jd)
        rec = q[0] + cum_s + cum_m
        sl = slice(0, iend + 1)
        err = rec[sl] - q[sl]
        rng = np.ptp(q[sl])
        # restrict the comparison to the part with a non-trivially large excursion
        corr = float(np.corrcoef(rec[sl] - q[0], q[sl] - q[0])[0, 1])
        rise_obs = q[ipk] - q[0]; rise_rec = rec[ipk] - q[0]
        moon_share = float(cum_m[ipk] / (cum_m[ipk] + cum_s[ipk])) if (cum_m[ipk] + cum_s[ipk]) != 0 else np.nan
        res.append(dict(norad=n, n=len(jd), i=float(np.mean(inc)), e0=float(e[0]), a0=float(a[0]),
                        hp0=float(hp[0]), hp_peak=float(hp[ipk]), t_peak_yr=float(t_yr[ipk]),
                        t_end_yr=float(t_yr[-1]), lead_peak_to_last_yr=float(t_yr[-1] - t_yr[ipk]),
                        rise_obs=float(rise_obs), rise_rec=float(rise_rec), corr=corr,
                        nrmse=float(np.sqrt(np.mean(err ** 2)) / rng), moon_share=moon_share,
                        w2m_start=float(np.degrees(np.angle(np.exp(1j * np.radians(2 * wpm[0]))))),
                        w2m_peak=float(np.degrees(np.angle(np.exp(1j * np.radians(2 * wpm[ipk]))))),
                        w2s_start=float(np.degrees(np.angle(np.exp(1j * np.radians(2 * wps[0]))))),
                        w2s_peak=float(np.degrees(np.angle(np.exp(1j * np.radians(2 * wps[ipk]))))),
                        ))
        series[n] = dict(t=t_yr.tolist(), hp=hp.tolist(), ha=ha.tolist(), rec_hp=(rec - RE).tolist(),
                         ds=ds.tolist(), dm=dm.tolist(), w2m=(2 * wpm).tolist(), w2s=(2 * wps).tolist(),
                         ipk=ipk, iend=iend, jd=jd.tolist())
    json.dump(series, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/molniya_series.json', 'w'))
    json.dump(res, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/molniya_cycle.json', 'w'), indent=1)
    print(f"{'NORAD':>6} {'i':>5} {'hp0':>6} {'hp_pk':>6} {'t_pk(yr)':>8} {'pk->last':>8} {'rise_obs':>8} {'rise_rec':>8} {'corr':>6} {'nRMSE':>6} {'Moon%':>6} {'2w_M@0':>7} {'2w_M@pk':>8}")
    for r in res:
        print(f"{r['norad']:>6} {r['i']:5.1f} {r['hp0']:6.0f} {r['hp_peak']:6.0f} {r['t_peak_yr']:8.2f} {r['lead_peak_to_last_yr']:8.2f} "
              f"{r['rise_obs']:8.0f} {r['rise_rec']:8.0f} {r['corr']:6.3f} {r['nrmse']:6.3f} {100*r['moon_share']:6.0f} {r['w2m_start']:7.0f} {r['w2m_peak']:8.0f}")
    c = np.array([r['corr'] for r in res]); nm = np.array([r['nrmse'] for r in res])
    print(f"\nmedian corr {np.median(c):.3f}  median nRMSE {np.median(nm):.3f}  median Moon share {np.nanmedian([r['moon_share'] for r in res]):.2f}")
    ro = np.array([r['rise_obs'] for r in res]); rr = np.array([r['rise_rec'] for r in res])
    print(f"rise (obs vs reconstructed) corr {np.corrcoef(ro, rr)[0,1]:.3f}  median ratio rec/obs {np.median(rr/ro):.2f}")
    lead = np.array([r['lead_peak_to_last_yr'] for r in res])
    print(f"peak -> last TLE lead: median {np.median(lead):.2f} yr  IQR [{np.percentile(lead,25):.2f}, {np.percentile(lead,75):.2f}]  range [{lead.min():.2f}, {lead.max():.2f}]")


if __name__ == '__main__':
    main()
