"""Trend analysis of the 28 objects in OREM_Validation_Report.pdf (evolution plots): perigee/apogee shape and the solar
azimuth psi = (Omega + omega - L_sun) mod 360 (the report's 'resonance angle'), from the same TLE data.

Per object: inclination, first-TLE e/hp, smoothed perigee peak (epoch, rise), time peak->last TLE, solar-azimuth
circulation period over the whole record and over the last 15 % (slow-down), minimum |psi_dot| window relative to the end,
and whether psi_dot changes sign in the last 30 % of the record (U-turn).
"""
import sys, json
import numpy as np
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import load_full_record, L_sun_deg   # noqa: E402
from validate_drive import clean_rows, RE                              # noqa: E402
from molniya_cycle import smooth                                       # noqa: E402

ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
OBJ = {12907: 'COSMOS 1030 DEB', 35009: 'FREGAT R/B', 8833: 'MOLNIYA 3-5', 37398: 'MERIDIAN 4', 9574: 'MOLNIYA 2-16',
       27902: 'COSMOS 1261 DEB', 15030: 'SL-6 R/B', 37151: 'CZ-3B R/B', 14297: 'COSMOS 1456 DEB', 11073: 'SL-6 R/B',
       10802: 'SL-6 R/B', 27526: 'PSLV R/B', 9411: 'SL-6 R/B', 9269: 'SL-6 R/B', 10605: 'MOLNIYA 3-9', 11602: 'SL-6 R/B',
       13107: 'MOLNIYA 3-18', 40943: 'ARIANE 5 DEB', 8844: 'SL-6 R/B', 7653: 'SL-6 R/B', 35497: 'ARIANE 5 R/B',
       7480: 'MOLNIYA 1-28', 9506: 'SL-6 R/B', 7641: 'MOLNIYA 2-12', 7903: 'MOLNIYA 1-30', 8187: 'MOLNIYA 1-31',
       59347: 'FALCON 9 R/B', 32977: 'COSMOS 1109 DEB'}


def psi_rate(t_days, psi_unwrapped_deg, win=120.0):
    out_t, out_r = [], []
    t0 = t_days[0]
    while t0 + win <= t_days[-1]:
        sel = (t_days >= t0) & (t_days < t0 + win)
        if sel.sum() >= 6:
            out_t.append(np.mean(t_days[sel])); out_r.append(np.polyfit(t_days[sel], psi_unwrapped_deg[sel], 1)[0])
        t0 += 30.0
    return np.array(out_t), np.array(out_r)


def main():
    rows_out = []; series = {}
    for n, name in OBJ.items():
        try:
            rows = load_full_record(f"{ROOT}/{n}/{n}.tle.txt")
        except FileNotFoundError:
            print(n, 'missing'); continue
        arr = np.array(clean_rows(rows)); jd, a, e, inc, O, w = arr.T
        hp, ha = a * (1 - e) - RE, a * (1 + e) - RE
        t_yr = (jd - jd[0]) / 365.25
        psi = np.degrees(np.unwrap(np.radians((O + w - np.array([L_sun_deg(j) for j in jd])) % 360.0)))
        tt, rr = psi_rate(jd - jd[0], psi)
        s = smooth(hp); ipk = int(np.argmax(s))
        t_end = t_yr[-1]
        per_all = 360.0 / abs(np.median(rr)) / 365.25 if len(rr) else np.nan   # yr, from median rate
        last = tt >= 0.85 * tt[-1]
        rate_last = np.median(rr[last]) if last.sum() else np.nan
        tail = tt >= 0.7 * tt[-1]
        flips = bool(len(rr[tail]) > 2 and (np.nanmax(rr[tail]) > 0 > np.nanmin(rr[tail])))
        rows_out.append(dict(norad=n, name=name, i=float(np.mean(inc)), e0=float(e[0]), a0=float(a[0]), hp0=float(hp[0]),
                             ha0=float(ha[0]), hp_peak=float(hp[ipk]), t_peak=float(t_yr[ipk]), t_end=float(t_end),
                             rise=float(hp[ipk] - hp[0]), peak_to_end=float(t_end - t_yr[ipk]),
                             psi_rate_med=float(np.median(rr)), psi_rate_last=float(rate_last),
                             psi_period_yr=float(per_all), psi_flip_tail=flips,
                             hp_min_after=float(np.min(hp[ipk:])), n=len(jd)))
        series[n] = dict(t_end=t_yr[-1], t=t_yr.tolist(), hp=hp.tolist(), ha=ha.tolist(), tt=(tt / 365.25).tolist(), rr=rr.tolist())
    json.dump(rows_out, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/trend_all28.json', 'w'), indent=1)
    json.dump(series, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/trend_all28_series.json', 'w'))
    print(f"{'NORAD':>6} {'name':<16} {'i':>5} {'e0':>5} {'hp0':>5} {'hp_pk':>6} {'t_pk':>5} {'rise':>6} {'pk->end':>7} {'T_end':>5} {'psi rate med':>12} {'last15%':>8} {'flip':>5}")
    for r in sorted(rows_out, key=lambda r: r['i']):
        print(f"{r['norad']:>6} {r['name']:<16} {r['i']:5.1f} {r['e0']:5.2f} {r['hp0']:5.0f} {r['hp_peak']:6.0f} {r['t_peak']:5.1f} {r['rise']:6.0f} {r['peak_to_end']:7.1f} {r['t_end']:5.1f} {r['psi_rate_med']:+9.3f} d/d {r['psi_rate_last']:+8.3f} {str(r['psi_flip_tail']):>5}")


if __name__ == '__main__':
    main()
