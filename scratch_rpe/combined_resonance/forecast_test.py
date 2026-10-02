"""Forecast test: start the secular propagator from TLE mean elements at epoch t_start and predict the perigee
cycle: (i) epoch of the perigee peak, (ii) epoch when perigee altitude first falls below 300 km after the peak
(E300; drag takes over below ~300 km), compared with the same events in the observed TLE series.
Variants: J2+Sun+Moon, J2+Moon only, J2+Sun only (tests which perturber is required).
"""
import sys, json, math
import numpy as np
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import load_full_record     # noqa: E402
from validate_drive import clean_rows                       # noqa: E402
from secular_propagator import propagate, RE                # noqa: E402
from molniya_cycle import NORADS, ROOT, smooth             # noqa: E402

VARIANTS = {'J2+Sun+Moon': (True, True, True), 'J2+Moon': (True, False, True), 'J2+Sun': (True, True, False),
            'J2 only': (True, False, False)}
LEADS = [0.0, 2.0, 4.0]   # years before the observed perigee peak at which the forecast starts (0 = at first TLE is separate)


def events(jd, hp):
    s = smooth(hp)
    ipk = int(np.argmax(s))
    after = np.where((np.arange(len(hp)) > ipk) & (s < 300.0))[0]
    i300 = int(after[0]) if len(after) else None
    return ipk, i300


def main():
    rows_out = []
    for n in NORADS:
        if n == 12907:
            continue
        rows = load_full_record(f"{ROOT}/{n}/{n}.tle.txt")
        arr = np.array(clean_rows(rows)); jd, a, e, inc, O, w = arr.T
        hp = a * (1 - e) - RE
        ipk, i300 = events(jd, hp)
        if i300 is None:
            print(n, 'no E300 in record'); continue
        t_pk_obs, t_300_obs = jd[ipk], jd[i300]
        # start epochs: first TLE (full-life forecast) and 4/2 years before the observed peak
        starts = [('first TLE', 0)]
        for lead in (4.0, 2.0):
            k = int(np.argmin(np.abs(jd - (t_pk_obs - lead * 365.25))))
            if k > 5:
                starts.append((f'{lead:.0f} yr before peak', k))
        a0 = float(np.median(a[:20]))
        for lab, k in starts:
            k0 = max(0, k - 5); k1 = k + 6
            e0, i0, O0, w0 = float(np.median(e[k0:k1])), float(np.median(inc[k0:k1])), float(np.median(O[k0:k1])), float(np.median(w[k0:k1]))
            # wrap-safe median for angles
            O0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[k0:k1])))))) % 360
            w0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[k0:k1])))))) % 360
            aa = float(np.median(a[k0:k1]))
            for vname, use in VARIANTS.items():
                traj = propagate(e0, i0, O0, w0, aa, jd[k], years=22.0, dt_days=10.0, use=use)
                hpp = aa * (1 - traj[:, 1]) - RE
                s = smooth(hpp)
                ipp = int(np.argmax(s[: max(5, int(len(s) * 0.9))]))
                aft = np.where((np.arange(len(hpp)) > ipp) & (s < 300.0))[0]
                t_pk_p = traj[ipp, 0]
                t_300_p = traj[aft[0], 0] if len(aft) else float('nan')
                rows_out.append(dict(norad=n, start=lab, variant=vname, t_start=float(jd[k]),
                                     pk_err_yr=float((t_pk_p - t_pk_obs) / 365.25),
                                     e300_err_yr=float((t_300_p - t_300_obs) / 365.25) if not math.isnan(t_300_p) else None,
                                     horizon_yr=float((t_300_obs - jd[k]) / 365.25),
                                     hp_pk_pred=float(hpp[ipp]), hp_pk_obs=float(hp[ipk])))
        print(n, 'done', flush=True)
    json.dump(rows_out, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/forecast_results.json', 'w'), indent=1)
    import collections
    print("\nERROR OF FORECAST EVENTS (years, predicted - observed); |err| median [IQR] over objects")
    for lab in ('first TLE', '4 yr before peak', '2 yr before peak'):
        for v in VARIANTS:
            sel = [r for r in rows_out if r['start'] == lab and r['variant'] == v]
            if not sel:
                continue
            pe = np.array([r['pk_err_yr'] for r in sel])
            ee = np.array([r['e300_err_yr'] for r in sel if r['e300_err_yr'] is not None])
            hz = np.median([r['horizon_yr'] for r in sel])
            print(f"{lab:>18} | {v:<12} n={len(sel):2d} horizon~{hz:4.1f} yr | peak epoch err med |{np.median(np.abs(pe)):5.2f}| "
                  f"[{np.percentile(np.abs(pe),25):.2f},{np.percentile(np.abs(pe),75):.2f}]  E300 err med |{np.median(np.abs(ee)) if len(ee) else float('nan'):5.2f}| "
                  f"[{np.percentile(np.abs(ee),25) if len(ee) else float('nan'):.2f},{np.percentile(np.abs(ee),75) if len(ee) else float('nan'):.2f}]  bias {np.median(ee) if len(ee) else float('nan'):+.2f}")


if __name__ == '__main__':
    main()
