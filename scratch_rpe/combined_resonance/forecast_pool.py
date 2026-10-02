"""Generality: first-TLE secular forecast (J2+Sun+Moon, a=const, no drag) of the epoch E300 at which the perigee
altitude first falls below 300 km, for every i >= 46.4 deg object of the 253-object pool whose record reaches it."""
import sys, json, math
import numpy as np
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from validate_drive import clean_rows                                              # noqa: E402
from secular_propagator import propagate, RE                                       # noqa: E402
from molniya_cycle import smooth                                                   # noqa: E402


def main():
    out = []
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
        arr = np.array(clean_rows(rows)); jd, a, e, inc, O, w = arr.T
        if np.median(e[:20]) < 0.3:
            continue
        hp = a * (1 - e) - RE
        s = smooth(hp)
        idx = np.where(s < 300.0)[0]
        if len(idx) == 0:
            continue                                   # record does not reach E300
        i300 = int(idx[0])
        if i300 < 40 or (jd[i300] - jd[0]) < 2 * 365.25:
            continue                                   # too short a forecast horizon to be meaningful
        k0, k1 = 0, 11
        O0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[k0:k1])))))) % 360
        w0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[k0:k1])))))) % 360
        e0, i0, aa = float(np.median(e[k0:k1])), float(np.median(inc[k0:k1])), float(np.median(a[k0:k1]))
        hp0 = aa * (1 - e0) - RE
        if hp0 < 350:
            continue
        traj = propagate(e0, i0, O0, w0, aa, jd[0], years=45.0, dt_days=10.0, hp_stop=300.0)
        hpp = aa * (1 - traj[:, 1]) - RE
        hit = np.where(hpp < 300.0)[0]
        t_pred = (traj[hit[0], 0] - jd[0]) / 365.25 if len(hit) else float('nan')
        t_obs = (jd[i300] - jd[0]) / 365.25
        out.append(dict(norad=o['norad'], i=float(i0), e0=e0, a0=aa, hp0=float(hp0), t_obs=float(t_obs), t_pred=float(t_pred)))
        print(o['norad'], f"i={i0:5.1f} e={e0:.3f} obs={t_obs:6.2f} pred={t_pred:6.2f}", flush=True)
    json.dump(out, open('E:/GitHub/OREM/scratch_rpe/combined_resonance/forecast_pool.json', 'w'), indent=1)
    t_obs = np.array([r['t_obs'] for r in out]); t_pred = np.array([r['t_pred'] for r in out])
    ok = np.isfinite(t_pred)
    print(f"\nobjects {len(out)}; model predicts a drop below 300 km within 45 yr for {ok.sum()}")
    d = t_pred[ok] - t_obs[ok]
    print(f"error (pred-obs): median {np.median(d):+.2f} yr, median |err| {np.median(np.abs(d)):.2f} yr, MAE {np.mean(np.abs(d)):.2f}")
    base = np.abs(t_obs[ok] - np.median(t_obs[ok]))
    print(f"baseline (sample-median lifetime) MAE {np.mean(base):.2f}  -> skill {1 - np.mean(np.abs(d)) / np.mean(base):.2f}")
    rel = np.abs(d) / t_obs[ok]
    print(f"relative error median {np.median(rel):.2f};  within 10%: {np.mean(rel<0.1):.2f}  within 25%: {np.mean(rel<0.25):.2f}")
    print(f"corr(log pred, log obs) = {np.corrcoef(np.log(t_pred[ok]), np.log(t_obs[ok]))[0,1]:.3f}")
    inc = np.array([r['i'] for r in out])[ok]
    for lo, hi in ((46.4, 56.0), (56.0, 62.0), (62.0, 66.0), (66.0, 90.0), (90.0, 181.0)):
        m = (inc >= lo) & (inc < hi)
        if m.sum() >= 3:
            print(f"  {lo:.0f}-{hi:.0f} deg: n={m.sum():3d}  median rel err {np.median(np.abs(d[m]) / t_obs[ok][m]):.2f}  within 25%: {np.mean(np.abs(d[m]) / t_obs[ok][m] < 0.25):.2f}")


if __name__ == '__main__':
    main()
