"""Extra mid-life forecast tests on the two other pool objects whose perigee dips below 300 km, recovers, and the record continues
> 2 yr (found by scan_second_hump.py): 42908 and 6231.  Forecast only (OREM needs per-object parameters not cached for these)."""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from resonance_campaign_sweep import cal2jd                                         # noqa: E402
from validate_drive import clean_rows, RE                                           # noqa: E402
from molniya_cycle import smooth                                                    # noqa: E402
from secular_propagator import propagate                                            # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
gap = {thr: pd.read_csv(f'{D}/cross_check_orem_{thr}.csv')['gap'] for thr in (200, 150)}
EXTRA = [42908, 6231]
objs = {o['norad']: o for o in build_object_list()}
for n in EXTRA:
    o = objs[n]
    rec = np.array(clean_rows(load_full_record(f"{BASE}/{o['tle_file']}"))); jd, a, e, inc, O, w = rec.T
    hp = a * (1 - e) - RE; s = smooth(hp); ipk = int(np.argmax(s))
    k = int(np.where((np.arange(len(hp)) > ipk) & (s < 300.0))[0][0])
    decay = cal2jd(*o['decay_ymd'])
    lead = decay - jd[k]
    sl = slice(max(0, k - 10), k + 1)
    e0, i0, aa = float(np.median(e[sl])), float(np.median(inc[sl])), float(np.median(a[sl]))
    O0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360)
    w0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360)
    traj = propagate(e0, i0, O0, w0, aa, jd[k], years=15.0, dt_days=5.0)
    hpp = aa * (1 - traj[:, 1]) - RE
    print(f"NORAD {n}: i={i0:.1f} e={e0:.3f} cutoff (first dip) JD {jd[k]:.1f}  decay in {lead:.0f} d ({lead/365.25:.2f} yr)  perigee at cutoff {hp[k]:.0f} km, model max next 8 yr {hpp[:int(8*365.25/5)].max():.0f} km")
    for thr in (200, 150):
        hit = np.where(hpp < thr)[0]
        if len(hit):
            est = traj[hit[0], 0] - traj[0, 0] + float(np.median(gap[thr]))
            print(f"   {thr} km event at +{traj[hit[0], 0] - traj[0, 0]:.0f} d -> decay estimate +{est:.0f} d, error {est - lead:+.0f} d")
        else:
            print(f"   {thr} km: model perigee never falls below within 15 yr")
