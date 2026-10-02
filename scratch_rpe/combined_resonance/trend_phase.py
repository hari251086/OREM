"""Does the starting nodal/apsidal phase set the timing of the perigee hump? (23 objects with i_mean >= 55 from the report)"""
import sys, json
import numpy as np
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import load_full_record   # noqa: E402
from validate_drive import clean_rows                     # noqa: E402
from trend_all28 import ROOT                              # noqa: E402

R = [r for r in json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/trend_all28.json')) if r['i'] >= 55 and r['norad'] != 12907]
for r in R:
    rows = np.array(clean_rows(load_full_record(f"{ROOT}/{r['norad']}/{r['norad']}.tle.txt")))
    r['O0'] = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(rows[:5, 4]))))) % 360)
    r['w0'] = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(rows[:5, 5]))))) % 360)
    r['jd0'] = float(rows[0, 0])
print(f"{'NORAD':>6} {'Omega0':>7} {'omega0':>7} {'t_pk':>5} {'pk->end':>7} {'rise':>6} {'T_end':>6}")
for r in sorted(R, key=lambda r: r['O0']):
    print(f"{r['norad']:>6} {r['O0']:7.1f} {r['w0']:7.1f} {r['t_peak']:5.1f} {r['peak_to_end']:7.1f} {r['rise']:6.0f} {r['t_end']:6.1f}")
O0 = np.array([r['O0'] for r in R]); pe = np.array([r['peak_to_end'] for r in R]); tp = np.array([r['t_peak'] for r in R]); rise = np.array([r['rise'] for r in R])
# circular-linear association: regress on cos/sin of Omega0
def circ_r2(y, ang):
    X = np.column_stack([np.ones_like(ang), np.cos(np.radians(ang)), np.sin(np.radians(ang))])
    beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2), beta
for lab, y in (('peak->end (yr)', pe), ('time of peak (yr)', tp), ('rise (km)', rise)):
    r2, b = circ_r2(y, O0)
    print(f"{lab:<18} R^2 vs first-TLE Omega (cos/sin fit) = {r2:.2f}   n={len(y)}")
