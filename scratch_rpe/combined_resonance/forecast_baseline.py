"""Is the forecast skill trivial?  Compare the secular propagator (first-TLE start, J2+Sun+Moon) against naive
baselines for the epoch the perigee falls below 300 km (E300): (a) every object has the sample-median lifetime;
(b) lifetime scales with nothing (same).  Also lists per-object results."""
import json
import numpy as np

R = json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/forecast_results.json'))
sel = [r for r in R if r['start'] == 'first TLE' and r['variant'] == 'J2+Sun+Moon' and r['e300_err_yr'] is not None]
h = np.array([r['horizon_yr'] for r in sel])
err = np.array([r['e300_err_yr'] for r in sel])
pk = np.array([r['pk_err_yr'] for r in sel])
print(f"objects: {len(sel)}")
print(f"observed horizon first-TLE -> E300: min {h.min():.2f}  median {np.median(h):.2f}  max {h.max():.2f} yr  (std {h.std():.2f})")
print(f"baseline (sample-median lifetime):   median |err| = {np.median(np.abs(h - np.median(h))):.2f} yr   MAE = {np.mean(np.abs(h - np.median(h))):.2f}")
print(f"secular propagator:                  median |err| = {np.median(np.abs(err)):.2f} yr   MAE = {np.mean(np.abs(err)):.2f}   max |err| = {np.max(np.abs(err)):.2f}")
print(f"skill (1 - MAE_model/MAE_baseline) = {1 - np.mean(np.abs(err)) / np.mean(np.abs(h - np.median(h))):.2f}")
print(f"corr(predicted horizon, observed horizon) = {np.corrcoef(h + err, h)[0,1]:.3f}")
print("\nper object (first-TLE forecast): NORAD  observed_horizon  E300_err(yr)  peak_err(yr)")
for r in sorted(sel, key=lambda r: r['horizon_yr']):
    print(f"   {r['norad']:>6}  {r['horizon_yr']:6.2f}  {r['e300_err_yr']:+6.2f}  {r['pk_err_yr']:+6.2f}")
