import json
import numpy as np
from scipy.stats import spearmanr

recs = json.load(open('E:/GitHub/OREM/scratch_rpe/combined_resonance/drive_windows.json'))
rng = np.random.default_rng(3)

rows = []
for r in recs:
    obs = np.array(r['obs']); pred = np.array(r['ps']) + np.array(r['pm'])
    if len(obs) < 8 or np.std(obs) == 0 or np.std(pred) == 0:
        continue
    rho = spearmanr(obs, pred)[0]
    c = np.corrcoef(obs, pred)[0, 1]
    beta = float(np.sum(obs * pred) / np.sum(pred * pred))
    # signal level: windows where either series is non-trivial
    amp = np.median(np.abs(obs))
    big = np.abs(obs) > amp
    sg = float(np.mean(np.sign(obs[big]) == np.sign(pred[big]))) if big.sum() else np.nan
    # lag-shift null: circularly shift predicted by a random lag (preserves its autocorrelation)
    nulls = []
    for _ in range(200):
        k = rng.integers(3, len(pred) - 2) if len(pred) > 6 else 1
        nulls.append(np.corrcoef(obs, np.roll(pred, k))[0, 1])
    rows.append((r['norad'], r['i'], len(obs), c, rho, beta, sg, np.mean(np.abs(np.array(nulls)) >= abs(c)), amp))
a = np.array(rows, dtype=float)
print(f"objects: {len(a)}")
print(f"per-object Pearson   : median {np.median(a[:,3]):+.3f}  IQR [{np.percentile(a[:,3],25):+.3f}, {np.percentile(a[:,3],75):+.3f}]  frac>0.5: {np.mean(a[:,3]>0.5):.2f}")
print(f"per-object Spearman  : median {np.median(a[:,4]):+.3f}  IQR [{np.percentile(a[:,4],25):+.3f}, {np.percentile(a[:,4],75):+.3f}]")
print(f"per-object slope     : median {np.median(a[:,5]):+.2f}  IQR [{np.percentile(a[:,5],25):+.2f}, {np.percentile(a[:,5],75):+.2f}]")
print(f"sign agreement (|obs|>median): median {np.nanmedian(a[:,6]):.2f}")
print(f"lag-shift null: fraction of objects with p<0.05 = {np.mean(a[:,7]<0.05):.2f}")
allo = np.concatenate([np.array(r['obs']) for r in recs]); allp = np.concatenate([np.array(r['ps']) + np.array(r['pm']) for r in recs])
print(f"pooled Spearman      : {spearmanr(allo, allp)[0]:+.3f}")
print(f"median |obs| dq/dt = {np.median(np.abs(allo)):.3f} km/day;  95th pct = {np.percentile(np.abs(allo),95):.2f} km/day")
print("amplitude concentration: top-5 objects' share of sum(obs^2):",
      f"{np.sort(np.array([np.sum(np.array(r['obs'])**2) for r in recs]))[-5:].sum()/np.sum(allo**2):.2f}")
