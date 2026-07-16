"""Offline evaluation of ensemble weighting schemes on the 7-object campaign
(issue #16, v1.20). Uses scratch_rpe/rpe_campaign.csv; recovers zone epochs
from the per-zone RPE definition. Result (|ensemble RPE| across objects):

  scheme     median    mean     max
  uniform    15.28%   45.29%  238.7%   (current t_mean)
  latest      8.18%    7.60%   14.4%   <-- winner: ALL objects within 15%
  w=idx      12.12%   31.73%  150.5%
  w=1/life    7.84%   12.85%   39.3%
  median      8.73%   38.31%  216.4%

The latest zone is closest to re-entry (shortest extrapolation) and its
fitted BN reflects the most recent attitude/altitude regime -- including
for 35497, where it lands -0.2% while earlier zones run +170..+520%.
"""
import csv, statistics

def cal2jd(y, m, d):
    if m <= 2: y -= 1; m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25*(y+4716)) + int(30.6001*(m+1)) + d + b - 1524.5

obs = {42928: cal2jd(2019,3,3), 35497: cal2jd(2016,10,31), 37151: cal2jd(2015,12,3),
       39615: cal2jd(2017,9,15), 27526: cal2jd(2012,5,9), 32007: cal2jd(2010,6,6),
       37819: cal2jd(2013,9,12)}

rows = {}
with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
    for r in csv.DictReader(f):
        n = int(r['norad']); rows.setdefault(n, []).append(r)

def rpe(pred, n, zep1):
    return (pred - obs[n]) / (obs[n] - zep1) * 100.0

agg = {k: [] for k in ['uniform','latest','idx','invlife','median']}
print(f"{'obj':>6} {'uniform':>9} {'latest':>9} {'w=idx':>9} {'w=1/life':>10} {'median':>9}")
for n, rs in rows.items():
    zs = []
    for r in rs:
        pred = float(r['reentry_jd']); rp = float(r['rpe_pct']); iz = int(r['zone'])
        if pred > 0 and abs(rp) > 1e-9:
            zep = obs[n] - (pred - obs[n]) / (rp/100.0)
            zs.append((iz, pred, zep))
    if not zs: continue
    zep1 = min(z[2] for z in zs)
    preds = [z[1] for z in zs]; idxs = [z[0] for z in zs]; zeps = [z[2] for z in zs]
    vals = dict(
        uniform=sum(preds)/len(preds),
        latest=max(zip(idxs,preds))[1],
        idx=sum(i*p for i,p in zip(idxs,preds))/sum(idxs),
        invlife=sum(p/(p-z) for p,z in zip(preds,zeps))/sum(1/(p-z) for p,z in zip(preds,zeps)),
        median=statistics.median(preds))
    out = [f"{n:6d}"]
    for k in agg:
        e = rpe(vals[k], n, zep1); agg[k].append(abs(e)); out.append(f"{e:9.1f}")
    print(" ".join(out))
print("\n|ens RPE| summary across objects:")
for k in agg:
    print(f"  {k:9s} median {statistics.median(agg[k]):7.2f}%   mean {statistics.mean(agg[k]):7.2f}%   max {max(agg[k]):7.1f}%")
