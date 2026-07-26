"""Issue #29/#32 follow-up (61734 zone-1 investigation): re-run the v1.20
ensemble-scheme comparison (scratch_rpe/ensemble_eval.py, issue #16) against
the CURRENT 50-object campaign instead of the original 7. So much has
changed since v1.20 (G3 BSTAR prior, zone_select recency fix, diurnal bulge,
max_zone_days 10->20, BN-carryover boundary-recenter fix, #31 noise-matched
apobs) that the original "latest zone wins decisively, median beats uniform
mean" finding is worth re-checking rather than assumed to still hold.

Unlike the original script (which hardcoded 7 real observed dates), this one
recovers each object's true observed re-entry date algebraically from the
CSV's own (reentry_jd, rpe_pct, zepoch) triple per zone -- rpe_pct =
(reentry_jd - t_obs_jd)/(t_obs_jd - zepoch)*100, solved for t_obs_jd. Uses
whichever zone has the largest |rpe_pct| (best-conditioned solve) as the
source, then verifies other zones of the same object agree (a real physical
date should recover consistently regardless of which zone's equation is used
-- this is also a free internal-consistency check on the CSV itself).
"""
import csv
import statistics
from collections import defaultdict

def recover_t_obs(reentry_jd, rpe_pct, zepoch):
    # rpe_pct/100 = (reentry_jd - t_obs)/(t_obs - zepoch)
    # => (t_obs - zepoch)*rpe_pct/100 = reentry_jd - t_obs
    # => t_obs*(rpe_pct/100 + 1) = reentry_jd + zepoch*rpe_pct/100
    r = rpe_pct / 100.0
    return (reentry_jd + zepoch * r) / (1.0 + r)

rows = defaultdict(list)
with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
    for r in csv.DictReader(f):
        if r['zone'] == 'ERR' or not r.get('zepoch', '').strip():
            continue
        rows[int(r['norad'])].append(r)

agg = {k: [] for k in ['uniform', 'median', 'latest']}
n_recovered = 0
n_inconsistent = 0

print(f"{'norad':>7} {'nzones':>7} {'t_obs_recovered':>16} "
      f"{'uniform':>9} {'median':>9} {'latest':>9}")

for norad, rs in sorted(rows.items()):
    valid = [r for r in rs if float(r['reentry_jd']) > 0
             and abs(float(r['rpe_pct'])) > 1e-6]
    if not valid:
        continue

    # recover t_obs from every valid zone, check consistency
    recovered = []
    for r in valid:
        t = recover_t_obs(float(r['reentry_jd']), float(r['rpe_pct']),
                           float(r['zepoch']))
        recovered.append(t)
    t_obs = statistics.median(recovered)
    spread = max(recovered) - min(recovered)
    if spread > 0.01:
        n_inconsistent += 1
    n_recovered += 1

    preds = [float(r['reentry_jd']) for r in valid]
    idxs = [int(r['zone']) for r in valid]
    zep1 = min(float(r['zepoch']) for r in valid)
    horiz = t_obs - zep1
    if abs(horiz) < 0.1:
        continue

    def rpe_of(pred):
        return (pred - t_obs) / horiz * 100.0

    vals = dict(
        uniform=sum(preds) / len(preds),
        median=statistics.median(preds),
        latest=preds[idxs.index(max(idxs))],
    )
    out = [f"{norad:7d}", f"{len(valid):7d}", f"{t_obs:16.4f}"]
    for k in ('uniform', 'median', 'latest'):
        e = rpe_of(vals[k])
        agg[k].append(abs(e))
        out.append(f"{e:9.1f}")
    print(" ".join(out))

print(f"\n(t_obs recovery: {n_recovered} objects, "
      f"{n_inconsistent} with >0.01-day inter-zone inconsistency)")
print("\n|ensemble RPE| summary across objects (current 50-object campaign):")
for k in ('uniform', 'median', 'latest'):
    vs = agg[k]
    print(f"  {k:9s} median {statistics.median(vs):7.2f}%   "
          f"mean {statistics.mean(vs):7.2f}%   max {max(vs):7.1f}%   n={len(vs)}")
