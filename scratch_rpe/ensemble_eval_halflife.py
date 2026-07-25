"""Offline evaluation of recency-weighted (half-life) ensemble estimator vs
the shipped 'latest zone' primary estimate and v1.20's other schemes,
against the CURRENT scratch_rpe/rpe_campaign.csv (50-object set, reflects
today's shipped pipeline incl. G3/trust-gate/noise-matched-apobs).

Unlike v1.20's ensemble_eval.py (7-object hardcoded ground-truth dict),
this derives t_obs per object directly from the CSV's own zepoch/rpe_pct
columns -- works for all 50 objects, no hardcoding, no staleness risk.

Weight bases (per issue #32 discussion, 2026-07-25):
  invlife    : w = 1/(reentry_jd - zepoch)          -- v1.20's w=1/life
  halfNN     : w = 2**(-(reentry_jd-zepoch)/NN)     -- new, smooth half-life
               decay in days, NN in {45,60,90,120,180}
"""
import csv, statistics

CURATED7 = {42928, 35497, 37151, 39615, 27526, 32007, 37819}

rows = {}
with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
    for r in csv.DictReader(f):
        n = int(r['norad'])
        rows.setdefault(n, []).append(r)

def derive_tobs(rs):
    """Back out true re-entry JD from any row with a real prediction+RPE."""
    vals = []
    for r in rs:
        if not r['reentry_jd'].strip() or not r['rpe_pct'].strip():
            continue
        pred = float(r['reentry_jd']); rp = float(r['rpe_pct'])
        zep = float(r['zepoch'])
        if pred > 0 and abs(rp) > 1e-9:
            vals.append((pred + zep*(rp/100.0)) / (1.0 + rp/100.0))
    if not vals:
        return None
    return statistics.median(vals)

halflives = [45, 60, 90, 120, 180]
schemes = ['uniform', 'latest', 'invlife'] + [f'half{h}' for h in halflives] + ['median']

agg_all = {k: [] for k in schemes}
agg_7   = {k: [] for k in schemes}
per_obj = []

for n, rs in rows.items():
    tobs = derive_tobs(rs)
    if tobs is None:
        continue
    zs = []
    for r in rs:
        if not r['reentry_jd'].strip():
            continue
        pred = float(r['reentry_jd']); iz = int(r['zone']); zep = float(r['zepoch'])
        if pred > 0:
            zs.append((iz, pred, zep))
    if not zs:
        continue
    zep1 = min(z[2] for z in zs)
    idxs  = [z[0] for z in zs]
    preds = [z[1] for z in zs]
    zeps  = [z[2] for z in zs]
    life  = [p - z for p, z in zip(preds, zeps)]

    vals = {}
    vals['uniform'] = sum(preds) / len(preds)
    vals['latest']  = max(zip(idxs, preds))[1]
    vals['median']  = statistics.median(preds)
    vals['invlife'] = (sum(p / l for p, l in zip(preds, life))
                        / sum(1.0 / l for l in life))
    for h in halflives:
        w = [2.0 ** (-l / h) for l in life]
        vals[f'half{h}'] = sum(p * wi for p, wi in zip(preds, w)) / sum(w)

    row_out = {'norad': n, 'nzones': len(zs)}
    for k in schemes:
        rpe = (vals[k] - tobs) / (tobs - zep1) * 100.0
        row_out[k] = rpe
        agg_all[k].append(abs(rpe))
        if n in CURATED7:
            agg_7[k].append(abs(rpe))
    per_obj.append(row_out)

def summarize(label, agg):
    print(f"\n--- {label} ---")
    print(f"{'scheme':>9} {'median':>8} {'mean':>8} {'max':>8}  n")
    for k in schemes:
        v = agg[k]
        if not v:
            continue
        print(f"{k:>9} {statistics.median(v):8.2f} {statistics.mean(v):8.2f} "
              f"{max(v):8.2f}  {len(v)}")

print(f"{'norad':>6} {'nz':>3} " + " ".join(f"{k:>9}" for k in schemes))
for r in sorted(per_obj, key=lambda x: x['norad']):
    tag = '*' if r['norad'] in CURATED7 else ' '
    print(f"{r['norad']:6d}{tag}{r['nzones']:3d} " +
          " ".join(f"{r[k]:9.1f}" for k in schemes))

summarize(f"ALL {len(per_obj)} objects", agg_all)
summarize(f"CURATED 7 only", agg_7)
