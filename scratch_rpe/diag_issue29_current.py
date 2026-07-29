"""Issue #29 re-check against the current (v1.46) stack.

Every prior #29 measurement in this investigation's history predates the
v1.46 revert of #31's apobs fix -- re-verify the generalization picture
(predict rate, curated-7 vs. broader-set accuracy) fresh rather than citing
a stale figure, per this project's own "always re-verify before citing an
old number" discipline (see project_orem_rpe_investigation_plan.md).

Reads scratch_rpe/rpe_campaign.csv directly (the 50-object campaign,
regenerated same session as the v1.46 revert). 'latest zone' RPE = rpe_pct
of the chronologically-last zone IF it predicts (reentry_jd>0) -- objects
whose latest zone doesn't predict are excluded from the latest-zone stat
(matches production's own primary-estimator semantics: it doesn't fall
back to an earlier zone). Ensemble RPE = the ens_rpe_pct column, already
computed by production (median-of-valid-zones formula, v1.44+).
"""
import csv
import statistics
from collections import defaultdict

CURATED_7 = {27526, 32007, 35497, 37151, 37819, 39615, 42928}

rows = defaultdict(list)
err_objs = []
with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
    for r in csv.DictReader(f):
        norad = int(r['norad'])
        if r['zone'] == 'ERR':
            err_objs.append(norad)
            continue
        rows[norad].append(r)

all_objs = set(rows) | set(err_objs)
print(f"Total objects: {len(all_objs)} "
      f"({len(err_objs)} ERR/zero-valid-zone, {len(rows)} with >=1 zone)")
print(f"ERR objects: {sorted(err_objs)}\n")


def classify(norad, rs):
    zones = sorted(rs, key=lambda r: int(r['zone']))
    latest = zones[-1]
    predicts = float(latest['reentry_jd']) > 0
    latest_rpe = abs(float(latest['rpe_pct'])) if predicts else None
    ens_rpe = abs(float(latest['ens_rpe_pct']))
    has_any_prediction = any(float(z['reentry_jd']) > 0 for z in zones)
    return predicts, latest_rpe, ens_rpe, has_any_prediction


def summarize(label, norads):
    total = len(norads)
    predict_latest = []
    predict_any = 0
    ens_vals = []
    zero_zone_no_predict = 0
    for n in norads:
        if n in err_objs:
            zero_zone_no_predict += 1
            continue
        predicts, latest_rpe, ens_rpe, has_any = classify(n, rows[n])
        if has_any:
            predict_any += 1
        if predicts:
            predict_latest.append(latest_rpe)
        ens_vals.append(ens_rpe)  # rpe_ens computed even if latest doesn't predict, as long as >=1 zone predicts; 0 if none

    print(f"--- {label} (n={total}) ---")
    print(f"  zero-valid-zone (ERR): {zero_zone_no_predict}/{total}")
    print(f"  >=1 zone predicts:     {predict_any}/{total}")
    print(f"  latest-zone predicts:  {len(predict_latest)}/{total} "
          f"(gate production actually uses as primary estimator)")
    if predict_latest:
        print(f"  latest-zone |RPE|: mean {statistics.mean(predict_latest):6.2f}%  "
              f"median {statistics.median(predict_latest):6.2f}%  "
              f"max {max(predict_latest):7.1f}%")
    nonzero_ens = [v for v in ens_vals if v > 0]
    if nonzero_ens:
        print(f"  ensemble |RPE| (obj w/ >=1 predicting zone): "
              f"mean {statistics.mean(nonzero_ens):6.2f}%  "
              f"median {statistics.median(nonzero_ens):6.2f}%  "
              f"max {max(nonzero_ens):7.1f}%  n={len(nonzero_ens)}")
    print()


rest = sorted(all_objs - CURATED_7)
summarize("Curated-7", sorted(CURATED_7))
summarize("Other 43", rest)
summarize("Full 50", sorted(all_objs))

# Per-object listing for the "other 43" set, to spot patterns
print("--- Other-43 per-object detail (norad: zero-zone / no-predict / latest_rpe / ens_rpe) ---")
for n in rest:
    if n in err_objs:
        print(f"  {n}: ERR (zero valid zones)")
        continue
    predicts, latest_rpe, ens_rpe, has_any = classify(n, rows[n])
    tag = "predicts" if predicts else ("has-earlier-pred" if has_any else "no-predict")
    lr = f"{latest_rpe:6.1f}%" if latest_rpe is not None else "   n/a"
    print(f"  {n}: {tag:17s} latest={lr}  ens={ens_rpe:6.1f}%")
