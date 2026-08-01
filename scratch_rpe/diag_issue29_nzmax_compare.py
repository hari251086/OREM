"""Issue #29 nzones_max sensitivity check -- same summary logic as
diag_issue29_current.py, parameterized by CSV path so it can be pointed at
the nz20/nz50 test campaign outputs as well as the shipped nzones_max=8
baseline (scratch_rpe/rpe_campaign.csv), for a direct before/after
comparison against today's full stack (G3, trust-gate, median ensemble,
zone_select recency, max_zone_days=20, Falcon9 swap, v1.47 perigee
indicator).

Usage: python scratch_rpe/diag_issue29_nzmax_compare.py <csv_path>
"""
import csv
import statistics
import sys
from collections import defaultdict

CURATED_7 = {27526, 32007, 35497, 37151, 37819, 39615, 42928}

csv_path = sys.argv[1] if len(sys.argv) > 1 else 'scratch_rpe/rpe_campaign.csv'

rows = defaultdict(list)
err_objs = []
with open(csv_path, newline='') as f:
    for r in csv.DictReader(f):
        norad = int(r['norad'])
        if r['zone'] == 'ERR':
            err_objs.append(norad)
            continue
        rows[norad].append(r)

all_objs = set(rows) | set(err_objs)
print(f"=== {csv_path} ===")
print(f"Total objects: {len(all_objs)} "
      f"({len(err_objs)} ERR/zero-valid-zone, {len(rows)} with >=1 zone)")


def classify(norad, rs):
    # sort by zepoch (chronological), not the 'zone' index column -- at
    # nzones_max>=10 the disposable test harness's I1 zone-index format
    # overflows to '*' for two-digit zone numbers (cosmetic only, doesn't
    # shift other fields; zepoch itself is unaffected and monotonic)
    zones = sorted(rs, key=lambda r: float(r['zepoch']))
    latest = zones[-1]
    predicts = float(latest['reentry_jd']) > 0
    latest_rpe = abs(float(latest['rpe_pct'])) if predicts else None
    ens_rpe = abs(float(latest['ens_rpe_pct']))
    return predicts, latest_rpe, ens_rpe


def summarize(label, norads):
    total = len(norads)
    predict_latest = []
    ens_vals = []
    zero_zone = 0
    for n in norads:
        if n in err_objs:
            zero_zone += 1
            continue
        predicts, latest_rpe, ens_rpe = classify(n, rows[n])
        if predicts:
            predict_latest.append(latest_rpe)
        ens_vals.append(ens_rpe)

    print(f"--- {label} (n={total}) ---")
    print(f"  zero-valid-zone (ERR): {zero_zone}/{total}")
    print(f"  latest-zone predicts:  {len(predict_latest)}/{total}")
    if predict_latest:
        print(f"  latest-zone |RPE|: mean {statistics.mean(predict_latest):6.2f}%  "
              f"median {statistics.median(predict_latest):6.2f}%  "
              f"max {max(predict_latest):7.1f}%")
    nonzero_ens = [v for v in ens_vals if v > 0]
    if nonzero_ens:
        print(f"  ensemble |RPE|: mean {statistics.mean(nonzero_ens):6.2f}%  "
              f"median {statistics.median(nonzero_ens):6.2f}%  "
              f"max {max(nonzero_ens):7.1f}%  n={len(nonzero_ens)}")
    print()


rest = sorted(all_objs - CURATED_7)
summarize("Curated-7", sorted(CURATED_7))
summarize("Other 43", rest)
summarize("Full 50", sorted(all_objs))
