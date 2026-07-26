"""Issue #35, Phase 1: compare the trust-gated BN-carryover chain (baseline,
src/orem.F as shipped) against removing it entirely (test variant,
scratch_rpe/orem_notrustgate_test.F), across all three validation campaigns,
under TODAY's codebase (post v1.44 median metric, post 4952a34 Falcon9 swap).

Uses the proper reentry_jd>0 predict gate (not the ens_rpe_pct==0 sentinel,
a methodology bug this project has already been bitten by once) and reports
mean/median |ensemble RPE|, mean |latest-zone RPE| (primary metric), and
predict rate for each campaign, plus a per-object table so individual
blowups/improvements are visible, not just aggregates.
"""
import csv
import statistics
from collections import defaultdict


def load_rows(path):
    rows = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            if row['zone'] == 'ERR':
                continue
            rows[row['norad']].append(row)
    return rows


def metrics(rows):
    """Return dict: norad -> (ens_rpe, latest_rpe) for objects with >=1
    predicting zone (reentry_jd>0)."""
    out = {}
    for n, rs in rows.items():
        valid = [r for r in rs if float(r['reentry_jd']) > 0]
        if not valid:
            continue
        ens = None
        for r in rs:
            if r.get('ens_rpe_pct', '').strip():
                ens = float(r['ens_rpe_pct'])
                break
        if ens is None:
            continue
        last_zone = max(valid, key=lambda r: int(r['zone']))
        latest = float(last_zone['rpe_pct'])
        out[n] = (ens, latest)
    return out


def summarize(name, baseline_path, test_path):
    base_rows = load_rows(baseline_path)
    test_rows = load_rows(test_path)
    base_m = metrics(base_rows)
    test_m = metrics(test_rows)

    print(f"\n{'='*70}\n{name}\n{'='*70}")
    print(f"predict rate: baseline={len(base_m)}  notrustgate={len(test_m)}")

    common = sorted(set(base_m) & set(test_m), key=int)
    only_base = sorted(set(base_m) - set(test_m), key=int)
    only_test = sorted(set(test_m) - set(base_m), key=int)
    if only_base:
        print(f"  predicted in baseline only: {only_base}")
    if only_test:
        print(f"  predicted in notrustgate only: {only_test}")

    print(f"\n{'norad':>7} {'ens_base':>10} {'ens_notg':>10} {'d_ens':>8}  "
          f"{'lat_base':>10} {'lat_notg':>10} {'d_lat':>8}")
    rows_out = []
    for n in common:
        eb, lb = base_m[n]
        et, lt = test_m[n]
        rows_out.append((n, eb, et, et - eb, lb, lt, lt - lb))
    rows_out.sort(key=lambda r: abs(r[3]), reverse=True)
    for n, eb, et, de, lb, lt, dl in rows_out:
        flag = "  <<<" if abs(de) > 20 or abs(dl) > 20 else ""
        print(f"{n:>7} {eb:>10.2f} {et:>10.2f} {de:>8.2f}  "
              f"{lb:>10.2f} {lt:>10.2f} {dl:>8.2f}{flag}")

    base_ens_abs = [abs(v[0]) for v in base_m.values()]
    test_ens_abs = [abs(v[0]) for v in test_m.values()]
    base_lat_abs = [abs(v[1]) for v in base_m.values()]
    test_lat_abs = [abs(v[1]) for v in test_m.values()]

    print(f"\nmean|ensemble RPE|   baseline={statistics.mean(base_ens_abs):.2f}%  "
          f"notrustgate={statistics.mean(test_ens_abs):.2f}%")
    print(f"median|ensemble RPE| baseline={statistics.median(base_ens_abs):.2f}%  "
          f"notrustgate={statistics.median(test_ens_abs):.2f}%")
    print(f"max|ensemble RPE|    baseline={max(base_ens_abs):.2f}%  "
          f"notrustgate={max(test_ens_abs):.2f}%")
    print(f"mean|latest RPE|     baseline={statistics.mean(base_lat_abs):.2f}%  "
          f"notrustgate={statistics.mean(test_lat_abs):.2f}%")
    print(f"median|latest RPE|   baseline={statistics.median(base_lat_abs):.2f}%  "
          f"notrustgate={statistics.median(test_lat_abs):.2f}%")


if __name__ == "__main__":
    summarize(
        "Curated-7 static (rpe_campaign_7obj_issue31)",
        "scratch_rpe/rpe_campaign_7obj_issue31_pretrustgateremoval_backup.csv",
        "scratch_rpe/rpe_campaign_7obj_issue31_notrustgate_result.csv",
    )
    summarize(
        "Curated-7 weather+bulge (rpe_campaign_weather)",
        "scratch_rpe/rpe_campaign_weather_pretrustgateremoval_backup.csv",
        "scratch_rpe/rpe_campaign_weather_notrustgate_result.csv",
    )
    summarize(
        "All-50 (rpe_campaign)",
        "scratch_rpe/rpe_campaign_pretrustgateremoval_backup.csv",
        "scratch_rpe/rpe_campaign_notrustgate_result.csv",
    )
