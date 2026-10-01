"""Validate resonance_detector.detect() against the 2026-09-25 manual review
(issue #56) over the full 253-object pool. Ground truth is the user's own
visual verdicts, not another detector run.

Policy-excluded objects (3 Falcon 9 R/Bs, 13913 bad decay date) are
reported but not scored: their exclusion is independent of the detector.
"""
import math
import sys

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
from prepare_plot_data_fullpool import (  # noqa: E402
    BASE, build_object_list, load_full_record,
)
from resonance_detector import detect  # noqa: E402

GENUINE = {40943, 59347, 35497, 10724, 37151, 27526, 18954}
POLICY_EXCLUDED = {39501, 44187, 42985, 13913}
# 19881 (i=59.5) was "uncertain" in the manual review; it is out of domain
# for the resonance, so it is neither a genuine positive nor scored.


def main():
    objs = build_object_list()
    res = {}
    for o in objs:
        try:
            rows = load_full_record(f"{BASE}/{o['tle_file']}")
        except FileNotFoundError:
            continue
        if len(rows) < 3:
            continue
        res[o['norad']] = detect([(r[0], r[1], r[2], r[3]) for r in rows])

    from collections import Counter
    print('status counts:', dict(Counter(v['status'] for v in res.values())))

    flagged = {n for n, v in res.items() if v['status'] == 'RESONANT'}
    scored_flagged = flagged - POLICY_EXCLUDED
    tp = scored_flagged & GENUINE
    fp = scored_flagged - GENUINE
    fn = GENUINE - flagged
    print(f'\nscored objects: {len(res)} (policy-excluded not scored: '
          f'{sorted(POLICY_EXCLUDED & set(res))})')
    print(f'genuine (manual): {len(GENUINE)}   detected: {len(tp)} '
          f'-> recall {len(tp)}/{len(GENUINE)}')
    print(f'false positives: {len(fp)} {sorted(fp)}')
    print(f'missed genuine : {len(fn)} '
          f'{[(n, res[n]["status"]) for n in sorted(fn)]}')
    print(f'precision: {len(tp)}/{len(scored_flagged)}')
    print('\nflagged policy-excluded (informational):',
          sorted(flagged & POLICY_EXCLUDED))
    print('\nper-genuine detail:')
    for n in sorted(GENUINE):
        v = res.get(n, {})
        print(' ', n, v.get('status'), 'i=%.1f' % v.get('i_mean', float('nan')),
              'n=', v.get('n_used'))


if __name__ == '__main__':
    main()
