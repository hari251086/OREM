"""Diagnostic (user request): 12 objects came back from the batched
gp_history fetch (fetch_new_candidates.py, 40-NORAD-ID-per-request
chunks) with only 1-2 TLEs each -- suspiciously low for objects that
aren't brand-new launches. Re-fetch each ONE AT A TIME (no batching,
generous pause between requests) to check whether the batched query
silently truncated/dropped rows for some objects in the chunk, or
whether Space-Track genuinely only has that few TLEs on file.
"""
import sys
import time

sys.path.insert(0, r"C:\Users\hari2\OneDrive\Documents\GitHub\OREM-Watchlist\src")
from heowatch.spacetrack_client import get_client  # noqa: E402

# (norad, name, n_tle from the batched fetch)
SUSPECT = [
    (17, 'THOR ABLE R/B', 1),
    (432, 'EXPLORER 14 (SERB 53A)', 2),
    (1460, 'ERS 17', 2),
    (6068, 'PROGNOZ 2', 1),
    (11700, 'COSMOS 1164', 1),
    (28014, 'SL-12 DEB', 1),
    (28015, 'SL-12 DEB', 2),
    (28037, 'SL-12 DEB', 2),
    (31801, 'CZ-3B R/B', 1),
    (40006, 'DELTA 4 DEB', 1),
    (68539, 'TACHELES', 1),
    (68541, 'K-RADCUBE', 1),
]

WAIT_SECONDS = 5  # generous pause between individual requests


def main():
    st = get_client()
    print(f'{"NORAD":>7} {"name":<26} {"batched":>8} {"solo":>6}  verdict')
    results = []
    for norad, name, batched_n in SUSPECT:
        tle_text = st.gp_history(norad_cat_id=norad, orderby='epoch', format='tle')
        lines = [l for l in tle_text.splitlines() if l.startswith('1 ') or l.startswith('2 ')]
        n_solo = len(lines) // 2
        verdict = 'MATCH (real sparse data)' if n_solo <= batched_n else \
                  'BATCH TRUNCATED -- solo fetch got more!'
        results.append((norad, name, batched_n, n_solo, verdict))
        print(f'{norad:>7} {name:<26} {batched_n:>8} {n_solo:>6}  {verdict}')
        time.sleep(WAIT_SECONDS)

    n_bug = sum(1 for r in results if r[3] > r[2])
    print(f'\n{n_bug} of {len(results)} objects got MORE TLEs when fetched solo '
          f'-- {"batching bug confirmed" if n_bug else "no evidence of a batching bug"}.')

    import csv
    with open('scratch_rpe/diag_low_tle_refetch_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['norad', 'name', 'n_tle_batched', 'n_tle_solo', 'verdict'])
        w.writerows(results)


if __name__ == '__main__':
    main()
