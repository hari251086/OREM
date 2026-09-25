"""Build post-resonance-only TLE input files for the 16 objects flagged
RESONANT by resonance_campaign_sweep.py, implementing the hard rule:
"run re-entry predictions only starting from the u-turn resonance period,
not before." Truncates each object's real TLE history to only the
records at/after its detected resonance epoch (last lambda/lambda_crit
crossing) -- OREM's zone_select/GA pipeline then only ever sees
post-resonance data, never anything from the pre-resonance regime.
"""
import csv

from resonance_campaign_sweep import load_record, parse_epoch_jd  # noqa: F401

BASE = 'E:/GitHub/OREM'


def main():
    with open(BASE + '/scratch_rpe/resonance_sweep_97obj.csv') as f:
        rows = list(csv.DictReader(f))

    summary = []
    for r in rows:
        if r['status'] != 'RESONANT':
            continue
        norad = r['norad']
        src = BASE + '/' + r['tle_file']
        resonance_jd = float(r['resonance_jd'])
        with open(src) as f:
            lines = [l.rstrip('\n') for l in f]

        kept_pairs = []
        for k in range(0, len(lines) - 1, 2):
            l1, l2 = lines[k], lines[k + 1]
            if not (l1.startswith('1 ') and l2.startswith('2 ')):
                continue
            jd = parse_epoch_jd(l1)
            if jd >= resonance_jd:
                kept_pairs.append((jd, l1, l2))
        kept_pairs.sort()

        dst = f'{BASE}/input/example_{norad}_postresonance.tle.txt'
        with open(dst, 'w') as f:
            for jd, l1, l2 in kept_pairs:
                f.write(l1 + '\n')
                f.write(l2 + '\n')

        summary.append((norad, r['name'], len(kept_pairs), dst))
        print(f'{norad:>7}  {r["name"][:28]:<28}  '
              f'{len(kept_pairs):4d} post-resonance TLE(s) -> {dst}')

    n_runnable = sum(1 for s in summary if s[2] >= 8)
    n_thin = sum(1 for s in summary if 0 < s[2] < 8)
    n_empty = sum(1 for s in summary if s[2] == 0)
    print(f'\n{len(summary)} resonant objects total: '
          f'{n_runnable} have >=8 post-resonance TLEs (zone_select needs '
          f'min_pts=8), {n_thin} have 1-7 (too few for a zone), '
          f'{n_empty} have ZERO (resonance detected at/after the last '
          f'tracked TLE -- no post-resonance data exists at all; the hard '
          f'rule means these cannot be run).')


if __name__ == '__main__':
    main()
