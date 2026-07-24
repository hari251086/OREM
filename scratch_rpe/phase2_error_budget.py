"""Phase 2 of the RPE investigation (issue #32): cheap, correlational error-
budget decomposition reusing the 30-object campaign CSV (rpe_campaign.csv,
regenerated with zepoch/rms_fit columns added to rpe_campaign.F -- both
already computed by orem_run, no new propagation) plus the on-disk TLE
files' own BSTAR field.

Three sub-analyses:
  A. BC/attitude variability: does the fitted BN track the TLE's own
     published BSTAR trend within an object (supports/refutes the
     literature finding that BC should not be treated as constant)?
  B. RPE-metric validity: does |RPE| inflate for short-horizon zones
     (a metric artifact) independent of real fit quality?
  C. Short-window identifiability floor: does the GA fit RMS predict
     prediction accuracy, and how does its distribution compare to the
     TLE noise floor already established (issue #31 investigation:
     independent per-TLE SGP4 conversion sets a ~5-6 km noise floor,
     which manifested as post-fix RMS of 0.17-0.51 km on the curated 7 --
     i.e. RMS there is NOT hitting the raw per-point noise floor, it's a
     fitted regression average over ~8-10 points).
"""
import csv, math

NORAD_FILE = {
    42928: 'input/example_42928.tle.txt', 35497: 'input/example_35497.tle.txt',
    37151: 'input/example_37151.tle.txt', 39615: 'input/example_39615.tle.txt',
    27526: 'input/example_27526.tle.txt', 32007: 'input/example_32007.tle.txt',
    37819: 'input/example_37819.tle.txt', 11550: 'input/example_11550.tle.txt',
    59347: 'input/example_59347.tle.txt', 40943: 'input/example_40943.tle.txt',
    66587: 'input/example_66587.tle.txt', 60328: 'input/example_60328.tle.txt',
    61734: 'input/example_61734.tle.txt', 57804: 'input/example_57804.tle.txt',
    56758: 'input/example_56758.tle.txt', 30799: 'input/example_30799.tle.txt',
    44187: 'input/example_44187.tle.txt', 35009: 'input/example_35009.tle.txt',
    41553: 'input/example_41553.tle.txt', 48259: 'input/example_48259.tle.txt',
    27906: 'input/example_27906.tle.txt', 27882: 'input/example_27882.tle.txt',
    39802: 'input/example_39802.tle.txt', 28572: 'input/example_28572.tle.txt',
    46429: 'input/example_46429.tle.txt', 44591: 'input/example_44591.tle.txt',
    41695: 'input/example_41695.tle.txt', 52205: 'input/example_52205.tle.txt',
    45349: 'input/example_45349.tle.txt', 23647: 'input/example_23647.tle.txt',
}

OBS_DECAY = {
    42928: (2019, 3, 3), 35497: (2016, 10, 31), 37151: (2015, 12, 3),
    39615: (2017, 9, 15), 27526: (2012, 5, 9), 32007: (2010, 6, 6),
    37819: (2013, 9, 12), 11550: (2025, 4, 10), 59347: (2026, 6, 8),
    40943: (2020, 2, 20), 66587: (2026, 2, 10), 60328: (2025, 4, 14),
    61734: (2025, 3, 22), 57804: (2025, 1, 3), 56758: (2024, 8, 25),
    30799: (2024, 5, 12), 44187: (2021, 6, 15), 35009: (2020, 8, 31),
    41553: (2019, 1, 6), 48259: (2026, 6, 18), 27906: (2018, 11, 19),
    27882: (2019, 4, 7), 39802: (2019, 11, 5), 28572: (2020, 10, 4),
    46429: (2021, 10, 18), 44591: (2022, 3, 3), 41695: (2022, 5, 15),
    52205: (2022, 6, 26), 45349: (2022, 10, 9), 23647: (2024, 1, 28),
}


def cal2jd(y, m, d):
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5


def parse_bstar(line1):
    s = line1[53:61]
    mant = s[:-2].strip()
    exp = s[-2:]
    if mant in ('', '-', '+'):
        return 0.0
    sign = 1.0
    if mant[0] == '-':
        sign = -1.0
        mant = mant[1:]
    elif mant[0] == '+':
        mant = mant[1:]
    return sign * float('0.' + mant) * (10.0 ** int(exp))


def parse_epoch_jd(line1):
    field = line1[18:32]
    yy = int(field[:2])
    yr = 2000 + yy if yy < 57 else 1900 + yy
    doy_frac = float(field[2:])
    jd0 = cal2jd(yr, 1, 1) - 1.0
    return jd0 + doy_frac


def load_tle_series(path):
    series = []
    with open(path) as f:
        for line in f:
            if line.startswith('1 '):
                jd = parse_epoch_jd(line)
                bstar = parse_bstar(line)
                series.append((jd, bstar))
    series.sort()
    return series


def nearest_bstar(series, jd_target):
    best = min(series, key=lambda p: abs(p[0] - jd_target))
    return best[1]


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def main():
    rows = []
    with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
        for r in csv.DictReader(f):
            if r['zone'] == 'ERR' or not r['zepoch'].strip():
                continue
            rows.append({
                'norad': int(r['norad']),
                'zone': int(r['zone']),
                'e_opt': float(r['e_opt']),
                'bn_opt': float(r['bn_opt']),
                'reentry_jd': float(r['reentry_jd']),
                'rpe_pct': float(r['rpe_pct']),
                'zstat': int(r['zstat']),
                'zepoch': float(r['zepoch']),
                'rms_fit': float(r['rms_fit']),
            })

    by_obj = {}
    for r in rows:
        by_obj.setdefault(r['norad'], []).append(r)
    for norad in by_obj:
        by_obj[norad].sort(key=lambda r: r['zone'])

    tle_cache = {}
    for norad, path in NORAD_FILE.items():
        try:
            tle_cache[norad] = load_tle_series(path)
        except FileNotFoundError:
            tle_cache[norad] = []

    # ---------------- A: BSTAR vs fitted BN ----------------
    print('=' * 70)
    print('A. BC/attitude variability: BSTAR (TLE, per zone epoch) vs fitted BN')
    print('=' * 70)
    all_r = []
    all_rd = []
    for norad, zrows in by_obj.items():
        series = tle_cache.get(norad, [])
        if not series or len(zrows) < 3:
            continue
        bstars = [nearest_bstar(series, r['zepoch']) for r in zrows]
        bns = [r['bn_opt'] for r in zrows]
        r = pearson(bstars, bns)
        dbstar = [b2 - b1 for b1, b2 in zip(bstars, bstars[1:])]
        dbn = [b2 - b1 for b1, b2 in zip(bns, bns[1:])]
        rd = pearson(dbstar, dbn)
        if r is not None:
            all_r.append((norad, r, len(zrows)))
            rd_str = f'{rd:+.3f}' if rd is not None else '  n/a'
            print(f'  {norad:6d}  n_zones={len(zrows):2d}  '
                  f'r(level)={r:+.3f}  r(zone-to-zone delta)={rd_str}')
            if rd is not None:
                all_rd.append(rd)
    if all_r:
        rs = [x[1] for x in all_r]
        print(f'\n  objects with >=3 zones: {len(all_r)}')
        print(f'  levels:  mean r = {sum(rs)/len(rs):+.3f}   '
              f'median r = {sorted(rs)[len(rs)//2]:+.3f}')
        pos = sum(1 for r in rs if r > 0.2)
        neg = sum(1 for r in rs if r < -0.2)
        print(f'  |r|>0.2: {pos} positive, {neg} negative, '
              f'{len(rs)-pos-neg} near-zero (of {len(rs)})')
    if all_rd:
        print(f'  deltas:  mean r = {sum(all_rd)/len(all_rd):+.3f}   '
              f'median r = {sorted(all_rd)[len(all_rd)//2]:+.3f}  '
              f'(n={len(all_rd)}, guards against shared-trend confound)')
        posd = sum(1 for r in all_rd if r > 0.2)
        negd = sum(1 for r in all_rd if r < -0.2)
        print(f'  |r|>0.2: {posd} positive, {negd} negative, '
              f'{len(all_rd)-posd-negd} near-zero (of {len(all_rd)})')

    # ---------------- B: RPE-metric validity ----------------
    print()
    print('=' * 70)
    print('B. RPE-metric validity: |RPE| vs prediction horizon length')
    print('=' * 70)
    pts = []
    for norad, zrows in by_obj.items():
        if norad not in OBS_DECAY:
            continue
        obs_jd = cal2jd(*OBS_DECAY[norad])
        for r in zrows:
            if r['reentry_jd'] <= 0:
                continue
            horizon = obs_jd - r['zepoch']
            if horizon <= 0:
                continue
            pts.append((horizon, abs(r['rpe_pct']), r['rms_fit'],
                        r['zstat'], norad, r['zone']))
    pts.sort()
    print(f'  {len(pts)} zone predictions with computable horizon\n')
    print(f'  {"horizon(d)":>10}  {"|RPE|%":>8}  {"rms":>8}  {"zstat":>5}  norad/zone')
    for h, rpe, rms, zstat, norad, zone in pts:
        print(f'  {h:10.1f}  {rpe:8.2f}  {rms:8.3f}  {zstat:5d}  {norad}/{zone}')
    short = [p for p in pts if p[0] < 100]
    long_ = [p for p in pts if p[0] >= 100]
    if short and long_:
        ms = sum(p[1] for p in short) / len(short)
        ml = sum(p[1] for p in long_) / len(long_)
        print(f'\n  mean |RPE|  horizon<100d (n={len(short)}): {ms:.1f}%')
        print(f'  mean |RPE|  horizon>=100d (n={len(long_)}): {ml:.1f}%')
    r_h = pearson([p[0] for p in pts], [p[1] for p in pts])
    r_invh = pearson([1.0 / p[0] for p in pts], [p[1] for p in pts])
    print(f'\n  r(horizon, |RPE|)     = {r_h:+.3f}' if r_h is not None else '')
    print(f'  r(1/horizon, |RPE|)   = {r_invh:+.3f}' if r_invh is not None else '')

    # ---------------- C: identifiability floor ----------------
    print()
    print('=' * 70)
    print('C. Short-window identifiability floor: fit RMS distribution')
    print('=' * 70)
    all_rms = [r['rms_fit'] for r in rows]
    all_rms.sort()
    n = len(all_rms)
    print(f'  n={n} zones total (all objects, all zstat)')
    print(f'  rms_fit: min={all_rms[0]:.4f}  '
          f'p25={all_rms[n//4]:.4f}  median={all_rms[n//2]:.4f}  '
          f'p75={all_rms[3*n//4]:.4f}  max={all_rms[-1]:.4f}')

    trusted = [r['rms_fit'] for r in rows if r['zstat'] == 0]
    boundary = [r['rms_fit'] for r in rows if r['zstat'] == 2]
    for label, arr in (('zstat=0 (ok/trusted)', trusted),
                        ('zstat=2 (boundary-pinned)', boundary)):
        if arr:
            arr2 = sorted(arr)
            m = len(arr2)
            print(f'  {label}: n={m}  median={arr2[m//2]:.4f}  '
                  f'max={arr2[-1]:.4f}')

    rms_all = [r['rms_fit'] for r in rows if r['reentry_jd'] > 0]
    rpe_all = [abs(r['rpe_pct']) for r in rows if r['reentry_jd'] > 0]
    r_rms_rpe = pearson(rms_all, rpe_all)
    print(f'\n  r(rms_fit, |RPE|) across all predicting zones (n={len(rms_all)}) '
          f'= {r_rms_rpe:+.3f}' if r_rms_rpe is not None else '')


if __name__ == '__main__':
    main()
