"""Issue #32 follow-up: calibrate a pooled, cross-object BN~BSTAR empirical
prior (log10(BN) = slope*log10(BSTAR) + intercept) from the 30-object
campaign, for narrowing the GA's BN search range per zone using that zone's
own TLE-published BSTAR -- the "not yet built" candidate flagged by Phase 2
finding 1 (BN anti-correlates with BSTAR, median r=-0.70, expected by
construction since BSTAR~Cd*A/m and BN~m/(Cd*A) are approximate
reciprocals -- log-log linear is the natural functional form for an
approximate power-law/reciprocal relationship).

Reuses phase2_error_budget.py's data loading (same CSV, same per-zone
nearest-BSTAR-by-epoch lookup) so the calibration is grounded in the same
already-validated correlation, not a fresh assumption.
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
                series.append((parse_epoch_jd(line), parse_bstar(line)))
    series.sort()
    return series


def nearest_bstar(series, jd_target):
    return min(series, key=lambda p: abs(p[0] - jd_target))[1]


def linreg(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    slope = sxy / sxx
    intercept = my - slope * mx
    resid = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
    rmse = math.sqrt(sum(r * r for r in resid) / n)
    r2 = 1.0 - sum(r * r for r in resid) / syy
    return slope, intercept, rmse, r2


def main():
    rows = []
    with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
        for r in csv.DictReader(f):
            if r['zone'] == 'ERR' or not r['zepoch'].strip():
                continue
            bn = float(r['bn_opt'])
            if bn <= 0:
                continue
            rows.append({
                'norad': int(r['norad']), 'bn_opt': bn,
                'zstat': int(r['zstat']), 'zepoch': float(r['zepoch']),
            })

    tle_cache = {}
    for norad, path in NORAD_FILE.items():
        try:
            tle_cache[norad] = load_tle_series(path)
        except FileNotFoundError:
            tle_cache[norad] = []

    log_bstar, log_bn = [], []
    log_bstar_trusted, log_bn_trusted = [], []
    for r in rows:
        series = tle_cache.get(r['norad'], [])
        if not series:
            continue
        b = nearest_bstar(series, r['zepoch'])
        if b <= 0:
            continue
        log_bstar.append(math.log10(b))
        log_bn.append(math.log10(r['bn_opt']))
        if r['zstat'] == 0:
            log_bstar_trusted.append(math.log10(b))
            log_bn_trusted.append(math.log10(r['bn_opt']))

    print(f'All zones (n={len(log_bstar)}):')
    slope, intercept, rmse, r2 = linreg(log_bstar, log_bn)
    print(f'  log10(BN) = {slope:.4f}*log10(BSTAR) + {intercept:.4f}')
    print(f'  rmse(log10 BN) = {rmse:.4f}  (factor of {10**rmse:.2f}x)  R2={r2:.3f}')

    print(f'\nTrusted zones only, zstat=0 (n={len(log_bstar_trusted)}):')
    slope_t, intercept_t, rmse_t, r2_t = linreg(log_bstar_trusted, log_bn_trusted)
    print(f'  log10(BN) = {slope_t:.4f}*log10(BSTAR) + {intercept_t:.4f}')
    print(f'  rmse(log10 BN) = {rmse_t:.4f}  (factor of {10**rmse_t:.2f}x)  R2={r2_t:.3f}')


if __name__ == '__main__':
    main()
