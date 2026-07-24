"""Phase 3 of the RPE investigation (issue #32): literature-grounded test of
Wang & Gurfil's solar-apsidal-resonance criterion against object 33587
(issue #27).

Background (see project_orem_status / reference_orem_reentry_literature
memory): 33587's last 131 days of real TLE data (2022-07-25 to 2022-12-03)
show hp collapsing 616->341 km with simultaneous circularization
(e: 0.37->0.27) -- textbook drag signature. But two leading hypotheses were
already ruled out: (1) GMAT full-force (near-exact point-mass Sun/Moon)
predicts flat hp(t), matching OREM's own truncated model -- third-body
truncation is NOT the gap; (2) a direct BN refit against the whole 131-day
window (bn_opt=126, not boundary-pinned) gives RMS~4.68 (~468 km residual,
10-100x worse than good fits) -- no physically plausible constant BN
explains the magnitude via King-Hele drag alone starting from a 616 km
perigee, where drag should be near-negligible.

Wang & Gurfil (2017, "The Role of Solar Apsidal Resonance...", Adv. Space
Res. 59, Eq. 6-12) give an exact, implementable resonance-detection
criterion: the GTO's combined RAAN+perigee secular drift rate (from J2)
becomes commensurate with the Sun's apparent mean motion when

    a(1-e^2)^(4/7) = [3 J2 sqrt(mu) RE^2 (5cos^2i - 2cosi - 1)
                      / (4 dM_solar/dt)]^(2/7)                 (their Eq. 12)

At that crossing, the previously-periodic (~180-day) perigee-height
oscillation turns into a MONOTONIC one-way increase or decrease -- exactly
matching 33587's own one-way hp(t) collapse, and mechanistically distinct
from both hypotheses already ruled out (this is neither third-body
truncation nor a pure-drag/BN misfit -- it's a resonance-driven secular
eccentricity change that then activates drag as perigee drops).

This script computes a(t), e(t), i(t) directly from every raw TLE mean
element in 33587's real 2009-2022 record (no propagation, no SGP4 --
just Kepler's third law on the TLE's own mean motion, the same basis
the earlier "textbook drag signature" diagnostic used) and checks
whether/when the record crosses its own resonance threshold.
"""
import math

MU = 398600.4415       # km^3/s^2 (input/const_new.dat)
RE = 6378.1363          # km (input/const_new.dat)
J2 = 1.08263e-3         # standard Earth J2 (literature value)
DM_SOLAR_DT = 2.0 * math.pi / (365.25 * 86400.0)   # rad/s, Sun's apparent mean motion


def cal2jd(y, m, d):
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5


def parse_epoch_jd(line1):
    field = line1[18:32]
    yy = int(field[:2])
    yr = 2000 + yy if yy < 57 else 1900 + yy
    doy_frac = float(field[2:])
    jd0 = cal2jd(yr, 1, 1) - 1.0
    return jd0 + doy_frac


def resonance_lambda_crit(i_rad):
    # (5cos^2(i) - 2cos(i) - 1) changes sign at i ~ 46.4 deg and ~106.8 deg
    # (roots of 5x^2-2x-1=0). Eq. 12's RHS is raised to the EVEN power 2/7
    # (= squaring an odd, real, 1/7-root branch), so it's real and positive
    # regardless of that sign -- take abs() before the fractional power to
    # get Python's real ** to agree with the real odd-root branch, instead
    # of the principal complex branch it defaults to for a negative base.
    ci = math.cos(i_rad)
    num = 3.0 * J2 * math.sqrt(MU) * RE * RE * (5.0 * ci * ci - 2.0 * ci - 1.0)
    return abs(num / (4.0 * DM_SOLAR_DT)) ** (2.0 / 7.0)


def load_record(path):
    rows = []
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]
    for k in range(0, len(lines) - 1, 2):
        l1, l2 = lines[k], lines[k + 1]
        if not (l1.startswith('1 ') and l2.startswith('2 ')):
            continue
        jd = parse_epoch_jd(l1)
        inc_deg = float(l2[8:16])
        ecc = float('0.' + l2[26:33].strip())
        n_rev_day = float(l2[52:63])
        n_rad_s = n_rev_day * 2.0 * math.pi / 86400.0
        a = (MU / (n_rad_s * n_rad_s)) ** (1.0 / 3.0)
        rows.append((jd, a, ecc, inc_deg))
    rows.sort()
    return rows


def jd_to_ymd(jd):
    jd = jd + 0.5
    z = int(jd)
    f = jd - z
    if z < 2299161:
        aa = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        aa = z + 1 + alpha - alpha // 4
    b = aa + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    return year, month, day


def main():
    rows = load_record('input/example_33587.tle.txt')
    print(f'33587: {len(rows)} real TLEs, '
          f'{jd_to_ymd(rows[0][0])[:2]} -> {jd_to_ymd(rows[-1][0])[:2]}')

    print('\n  date        a(km)     e       i(deg)   lambda    lambda_crit   '
          'ratio  regime')
    prev_regime = None
    crossings = []
    for jd, a, e, i_deg in rows:
        i_rad = math.radians(i_deg)
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        lam_crit = resonance_lambda_crit(i_rad)
        ratio = lam / lam_crit
        regime = 'PRE (above crit)' if lam > lam_crit else 'POST (below crit)'
        if prev_regime is not None and regime != prev_regime:
            crossings.append((jd, a, e, i_deg, lam, lam_crit))
        prev_regime = regime

    # Print a decimated view: every ~10th point plus anything near the end
    n = len(rows)
    step = max(1, n // 40)
    for idx in range(0, n, step):
        jd, a, e, i_deg = rows[idx]
        i_rad = math.radians(i_deg)
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        lam_crit = resonance_lambda_crit(i_rad)
        y, m, d = jd_to_ymd(jd)
        print(f'  {y}-{m:02d}-{int(d):02d}  {a:9.2f}  {e:.4f}  {i_deg:7.3f}  '
              f'{lam:8.2f}  {lam_crit:10.2f}  {lam/lam_crit:.4f}  '
              f'{"PRE" if lam > lam_crit else "POST"}')
    # Always show the tail densely (last 20 points -- where the collapse is)
    print('\n  --- last 20 TLEs (dense) ---')
    for jd, a, e, i_deg in rows[-20:]:
        i_rad = math.radians(i_deg)
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        lam_crit = resonance_lambda_crit(i_rad)
        y, m, d = jd_to_ymd(jd)
        print(f'  {y}-{m:02d}-{int(d):02d}  {a:9.2f}  {e:.4f}  {i_deg:7.3f}  '
              f'{lam:8.2f}  {lam_crit:10.2f}  {lam/lam_crit:.4f}  '
              f'{"PRE" if lam > lam_crit else "POST"}')

    print(f'\n  {len(crossings)} regime crossing(s) found in the full record:')
    for jd, a, e, i_deg, lam, lam_crit in crossings:
        y, m, d = jd_to_ymd(jd)
        print(f'    crossing at {y}-{m:02d}-{int(d):02d}  '
              f'a={a:.1f} e={e:.4f} i={i_deg:.3f}  '
              f'lambda={lam:.2f} vs crit={lam_crit:.2f}')


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


def campaign_sweep():
    print('\n' + '=' * 70)
    print('Generalization check: resonance crossings across all 30 campaign objects')
    print('=' * 70)
    for norad, path in NORAD_FILE.items():
        try:
            rows = load_record(path)
        except FileNotFoundError:
            continue
        if len(rows) < 3:
            continue
        prev_regime = None
        crossings = []
        for jd, a, e, i_deg in rows:
            lam = a * (1.0 - e * e) ** (4.0 / 7.0)
            lam_crit = resonance_lambda_crit(math.radians(i_deg))
            regime = 'PRE' if lam > lam_crit else 'POST'
            if prev_regime is not None and regime != prev_regime:
                crossings.append(jd)
            prev_regime = regime
        decay_jd = cal2jd(*OBS_DECAY[norad]) if norad in OBS_DECAY else None
        last_jd = rows[-1][0]
        i_mean = sum(r[3] for r in rows) / len(rows)
        if crossings:
            y, m, d = jd_to_ymd(crossings[-1])
            gap_to_decay = f'{decay_jd - crossings[-1]:.0f}d to decay' if decay_jd else 'n/a'
            gap_to_last = f'{last_jd - crossings[-1]:.0f}d to last TLE'
            print(f'  {norad:6d}  i~{i_mean:5.1f} deg  {len(crossings)} crossing(s), '
                  f'last at {y}-{m:02d}-{int(d):02d}  ({gap_to_last}, {gap_to_decay})')
        else:
            print(f'  {norad:6d}  i~{i_mean:5.1f} deg  no crossing in record')


if __name__ == '__main__':
    main()
    campaign_sweep()
