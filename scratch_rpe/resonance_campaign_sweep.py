"""Wang & Gurfil (2017) solar-apsidal-resonance ("U-turn") sweep across the
full OREM campaign (97 objects, `scratch_rpe/rpe_campaign.F`).

Generalizes `phase3_resonance_33587.py`'s 30-object sweep to the full
current campaign, and additionally reports the U-turn resonance EPOCH
(JD) per object -- the point needed to implement the hard rule "only run
re-entry prediction from the resonance period onward."

Resonance criterion (Wang & Gurfil Eq. 12): the GTO's J2-driven
(RAAN+perigee) secular drift rate becomes commensurate with the Sun's
apparent mean motion when

    a(1-e^2)^(4/7) = lambda_crit(i)
                   = [3 J2 sqrt(mu) RE^2 (5cos^2i - 2cosi - 1)
                      / (4 dM_solar/dt)]^(2/7)

Before the crossing, eccentricity/perigee-height undergo the ~180-day
short-period oscillation (the "PRE" regime, lambda > lambda_crit for the
typical GTO inclination regime where 5cos^2i-2cosi-1 > 0, i.e. i below
~46.4 deg -- OREM's whole campaign is GTO/HEO-heavy so this branch
applies almost everywhere). The crossing IS the U-turn: the apsidal line
locks onto the solar azimuth, the periodic oscillation turns into a
monotonic increase/decrease of eccentricity, and per the paper this is
"the prelude to the final re-entry of GTO."

Only real TLE data is used (a, e, i straight from each TLE's own mean
elements via Kepler's third law on the mean motion -- no propagation),
same basis as the validated 30-object sweep this extends.
"""
import math
import re

MU = 398600.4415
RE = 6378.1363
J2 = 1.08263e-3
DM_SOLAR_DT = 2.0 * math.pi / (365.25 * 86400.0)


def cal2jd(y, m, d):
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5


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
    return year, month, int(day)


def parse_epoch_jd(line1):
    field = line1[18:32]
    yy = int(field[:2])
    yr = 2000 + yy if yy < 57 else 1900 + yy
    doy_frac = float(field[2:])
    jd0 = cal2jd(yr, 1, 1) - 1.0
    return jd0 + doy_frac


def resonance_lambda_crit(i_rad):
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


def parse_campaign_fortran(path):
    """Extract norad / tle_file / obj_name / oyr,omo,ody straight out of
    rpe_campaign.F's own DATA/assignment statements -- regex, no
    re-derivation, same pattern as orem_ground_truth.py."""
    with open(path) as f:
        text = f.read()

    def grab_data_block(name):
        m = re.search(r'data\s+' + name + r'\s*/\s*(.*?)/', text,
                       re.IGNORECASE | re.DOTALL)
        body = m.group(1)
        body = re.sub(r'\n\s*&', ' ', body)
        toks = [re.sub(r'[dD]0$', '', t.strip()) for t in body.split(',')]
        return toks

    norad = [int(t) for t in grab_data_block('norad')]
    oyr = [int(float(t)) for t in grab_data_block('oyr')]
    omo = [int(float(t)) for t in grab_data_block('omo')]
    ody = [int(float(t)) for t in grab_data_block('ody')]

    tle_file = {}
    for m in re.finditer(r"tle_file\((\d+)\)\s*=\s*'([^']+)'", text):
        tle_file[int(m.group(1))] = m.group(2)
    obj_name = {}
    for m in re.finditer(r"obj_name\((\d+)\)\s*=\s*'([^']+)'", text):
        obj_name[int(m.group(1))] = m.group(2)

    n = len(norad)
    assert len(oyr) == n and len(omo) == n and len(ody) == n, \
        f'array length mismatch: norad={n} oyr={len(oyr)} omo={len(omo)} ody={len(ody)}'
    objs = []
    for idx in range(1, n + 1):
        objs.append({
            'norad': norad[idx - 1],
            'name': obj_name.get(idx, ''),
            'tle_file': tle_file.get(idx, ''),
            'decay_ymd': (oyr[idx - 1], omo[idx - 1], ody[idx - 1]),
        })
    return objs


def analyze_object(obj, base_dir):
    path = base_dir + '/' + obj['tle_file']
    try:
        rows = load_record(path)
    except FileNotFoundError:
        return {**obj, 'status': 'NO_TLE_FILE'}
    if len(rows) < 3:
        return {**obj, 'status': 'TOO_FEW_TLES', 'n_tle': len(rows)}

    prev_regime = None
    crossings = []
    for jd, a, e, i_deg in rows:
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        lam_crit = resonance_lambda_crit(math.radians(i_deg))
        regime = 'PRE' if lam > lam_crit else 'POST'
        if prev_regime is not None and regime != prev_regime:
            crossings.append(jd)
        prev_regime = regime

    decay_jd = cal2jd(*obj['decay_ymd'])
    last_jd = rows[-1][0]
    first_jd = rows[0][0]
    i_mean = sum(r[3] for r in rows) / len(rows)

    result = {
        **obj,
        'n_tle': len(rows),
        'i_mean_deg': i_mean,
        'first_tle_jd': first_jd,
        'last_tle_jd': last_jd,
        'decay_jd': decay_jd,
        'n_crossings': len(crossings),
    }
    if crossings:
        result['status'] = 'RESONANT'
        result['resonance_jd'] = crossings[-1]
        result['resonance_ymd'] = jd_to_ymd(crossings[-1])
        result['days_resonance_to_decay'] = decay_jd - crossings[-1]
        result['days_resonance_to_last_tle'] = last_jd - crossings[-1]
        result['frac_record_pre_resonance'] = (
            (crossings[-1] - first_jd) / (last_jd - first_jd)
            if last_jd > first_jd else None)
    else:
        result['status'] = 'NOT_RESONANT'
    return result


def main():
    base = 'E:/GitHub/OREM'
    objs = parse_campaign_fortran(base + '/scratch_rpe/rpe_campaign.F')
    print(f'Parsed {len(objs)} objects from rpe_campaign.F\n')

    results = [analyze_object(o, base) for o in objs]

    resonant = [r for r in results if r['status'] == 'RESONANT']
    not_resonant = [r for r in results if r['status'] == 'NOT_RESONANT']
    other = [r for r in results if r['status'] not in ('RESONANT', 'NOT_RESONANT')]

    print(f'RESONANT (real U-turn crossing found): {len(resonant)}/{len(results)}')
    print(f'NOT_RESONANT (no crossing in tracked record): {len(not_resonant)}')
    print(f'Other (no/short TLE data): {len(other)}\n')

    print(f'{"NORAD":>7} {"name":<28} {"i(deg)":>7} {"#TLE":>5} {"#X":>3} '
          f'{"resonance date":>14} {"->decay(d)":>10} {"pre-frac":>8}')
    for r in sorted(resonant, key=lambda x: x['days_resonance_to_decay']):
        y, m, d = r['resonance_ymd']
        frac = r['frac_record_pre_resonance']
        print(f'{r["norad"]:7d} {r["name"][:28]:<28} {r["i_mean_deg"]:7.2f} '
              f'{r["n_tle"]:5d} {r["n_crossings"]:3d} '
              f'{y}-{m:02d}-{d:02d}   {r["days_resonance_to_decay"]:10.0f} '
              f'{frac:8.3f}' if frac is not None else '    n/a')

    print('\n--- NOT_RESONANT objects ---')
    for r in not_resonant:
        print(f'{r["norad"]:7d} {r["name"][:28]:<28} i~{r["i_mean_deg"]:.1f}  '
              f'{r["n_tle"]} TLEs')

    print('\n--- other (missing/short TLE) ---')
    for r in other:
        print(f'{r["norad"]:7d} {r["name"][:28]:<28} status={r["status"]}')

    import csv
    with open(base + '/scratch_rpe/resonance_sweep_97obj.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['norad', 'name', 'tle_file', 'status', 'i_mean_deg', 'n_tle',
                    'n_crossings', 'resonance_jd', 'resonance_date',
                    'decay_jd', 'decay_date', 'days_resonance_to_decay',
                    'days_resonance_to_last_tle', 'frac_record_pre_resonance'])
        for r in results:
            if r['status'] == 'RESONANT':
                y, m, d = r['resonance_ymd']
                dy, dm, dd = jd_to_ymd(r['decay_jd'])
                w.writerow([r['norad'], r['name'], r['tle_file'], r['status'],
                            f'{r["i_mean_deg"]:.3f}', r['n_tle'], r['n_crossings'],
                            f'{r["resonance_jd"]:.4f}', f'{y}-{m:02d}-{d:02d}',
                            f'{r["decay_jd"]:.4f}', f'{dy}-{dm:02d}-{dd:02d}',
                            f'{r["days_resonance_to_decay"]:.1f}',
                            f'{r["days_resonance_to_last_tle"]:.1f}',
                            f'{r["frac_record_pre_resonance"]:.4f}'
                            if r['frac_record_pre_resonance'] is not None else ''])
            else:
                w.writerow([r['norad'], r['name'], r['tle_file'], r['status'],
                            f'{r.get("i_mean_deg", ""):.3f}' if 'i_mean_deg' in r else '',
                            r.get('n_tle', ''), r.get('n_crossings', ''),
                            '', '', '', '', '', '', ''])
    print(f'\nWritten: scratch_rpe/resonance_sweep_97obj.csv')


if __name__ == '__main__':
    main()
