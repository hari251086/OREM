"""Wang & Gurfil (2017) U-turn resonance sweep applied to the 180 new
GTO/HEO candidates found beyond OREM's original 97-object campaign
(screen_new_candidates.py + fetch_new_candidates.py). Reuses the exact
same detection logic as resonance_campaign_sweep.py (its
`resonance_lambda_crit`/`load_record`/`cal2jd`/`jd_to_ymd` functions),
just against the new candidate list + their real satcat DECAY dates
instead of rpe_campaign.F's DATA arrays.
"""
import csv
import json
import math

from resonance_campaign_sweep import (
    cal2jd, jd_to_ymd, load_record, resonance_lambda_crit,
)

BASE = 'E:/GitHub/OREM'


def parse_decay_ymd(decay_str):
    # satcat DECAY column: 'YYYY-MM-DD'
    y, m, d = decay_str.split('-')
    return int(y), int(m), int(d)


def analyze(norad, name, decay_str):
    path = f'{BASE}/input/example_{norad}.tle.txt'
    try:
        rows = load_record(path)
    except FileNotFoundError:
        return {'norad': norad, 'name': name, 'status': 'NO_TLE_FILE'}
    if len(rows) < 3:
        return {'norad': norad, 'name': name, 'status': 'TOO_FEW_TLES',
                'n_tle': len(rows)}

    prev_regime = None
    crossings = []
    for jd, a, e, i_deg in rows:
        lam = a * (1.0 - e * e) ** (4.0 / 7.0)
        lam_crit = resonance_lambda_crit(math.radians(i_deg))
        regime = 'PRE' if lam > lam_crit else 'POST'
        if prev_regime is not None and regime != prev_regime:
            crossings.append(jd)
        prev_regime = regime

    decay_jd = cal2jd(*parse_decay_ymd(decay_str))
    last_jd = rows[-1][0]
    first_jd = rows[0][0]
    i_mean = sum(r[3] for r in rows) / len(rows)

    result = {
        'norad': norad, 'name': name, 'n_tle': len(rows),
        'i_mean_deg': i_mean, 'first_tle_jd': first_jd,
        'last_tle_jd': last_jd, 'decay_jd': decay_jd,
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
    with open(f'{BASE}/scratch_rpe/new_candidates.json') as f:
        data = json.load(f)
    candidates = data['have_cache'] + data['need_fetch']
    print(f'{len(candidates)} new candidates to sweep\n')

    results = [analyze(c['norad'], c['name'], c['decay']) for c in candidates]

    resonant = [r for r in results if r['status'] == 'RESONANT']
    not_resonant = [r for r in results if r['status'] == 'NOT_RESONANT']
    other = [r for r in results if r['status'] not in ('RESONANT', 'NOT_RESONANT')]

    print(f'RESONANT: {len(resonant)}/{len(results)}')
    print(f'NOT_RESONANT: {len(not_resonant)}')
    print(f'Other (no/too-few TLE): {len(other)}\n')

    print(f'{"NORAD":>7} {"name":<26} {"i(deg)":>7} {"#TLE":>5} {"#X":>3} '
          f'{"resonance date":>14} {"->decay(d)":>10} {"pre-frac":>8}')
    for r in sorted(resonant, key=lambda x: x['days_resonance_to_decay']):
        y, m, d = r['resonance_ymd']
        frac = r['frac_record_pre_resonance']
        frac_s = f'{frac:8.3f}' if frac is not None else '     n/a'
        print(f'{r["norad"]:7d} {r["name"][:26]:<26} {r["i_mean_deg"]:7.2f} '
              f'{r["n_tle"]:5d} {r["n_crossings"]:3d} '
              f'{y}-{m:02d}-{d:02d}   {r["days_resonance_to_decay"]:10.0f} '
              f'{frac_s}')

    with open(f'{BASE}/scratch_rpe/resonance_sweep_newcand.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['norad', 'name', 'status', 'i_mean_deg', 'n_tle',
                    'n_crossings', 'resonance_jd', 'resonance_date',
                    'decay_jd', 'decay_date', 'days_resonance_to_decay',
                    'days_resonance_to_last_tle', 'frac_record_pre_resonance'])
        for r in results:
            if r['status'] == 'RESONANT':
                y, m, d = r['resonance_ymd']
                dy, dm, dd = jd_to_ymd(r['decay_jd'])
                w.writerow([r['norad'], r['name'], r['status'],
                            f'{r["i_mean_deg"]:.3f}', r['n_tle'], r['n_crossings'],
                            f'{r["resonance_jd"]:.4f}', f'{y}-{m:02d}-{d:02d}',
                            f'{r["decay_jd"]:.4f}', f'{dy}-{dm:02d}-{dd:02d}',
                            f'{r["days_resonance_to_decay"]:.1f}',
                            f'{r["days_resonance_to_last_tle"]:.1f}',
                            f'{r["frac_record_pre_resonance"]:.4f}'
                            if r['frac_record_pre_resonance'] is not None else ''])
            else:
                w.writerow([r['norad'], r['name'], r['status'],
                            f'{r.get("i_mean_deg", ""):.3f}' if 'i_mean_deg' in r else '',
                            r.get('n_tle', ''), r.get('n_crossings', ''),
                            '', '', '', '', '', '', ''])
    print(f'\nWritten: scratch_rpe/resonance_sweep_newcand.csv')


if __name__ == '__main__':
    main()
