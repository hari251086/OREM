"""Data prep for the FULL GTO/HEO candidate pool (all 277 objects
screened for issue #56: the original 97-object RPE campaign +
the 180 found by expanding the pool), not just the 23 my own Eq.12
crossing-detector flagged RESONANT. User wants to eyeball every object
and pick the real U-turn cases manually, rather than trust the
automated detector alone (some flagged crossings turned out ambiguous
on inspection, and the detector may also be missing real ones with a
messier/multi-crossing signature).

Recomputes resonance status directly per object (single source of
truth, rather than joining the two earlier sweep CSVs) so this script
also re-includes the 4 objects previously excluded by policy (3 Falcon
9 R/Bs, 1 inconsistent-decay-date object) -- flagged with a note in the
info panel instead of being silently dropped, since the point now is
manual review, not an automated-pipeline decision.

Time series are decimated to <=800 points (uniform stride, always
keeping the first/last point) before writing -- several objects in the
expanded pool have 5,000-11,000+ real TLEs, and full-resolution vector
plots at that density would make Scilab slow and the final PDF huge.
800 points is far more than needed to see a ~180-day short-period
oscillation or a resonance U-turn by eye.
"""
import csv
import json
import math
import re
import sys

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
from resonance_campaign_sweep import (  # noqa: E402
    MU, RE, cal2jd, resonance_lambda_crit,
)

BASE = 'E:/GitHub/OREM'
OUTDIR = f'{BASE}/scratch_rpe/resonance_plots/data_fullpool'
MAX_PTS = 800


def parse_campaign_fortran(path):
    with open(path) as f:
        text = f.read()

    def grab_data_block(name):
        m = re.search(r'data\s+' + name + r'\s*/\s*(.*?)/', text,
                       re.IGNORECASE | re.DOTALL)
        body = re.sub(r'\n\s*&', ' ', m.group(1))
        return [re.sub(r'[dD]0$', '', t.strip()) for t in body.split(',')]

    norad = [int(t) for t in grab_data_block('norad')]
    oyr = [int(float(t)) for t in grab_data_block('oyr')]
    omo = [int(float(t)) for t in grab_data_block('omo')]
    ody = [int(float(t)) for t in grab_data_block('ody')]

    obj_name = {}
    for m in re.finditer(r"obj_name\((\d+)\)\s*=\s*'([^']+)'", text):
        obj_name[int(m.group(1))] = m.group(2)

    objs = []
    for idx in range(1, len(norad) + 1):
        objs.append({
            'norad': norad[idx - 1],
            'name': obj_name.get(idx, str(norad[idx - 1])),
            'decay_ymd': (oyr[idx - 1], omo[idx - 1], ody[idx - 1]),
        })
    return objs


def parse_epoch_jd(line1):
    field = line1[18:32]
    yy = int(field[:2])
    yr = 2000 + yy if yy < 57 else 1900 + yy
    doy_frac = float(field[2:])
    jd0 = cal2jd(yr, 1, 1) - 1.0
    return jd0 + doy_frac


def load_full_record(path):
    rows = []
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]
    for k in range(0, len(lines) - 1, 2):
        l1, l2 = lines[k], lines[k + 1]
        if not (l1.startswith('1 ') and l2.startswith('2 ')):
            continue
        jd = parse_epoch_jd(l1)
        inc_deg = float(l2[8:16])
        raan_deg = float(l2[17:25])
        ecc = float('0.' + l2[26:33].strip())
        argp_deg = float(l2[34:42])
        n_rev_day = float(l2[52:63])
        n_rad_s = n_rev_day * 2.0 * math.pi / 86400.0
        a = (MU / (n_rad_s * n_rad_s)) ** (1.0 / 3.0)
        rows.append((jd, a, ecc, inc_deg, raan_deg, argp_deg))
    rows.sort()
    return rows


def jd_to_decimal_year(jd):
    return 2000.0 + (jd - 2451545.0) / 365.25


def L_sun_deg(jd):
    d = jd - 2451545.0
    return (280.46 + (360.0 / 365.25) * d) % 360.0


def decimate(rows, max_pts):
    n = len(rows)
    if n <= max_pts:
        return rows
    stride = n / float(max_pts)
    idx = sorted(set(int(round(i * stride)) for i in range(max_pts)))
    idx[-1] = n - 1
    idx[0] = 0
    return [rows[i] for i in idx]


def build_object_list():
    objs = []
    campaign_objs = parse_campaign_fortran(f'{BASE}/scratch_rpe/rpe_campaign.F')
    for o in campaign_objs:
        objs.append({
            'norad': o['norad'],
            'name': o['name'],
            'tle_file': f"input/example_{o['norad']}.tle.txt",
            'decay_ymd': o['decay_ymd'],
            'pool': 'original97',
        })

    with open(f'{BASE}/scratch_rpe/new_candidates.json') as f:
        newcand = json.load(f)
    for c in newcand['have_cache'] + newcand['need_fetch']:
        y, m, d = (int(x) for x in c['decay'].split('-'))
        objs.append({
            'norad': c['norad'],
            'name': c['name'],
            'tle_file': f"input/example_{c['norad']}.tle.txt",
            'decay_ymd': (y, m, d),
            'pool': 'expanded180',
        })
    return objs


def main():
    import os
    os.makedirs(OUTDIR, exist_ok=True)

    objs = build_object_list()
    print(f'{len(objs)} total objects in the GTO/HEO candidate pool\n')

    manifest = []
    n_skipped = 0
    for o in objs:
        path = f"{BASE}/{o['tle_file']}"
        try:
            rows = load_full_record(path)
        except FileNotFoundError:
            n_skipped += 1
            continue
        if len(rows) < 3:
            n_skipped += 1
            continue

        i_mean = sum(r[3] for r in rows) / len(rows)
        lam_crit = resonance_lambda_crit(math.radians(i_mean))

        prev_regime = None
        crossings = []
        for jd, a, e, i_deg, raan, argp in rows:
            lam = a * (1.0 - e * e) ** (4.0 / 7.0)
            regime = 'PRE' if lam > lam_crit else 'POST'
            if prev_regime is not None and regime != prev_regime:
                crossings.append(jd)
            prev_regime = regime
        is_resonant = len(crossings) > 0
        resonance_jd = crossings[-1] if crossings else None

        decay_jd = cal2jd(*o['decay_ymd'])
        note = ''
        if 'FALCON 9' in o['name'].upper():
            note = 'Falcon 9 -- active deorbit/passivation burns possible (drag-only model may not apply)'
        if decay_jd < rows[-1][0]:
            note = (note + '; ' if note else '') + \
                'inconsistent SATCAT decay date (< last tracked TLE epoch)'

        plot_rows = decimate(rows, MAX_PTS)
        csv_path = f'{OUTDIR}/obj_{o["norad"]}.csv'
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['jd', 'year', 'a_km', 'e', 'i_deg', 'raan_deg',
                        'argp_deg', 'hp_km', 'azimuth_deg', 'lambda'])
            for jd, a, e, i_deg, raan, argp in plot_rows:
                year = jd_to_decimal_year(jd)
                hp = a * (1.0 - e) - RE
                lam = a * (1.0 - e * e) ** (4.0 / 7.0)
                lsun = L_sun_deg(jd)
                azimuth = (raan + argp - lsun) % 360.0
                w.writerow([jd, year, a, e, i_deg, raan, argp, hp, azimuth, lam])

        manifest.append({
            'norad': o['norad'], 'name': o['name'], 'pool': o['pool'],
            'n_tle': len(rows), 'n_plotted': len(plot_rows),
            'i_mean_deg': i_mean, 'lambda_crit': lam_crit,
            'is_resonant': is_resonant,
            'resonance_year': jd_to_decimal_year(resonance_jd) if resonance_jd else None,
            'n_crossings': len(crossings),
            'decay_year': jd_to_decimal_year(decay_jd),
            'first_year': jd_to_decimal_year(rows[0][0]),
            'last_year': jd_to_decimal_year(rows[-1][0]),
            'note': note,
            'csv_path': f'data_fullpool/obj_{o["norad"]}.csv',
        })
        print(f"{o['norad']:>7} {o['name'][:28]:<28} {len(rows):5d} pts "
              f"({len(plot_rows):4d} plotted)  i={i_mean:6.2f}  "
              f"{'RESONANT' if is_resonant else 'not resonant':<13} "
              f"{note}")

    manifest_path = f'{BASE}/scratch_rpe/resonance_plots/manifest_fullpool.csv'
    with open(manifest_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['norad', 'name', 'pool', 'n_tle', 'n_plotted',
                    'i_mean_deg', 'lambda_crit', 'is_resonant',
                    'resonance_year', 'n_crossings', 'decay_year',
                    'first_year', 'last_year', 'note', 'csv_path'])
        for m in manifest:
            w.writerow([m['norad'], m['name'], m['pool'], m['n_tle'],
                        m['n_plotted'], f"{m['i_mean_deg']:.4f}",
                        f"{m['lambda_crit']:.4f}", int(m['is_resonant']),
                        f"{m['resonance_year']:.4f}" if m['resonance_year'] else '',
                        m['n_crossings'], f"{m['decay_year']:.4f}",
                        f"{m['first_year']:.4f}", f"{m['last_year']:.4f}",
                        m['note'], m['csv_path']])

    n_resonant = sum(1 for m in manifest if m['is_resonant'])
    print(f"\n{len(manifest)} objects written ({n_resonant} resonant, "
          f"{len(manifest) - n_resonant} not), {n_skipped} skipped "
          f"(no/too-few TLE data). Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
