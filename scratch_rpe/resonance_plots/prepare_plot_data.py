"""Data prep for the Wang & Gurfil (2017)-style resonance plot PDF
(user request, follow-on to issue #56). Scilab does the plotting; this
script does all the orbital-mechanics work Scilab shouldn't have to
duplicate -- TLE parsing (via the already-validated
resonance_campaign_sweep.load_record), Sun mean-longitude, azimuth
angle, perigee height, and the resonance-condition/decay-condition
curve parameters -- and writes one plain CSV per object plus a manifest
CSV Scilab reads to drive the batch.

Azimuth definition: the paper's own solar-azimuth-vs-apsidal-line angle
isn't literally reproducible from TLEs alone (it needs a full 3D Sun
position projected into the orbit plane), but the paper's own resonance
condition (Eq. 9) is exactly d(RAAN+AOP)/dt ~= dM_solar/dt -- so the
same U-turn signature (the angle increasing, then reversing, then
decreasing) is captured by the resonance angle
    psi(t) = mod(RAAN(t) + AOP(t) - L_sun(t), 360)
where L_sun is the Sun's mean ecliptic longitude. This is a standard
substitute for "local time of perigee"/solar-azimuth-vs-apsidal-line in
this exact literature (the paper's own introduction notes the two are
directly linked). RAAN/AOP come straight from each TLE's own mean
elements (columns 18-25 / 35-42 of line 2); L_sun uses the same
365.25-day solar year already baked into resonance_lambda_crit's
DM_SOLAR_DT, for self-consistency with the resonance-crossing detection
already done in resonance_campaign_sweep.py.
"""
import csv
import math
import sys

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe')
from resonance_campaign_sweep import (  # noqa: E402
    MU, RE, cal2jd, resonance_lambda_crit,
)

BASE = 'E:/GitHub/OREM'
OUTDIR = f'{BASE}/scratch_rpe/resonance_plots/data'

# (norad, tle_file, name, resonance_jd, decay_jd) -- resonance_jd/decay_jd
# pulled from the two sweep CSVs already produced for issue #56.
OBJECTS = [
    # --- original 16 (scratch_rpe/resonance_sweep_97obj.csv) ---
    (35497, 'input/example_35497.tle.txt', 'Ariane 5 ESC-A'),
    (37151, 'input/example_37151.tle.txt', 'Long March 3B'),
    (27526, 'input/example_27526.tle.txt', 'Ariane 5 R/B'),
    (37819, 'input/example_37819.tle.txt', 'Proton-M R/B'),
    (59347, 'input/example_59347.tle.txt', 'GTO R/B'),
    (40943, 'input/example_40943.tle.txt', 'GTO debris'),
    (60328, 'input/example_60328.tle.txt', 'CZ-3B R/B'),
    (18983, 'input/example_18983.tle.txt', 'SL-6 R/B(2)'),
    (13112, 'input/example_13112.tle.txt', 'SL-6 R/B(2)'),
    (28141, 'input/example_28141.tle.txt', 'CZ-2C R/B'),
    (19881, 'input/example_19881.tle.txt', 'Cosmos931 deb'),
    (5980, 'input/example_5980.tle.txt', 'Titan3C TS deb'),
    (4570, 'input/example_4570.tle.txt', 'SL-6 R/B(2)'),
    (10724, 'input/example_10724.tle.txt', 'AtlasCentaur deb'),
    (28140, 'input/example_28140.tle.txt', 'Doublestar TC-1'),
    (44911, 'input/example_44911.tle.txt', 'CZ-5 R/B'),
    # --- new 7 (scratch_rpe/resonance_sweep_newcand.csv) ---
    (7096, 'input/example_7096.tle.txt', 'SKYNET 2A'),
    (68538, 'input/example_68538.tle.txt', 'ARTEMIS 2 (INTEGRITY)'),
    (6294, 'input/example_6294.tle.txt', 'MOLNIYA 1-22'),
    (18954, 'input/example_18954.tle.txt', 'ARIANE 3 DEB (SYLDA)'),
    (29462, 'input/example_29462.tle.txt', 'SL-12 DEB'),
    (5941, 'input/example_5941.tle.txt', 'PROGNOZ 1'),
    (6231, 'input/example_6231.tle.txt', 'MOLNIYA 1-21'),
]


def parse_epoch_jd(line1):
    field = line1[18:32]
    yy = int(field[:2])
    yr = 2000 + yy if yy < 57 else 1900 + yy
    doy_frac = float(field[2:])
    jd0 = cal2jd(yr, 1, 1) - 1.0
    return jd0 + doy_frac


def load_full_record(path):
    """Like resonance_campaign_sweep.load_record, but also keeps RAAN
    and argument of perigee (needed for the azimuth angle, which the
    campaign sweep never needed)."""
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
    # good enough for axis labeling (not calendar-exact at the day level)
    return 2000.0 + (jd - 2451545.0) / 365.25


def L_sun_deg(jd):
    d = jd - 2451545.0
    return (280.46 + (360.0 / 365.25) * d) % 360.0


def load_resonance_epoch(norad):
    import csv as _csv
    for fname in ('resonance_sweep_97obj.csv', 'resonance_sweep_newcand.csv'):
        path = f'{BASE}/scratch_rpe/{fname}'
        with open(path) as f:
            for row in _csv.DictReader(f):
                if int(row['norad']) == norad and row['status'] == 'RESONANT':
                    return float(row['resonance_jd']), row['decay_jd']
    raise KeyError(f'{norad} not found as RESONANT in either sweep CSV')


def main():
    manifest = []
    for norad, tle_file, name in OBJECTS:
        rows = load_full_record(f'{BASE}/{tle_file}')
        resonance_jd, decay_jd_s = load_resonance_epoch(norad)
        i_mean = sum(r[3] for r in rows) / len(rows)
        lam_crit = resonance_lambda_crit(math.radians(i_mean))

        out_rows = []
        for jd, a, e, i_deg, raan, argp in rows:
            year = jd_to_decimal_year(jd)
            hp = a * (1.0 - e) - RE
            lam = a * (1.0 - e * e) ** (4.0 / 7.0)
            lsun = L_sun_deg(jd)
            azimuth = (raan + argp - lsun) % 360.0
            out_rows.append((jd, year, a, e, i_deg, raan, argp, hp,
                              azimuth, lam))

        csv_path = f'{OUTDIR}/obj_{norad}.csv'
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['jd', 'year', 'a_km', 'e', 'i_deg', 'raan_deg',
                        'argp_deg', 'hp_km', 'azimuth_deg', 'lambda'])
            w.writerows(out_rows)

        resonance_year = jd_to_decimal_year(resonance_jd)
        decay_year = jd_to_decimal_year(float(decay_jd_s)) if decay_jd_s else None

        manifest.append({
            'norad': norad, 'name': name, 'n_tle': len(rows),
            'i_mean_deg': i_mean, 'lambda_crit': lam_crit,
            'resonance_jd': resonance_jd, 'resonance_year': resonance_year,
            'decay_year': decay_year,
            'first_year': jd_to_decimal_year(rows[0][0]),
            'last_year': jd_to_decimal_year(rows[-1][0]),
            'csv_path': f'data/obj_{norad}.csv',
        })
        print(f'{norad:>7} {name[:28]:<28} {len(rows):5d} pts  '
              f'i={i_mean:6.2f}  lam_crit={lam_crit:8.2f}  '
              f'resonance_year={resonance_year:8.2f}')

    manifest_path = f'{BASE}/scratch_rpe/resonance_plots/manifest.csv'
    with open(manifest_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['norad', 'name', 'n_tle', 'i_mean_deg', 'lambda_crit',
                    'resonance_jd', 'resonance_year', 'decay_year',
                    'first_year', 'last_year', 'csv_path'])
        for m in manifest:
            w.writerow([m['norad'], m['name'], m['n_tle'],
                        f"{m['i_mean_deg']:.4f}", f"{m['lambda_crit']:.4f}",
                        f"{m['resonance_jd']:.4f}", f"{m['resonance_year']:.4f}",
                        f"{m['decay_year']:.4f}" if m['decay_year'] else '',
                        f"{m['first_year']:.4f}", f"{m['last_year']:.4f}",
                        m['csv_path']])
    print(f'\n{len(manifest)} objects written. Manifest: {manifest_path}')


if __name__ == '__main__':
    main()
