"""Validation RPE campaign -- writes per-zone TLE subset files (2-line
TLE format, matching the COSPAR ASR reference's own per-zone TLE-split
convention: 35497/35497-0.txt, -2.txt, -3.txt, -4.txt) plus the zone
table CSVs, and generates the Fortran driver's object/zone DATA list so
it isn't hand-transcribed.

Reuses validation_rpe_zones.build_zones (the already-decided zone
selection, forward-only cursor, 3-29 TLEs/zone, first zone >=6mo before
resonance, spread to last TLE) -- no new zone-selection logic here.
"""
import csv
import json
import os
from validation_rpe_zones import OBJECTS, build_zones, BASE

OUT_ROOT = f'{BASE}/scratch_rpe/validation_RPE'
INPUT_DIR = f'{BASE}/input'


def read_tle_lines(path):
    with open(path) as f:
        return [l.rstrip('\n') for l in f]


def main():
    driver_entries = []  # (norad, zone_idx, label, tle_filename)
    for norad, cfg in OBJECTS.items():
        rows, zones, notes = build_zones(norad, cfg['tle'], cfg['resonance_jd'])
        full_lines = read_tle_lines(f'{BASE}/{cfg["tle"]}')
        # full_lines is 2 lines per TLE, same order as load_record's
        # sorted rows -- rows are already sorted by load_record, and
        # the raw file is written in epoch order too (verified: this
        # campaign's TLE files come straight from Space-Track history
        # fetches, already chronological), so index i in `rows` maps to
        # lines[2*i], lines[2*i+1] directly.
        obj_dir = f'{OUT_ROOT}/{norad}'
        os.makedirs(obj_dir, exist_ok=True)

        zone_table_path = f'{obj_dir}/zones.csv'
        with open(zone_table_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['zone', 'label', 'start_epoch', 'end_epoch',
                        'apogee_km_min', 'apogee_km_max', 'n_tle',
                        'days_before_last_tle'])
            for zi, (label, s, e) in enumerate(zones, start=1):
                sub_rows = rows[s:e + 1]
                apo = [a * (1.0 + ecc) - 6378.1363 for (_, a, ecc, _) in sub_rows]
                from resonance_campaign_sweep import jd_to_ymd
                y1, m1, d1 = jd_to_ymd(sub_rows[0][0])
                y2, m2, d2 = jd_to_ymd(sub_rows[-1][0])
                w.writerow([
                    zi, label, f'{y1}-{m1:02d}-{d1:02d}', f'{y2}-{m2:02d}-{d2:02d}',
                    f'{min(apo):.1f}', f'{max(apo):.1f}', len(sub_rows),
                    f'{rows[-1][0] - sub_rows[-1][0]:.1f}',
                ])

                tle_name = f'example_{norad}_valzone{zi}.tle.txt'
                tle_path = f'{INPUT_DIR}/{tle_name}'
                with open(tle_path, 'w') as tf:
                    for i in range(s, e + 1):
                        tf.write(full_lines[2 * i] + '\n')
                        tf.write(full_lines[2 * i + 1] + '\n')
                driver_entries.append((norad, zi, label, tle_name))

        with open(f'{obj_dir}/notes.txt', 'w') as f:
            f.write('\n'.join(notes) + ('\n' if notes else ''))

        print(f'{norad}: {len(zones)} zones written -> {obj_dir}, notes={len(notes)}')

    # Fortran driver DATA block generator -- fixed-form 72-col limit
    # (fpm.toml source-form=fixed truncates past col 72 silently, so
    # values must be wrapped across continuation lines, not one long
    # line like the hand-written rpe_campaign_resonance16.F precedent
    # got away with for only 16 entries).
    def wrap_data(name, values, per_line=6):
        lines = [f'      data {name} /']
        for i in range(0, len(values), per_line):
            chunk = ', '.join(values[i:i + per_line])
            sep = ',' if i + per_line < len(values) else ''
            lines.append(f'     &   {chunk}{sep}')
        lines[-1] += ' /'
        return '\n'.join(lines) + '\n'

    # DATA block (declaration section) and string-assignment block
    # (executable statements -- must follow all declarations, so they
    # cannot share one Fortran INCLUDE placed in the declaration area).
    # Spliced directly into the driver template below rather than left
    # as separate .inc files, since a moved/renamed source file would
    # otherwise break the include's relative path silently.
    data_block = (
        wrap_data('norad', [str(n) for n, _, _, _ in driver_entries])
        + wrap_data('oyr', [f'{OBJECTS[n]["decay_ymd"][0]}.d0' for n, _, _, _ in driver_entries])
        + wrap_data('omo', [f'{OBJECTS[n]["decay_ymd"][1]}.d0' for n, _, _, _ in driver_entries])
        + wrap_data('ody', [f'{OBJECTS[n]["decay_ymd"][2]}.d0' for n, _, _, _ in driver_entries])
    )
    assign_lines = []
    for i, (norad, zi, label, tle_name) in enumerate(driver_entries, start=1):
        assign_lines.append(f"      tle_file({i}) = 'input/{tle_name}'")
    for i, (norad, zi, label, tle_name) in enumerate(driver_entries, start=1):
        short = f'{norad} z{zi} {label}'[:24]
        assign_lines.append(f"      obj_name({i}) = '{short}'")
    assign_block = '\n'.join(assign_lines) + '\n'

    out_dir = f'{BASE}/scratch_rpe/validation_rpe'
    os.makedirs(out_dir, exist_ok=True)
    from driver_template import DRIVER_TEMPLATE
    with open(f'{out_dir}/validation_rpe_campaign.F', 'w') as f:
        f.write(DRIVER_TEMPLATE.format(
            nobj=len(driver_entries), data_block=data_block, assign_block=assign_block,
        ))

    with open(f'{BASE}/scratch_rpe/validation_rpe_zone_map.json', 'w') as f:
        json.dump(
            [{'order': i, 'norad': n, 'zone': zi, 'label': lb, 'tle': tn}
             for i, (n, zi, lb, tn) in enumerate(driver_entries)],
            f, indent=2,
        )

    print(f'\n{len(driver_entries)} total zones -> {out_dir}/validation_rpe_campaign.F '
          f'(+ validation_rpe_zone_map.json for result joining)')


if __name__ == '__main__':
    main()
