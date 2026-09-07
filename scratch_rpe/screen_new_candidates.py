"""Screen Space-Track's SATCAT (local cached snapshot, reused per the
Space-Track usage policy's "store and reuse" rule -- OREM-Watchlist's
`data_root/catalog/latest/satcat.csv`, 2026-08-30, ~1 week old, well
within the established weekly-cadence practice) for GTO/HEO candidates
beyond the 97 already in OREM's `scratch_rpe/rpe_campaign.F`, using the
exact same established filter this campaign has used for every prior
expansion (see rpe_campaign.F's own header comments): APOGEE>8000km,
PERIGEE<3000km, APOGEE<=280,000km (excludes TLI/lunar-regime stages),
OBJECT_TYPE in (ROCKET BODY, DEBRIS, PAYLOAD), DECAY not null (real
observed decay needed for RPE ground truth).

Also cross-references OREM-Watchlist's already-fetched per-object TLE
cache (`data_root/objects/<norad>/<norad>.tle.txt`, 729 objects as of
this run) so candidates that already have local history need no new
Space-Track query at all -- per the usage policy's "store and reuse
local data" rule (SPACE_TRACK_USAGE_POLICY.md �4).
"""
import csv
import re
from pathlib import Path

DATA_ROOT = Path("E:/Research/1. R&D/Re-entry/2026/Data")
SATCAT = DATA_ROOT / "catalog/latest/satcat.csv"
OBJECTS_DIR = DATA_ROOT / "objects"
CAMPAIGN_F = Path("E:/GitHub/OREM/scratch_rpe/rpe_campaign.F")


def existing_norads():
    text = CAMPAIGN_F.read_text()
    m = re.search(r'data\s+norad\s*/\s*(.*?)/', text, re.IGNORECASE | re.DOTALL)
    body = re.sub(r'\n\s*&', ' ', m.group(1))
    return set(int(t.strip()) for t in body.split(','))


def main():
    existing = existing_norads()
    print(f'{len(existing)} objects already in the campaign (excluded)')

    candidates = []
    with open(SATCAT, encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f):
            try:
                norad = int(row['NORAD_CAT_ID'])
                apogee = float(row['APOGEE']) if row['APOGEE'] else None
                perigee = float(row['PERIGEE']) if row['PERIGEE'] else None
            except (ValueError, KeyError):
                continue
            if norad in existing:
                continue
            if apogee is None or perigee is None:
                continue
            if not (apogee > 8000 and perigee < 3000 and apogee <= 280000):
                continue
            if row['OBJECT_TYPE'] not in ('ROCKET BODY', 'DEBRIS', 'PAYLOAD'):
                continue
            if not row['DECAY']:
                continue
            candidates.append({
                'norad': norad, 'name': row['SATNAME'],
                'type': row['OBJECT_TYPE'], 'apogee': apogee,
                'perigee': perigee, 'inc': row['INCLINATION'],
                'decay': row['DECAY'],
            })

    print(f'{len(candidates)} candidates pass the established filter '
          f'(APOGEE>8000/PERIGEE<3000/APOGEE<=280000/ROCKET BODY|DEBRIS|'
          f'PAYLOAD/DECAY not null), excluding the 97 already used')

    have_cache = []
    need_fetch = []
    for c in candidates:
        local = OBJECTS_DIR / str(c['norad']) / f"{c['norad']}.tle.txt"
        if local.exists() and local.stat().st_size > 200:
            have_cache.append(c)
        else:
            need_fetch.append(c)

    print(f'{len(have_cache)} already have cached TLE history locally '
          f'(OREM-Watchlist objects/ dir, zero new Space-Track requests)')
    print(f'{len(need_fetch)} need a fresh Space-Track gp_history fetch')

    import json
    out = {
        'have_cache': have_cache,
        'need_fetch': need_fetch,
    }
    with open('scratch_rpe/new_candidates.json', 'w') as f:
        json.dump(out, f, indent=2)
    print('\nWritten: scratch_rpe/new_candidates.json')

    print('\n--- sample of need_fetch (first 20) ---')
    for c in need_fetch[:20]:
        print(f"  {c['norad']:>7} {c['name'][:30]:<30} {c['type']:<12} "
              f"apo={c['apogee']:.0f} per={c['perigee']:.0f} inc={c['inc']} "
              f"decay={c['decay']}")


if __name__ == '__main__':
    main()
