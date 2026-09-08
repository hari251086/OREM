"""Re-screen the freshly refreshed SATCAT snapshot (2026-09-08, vs. the
2026-08-30 one screen_new_candidates.py used) against ALL 277 objects
already in the pool (97 original + 180 already screened/fetched), to
answer "any more objects?" -- new decays or SATCAT corrections in the
9 days since the last screen could have made new candidates eligible.
Same established filter as every prior round.
"""
import csv
import json
import re
from pathlib import Path

DATA_ROOT = Path("E:/Research/1. R&D/Re-entry/2026/Data")
SATCAT = DATA_ROOT / "catalog/latest/satcat.csv"
CAMPAIGN_F = Path("E:/GitHub/OREM/scratch_rpe/rpe_campaign.F")
NEW_CANDIDATES_JSON = Path("E:/GitHub/OREM/scratch_rpe/new_candidates.json")


def existing_norads():
    text = CAMPAIGN_F.read_text()
    m = re.search(r'data\s+norad\s*/\s*(.*?)/', text, re.IGNORECASE | re.DOTALL)
    body = re.sub(r'\n\s*&', ' ', m.group(1))
    original97 = set(int(t.strip()) for t in body.split(','))

    with open(NEW_CANDIDATES_JSON) as f:
        data = json.load(f)
    expanded180 = set(c['norad'] for c in data['have_cache'] + data['need_fetch'])

    return original97 | expanded180


def main():
    excluded = existing_norads()
    print(f'{len(excluded)} objects already in the pool (excluded)')

    candidates = []
    with open(SATCAT, encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f):
            try:
                norad = int(row['NORAD_CAT_ID'])
                apogee = float(row['APOGEE']) if row['APOGEE'] else None
                perigee = float(row['PERIGEE']) if row['PERIGEE'] else None
            except (ValueError, KeyError):
                continue
            if norad in excluded:
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

    print(f'{len(candidates)} NEW candidates found (same filter, fresh SATCAT, '
          f'excluding all 277 already screened)\n')
    for c in candidates:
        print(f"  {c['norad']:>7} {c['name'][:30]:<30} {c['type']:<12} "
              f"apo={c['apogee']:.0f} per={c['perigee']:.0f} inc={c['inc']} "
              f"decay={c['decay']}")

    with open('scratch_rpe/more_candidates_20260908.json', 'w') as f:
        json.dump(candidates, f, indent=2)
    print(f"\nWritten: scratch_rpe/more_candidates_20260908.json")


if __name__ == '__main__':
    main()
