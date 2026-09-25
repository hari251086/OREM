"""Validation RPE campaign -- zone selection for 7 confirmed/plausible
resonant objects (issue #56 manual review). Reuses
resonance_campaign_sweep.load_record for TLE parsing (no re-derivation
of epoch/a/e/i math).

Zone criteria (per task spec, not OREM's own decay-linearity
zone_select.F -- that algorithm only accepts quasi-linear apogee-decay
windows and would reject most pre-resonance data by design, since
pre-resonance the orbit is still in its ~180-day oscillation, not
monotonic decay; the point of this campaign is to explicitly cover
that regime too):
  1. First zone must start >=6 months (182.5d) before the object's own
     U-turn resonance epoch (from resonance_sweep_97obj.csv /
     resonance_sweep_newcand.csv -- already computed via the existing
     Wang & Gurfil Eq.12 detector). Zones continue through last TLE.
  2. 3-29 TLEs per zone.
  3. Spread across the full observed range with gaps, matching the
     COSPAR ASR reference style (E:\\Research\\1. R&D\\Re-entry\\COSPAR ASR
     35497/Zone-0,2,3,4: 4 non-contiguous zones, ~6-11 day real-time
     width, TLE counts 8/12/14/33 -- increasing density approaching
     decay). Target 4 anchors: early baseline, pre-resonance approach,
     post-resonance, terminal (near last TLE).
"""
import sys
sys.path.insert(0, '.')
from resonance_campaign_sweep import load_record, jd_to_ymd, RE

BASE = 'E:/GitHub/OREM'
MIN_PTS = 3
MAX_PTS = 29

OBJECTS = {
    # decay_ymd from resonance_sweep_97obj.csv / resonance_sweep_newcand.csv
    # (decay_date column) -- the real, catalogued re-entry date, used as
    # t_obs_cal so orem_run computes a real RPE against ground truth.
    35497: dict(tle='input/example_35497.tle.txt', resonance_jd=2457517.5421, decay_ymd=(2016, 10, 31)),
    37151: dict(tle='input/example_37151.tle.txt', resonance_jd=2456911.2870, decay_ymd=(2015, 12, 3)),
    27526: dict(tle='input/example_27526.tle.txt', resonance_jd=2455566.9414, decay_ymd=(2012, 5, 9)),
    59347: dict(tle='input/example_59347.tle.txt', resonance_jd=2461169.8442, decay_ymd=(2026, 6, 8)),
    40943: dict(tle='input/example_40943.tle.txt', resonance_jd=2458756.6180, decay_ymd=(2020, 2, 20)),
    10724: dict(tle='input/example_10724.tle.txt', resonance_jd=2444483.9536, decay_ymd=(1982, 8, 31)),
    18954: dict(tle='input/example_18954.tle.txt', resonance_jd=2447541.1583, decay_ymd=(1989, 2, 1)),
}


def apogee_km(a, e):
    return a * (1.0 + e) - RE


def build_zones(norad, tle_path, resonance_jd, min_pts=MIN_PTS, max_pts=MAX_PTS):
    """Full contiguous tiling of the ENTIRE TLE record, unlimited zone
    count -- per explicit user correction: the earlier 4-fixed-anchor
    design (baseline/approach/post-resonance/terminal, with large
    intentional gaps between them) was wrong. The only real constraint
    is the per-zone TLE-count range (3-29); zone count itself is
    whatever it takes to tile the full record with no gaps.

    Balanced partition (ceil(n/max_pts) zones, sizes within 1 of each
    other) rather than greedy-fill-then-remainder, so every zone lands
    inside [min_pts, max_pts] by construction -- no special-casing a
    short tail zone.

    Criterion 1 (first zone >=6mo before the object's own U-turn) is
    satisfied automatically: the first zone always starts at the
    object's very first tracked TLE, and every one of these 7 objects'
    records already starts well over 6 months before its own resonance
    epoch (checked below, reported as a note if ever violated). Each
    zone is labelled pre-/post-resonance purely for reporting, based on
    whether its epoch range falls before or after resonance_jd -- this
    has no effect on how the tiling itself is built."""
    rows = load_record(f'{BASE}/{tle_path}')
    n = len(rows)
    first_jd, last_jd = rows[0][0], rows[-1][0]
    gap_days = resonance_jd - first_jd

    notes = []
    if n < min_pts:
        notes.append(f'only {n} TLEs total (<{min_pts}) -- no zone possible')
        return rows, [], notes

    n_zones = max(1, -(-n // max_pts))  # ceil(n / max_pts)
    base = n // n_zones
    rem = n % n_zones

    zones = []
    start = 0
    for i in range(n_zones):
        size = base + (1 if i < rem else 0)
        end = start + size - 1
        mid_jd = rows[(start + end) // 2][0]
        label = 'pre-resonance' if mid_jd < resonance_jd else 'post-resonance'
        zones.append((f'{label} (tile {i + 1}/{n_zones})', start, end))
        start = end + 1

    if gap_days < 182.5:
        notes.append(
            f'first zone starts only {gap_days:.0f}d before resonance '
            f'(<6mo requirement) -- record itself starts late'
        )

    return rows, zones, notes


def main():
    print(f'{"norad":>6} {"zone":<26} {"start":>10} {"end":>10} '
          f'{"n_tle":>5} {"apogee_km_range":>22} {"days_to_last":>12}')
    all_out = {}
    for norad, cfg in OBJECTS.items():
        rows, zones, notes = build_zones(norad, cfg['tle'], cfg['resonance_jd'])
        last_jd = rows[-1][0]
        out_zones = []
        for label, s, e in zones:
            sub = rows[s:e + 1]
            apo = [apogee_km(a, ecc) for (_, a, ecc, _) in sub]
            y1, m1, d1 = jd_to_ymd(sub[0][0])
            y2, m2, d2 = jd_to_ymd(sub[-1][0])
            print(f'{norad:>6} {label:<26} {y1}-{m1:02d}-{d1:02d} '
                  f'{y2}-{m2:02d}-{d2:02d} {len(sub):>5} '
                  f'{min(apo):>9.0f}-{max(apo):<9.0f} '
                  f'{last_jd - sub[-1][0]:>12.1f}')
            out_zones.append(dict(
                label=label, start_jd=sub[0][0], end_jd=sub[-1][0],
                n_tle=len(sub), apo_min=min(apo), apo_max=max(apo),
                days_to_last=last_jd - sub[-1][0],
            ))
        if notes:
            for note in notes:
                print(f'{norad:>6}   NOTE: {note}')
        all_out[norad] = dict(rows=rows, zones=out_zones, notes=notes,
                               tle_path=cfg['tle'])
    return all_out


if __name__ == '__main__':
    main()
