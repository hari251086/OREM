"""npoe_compare.py -- compares orem_trajectory.csv (from npoe_xval.F)
against the NPOE reference propagator's own output for object 42928,
epoch 2017-10-03 (see npoe_xval.F's header for the full context).

Reference data lives outside this repo -- it's the original NPOE run
under E:\\Research\\1. R&D\\Re-entry\\COSPAR ASR\\42928\\Zone - 2\\,
not something this script can regenerate. Kept here as the analysis
that produced the numbers in issue #38 / README.md Version History,
not as a rerunnable CI-style check.
"""
import csv
import datetime
import math

R_EARTH = 6378.1363

REFERENCE_DIR = (
    r"E:\Research\1. R&D\Re-entry\COSPAR ASR\42928\Zone - 2"
)


def parse_npoe(path):
    """Parse NPOE osculating/mean elements output: returns list of
    (day, apogee_alt_km, sma_km, ecc)."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = [p.strip() for p in line.split(',')]
            day = float(parts[2])
            sma = float(parts[3])
            ecc = float(parts[4])
            apo = sma * (1.0 + ecc) - R_EARTH
            rows.append((day, apo, sma, ecc))
    return rows


def parse_orem(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append((float(row['day']), float(row['apogee_alt_km']),
                         float(row['sma_km']), float(row['ecc'])))
    return rows


def interp(rows, day):
    """Linear interpolation of OREM's (sparse, ~1/rev) trajectory at
    an arbitrary day, to compare against NPOE's dense day marks."""
    if day <= rows[0][0]:
        return rows[0][1]
    if day >= rows[-1][0]:
        return rows[-1][1]
    lo, hi = 0, len(rows) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rows[mid][0] <= day:
            lo = mid
        else:
            hi = mid
    d0, a0 = rows[lo][0], rows[lo][1]
    d1, a1 = rows[hi][0], rows[hi][1]
    if d1 == d0:
        return a0
    frac = (day - d0) / (d1 - d0)
    return a0 + frac * (a1 - a0)


def find_crossing(rows, threshold=80.0):
    """First day the apogee altitude crosses below threshold (OREM's
    own re-entry criterion), linearly interpolated."""
    for i in range(1, len(rows)):
        if rows[i - 1][1] >= threshold and rows[i][1] < threshold:
            d0, a0 = rows[i - 1][0], rows[i - 1][1]
            d1, a1 = rows[i][0], rows[i][1]
            frac = (threshold - a0) / (a1 - a0)
            return d0 + frac * (d1 - d0)
    return None


def main():
    npoe_osc = parse_npoe(REFERENCE_DIR + r"\Z2osc.out")
    npoe_mean = parse_npoe(REFERENCE_DIR + r"\Z2.out")
    orem = parse_orem("orem_trajectory.csv")

    print(f"NPOE osc points: {len(npoe_osc)}, span "
          f"{npoe_osc[0][0]:.2f}-{npoe_osc[-1][0]:.2f} days")
    print(f"NPOE mean points: {len(npoe_mean)}, span "
          f"{npoe_mean[0][0]:.2f}-{npoe_mean[-1][0]:.2f} days")
    print(f"OREM points: {len(orem)}, span "
          f"{orem[0][0]:.2f}-{orem[-1][0]:.2f} days")
    print()

    checkpoints = [0, 10, 30, 60, 100, 150, 200, 250, 300, 350, 400,
                   430, 440, 443]
    print(f"{'day':>6s} {'NPOE_osc_apo':>13s} {'NPOE_mean_apo':>14s} "
          f"{'OREM_apo':>10s} {'OREM-NPOEosc':>13s} {'pct_err':>9s}")
    for day in checkpoints:
        if day > orem[-1][0] or day > npoe_osc[-1][0]:
            continue
        no = interp(npoe_osc, day)
        nm = interp(npoe_mean, day)
        oo = interp(orem, day)
        diff = oo - no
        pct = diff / no * 100 if no != 0 else float('nan')
        print(f"{day:6.1f} {no:13.3f} {nm:14.3f} {oo:10.3f} "
              f"{diff:13.3f} {pct:9.3f}")

    print()
    common_max_day = min(orem[-1][0], npoe_osc[-1][0])
    errs = []
    for (day, apo, sma, ecc) in npoe_osc:
        if day > common_max_day:
            break
        oo = interp(orem, day)
        errs.append(oo - apo)

    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    mean_err = sum(errs) / len(errs)
    print(f"Full-trajectory comparison (n={len(errs)} NPOE sample "
          f"points, day 0-{common_max_day:.1f}):")
    print(f"  RMS(OREM - NPOE_osc apogee alt) = {rms:.3f} km")
    print(f"  mean signed error = {mean_err:.3f} km")
    print(f"  max |error| = {max(abs(e) for e in errs):.3f} km")

    orem_cross = find_crossing(orem)
    npoe_osc_cross = find_crossing(npoe_osc)
    npoe_mean_cross = find_crossing(npoe_mean)
    print()
    print("Day of 80km apogee-altitude crossing (re-entry proxy):")
    print(f"  OREM:      {orem_cross}")
    print(f"  NPOE osc:  {npoe_osc_cross}")
    print(f"  NPOE mean: {npoe_mean_cross}")
    print(f"  OREM last day: {orem[-1][0]:.2f}, last apogee alt: "
          f"{orem[-1][1]:.2f} km")
    print(f"  NPOE osc last day: {npoe_osc[-1][0]:.2f}, last apogee "
          f"alt: {npoe_osc[-1][1]:.2f} km")

    epoch0 = datetime.datetime(2017, 10, 3, 14, 35, 56)
    if orem_cross:
        print(f"  OREM re-entry date:     "
              f"{epoch0 + datetime.timedelta(days=orem_cross)}")
    if npoe_osc_cross:
        print(f"  NPOE osc re-entry date: "
              f"{epoch0 + datetime.timedelta(days=npoe_osc_cross)}")
    print("  Real observed re-entry: 2019-03-03")


if __name__ == "__main__":
    main()
