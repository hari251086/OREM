#!/usr/bin/env python3
"""
Diagnostic for issue #32 follow-up: does the Sharma et al. (2009) KSGEN
technique of splitting the ballistic coefficient by radial-distance regime
(BC above vs below r=6500 km, i.e. ~122 km altitude) even apply to OREM's
own object population, before building the full RSM/GA/orem.F architecture
needed to search two independently-fitted BNs per zone?

Sharma's split point was fit against near-circular LEO decay objects
(SROSS-C2, SL-12 R/B tracked via CIRA 1972 lineage) where perigee radius
falls below the threshold only in the terminal weeks of decay, uniformly
across the whole orbit (near-circular => perigee ~ apogee ~ r). OREM's
curated-7 objects are GTO-insertion-lineage debris with eccentricity often
0.3-0.6 -- for those, perigee radius can already sit near/below the 6500 km
threshold for YEARS while apogee remains at GTO altitude (~36000 km), which
is a structurally different regime than what Sharma's technique targets.

This script parses each curated-7 object's full TLE history, computes
perigee/apogee radius per TLE (from mean motion + eccentricity, standard
two-line element field defs, no external library), and reports:
  - the eccentricity and perigee-altitude range across the tracked history
  - the fraction of TLEs where perigee is already below the 122 km split
  - the day of first crossing relative to total tracked span (does it only
    happen in the terminal few percent of history, or much earlier?)
  - persistence after first crossing: of all TLEs after that point, what
    fraction stay below the split? Distinguishes genuine sustained terminal
    decay (persist ~100%, Sharma's assumed regime) from a transient dip or
    perigee-precession oscillation around the boundary (persist low, object
    recovers above the split for a long stretch afterward)

No changes to any production code. Read-only analysis of input/*.tle.txt.
"""
import math
import glob
import os

MU = 398600.4418  # km^3/s^2, matches OREM's amue (Earth GM)
R_EARTH = 6378.137  # km
R_SPLIT = 6500.0  # km, Sharma et al. 2009's regime boundary

CURATED7 = ["42928", "35497", "37151", "39615", "27526", "32007", "37819"]


def parse_tle_file(path):
    """Return list of (epoch_days_since_first, ecc, sma_km, perigee_km, apogee_km)."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]

    recs = []
    epoch0 = None
    i = 0
    while i + 1 < len(lines):
        l1, l2 = lines[i], lines[i + 1]
        if not (l1.startswith("1 ") and l2.startswith("2 ")):
            i += 1
            continue
        try:
            epoch_yr = int(l1[18:20])
            epoch_day = float(l1[20:32])
            ecc = float("0." + l2[26:33].strip())
            mean_motion = float(l2[52:63])  # rev/day
        except ValueError:
            i += 2
            continue

        year = 2000 + epoch_yr if epoch_yr < 57 else 1900 + epoch_yr
        # days since an arbitrary common reference (Jan 1 of epoch year) is
        # enough for a *relative* time axis within one object's file
        epoch_abs_days = (year - 2000) * 365.25 + epoch_day

        n_rad_s = mean_motion * 2.0 * math.pi / 86400.0
        sma = (MU / (n_rad_s ** 2)) ** (1.0 / 3.0)
        perigee = sma * (1.0 - ecc) - R_EARTH
        apogee = sma * (1.0 + ecc) - R_EARTH

        recs.append((epoch_abs_days, ecc, sma, perigee, apogee))
        i += 2

    recs.sort(key=lambda r: r[0])
    if recs:
        t0 = recs[0][0]
        recs = [(t - t0, e, a, p, ap) for (t, e, a, p, ap) in recs]
    return recs


def main():
    alt_split = R_SPLIT - R_EARTH
    print(f"(split point r={R_SPLIT:.0f} km = altitude {alt_split:.1f} km)\n")
    print(f"{'norad':>7} {'n_tle':>6} {'e_min':>7} {'e_max':>7} "
          f"{'peri_min_alt':>12} {'peri_max_alt':>12} "
          f"{'frac_below_split':>17} {'first_below_day':>16} {'span_days':>10} "
          f"{'persist_after':>13} {'sustained?':>10}")

    for norad in CURATED7:
        candidates = glob.glob(f"input/example_{norad}.tle.txt")
        if not candidates:
            candidates = glob.glob(f"input/example_{norad}*.tle.txt")
        if not candidates:
            print(f"{norad:>7}  (no TLE file found)")
            continue
        # prefer the un-suffixed full-history file if present
        path = min(candidates, key=len)

        recs = parse_tle_file(path)
        if not recs:
            print(f"{norad:>7}  (no parsable TLEs in {path})")
            continue

        peri_alt = [r[3] for r in recs]  # perigee altitude, km
        ecc = [r[1] for r in recs]
        span = recs[-1][0]

        n_below = sum(1 for pa in peri_alt if pa < alt_split)
        frac_below = n_below / len(recs)

        first_below_day = None
        for (t, e, a, p, ap) in recs:
            if p < alt_split:
                first_below_day = t
                break

        # Persistence: of all TLEs *after* the first crossing, what fraction
        # stay below the split? Distinguishes a genuine sustained terminal
        # decay (Sharma's assumed regime, persist~100%) from a transient
        # dip/perigee-precession oscillation around the boundary (persist
        # low, object recovers above the split for a long stretch after).
        persist_pct = None
        if first_below_day is not None:
            after = [p for (t, e, a, p, ap) in recs if t >= first_below_day]
            n_after_below = sum(1 for p in after if p < alt_split)
            persist_pct = n_after_below / len(after)

        sustained = "yes" if (persist_pct is not None and persist_pct > 0.5) \
            else ("no" if persist_pct is not None else "n/a")

        print(f"{norad:>7} {len(recs):>6} {min(ecc):>7.3f} {max(ecc):>7.3f} "
              f"{min(peri_alt):>12.1f} {max(peri_alt):>12.1f} "
              f"{frac_below:>16.1%} "
              f"{('%.0f' % first_below_day) if first_below_day is not None else 'never':>16} "
              f"{span:>10.0f} "
              f"{(f'{persist_pct:.1%}') if persist_pct is not None else 'n/a':>13} "
              f"{sustained:>10}")


if __name__ == "__main__":
    main()
