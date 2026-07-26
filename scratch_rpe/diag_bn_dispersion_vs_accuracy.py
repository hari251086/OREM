#!/usr/bin/env python3
"""
Diagnostic for issue #32 follow-up: does Sharma et al. (2009) KSGEN's GA
objective -- minimize the DISPERSION of independently-predicted re-entry
times across TLEs -- have a shot at improving OREM, before building it?

Sharma's GA doesn't fit a single BC to observed apogee data within one
window (OREM's current ga_fitness, in ga.F). It instead searches for the
ONE ballistic coefficient value that makes MANY independent single-TLE
propagations -- each run forward to its own predicted re-entry time -- all
agree with each other. Cross-prediction self-consistency, not curve-fit
goodness.

OREM never explicitly optimizes for this, but it already produces exactly
the data needed to check whether the premise holds: each zone independently
fits its own bn_opt via the current RMS-based GA, and the campaign harness
already reports both bn_opt per zone and the final ensemble accuracy
(ens_rpe_pct) against the real observed re-entry date. If Sharma's premise
holds for OREM's object population, objects whose independently-fitted
bn_opt already agrees closely across zones (low dispersion) should show
better final accuracy than objects where it swings wildly -- i.e. cross-zone
BN agreement should predict accuracy.

This is a correlational check using data already committed in
scratch_rpe/rpe_campaign.csv -- no new propagation runs, no Fortran, no
production code touched. It is NOT a direct test of "would explicitly
optimizing for cross-TLE self-consistency improve accuracy" (that would
require actually running the alternative GA objective) -- it tests whether
the WEAKER premise the objective relies on is even true of OREM's existing
per-zone fits.
"""
import csv
import math
from collections import defaultdict

CSV_PATH = "scratch_rpe/rpe_campaign.csv"


def to_f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def main():
    rows = defaultdict(list)
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            try:
                norad = int(row["norad"])
            except (ValueError, KeyError):
                continue
            rows[norad].append(row)

    stats = []
    for norad, zrows in rows.items():
        valid = [z for z in zrows if to_f(z["reentry_jd"]) not in (None, 0.0)]
        if len(valid) < 2:
            continue
        bns = [to_f(z["bn_opt"]) for z in valid]
        ens_rpe = to_f(valid[0]["ens_rpe_pct"])
        if ens_rpe is None:
            continue
        mean_bn = sum(bns) / len(bns)
        std_bn = math.sqrt(sum((b - mean_bn) ** 2 for b in bns) / len(bns))
        cv_bn = std_bn / mean_bn if mean_bn else None
        stats.append((norad, len(valid), mean_bn, std_bn, cv_bn, abs(ens_rpe)))

    print(f"{'norad':>7} {'nvalid':>6} {'mean_bn':>9} {'std_bn':>8} "
          f"{'cv_bn':>7} {'abs_ens_rpe':>12}")
    for s in sorted(stats, key=lambda x: x[4] if x[4] is not None else 0):
        norad, n, mb, sb, cv, ar = s
        print(f"{norad:>7} {n:>6} {mb:>9.2f} {sb:>8.2f} {cv:>7.3f} {ar:>12.2f}")

    xs = [s[4] for s in stats if s[4] is not None]
    ys = [s[5] for s in stats if s[4] is not None]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
    r = cov / (sx * sy) if sx * sy else float("nan")
    print(f"\nn objects with >=2 valid zones: {n}")
    print(f"Pearson r(cv_bn, |ens_rpe_pct|) = {r:.3f}  (R^2 = {r*r:.3f})")


if __name__ == "__main__":
    main()
