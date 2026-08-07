"""Compare OREM's own re-entry-date predictions against DRAMA/OSCAR's,
using the results run_oscar_campaign.py already wrote to
drama/output/oscar_campaign_results.csv.

Prints an aggregate summary, a per-object table of the biggest OREM/DRAMA
disagreements, and appends a short markdown summary to drama/README.md.
"""

import csv
import statistics
from pathlib import Path

RESULTS_CSV = Path(__file__).resolve().parent.parent / "output" / "oscar_campaign_results.csv"
README = Path(__file__).resolve().parent.parent / "README.md"


EXTREME_OUTLIER_DELTA = 500  # |orem_rpe - oscar_rpe| threshold, see analyze_outliers
MOLNIYA_INC = 63.4  # critical inclination, deg -- HEO Molniya/GLONASS-family debris orbits


def load_comparable_rows():
    with open(RESULTS_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    comparable = []
    for r in rows:
        if r["orem_rpe_pct"] and r["oscar_rpe_pct"]:
            r["orem_rpe_pct"] = float(r["orem_rpe_pct"])
            r["oscar_rpe_pct"] = float(r["oscar_rpe_pct"])
            r["norad"] = int(r["norad"])
            if r.get("inc_deg"):
                r["inc_deg"] = float(r["inc_deg"])
                r["perigee_alt_km"] = float(r["perigee_alt_km"])
            comparable.append(r)
    return rows, comparable


def analyze_outliers(comparable):
    """Check whether the most extreme OREM/DRAMA disagreements share an orbital
    regime (inclination, perigee altitude) -- distinguishes "these are all the
    same kind of orbit" from "scattered, unrelated mispredictions."
    """
    extreme = [
        r for r in comparable
        if abs(r["orem_rpe_pct"] - r["oscar_rpe_pct"]) > EXTREME_OUTLIER_DELTA
        and "inc_deg" in r
    ]
    if not extreme:
        return None
    near_molniya = [r for r in extreme if abs(r["inc_deg"] - MOLNIYA_INC) < 3.0]
    return {
        "n_extreme": len(extreme),
        "n_near_molniya": len(near_molniya),
        "extreme_norads": sorted(r["norad"] for r in extreme),
        "molniya_norads": sorted(r["norad"] for r in near_molniya),
        "inc_range": (min(r["inc_deg"] for r in extreme), max(r["inc_deg"] for r in extreme)),
        "perigee_range": (
            min(r["perigee_alt_km"] for r in extreme),
            max(r["perigee_alt_km"] for r in extreme),
        ),
    }


def summarize(rows, comparable):
    n_total = len(rows)
    n_no_orem = sum(1 for r in rows if not r["orem_rpe_pct"])
    n_no_reentry = sum(1 for r in rows if r["oscar_status"] == "no_reentry_in_100yr")
    n_error = sum(1 for r in rows if r["oscar_status"] == "run_error")

    orem_abs = [abs(r["orem_rpe_pct"]) for r in comparable]
    oscar_abs = [abs(r["oscar_rpe_pct"]) for r in comparable]

    summary = {
        "n_total": n_total,
        "n_comparable": len(comparable),
        "n_no_orem_primary": n_no_orem,
        "n_oscar_no_reentry": n_no_reentry,
        "n_oscar_error": n_error,
        "orem_mean_abs_rpe": statistics.mean(orem_abs),
        "orem_median_abs_rpe": statistics.median(orem_abs),
        "oscar_mean_abs_rpe": statistics.mean(oscar_abs),
        "oscar_median_abs_rpe": statistics.median(oscar_abs),
        "n_drama_closer": sum(
            1 for r in comparable if abs(r["oscar_rpe_pct"]) < abs(r["orem_rpe_pct"])
        ),
    }
    return summary


def print_report(summary, comparable):
    print(f"Campaign objects total:        {summary['n_total']}")
    print(f"Comparable (both predictions): {summary['n_comparable']}")
    print(f"  OREM had no primary estimate:      {summary['n_no_orem_primary']}")
    print(f"  OSCAR predicted no re-entry (100yr): {summary['n_oscar_no_reentry']}")
    print(f"  OSCAR run errors:                   {summary['n_oscar_error']}")
    print()
    print(f"OREM  |RPE%|  median={summary['orem_median_abs_rpe']:.2f}  "
          f"mean={summary['orem_mean_abs_rpe']:.2f}")
    print(f"OSCAR |RPE%|  median={summary['oscar_median_abs_rpe']:.2f}  "
          f"mean={summary['oscar_mean_abs_rpe']:.2f}")
    print("(mean is outlier-dominated for OSCAR -- see the top-10 table below; "
          "median is the fairer summary)")
    print(f"DRAMA/OSCAR closer to truth than OREM: "
          f"{summary['n_drama_closer']}/{summary['n_comparable']}")
    print()

    print("Biggest OREM vs. DRAMA disagreements (top 10):")
    print(f"{'norad':>7} {'orem_rpe%':>10} {'oscar_rpe%':>11} {'|delta|':>8}")
    ranked = sorted(
        comparable,
        key=lambda r: abs(r["orem_rpe_pct"] - r["oscar_rpe_pct"]),
        reverse=True,
    )
    for r in ranked[:10]:
        delta = abs(r["orem_rpe_pct"] - r["oscar_rpe_pct"])
        print(f"{r['norad']:>7} {r['orem_rpe_pct']:>10.2f} {r['oscar_rpe_pct']:>11.2f} {delta:>8.2f}")

    outliers = analyze_outliers(comparable)
    if outliers:
        print(f"\nExtreme outliers (|delta| > {EXTREME_OUTLIER_DELTA}%): "
              f"{outliers['n_extreme']} objects {outliers['extreme_norads']}")
        print(f"  inclination range: {outliers['inc_range'][0]:.1f}-{outliers['inc_range'][1]:.1f} deg")
        print(f"  perigee altitude range: {outliers['perigee_range'][0]:.0f}-{outliers['perigee_range'][1]:.0f} km")
        print(f"  within 3 deg of Molniya critical inclination (63.4 deg): "
              f"{outliers['n_near_molniya']}/{outliers['n_extreme']}  {outliers['molniya_norads']}")


def update_readme(summary, outliers):
    if outliers and outliers["n_near_molniya"] >= outliers["n_extreme"] - 1:
        outlier_para = f"""Median is the fairer headline comparison -- OSCAR's mean is dragged out by
{outliers['n_extreme']} extreme outliers (NORAD {outliers['extreme_norads']}, all with
\\|orem_rpe - oscar_rpe\\| > {EXTREME_OUTLIER_DELTA}%: OSCAR predicts a multi-decade lifetime
where the object actually decayed within ~1-1.5 years). Two hypotheses were tested and
ruled out: an eccentricity mismatch (OREM's own fitted eccentricity for the same zone is
close to the raw TLE eccentricity used here) and OSCAR's solar/geomagnetic activity
scenario choice (tested directly across all 4 of OSCAR's scenario options for NORAD
39802 -- lifetime only varies 27-37 years across all of them, nowhere near enough to
explain a 30x gap from the true ~1.35-year lifetime). What the outliers actually share:
**{outliers['n_near_molniya']}/{outliers['n_extreme']}** sit within 3&deg; of the Molniya
critical inclination (63.4&deg;, inclination range {outliers['inc_range'][0]:.1f}-{outliers['inc_range'][1]:.1f}&deg;),
perigee altitude {outliers['perigee_range'][0]:.0f}-{outliers['perigee_range'][1]:.0f} km --
and their NORAD IDs cluster tightly (41679/41686/41695), consistent with fragments of the
same debris-generating breakup event. The most likely explanation is that OREM's own
ballistic-number fit for these small-sample, late-zone objects doesn't transfer as a
physically portable quantity into OSCAR's own density/decay model for this specific
regime -- not root-caused further than that here (would need per-fragment attitude/shape
modeling or a direct atmosphere-model diff to close out)."""
    else:
        outlier_para = ("OSCAR's mean is dragged out by a handful of extreme outliers "
                         "-- see `oscar_campaign_results.csv` for the worst-disagreement "
                         "objects (inc_deg/perigee_alt_km columns included for follow-up).")

    section = f"""
## Campaign results (OREM vs. DRAMA/OSCAR, {summary['n_comparable']} comparable objects)

Full run of `run_oscar_campaign.py` + `compare_orem_drama.py` against all
{summary['n_total']} objects in `scratch_rpe/rpe_campaign.csv`:

| metric | OREM | DRAMA/OSCAR |
|---|---|---|
| median \\|RPE%\\| | {summary['orem_median_abs_rpe']:.2f} | {summary['oscar_median_abs_rpe']:.2f} |
| mean \\|RPE%\\| | {summary['orem_mean_abs_rpe']:.2f} | {summary['oscar_mean_abs_rpe']:.2f} |

{outlier_para}

DRAMA/OSCAR's independent prediction landed closer to the true observed
re-entry date than OREM's own primary-zone estimate for
**{summary['n_drama_closer']}/{summary['n_comparable']}** objects.

{summary['n_no_orem_primary']} objects had no OREM primary-zone estimate to
compare against (DRAMA/OSCAR was still run on these -- see
`oscar_campaign_results.csv`); {summary['n_oscar_no_reentry']} objects showed
no OSCAR re-entry within the 100-year propagation window;
{summary['n_oscar_error']} OSCAR runs errored.

This campaign used the "simpler" orbit-state method (raw TLE semi-major
axis/eccentricity; only OREM's fitted ballistic number feeds in, via mass)
rather than perigee-preserving state matching -- see `drama/README.md`'s
setup notes. That choice trades some physical consistency for simplicity
and is a plausible secondary contributor to the outliers above.

Full per-object results: `output/oscar_campaign_results.csv`.
"""
    text = README.read_text()
    marker = "## Campaign results"
    if marker in text:
        text = text[: text.index(marker)]
    README.write_text(text.rstrip() + "\n" + section)


def main():
    rows, comparable = load_comparable_rows()
    summary = summarize(rows, comparable)
    outliers = analyze_outliers(comparable)
    print_report(summary, comparable)
    update_readme(summary, outliers)
    print(f"\nUpdated {README}")


if __name__ == "__main__":
    main()
