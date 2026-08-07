"""Batch OSCAR run: DRAMA's independent re-entry-date prediction for every
object in OREM's own scratch_rpe/rpe_campaign.csv, for later comparison
against OREM's own predictions (see compare_orem_drama.py).

For each object, uses OREM's own PRIMARY-zone selection (the latest zone
with a valid prediction -- mirrors src/report.F's own logic) to pick which
zone's fitted ballistic number (bn_opt) and epoch (zepoch) to build an OSCAR
run from. The orbital state itself comes straight from the TLE nearest that
epoch (semiMajorAxis/eccentricity/inclination/RAAN/argp/meanAnomaly) --
only bn_opt feeds in from OREM's own fit, via spacecraft mass, exactly as
validated for NORAD 21670 in run_oscar_reentry.py.

Objects where OREM itself has no valid primary-zone prediction anywhere
(e.g. every zone errored or hit a GA boundary with reentry_jd==0) are still
run through OSCAR and recorded -- DRAMA predicting a re-entry where OREM's
own pipeline didn't is itself a useful finding, just not something with an
orem_rpe_pct to compare against.

Usage:
    .venv\\Scripts\\python run_oscar_campaign.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from drama import oscar
from orem_ground_truth import load_ground_truth
from tle_utils import read_tle_near_epoch, tle_to_oscar_elements

REPO_ROOT = Path(__file__).resolve().parents[2]
RPE_CSV = REPO_ROOT / "scratch_rpe" / "rpe_campaign.csv"
OUTPUT_CSV = Path(__file__).resolve().parent.parent / "output" / "oscar_campaign_results.csv"

DRAG_COEFFICIENT = 2.2
CROSS_SECTION_AREA_M2 = 10.0
NCPUS = 4  # repo-wide 4-core cap (GitHub\CLAUDE.md SS1) -- never leave this at the module default


def load_primary_zones():
    """Return (run_zones, primary_norads, by_norad).

    run_zones: {norad: row_dict} -- the zone to actually run through OSCAR for
        every object that has at least one real (non-ERR) row. Preference is
        OREM's own PRIMARY zone (latest zone with reentry_jd > 0); objects
        with no such zone still get a run, using their own latest zone
        regardless of validity (e.g. NORAD 21670: only zone has zstat=2,
        reentry_jd==0 -- OREM has no prediction, but DRAMA/OSCAR can still be
        asked what it thinks, using that zone's own bn_opt/zepoch).
    primary_norads: the subset of run_zones keys that have a real OREM
        prediction to compare against (i.e. an orem_rpe_pct exists).
    """
    by_norad = defaultdict(list)
    with open(RPE_CSV, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < len(header):
                continue  # ERR rows (pipeline error, no zone data)
            rec = dict(zip(header, row))
            try:
                rec["norad"] = int(rec["norad"])
                rec["zone"] = int(rec["zone"])
                rec["reentry_jd"] = float(rec["reentry_jd"])
                rec["rpe_pct"] = float(rec["rpe_pct"])
                rec["bn_opt"] = float(rec["bn_opt"])
                rec["zepoch"] = float(rec["zepoch"])
            except (KeyError, ValueError):
                continue
            by_norad[rec["norad"]].append(rec)

    run_zones = {}
    primary_norads = set()
    for norad, rows in by_norad.items():
        valid = [r for r in rows if r["reentry_jd"] > 0]
        if valid:
            run_zones[norad] = max(valid, key=lambda r: r["zone"])
            primary_norads.add(norad)
        else:
            run_zones[norad] = max(rows, key=lambda r: r["zone"])
    return run_zones, primary_norads, by_norad


def build_configs(run_zones, ground_truth):
    configs = []
    meta = []  # parallel list: (norad, obs_jd, zone_row) per config
    skipped_no_tle = []
    for norad, row in run_zones.items():
        tle_path = REPO_ROOT / "input" / f"example_{norad}.tle.txt"
        if not tle_path.exists():
            skipped_no_tle.append(norad)
            continue
        line1, line2 = read_tle_near_epoch(tle_path, row["zepoch"])
        elements = tle_to_oscar_elements(line1, line2)
        mass_kg = row["bn_opt"] * DRAG_COEFFICIENT * CROSS_SECTION_AREA_M2
        config = {
            **elements,
            "runId": str(norad),
            "spacecraftMass": mass_kg,
            "spacecraftCrossSectionArea": CROSS_SECTION_AREA_M2,
            "dragCoefficient": DRAG_COEFFICIENT,
            "disposalOption": 6,
            "propagationTime": 100.0,
        }
        configs.append(config)
        meta.append((norad, ground_truth.get(norad), row))
    return configs, meta, skipped_no_tle


def main():
    run_zones, primary_norads, by_norad = load_primary_zones()
    ground_truth = load_ground_truth()
    print(f"{len(by_norad)} objects in campaign CSV; "
          f"{len(primary_norads)} have an OREM primary-zone prediction; "
          f"{len(run_zones)} total will be run through OSCAR.")

    configs, meta, skipped_no_tle = build_configs(run_zones, ground_truth)
    if skipped_no_tle:
        print(f"Skipped {len(skipped_no_tle)} objects with no TLE file: {skipped_no_tle}")

    print(f"Running OSCAR for {len(configs)} objects (parallel, ncpus={NCPUS})...")
    results = oscar.run(config=configs, parallel=True, ncpus=NCPUS)

    # oscar.run's 'results' + 'errors' lists don't preserve input order/count 1:1
    # with 'config' when some runs error, and runId is silently ignored by
    # pydrama (always echoed back as "PyOscar" -- confirmed via _test_batch.py),
    # so match back up via the run's own echoed orbital elements instead --
    # practically guaranteed unique across 97 distinct real objects.
    def config_key(cfg):
        return (
            round(float(cfg["semiMajorAxis"]), 3),
            round(float(cfg["eccentricity"]), 6),
            round(float(cfg["inclination"]), 4),
        )

    ok_by_key = {config_key(r["config"]): r for r in results["results"]}
    err_by_key = {config_key(e["config"]): e for e in results["errors"]}

    OUTPUT_CSV.parent.mkdir(exist_ok=True)
    n_ok = n_no_reentry = n_error = n_no_orem_primary = 0
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "norad", "primary_zone", "e_tle", "inc_deg", "perigee_alt_km",
            "bn_opt", "zepoch_jd", "orem_reentry_jd", "orem_rpe_pct",
            "oscar_reentry_jd", "oscar_rpe_pct", "obs_jd", "oscar_status",
        ])
        for cfg, (norad, obs_jd, row) in zip(configs, meta):
            key = config_key(cfg)
            zepoch = row["zepoch"]
            has_orem = norad in primary_norads
            if not has_orem:
                n_no_orem_primary += 1
            perigee_alt = cfg["semiMajorAxis"] * (1 - cfg["eccentricity"]) - 6378.14

            if key in err_by_key:
                n_error += 1
                writer.writerow([
                    norad, row["zone"], cfg["eccentricity"], cfg["inclination"],
                    perigee_alt, row["bn_opt"], zepoch,
                    row["reentry_jd"] if has_orem else "",
                    row["rpe_pct"] if has_orem else "", "", "",
                    obs_jd if obs_jd else "", "run_error",
                ])
                continue

            result = ok_by_key[key]
            if not result["reentry"]:
                n_no_reentry += 1
                oscar_jd = ""
                oscar_rpe = ""
                status = "no_reentry_in_100yr"
            else:
                from tle_utils import ymd_to_jd
                fs_epoch = result["final_state"]["epoch"]
                oscar_jd = ymd_to_jd(fs_epoch.year, fs_epoch.month, fs_epoch.day) + (
                    fs_epoch.hour / 24 + fs_epoch.minute / 1440 + fs_epoch.second / 86400
                )
                if obs_jd:
                    horizon = obs_jd - zepoch
                    oscar_rpe = (oscar_jd - obs_jd) / horizon * 100.0 if horizon > 0 else ""
                else:
                    oscar_rpe = ""
                status = "ok"
                n_ok += 1

            writer.writerow([
                norad, row["zone"], cfg["eccentricity"], cfg["inclination"],
                perigee_alt, row["bn_opt"], zepoch,
                row["reentry_jd"] if has_orem else "",
                row["rpe_pct"] if has_orem else "",
                oscar_jd, oscar_rpe, obs_jd if obs_jd else "", status,
            ])

    print(f"\nWrote {OUTPUT_CSV}")
    print(f"  ok (OSCAR predicted re-entry): {n_ok}")
    print(f"  no re-entry within 100yr:      {n_no_reentry}")
    print(f"  OSCAR run errors:              {n_error}")
    print(f"  (of which, no OREM primary):   {n_no_orem_primary}")


if __name__ == "__main__":
    main()
