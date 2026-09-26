"""Manuscript (issue #45/#48): lead-time hindcast campaign for the 23
solar-apsidal-resonance objects from issue #56 (16 from the original
97-object campaign + 7 from the expanded pool, policy exclusions
already removed).

For every object and every lead time L, truncates the object's locally
stored TLE history (input/example_<norad>.tle.txt -- no Space-Track
query) to epochs strictly before (SATCAT decay date - L days) and
writes the truncated file + a run manifest for manuscript23_leadtime.F.
A prediction made from a truncated file uses no data from after its
cutoff -- a genuine hindcast at that lead time.

Also writes runs_index.csv: per-run metadata (name/type/inclination
from SATCAT, how many TLEs survived the cut, the gap between the last
stored TLE and the catalogued decay) so data-availability failures can
be reported honestly rather than silently dropped.

Run from anywhere:
    python scratch_rpe/manuscript23/prepare_leadtime_runs.py
"""
import csv
import sys
from pathlib import Path

import pandas as pd

OREM = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(OREM / "scratch_rpe"))
from resonance_campaign_sweep import cal2jd, jd_to_ymd, parse_epoch_jd  # noqa: E402

SATCAT = Path("E:/Research/1. R&D/Re-entry/2026/Data/catalog/latest/satcat.csv")
OUT = OREM / "scratch_rpe" / "manuscript23"
LEADS = [30, 90, 180, 365]

# issue #56: 16 original-97 crossings + 7 expanded-pool crossings
# (3 Falcon 9 R/Bs 39501/44187/42985 and 13913 already excluded there).
OBJECTS = [35497, 37151, 27526, 37819, 59347, 40943, 60328, 18983, 13112,
           28141, 19881, 5980, 4570, 10724, 28140, 44911,
           7096, 68538, 6294, 18954, 29462, 5941, 6231]


def tle_pairs(path):
    lines = [l.rstrip("\n") for l in open(path)]
    pairs = []
    for k in range(len(lines) - 1):
        if lines[k].startswith("1 ") and lines[k + 1].startswith("2 "):
            pairs.append((parse_epoch_jd(lines[k]), lines[k], lines[k + 1]))
    pairs.sort(key=lambda p: p[0])
    return pairs


def ymd(jd):
    y, m, d = jd_to_ymd(jd)
    return f"{y:04d}-{m:02d}-{d:02d}"


def main():
    sat = pd.read_csv(SATCAT).set_index("NORAD_CAT_ID")
    (OUT / "tle").mkdir(parents=True, exist_ok=True)
    index_rows, manifest = [], []

    for norad in OBJECTS:
        s = sat.loc[norad]
        dy, dm, dd = (int(x) for x in str(s["DECAY"])[:10].split("-"))
        decay_jd = cal2jd(dy, dm, dd)
        pairs = tle_pairs(OREM / "input" / f"example_{norad}.tle.txt")
        last_jd = pairs[-1][0]

        for lead in LEADS:
            cutoff_jd = decay_jd - lead
            kept = [p for p in pairs if p[0] < cutoff_jd]
            row = {
                "norad": norad, "name": s["OBJECT_NAME"], "type": s["OBJECT_TYPE"],
                "incl_deg": s["INCLINATION"], "decay": f"{dy:04d}-{dm:02d}-{dd:02d}",
                "lead_days": lead, "cutoff": ymd(cutoff_jd),
                "n_tle_total": len(pairs), "n_tle_kept": len(kept),
                "first_tle": ymd(pairs[0][0]),
                "last_kept_tle": ymd(kept[-1][0]) if kept else "",
                "last_tle_to_decay_days": round(decay_jd - last_jd, 1),
                "run_idx": "",
            }
            if kept:
                rel = f"scratch_rpe/manuscript23/tle/{norad}_L{lead}.tle.txt"
                with open(OREM / rel, "w") as f:
                    for _, l1, l2 in kept:
                        f.write(l1 + "\n" + l2 + "\n")
                manifest.append((norad, dy, dm, dd, lead, rel))
                row["run_idx"] = len(manifest)
            index_rows.append(row)

    # Two lines per run: "norad yr mo dy lead" then the TLE path on its
    # own line (list-directed Fortran input treats '/' as end-of-record,
    # so the path must be read with an explicit (A) format).
    with open(OUT / "runs_manifest.txt", "w") as f:
        f.write(f"{len(manifest)}\n")
        for norad, y, m, d, lead, rel in manifest:
            f.write(f"{norad} {y} {m} {d} {lead}\n{rel}\n")

    with open(OUT / "runs_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(index_rows[0].keys()))
        w.writeheader()
        w.writerows(index_rows)

    no_data = [r for r in index_rows if r["n_tle_kept"] == 0]
    print(f"{len(manifest)} runs written ({len(OBJECTS)} objects x {len(LEADS)} lead times, "
          f"{len(no_data)} object/lead combinations have no TLEs before the cutoff)")
    for r in no_data:
        print(f"  no data: {r['norad']} L={r['lead_days']} (first TLE {r['first_tle']}, cutoff {r['cutoff']})")


if __name__ == "__main__":
    main()
