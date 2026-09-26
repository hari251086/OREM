"""Manuscript (issues #45/#48): collate manuscript23_leadtime.F's
4-part CSV output into the summary table for the 23-object lead-time
campaign. Flags the two policy cases (59347 Falcon 9 R/B, 68538
Artemis II crewed capsule) rather than dropping them -- included in
every number per the explicit "include in everything" decision, but
footnoted so a reader can see them and reproduce the stats with or
without them.

Usage: python scratch_rpe/manuscript23/build_manuscript_table.py
"""
import glob
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
POLICY_NOTE = {
    59347: "Falcon 9 R/B -- possible post-separation deorbit burn (no thrust model in OREM)",
    68538: "ARTEMIS 2 (INTEGRITY) -- crewed Orion capsule, controlled re-entry, not natural decay",
}


def load(pattern):
    parts = sorted(glob.glob(str(HERE / pattern)))
    if not parts:
        raise SystemExit(f"no files matched {pattern} -- has the campaign finished?")
    return pd.concat([pd.read_csv(p, skipinitialspace=True) for p in parts], ignore_index=True)


def main():
    summ = load("results_summary_part_*.csv")
    idx = pd.read_csv(HERE / "runs_index.csv")

    summ.columns = [c.strip() for c in summ.columns]
    merged = idx.merge(summ, on=["norad", "lead_days"], how="left", suffixes=("", "_r"))
    merged["predicted"] = merged["npred"].fillna(0) > 0
    # orem_run initializes rpe_ens/rpe_pri to 0.0 before checking
    # whether any zone actually predicted a re-entry -- without this,
    # a genuine "no zone predicted re-entry" (npred=0) is indistinguishable
    # from a genuine 0.00% RPE. Null them out explicitly.
    merged.loc[~merged["predicted"], ["ens_rpe_pct", "primary_rpe_pct"]] = pd.NA
    merged["policy_case"] = merged["norad"].map(POLICY_NOTE)

    out_cols = ["norad", "name", "type", "incl_deg", "lead_days", "n_tle_kept",
                "predicted", "primary_rpe_pct", "ens_rpe_pct", "policy_case"]
    merged[out_cols].to_csv(HERE / "manuscript_table_full.csv", index=False)

    print(f"{merged['norad'].nunique()} objects x {sorted(merged['lead_days'].unique())} lead times "
          f"= {len(merged)} rows -> manuscript_table_full.csv\n")

    for lead in sorted(merged["lead_days"].unique()):
        sub = merged[merged["lead_days"] == lead]
        pred = sub[sub["predicted"]]
        n_pred = len(pred)
        n_total = sub["n_tle_kept"].gt(0).sum()  # only counts where data existed at all
        print(f"lead={lead:>3}d  predicted {n_pred}/{n_total} (data-available) "
              f"of {len(sub)} objects total")
        if n_pred:
            e = pred["ens_rpe_pct"].abs()
            print(f"    ensemble |RPE|  median={e.median():.1f}%  mean={e.mean():.1f}%  "
                  f"max={e.max():.1f}%")
            p = pred["primary_rpe_pct"].abs()
            print(f"    primary  |RPE|  median={p.median():.1f}%  mean={p.mean():.1f}%  "
                  f"max={p.max():.1f}%")
        no_data = sub[sub["n_tle_kept"] == 0]
        if len(no_data):
            print(f"    no TLE data before this cutoff: "
                  f"{', '.join(str(n) for n in no_data['norad'])}")

    print("\nPolicy cases (included above; flagged for the manuscript to discuss "
          "separately or as a sensitivity check):")
    for norad, note in POLICY_NOTE.items():
        print(f"  {norad}: {note}")


if __name__ == "__main__":
    main()
