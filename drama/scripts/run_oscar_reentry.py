"""Run a DRAMA OSCAR natural-decay lifetime analysis for an OREM TLE input.

This is the first working link in OREM's DRAMA-based re-entry validation:
it takes one of OREM's existing curated `input/example_*.tle.txt` objects,
converts its first TLE entry into the classical elements OSCAR wants, runs
OSCAR with disposal option "none" (pure atmospheric-drag decay, matching
what OREM itself predicts), and writes the result to drama/output/.

Physical parameters (mass, cross-section area, drag coefficient) are not in
a TLE. Rather than guessing them, this script derives them from OREM's own
fitted ballistic number BN = mass / (Cd * Area) for the same object/epoch
(see OREM/scratch_rpe/rpe_campaign.csv), holding Cd and Area at OSCAR's own
defaults and solving for mass. That keeps the OSCAR run traceable to what
OREM already believes about this object, instead of an arbitrary guess.

Usage:
    .venv\\Scripts\\python run_oscar_reentry.py [norad_id]

Defaults to NORAD 21670 (H-1 R/B(2), the object this script was validated
against: OREM's own fit for this object/epoch is BN_opt=10.246 kg/m^2,
e_opt=0.739988 -- matching the TLE's own eccentricity of 0.7399072).
"""

import json
import sys
from datetime import datetime
from pathlib import Path

from drama import oscar
from tle_utils import read_first_tle, tle_to_oscar_elements

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"

# OREM's own fitted ballistic number for this object/epoch (rpe_campaign.csv,
# norad=21670, zone=1). BN = mass / (dragCoefficient * area).
OREM_BN_KG_PER_M2 = 10.246

# Held at OSCAR's own defaults so mass can be solved from OREM_BN_KG_PER_M2.
DRAG_COEFFICIENT = 2.2
CROSS_SECTION_AREA_M2 = 10.0


def main():
    norad_id = sys.argv[1] if len(sys.argv) > 1 else "21670"
    tle_path = REPO_ROOT / "input" / f"example_{norad_id}.tle.txt"
    line1, line2 = read_first_tle(tle_path)
    elements = tle_to_oscar_elements(line1, line2)

    mass_kg = OREM_BN_KG_PER_M2 * DRAG_COEFFICIENT * CROSS_SECTION_AREA_M2

    config = {
        **elements,
        "runId": f"orem_{norad_id}",
        "spacecraftMass": mass_kg,
        "spacecraftCrossSectionArea": CROSS_SECTION_AREA_M2,
        "dragCoefficient": DRAG_COEFFICIENT,
        "disposalOption": 6,  # none -- pure natural decay, matches OREM
        "propagationTime": 100.0,  # years; OSCAR default, plenty for this object
    }

    print(f"Running OSCAR for NORAD {norad_id} with config:")
    for k, v in config.items():
        print(f"  {k}: {v}")

    results = oscar.run(config=config, parallel=False)

    if results["errors"]:
        print("OSCAR run FAILED:")
        for err in results["errors"]:
            print(err["status"])
            print(err["output"])
        sys.exit(1)

    result = results["results"][0]
    print(f"\nlifetime (years): {result['lifetime']}")
    print(f"reentry within propagation span: {result['reentry']}")

    OUTPUT_DIR.mkdir(exist_ok=True)
    out_path = OUTPUT_DIR / f"oscar_{norad_id}_{datetime.now():%Y%m%d_%H%M%S}.json"
    with open(out_path, "w") as f:
        json.dump(
            {
                "config": {k: str(v) for k, v in config.items()},
                "lifetime_years": result["lifetime"],
                "reentry": result["reentry"],
                "final_state": result.get("final_state"),
            },
            f,
            indent=2,
            default=str,
        )
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
