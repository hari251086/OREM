"""User request: replace the 3 Falcon9 R/B objects (44187, 56758, 41553) in
the 50-object campaign -- suspected of active deorbit/passivation burns,
which OREM's drag-only physics (no thrust modeling) cannot represent -- with
non-maneuverable rocket-body/debris replacements, matching the established
sourcing criteria (APOGEE>8000km, PERIGEE<3000km, ROCKET BODY or DEBRIS,
decay after 2015), avoiding objects already excluded in prior rounds (the
2025-05-07/08 SL-12 debris cluster already represented once in the set;
38587/27900, previously rejected for a 12-16yr last-TLE-to-decay gap;
CZ-5 apogee>280,000km trans-lunar-injection-stage exclusions).

Uses OREM-Watchlist's existing, policy-compliant Space-Track client
(SPACE_TRACK_USAGE_POLICY.md requires routing every query through
heowatch.spacetrack_client.get_client(), never a second raw client) --
imported cross-repo via sys.path, not reimplemented. gp_history is a
"1/lifetime per range" data class; this is each candidate's first-ever
fetch (full backfill), same as every other object already in this project's
input/ directory.
"""
import sys
sys.path.insert(0, r"C:\Users\hari2\OneDrive\Documents\GitHub\OREM-Watchlist\src")

from heowatch.spacetrack_client import get_client
from spacetrack import operators as op

CANDIDATES = {
    60607: "CZ-7A R/B",
    36518: "Breeze-M deb (tank)",
    40777: "Ariane 1 deb",
}

st = get_client()

for norad, name in CANDIDATES.items():
    print(f"fetching {norad} ({name}) gp_history...")
    tle_text = st.gp_history(
        norad_cat_id=norad, orderby="epoch", format="tle"
    )
    out_path = f"input/example_{norad}.tle.txt"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(tle_text)
    n_lines = tle_text.count("\n")
    print(f"  wrote {out_path} ({n_lines} lines, ~{n_lines//2} TLEs)")
