"""Parse OREM's real observed re-entry dates out of scratch_rpe/rpe_campaign.F.

The campaign's ground truth is hardcoded there as three parallel Fortran DATA
arrays (norad / oyr / omo / ody), not exposed in any of the campaign's CSV
outputs. Parsing it directly (rather than re-deriving observed dates by
inverting OREM's own rpe_pct formula) keeps this tied to the actual source of
truth, immune to any rpe_pct-formula drift between the many campaign CSV
variants in scratch_rpe/.
"""

import re
from pathlib import Path

from tle_utils import ymd_to_jd

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_F = REPO_ROOT / "scratch_rpe" / "rpe_campaign.F"


def _parse_data_array(text, name):
    """Extract the numbers in a Fortran `data <name> / ... /` statement.

    Handles the standard fixed-form continuation convention used throughout
    this file: a line starting with `     &` continues the previous line's
    list, up to the closing `/`.
    """
    match = re.search(
        rf"data\s+{name}\s*/(.*?)/", text, re.IGNORECASE | re.DOTALL
    )
    if not match:
        raise ValueError(f"DATA array '{name}' not found in {CAMPAIGN_F}")
    body = match.group(1)
    # Strip Fortran double-precision exponent suffixes ("2019.d0" -> "2019.")
    # before scanning for numbers -- otherwise the '0' in 'd0' is picked up
    # as its own bogus token.
    body = re.sub(r"[dD]0", "", body)
    tokens = re.findall(r"-?\d+\.?\d*", body)
    return [float(t) for t in tokens]


def load_ground_truth():
    """Return {norad_id (int): observed_reentry_jd (float)} for all campaign objects."""
    text = CAMPAIGN_F.read_text()
    norad = _parse_data_array(text, "norad")
    oyr = _parse_data_array(text, "oyr")
    omo = _parse_data_array(text, "omo")
    ody = _parse_data_array(text, "ody")
    assert len(norad) == len(oyr) == len(omo) == len(ody), (
        f"DATA array length mismatch: norad={len(norad)} oyr={len(oyr)} "
        f"omo={len(omo)} ody={len(ody)}"
    )
    return {
        int(n): ymd_to_jd(int(y), int(m), d)
        for n, y, m, d in zip(norad, oyr, omo, ody)
    }


if __name__ == "__main__":
    gt = load_ground_truth()
    print(f"Parsed {len(gt)} objects' ground-truth re-entry dates.")
    for norad in (42928, 35497, 39615):
        print(f"  {norad}: JD {gt[norad]:.4f}")
