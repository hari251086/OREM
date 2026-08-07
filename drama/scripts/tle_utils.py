"""Shared helpers for turning OREM's TLE inputs into DRAMA/OSCAR run configs."""

from datetime import datetime, timedelta

from sgp4.api import Satrec

MU_EARTH = 398600.4418  # km^3/s^2


def read_first_tle(tle_path):
    """Return the (line1, line2) of the first TLE entry in an OREM input file.

    OREM's example_*.tle.txt files are a name line (optionally) followed by
    the standard two TLE lines, repeated for many epochs. This grabs the
    first entry only.
    """
    line1 = line2 = None
    with open(tle_path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("1 "):
                line1 = line
            elif line.startswith("2 ") and line1 is not None:
                line2 = line
                break
    if line1 is None or line2 is None:
        raise ValueError(f"No TLE entry found in {tle_path}")
    return line1, line2


def tle_epoch_to_datetime(line1):
    epoch_year = int(line1[18:20])
    epoch_day = float(line1[20:32])
    year = 1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year
    return datetime(year, 1, 1) + timedelta(days=epoch_day - 1)


def ymd_to_jd(year, month, day):
    """Julian Date at 0h UT for a Gregorian calendar date (standard algorithm)."""
    y, m = year, month
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + day + b - 1524.5


def tle_epoch_to_jd(line1):
    epoch_year = int(line1[18:20])
    epoch_day = float(line1[20:32])
    year = 1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year
    return ymd_to_jd(year, 1, 1) + (epoch_day - 1)


def read_all_tles(tle_path):
    """Return every (line1, line2) TLE entry in an OREM input file, in order."""
    entries = []
    line1 = None
    with open(tle_path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("1 "):
                line1 = line
            elif line.startswith("2 ") and line1 is not None:
                entries.append((line1, line))
                line1 = None
    if not entries:
        raise ValueError(f"No TLE entries found in {tle_path}")
    return entries


def read_tle_near_epoch(tle_path, target_jd):
    """Return the (line1, line2) TLE entry whose epoch is closest to target_jd."""
    entries = read_all_tles(tle_path)
    return min(entries, key=lambda e: abs(tle_epoch_to_jd(e[0]) - target_jd))


def tle_to_oscar_elements(line1, line2):
    """Convert a TLE pair into the classical elements OSCAR's config expects.

    Uses sgp4 to parse the fixed-width TLE fields (rather than hand-rolling
    column slicing) and Kepler's third law to turn TLE mean motion into
    semi-major axis.
    """
    sat = Satrec.twoline2rv(line1, line2)
    n_rad_per_s = sat.no_kozai / 60.0
    sma_km = (MU_EARTH / n_rad_per_s**2) ** (1 / 3)

    import math

    return {
        "beginDate": tle_epoch_to_datetime(line1),
        "semiMajorAxis": sma_km,
        "eccentricity": sat.ecco,
        "inclination": math.degrees(sat.inclo),
        "rightAscensionOfTheAscendingNode": math.degrees(sat.nodeo),
        "argumentOfPerigee": math.degrees(sat.argpo),
        "meanAnomaly": math.degrees(sat.mo),
    }
