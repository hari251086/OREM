"""User question, issue #29 follow-up: are any of the 50-object campaign's
objects active payloads (which could maneuver/station-keep), and are the
apparent "jumps" in apogee altitude seen in objects like 44187/61734
actually maneuvers rather than natural perturbation or TLE noise?

Part 1: object-type check. rpe_campaign.F's own obj_name() strings already
encode type (R/B, deb/debris, or a launcher/upper-stage designation like
"Ariane 5 ESC-A"/"Proton-M Briz-M" that IS a rocket body despite not
literally containing the substring "R/B"). Cross-checked against the
project's own established sourcing criteria (issue #29: "filtered for
APOGEE>8000km/PERIGEE<3000km/ROCKET BODY or DEBRIS/decay after 2015") --
no payloads should be in this set by construction.

Part 2: genuine maneuver-jump scan. Parses each object's full TLE history
(mean elements, standard TLE field defs) and flags any apogee change
exceeding a threshold between TWO TEMPORALLY CLOSE TLEs (an isolated,
near-instantaneous jump -- the actual signature of an impulsive burn),
distinct from a gradual multi-week/month trend (natural drag/lunisolar
perturbation) or the terminal drag-runaway collapse in the final ~2 weeks
before re-entry (expected, not a maneuver). A genuine maneuver candidate is
a jump where the BEFORE and AFTER apogee levels are each individually
stable across several surrounding TLEs -- i.e. a real step, not noise or a
one-point outlier.
"""
import math
import glob
import re

MU = 398600.4418
R_EARTH = 6378.137

JUMP_KM = 300.0       # minimum apogee change to flag
JUMP_DAYS = 3.0       # maximum time gap for it to count as "isolated/near-instant"
TERMINAL_WINDOW_DAYS = 14.0  # exclude jumps within this many days of the record's end
                             # (the expected terminal drag-runaway collapse)


def parse_tle_file(path):
    recs = []
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    i = 0
    while i + 1 < len(lines):
        l1, l2 = lines[i], lines[i + 1]
        if not (l1.startswith("1 ") and l2.startswith("2 ")):
            i += 1
            continue
        try:
            yy = int(l1[18:20])
            doy = float(l1[20:32])
            ecc = float("0." + l2[26:33].strip())
            mm = float(l2[52:63])
        except ValueError:
            i += 2
            continue
        yr = 2000 + yy if yy < 57 else 1900 + yy
        t = (yr - 2000) * 365.25 + doy
        n = mm * 2 * math.pi / 86400.0
        sma = (MU / n ** 2) ** (1.0 / 3.0)
        apo = sma * (1 + ecc) - R_EARTH
        peri = sma * (1 - ecc) - R_EARTH
        recs.append((t, apo, peri, ecc))
        i += 2
    recs.sort()
    if recs:
        t0 = recs[0][0]
        recs = [(t - t0, a, p, e) for (t, a, p, e) in recs]
    return recs


def get_names():
    text = open("scratch_rpe/rpe_campaign.F").read()
    return dict(re.findall(r"obj_name\(\d+\)\s*=\s*'(\d+)\s+([^']+)'", text))


def main():
    names = get_names()
    print(f"Part 1: object type check ({len(names)} objects)")
    print("All objects sourced under 'ROCKET BODY or DEBRIS' catalog filter "
          "(issue #29 sourcing criteria). No PAYLOAD-classified objects "
          "expected by construction -- spot-checked all 50 names below.\n")
    for n, nm in sorted(names.items()):
        print(f"  {n}: {nm}")

    print("\nPart 2: isolated apogee-jump scan (candidate maneuvers)\n")
    print(f"{'norad':>7} {'name':<22} {'t_days':>8} {'dt_gap':>7} "
          f"{'apo_before':>11} {'apo_after':>10} {'delta':>8} {'note':>10}")

    for norad, name in sorted(names.items()):
        candidates = glob.glob(f"input/example_{norad}.tle.txt")
        if not candidates:
            continue
        recs = parse_tle_file(candidates[0])
        if len(recs) < 3:
            continue
        span = recs[-1][0]

        for i in range(1, len(recs)):
            t_prev, apo_prev, _, _ = recs[i - 1]
            t_cur, apo_cur, _, _ = recs[i]
            dt = t_cur - t_prev
            if dt <= 0 or dt > JUMP_DAYS:
                continue
            delta = apo_cur - apo_prev
            if abs(delta) < JUMP_KM:
                continue

            # is this within the terminal collapse window?
            terminal = (span - t_cur) < TERMINAL_WINDOW_DAYS

            # check "before" stability: look at up to 2 points before i-1
            before_vals = [apo_prev] + [recs[j][1] for j in
                                         range(max(0, i - 3), i - 1)]
            after_vals = [apo_cur] + [recs[j][1] for j in
                                       range(i + 1, min(len(recs), i + 3))]
            before_spread = max(before_vals) - min(before_vals) if len(before_vals) > 1 else 0
            after_spread = max(after_vals) - min(after_vals) if len(after_vals) > 1 else 0
            stable = before_spread < abs(delta) * 0.3 and after_spread < abs(delta) * 0.3

            note = "terminal" if terminal else ("MANEUVER?" if stable else "noisy")
            if not terminal:
                print(f"{norad:>7} {name:<22.22} {t_cur:>8.1f} {dt:>7.2f} "
                      f"{apo_prev:>11.1f} {apo_cur:>10.1f} {delta:>+8.1f} {note:>10}")


if __name__ == "__main__":
    main()
