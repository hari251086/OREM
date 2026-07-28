"""Issue #29 follow-up: test the Phase-0 literature hypothesis (Gupta &
Anilkumar 2014's best-in-class accuracy is only achieved in the terminal/
drag-dominated decay phase -- see project_orem_rpe_investigation_plan.md,
"Phase 0 Literature Reading" section) directly against the current 50-
object campaign, rather than leaving it as an unconfirmed framing note.

Hypothesis: whether an object predicts at all (within the 5-year horizon,
src/orem.F:587 `nrev_re = int(5*365.25/period_d)+1`) and how accurate that
prediction is, is explained by how close the object's LAST tracked TLE
already is to its actual terminal decay phase (low perigee altitude) --
not a uniform algorithm identifiability floor. The curated-7 were all
selected because they have real, documented decay events with GTO-debris-
like final TLEs; the broader 43 are a much more heterogeneous mix of
tracked-catalog objects, many of which may not be anywhere near end-of-
life at their last tracked TLE.

Read-only: parses each object's LAST '2 ' TLE line directly (no orem_run
call, no production code touched), same class of diagnostic as
diag_altitude_regime.py earlier in this investigation.
"""
import csv
import statistics
from collections import defaultdict

MU = 398600.4418  # km^3/s^2, standard Earth GM (coarse diagnostic only)
R_EARTH = 6378.137  # km, equatorial

CURATED_7 = {27526, 32007, 35497, 37151, 37819, 39615, 42928}

# norad -> tle_file, exactly as registered in rpe_campaign.F
TLE_FILES = {}
with open('scratch_rpe/rpe_campaign.F') as f:
    for line in f:
        line = line.strip()
        if line.startswith('tle_file(') and '=' in line:
            path = line.split("'")[1]
            norad = int(path.split('_')[-1].split('.')[0])
            TLE_FILES[norad] = path


def last_tle_perigee_alt(path):
    last_line2 = None
    with open(path) as f:
        for line in f:
            if line.startswith('2 '):
                last_line2 = line
    if last_line2 is None:
        return None
    ecc = float('0.' + last_line2[26:33].strip())
    mean_motion = float(last_line2[52:63])  # rev/day
    n_rad_s = mean_motion * 2 * 3.141592653589793 / 86400.0
    a = (MU / n_rad_s**2) ** (1.0 / 3.0)
    return a * (1.0 - ecc) - R_EARTH


# Load current campaign results
rows = defaultdict(list)
err_objs = set()
with open('scratch_rpe/rpe_campaign.csv', newline='') as f:
    for r in csv.DictReader(f):
        norad = int(r['norad'])
        if r['zone'] == 'ERR':
            err_objs.add(norad)
        else:
            rows[norad].append(r)


def status(norad):
    if norad in err_objs:
        return 'ERR', None, None
    zones = sorted(rows[norad], key=lambda r: int(r['zone']))
    latest = zones[-1]
    predicts = float(latest['reentry_jd']) > 0
    latest_rpe = abs(float(latest['rpe_pct'])) if predicts else None
    ens_rpe = abs(float(latest['ens_rpe_pct']))
    return ('predicts' if predicts else 'no-predict'), latest_rpe, ens_rpe


print(f"{'norad':>7} {'set':>9} {'status':>11} {'perigee_alt_km':>15} "
      f"{'latest_rpe':>11} {'ens_rpe':>8}")

by_status_perigee = defaultdict(list)
predict_rpe_pairs = []

for norad in sorted(TLE_FILES):
    path = TLE_FILES[norad]
    perigee = last_tle_perigee_alt(path)
    st, latest_rpe, ens_rpe = status(norad)
    setname = 'curated7' if norad in CURATED_7 else 'other43'
    by_status_perigee[st].append(perigee)
    if latest_rpe is not None:
        predict_rpe_pairs.append((perigee, latest_rpe, ens_rpe))
    lr = f"{latest_rpe:11.1f}" if latest_rpe is not None else "        n/a"
    er = f"{ens_rpe:8.1f}" if ens_rpe is not None else "     n/a"
    print(f"{norad:7d} {setname:>9} {st:>11} {perigee:15.1f} {lr} {er}")

print("\n--- Perigee altitude by predict status ---")
for st in ('ERR', 'no-predict', 'predicts'):
    vs = by_status_perigee.get(st, [])
    if vs:
        print(f"  {st:11s} n={len(vs):3d}  mean={statistics.mean(vs):8.1f} km  "
              f"median={statistics.median(vs):8.1f} km  "
              f"min={min(vs):8.1f}  max={max(vs):8.1f}")

print("\n--- Correlation: last-TLE perigee altitude vs. accuracy "
      "(objects that predict) ---")
peri = [p for p, _, _ in predict_rpe_pairs]
lrpe = [r for _, r, _ in predict_rpe_pairs]
erpe = [e for _, _, e in predict_rpe_pairs]
n = len(peri)
mp, ml, me = statistics.mean(peri), statistics.mean(lrpe), statistics.mean(erpe)


def pearson(xs, ys, mx, my):
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = (sum((x - mx) ** 2 for x in xs)) ** 0.5
    dy = (sum((y - my) ** 2 for y in ys)) ** 0.5
    return num / (dx * dy) if dx > 0 and dy > 0 else float('nan')


print(f"  n={n}")
print(f"  r(perigee, latest-zone |RPE|) = {pearson(peri, lrpe, mp, ml):.3f}")
print(f"  r(perigee, ensemble |RPE|)    = {pearson(peri, erpe, mp, me):.3f}")
