"""Find objects in the 253-object pool whose perigee first falls below 300 km after its peak and then RECOVERS (second hump):
the record continues > 2 yr after that first dip and the perigee rises above 500 km again."""
import sys
import numpy as np
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record   # noqa: E402
from validate_drive import clean_rows, RE                                          # noqa: E402
from molniya_cycle import smooth                                                   # noqa: E402

out = []
for o in build_object_list():
    try:
        rows = load_full_record(f"{BASE}/{o['tle_file']}")
    except FileNotFoundError:
        continue
    if len(rows) < 300:
        continue
    inc = sum(r[3] for r in rows) / len(rows)
    arr = np.array(clean_rows(rows)); jd, a, e = arr[:, 0], arr[:, 1], arr[:, 2]
    if inc < 46.4 or np.median(e[:20]) < 0.3:
        continue
    hp = a * (1 - e) - RE
    s = smooth(hp); ipk = int(np.argmax(s))
    idx = np.where((np.arange(len(hp)) > ipk) & (s < 300.0))[0]
    if len(idx) == 0:
        continue
    k = int(idx[0])
    lead = (jd[-1] - jd[k]) / 365.25
    rec = float(np.max(s[k:])) if k < len(s) - 1 else 0.0
    out.append((o['norad'], inc, lead, rec, (jd[k] - jd[0]) / 365.25, len(jd)))
out.sort(key=lambda r: -r[2])
print(f"{'NORAD':>6} {'i':>5} {'yr from 1st dip to last TLE':>28} {'max perigee after dip':>22} {'t_dip (yr)':>10} {'n':>6}")
for r in out[:25]:
    print(f"{r[0]:>6} {r[1]:5.1f} {r[2]:28.2f} {r[3]:22.0f} {r[4]:10.1f} {r[5]:6d}")
print('objects with a post-peak <300 km event:', len(out), '| second-hump candidates (lead > 2 yr and recovery > 500 km):',
      sum(1 for r in out if r[2] > 2.0 and r[3] > 500))
