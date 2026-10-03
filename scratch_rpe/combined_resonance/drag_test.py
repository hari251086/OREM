"""Forecast WITH drag: event = decay itself (perigee altitude < 80 km), no gap conversion.
Modes: (a) mid-life cutoff = first dip < 300 km after the peak (as midlife_test.py); (b) from the first TLE (full life).
B = Cd A / m from the Watchlist object params (DISCOS where available; independent of OREM's fits), with B x0.5 / x1 / x2 sensitivity.
Compared with the earlier no-drag + leave-one-out-gap forecasts."""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                              # noqa: E402
from heowatch.object_info import get_object_params                   # noqa: E402
from prepare_plot_data_fullpool import load_full_record              # noqa: E402
from validate_drive import clean_rows                                 # noqa: E402
from secular_drag_prop import propagate_drag                          # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
cfg = load_config()
cut = pd.read_csv(f'{D}/midlife_cutoffs.csv').set_index('norad')
prev_mid = pd.read_csv(f'{D}/midlife_summary.csv').set_index('norad')
prev_full = pd.read_csv(f'{D}/cross_check_orem_300.csv').set_index('norad')
SCALES = (0.5, 1.0, 2.0)
rows = []
for n in [int(x) for x in cut.index]:
    p = get_object_params(n, cfg)
    B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6                      # km^2/kg
    rec = np.array(clean_rows(load_full_record(f"{ROOT}/{n}/{n}.tle.txt"))); jd, a, e, inc, O, w = rec.T
    decay = float(cut.loc[n, 'decay_jd']); c = float(cut.loc[n, 'cutoff_jd'])

    def init(sl):
        return (float(np.median(a[sl])), float(np.median(e[sl])), float(np.median(inc[sl])),
                float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360),
                float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360))
    i1 = int(np.searchsorted(jd, c, side='right'))
    st_mid = init(slice(max(0, i1 - 11), i1)); jd_mid = jd[i1 - 1]
    st_full = init(slice(0, 11)); jd_full = jd[0]
    r = dict(norad=n, BN=1.0 / (B0 * 1e6) * 1.0)
    r['BN_kg_m2'] = p.mass_kg / (p.cd * p.area_m2)
    for s in SCALES:
        dm, _ = propagate_drag(*st_mid, jd_mid, B0 * s, years=15.0)
        df_, _ = propagate_drag(*st_full, jd_full, B0 * s, years=30.0)
        r[f'mid_err_x{s}'] = dm - decay
        r[f'full_err_x{s}'] = df_ - decay
    r['mid_prev_150'] = prev_mid.loc[n, 'fc150_err_d']
    r['full_prev_300'] = prev_full.loc[n, 'err_fc_d'] if n in prev_full.index else np.nan
    r['lead_mid_d'] = decay - c
    rows.append(r)
    print(n, {k: (round(v) if isinstance(v, float) and not np.isnan(v) else v) for k, v in r.items() if k not in ('norad', 'BN')}, flush=True)
df = pd.DataFrame(rows)
df.to_csv(f'{D}/drag_test.csv', index=False)

def stats(col):
    allv = df[col]
    x = allv.dropna()
    return (f"decay date produced {len(x):2d}/{len(allv)} | of those: median |err| {x.abs().median():6.0f} d  mean {x.abs().mean():6.0f}  bias {x.median():+6.0f} | "
            f"within 250 d {int((x.abs() <= 250).sum())}/{len(allv)} (no date counted as a miss)  >1 yr off {int((x.abs() > 365).sum())}")
print('\n=== mid-life cutoff (first dip): decay-date error, days ===')
for s in SCALES:
    print(f"drag, B x{s:<4}: ", stats(f'mid_err_x{s}'))
print("no drag + gap, 150 km event:", stats('mid_prev_150'))
print('\n=== from the first TLE (full life, ~10-15 yr ahead) ===')
for s in SCALES:
    print(f"drag, B x{s:<4}: ", stats(f'full_err_x{s}'))
print("no drag + gap, 300 km event:", stats('full_prev_300'))
