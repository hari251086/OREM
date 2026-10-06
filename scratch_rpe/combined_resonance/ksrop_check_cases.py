"""Case table for the KSROP full-force check of the secular + drag forecast: every test case of THEORY sections 13-14 (mid-life cutoff, first TLE,
leads 0.25-4 yr on the 17 Molniya-type objects, leads 0.5/1/2 yr on the 40 out-of-sample pool objects).  Start elements are built exactly as in
drag_test.py / window_coverage_test.py / window_coverage_pool.py (median a, e, i and circular-mean RAAN, argp of the 11 TLEs ending at the start
epoch); mean anomaly is the last TLE's.  Also recomputes the secular forecast decay date for each case (same code) so both models can be compared
against the same real decay.  Writes ksrop_check_cases.csv."""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from heowatch.config import load_config                                              # noqa: E402
from heowatch.object_info import get_object_params                                   # noqa: E402
from prepare_plot_data_fullpool import BASE, build_object_list, load_full_record, parse_epoch_jd   # noqa: E402
from resonance_campaign_sweep import cal2jd                                           # noqa: E402
from validate_drive import clean_rows                                                 # noqa: E402
from secular_drag_prop import propagate_drag                                          # noqa: E402

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"
cfg = load_config()
cut = pd.read_csv(f'{D}/midlife_cutoffs.csv').set_index('norad')


def load_M(path):
    """jd -> TLE mean anomaly (deg), keyed by the same epoch the record loader uses."""
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]
    out = {}
    for k in range(0, len(lines) - 1, 2):
        l1, l2 = lines[k], lines[k + 1]
        if l1.startswith('1 ') and l2.startswith('2 '):
            out[round(parse_epoch_jd(l1), 6)] = float(l2[43:51])
    return out


def M_at(Mmap, jd0):
    keys = np.array(list(Mmap))
    j = int(np.argmin(np.abs(keys - jd0)))
    assert abs(keys[j] - jd0) < 1e-3, (jd0, keys[j])
    return Mmap[keys[j]]


def init(a, e, inc, O, w, sl):
    return (float(np.median(a[sl])), float(np.median(e[sl])), float(np.median(inc[sl])),
            float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360),
            float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360))


rows = []


def add(family, n, lead, st, jd0, M0, B0, decay, years, B_src):
    dj, _ = propagate_drag(*st, jd0, B0, years=years)
    rows.append(dict(case=f'{family}_{n}' + ('' if lead is None else f'_L{lead}'), family=family, norad=n, lead_yr=lead,
                     jd0=jd0, a=st[0], e=st[1], inc=st[2], raan=st[3], argp=st[4], M=M0, B_km2_kg=B0, BN_kg_m2=1.0 / (B0 * 1e6),
                     B_src=B_src, years_cap=years, real_decay_jd=decay, sec_decay_jd=dj))


# --- the 17 Molniya-type objects: mid-life, first TLE, and leads 0.25-4 yr -----------------------------------------------------------------
for n in [int(x) for x in cut.index]:
    p = get_object_params(n, cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6
    path = f"{ROOT}/{n}/{n}.tle.txt"
    rec = np.array(clean_rows(load_full_record(path))); jd, a, e, inc, O, w = rec.T
    Mmap = load_M(path)
    decay = float(cut.loc[n, 'decay_jd']); c = float(cut.loc[n, 'cutoff_jd'])
    i1 = int(np.searchsorted(jd, c, side='right'))
    add('mid', n, None, init(a, e, inc, O, w, slice(max(0, i1 - 11), i1)), jd[i1 - 1], M_at(Mmap, jd[i1 - 1]), B0, decay, 15.0, p.source)
    add('full', n, None, init(a, e, inc, O, w, slice(0, 11)), jd[0], M_at(Mmap, jd[0]), B0, decay, 30.0, p.source)
    for L in (0.25, 0.5, 1.0, 2.0, 4.0):
        t0 = decay - L * 365.25
        k = int(np.searchsorted(jd, t0, side='right'))
        if k < 12 or jd[k - 1] < t0 - 120:
            continue
        add('lead', n, L, init(a, e, inc, O, w, slice(k - 11, k)), jd[k - 1], M_at(Mmap, jd[k - 1]), B0, decay, 20.0, p.source)
    print(n, 'done', flush=True)

# --- the out-of-sample pool (as window_coverage_pool.py) ------------------------------------------------------------------------------
seen = set(int(x) for x in cut.index)
for o in build_object_list():
    n = o['norad']
    if n in seen:
        continue
    path = f"{BASE}/{o['tle_file']}"
    try:
        raw = load_full_record(path)
    except FileNotFoundError:
        continue
    if len(raw) < 300:
        continue
    inc_mean = sum(r[3] for r in raw) / len(raw)
    arr = np.array(clean_rows(raw)); jd, a, e, inc, O, w = arr.T
    if inc_mean < 46.4 or np.median(e[:20]) < 0.3:
        continue
    decay = cal2jd(*o['decay_ymd'])
    if decay < jd[0] or decay - jd[-1] > 400:
        continue
    try:
        p = get_object_params(n, cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6; src = p.source
    except Exception:
        B0 = 2.2 / 56.0 * 1e-6; src = 'default BN 56'
    Mmap = load_M(path)
    for L in (0.5, 1.0, 2.0):
        t0 = decay - L * 365.25
        k = int(np.searchsorted(jd, t0, side='right'))
        if k < 12 or jd[k - 1] < t0 - 120 or a[k - 11:k].min() <= 0:
            continue
        st = init(a, e, inc, O, w, slice(k - 11, k))
        if st[1] < 0.3 or st[2] < 46.4:
            continue
        add('pool', n, L, st, jd[k - 1], M_at(Mmap, jd[k - 1]), B0, decay, L + 4.0, src)
    print(n, 'done', flush=True)

df = pd.DataFrame(rows)
df['sec_err_d'] = df['sec_decay_jd'] - df['real_decay_jd']
df['real_horizon_d'] = df['real_decay_jd'] - df['jd0']
df.to_csv(f'{D}/ksrop_check_cases.csv', index=False)
print(f"\n{len(df)} cases: " + ', '.join(f"{k} {v}" for k, v in df['family'].value_counts().items()) + f"; {df['norad'].nunique()} objects")
