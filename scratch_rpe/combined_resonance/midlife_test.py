"""Mid-life cutoff test (OREM #56 follow-up).  Subcommands:
    prep            write TLE files truncated at each object's first perigee dip < 300 km after its peak
    orem NORAD      run production OREM (v1.48 exe, Watchlist cfg conventions) on the truncated history
    forecast        secular forecast started at the same epoch
    summary         compare both with the real decay
Nothing here touches the Watchlist database or the Object Detail cache; all files under E:/claude/p58/midlife."""
import sys, json, math, datetime as dt
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, 'E:/GitHub/OREM-Watchlist/src')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe'); sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/resonance_plots')
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')

WORK = Path('E:/claude/p58/midlife'); WORK.mkdir(parents=True, exist_ok=True)
D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
KEY = [8844, 9269, 11073]            # second-hump objects (forecast failed at 300 km)
CONTROL = [7480, 8833, 14297, 7903]  # ordinary objects: the dip is the final approach
_CUT = pd.read_csv(f"{D}/cross_check_orem_300.csv")        # the 17 Molniya-type objects with a post-peak <300 km event
NEW = [int(n) for n in _CUT["norad"] if int(n) not in KEY + CONTROL]   # not examined before the 150-vs-200 km choice
ALL = KEY + CONTROL + NEW


def group_of(n):
    return "second-hump" if n in KEY else "control" if n in CONTROL else "new (out-of-sample)"
ROOT = "E:/Research/1. R&D/Re-entry/2026/Data/objects"


def jd_to_date(jd):
    return dt.date(2000, 1, 1) + dt.timedelta(days=jd - 2451544.5)


def cutoffs():
    df = pd.read_csv(f'{D}/cross_check_orem_300.csv').set_index('norad')
    return df


def prep():
    from heowatch.lunisolar_drive import _cal_to_jd
    df = cutoffs()
    rows = []
    for n in ALL:
        c = float(df.loc[n, 'e300_obs'])
        lines = Path(f"{ROOT}/{n}/{n}.tle.txt").read_text(errors='ignore').splitlines()
        keep, k, kept = [], 0, 0
        while k < len(lines) - 1:
            if lines[k].startswith('1 ') and lines[k + 1].startswith('2 '):
                l1 = lines[k]; yy = int(l1[18:20]); year = 2000 + yy if yy < 57 else 1900 + yy
                jd = _cal_to_jd(year, float(l1[20:32]))
                if jd <= c:
                    if k > 0 and lines[k - 1].startswith('0 '):
                        keep.append(lines[k - 1])
                    keep += [lines[k], lines[k + 1]]; kept += 1
                k += 2
            else:
                k += 1
        out = WORK / f'{n}.tle.txt'
        out.write_text('\n'.join(keep) + '\n')
        decay = float(df.loc[n, 'decay'])
        rows.append(dict(norad=n, cutoff_jd=c, cutoff_date=str(jd_to_date(c)), decay_jd=decay, decay_date=str(jd_to_date(decay)),
                         lead_to_decay_d=decay - c, n_tle=kept, group=group_of(n)))
        print(n, 'cutoff', jd_to_date(c), 'decay', jd_to_date(decay), f'lead {decay - c:.0f} d', 'TLEs kept', kept)
    pd.DataFrame(rows).to_csv(f'{D}/midlife_cutoffs.csv', index=False)


def run_orem(n):
    from heowatch.config import load_config
    from heowatch.object_info import get_object_params
    from heowatch.orem_wrapper import bn_range_from_params, write_orem_cfg, run_and_parse, OremRunError
    cfg = load_config()
    exe = Path(cfg['orem_exe_path'])
    tle = WORK / f'{n}.tle.txt'
    params = get_object_params(n, cfg)
    lo, hi = bn_range_from_params(params, cfg)
    cfg_path = WORK / f'orem_{n}.cfg'
    write_orem_cfg(cfg_path, tle, n, observed_reentry=None, bn_lo=lo, bn_hi=hi)
    report = exe.parent / 'output' / f"OREM_{n:05d}_{dt.date.today().strftime('%Y%m%d')}.txt"
    res = {'norad': n}
    try:
        r = run_and_parse(exe, cfg_path, report, n)
        res.update(primary=r.primary_reentry_utc, ensemble_mean=r.ensemble_mean_utc, ensemble_std_days=r.ensemble_std_days,
                   n_zones=len(r.zones), no_reentry_predicted=bool(getattr(r, 'no_reentry_predicted', False)), status='ok')
        txt = report.read_text(encoding='utf-8', errors='ignore')
        res['report_excerpt'] = '\n'.join(l for l in txt.splitlines() if 'PRIMARY' in l or 'STALE' in l or 'Ensemble' in l)[:600]
    except OremRunError as exc:
        res.update(status='error', error=str(exc)[:300])
    json.dump(res, open(WORK / f'orem_result_{n}.json', 'w'), indent=1)
    print(json.dumps(res, indent=1)[:900])


def forecast():
    from validate_drive import clean_rows
    from prepare_plot_data_fullpool import load_full_record
    from secular_propagator import propagate, RE
    cut = pd.read_csv(f'{D}/midlife_cutoffs.csv').set_index('norad')
    # observed gap (decay - E150/E200) from the full-record cross-check, leave-one-out medians
    g = {thr: pd.read_csv(f'{D}/cross_check_orem_{thr}.csv').set_index('norad')['gap'] for thr in (200, 150)}
    out = []
    for n in ALL:
        c = float(cut.loc[n, 'cutoff_jd'])
        rec = np.array(clean_rows(load_full_record(f"{ROOT}/{n}/{n}.tle.txt"))); jd, a, e, inc, O, w = rec.T
        i1 = int(np.searchsorted(jd, c, side='right')); sl = slice(max(0, i1 - 11), i1)
        e0, i0, aa = float(np.median(e[sl])), float(np.median(inc[sl])), float(np.median(a[sl]))
        O0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(O[sl]))))) % 360)
        w0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(w[sl]))))) % 360)
        traj = propagate(e0, i0, O0, w0, aa, jd[i1 - 1], years=15.0, dt_days=5.0)
        hpp = aa * (1 - traj[:, 1]) - RE
        row = dict(norad=n, hp_at_cutoff=float(np.median(a[sl] * (1 - e[sl]) - RE)), hp_max_next_8yr=float(np.max(hpp[: int(8 * 365.25 / 5)])))
        for thr in (200, 150):
            hit = np.where(hpp < thr)[0]
            gap = float(np.median(g[thr].drop(index=n, errors='ignore'))) if n in g[thr].index or True else float('nan')
            if len(hit):
                row[f'event_{thr}_d_after_cutoff'] = float(traj[hit[0], 0] - traj[0, 0])
                row[f'fc_decay_jd_{thr}'] = float(traj[hit[0], 0] + gap)
            else:
                row[f'event_{thr}_d_after_cutoff'] = float('nan'); row[f'fc_decay_jd_{thr}'] = float('nan')
            row[f'gap_loo_{thr}'] = gap
        out.append(row)
        print(n, {k: round(v, 1) for k, v in row.items() if k != 'norad'})
    pd.DataFrame(out).to_csv(f'{D}/midlife_forecast.csv', index=False)


def summary():
    cut = pd.read_csv(f'{D}/midlife_cutoffs.csv').set_index('norad')
    fc = pd.read_csv(f'{D}/midlife_forecast.csv').set_index('norad')
    rows = []
    for n in ALL:
        decay = float(cut.loc[n, 'decay_jd'])
        try:
            o = json.load(open(WORK / f'orem_result_{n}.json'))
        except FileNotFoundError:
            o = {'status': 'not run'}
        def to_jd(s):
            return (dt.date.fromisoformat(s[:10]) - dt.date(1970, 1, 1)).days + 2440587.5 if s else float('nan')
        orem_jd = to_jd(o.get('primary')) if o.get('primary') else (to_jd(o.get('ensemble_mean')) if o.get('ensemble_mean') else float('nan'))
        r = dict(norad=n, group=cut.loc[n, 'group'], lead_d=round(float(cut.loc[n, 'lead_to_decay_d'])), orem_status=o.get('status'),
                 orem_primary=o.get('primary'), orem_noreentry=o.get('no_reentry_predicted'),
                 orem_err_d=round(orem_jd - decay) if not math.isnan(orem_jd) else None,
                 hp_cut=round(fc.loc[n, 'hp_at_cutoff']), hp_max8=round(fc.loc[n, 'hp_max_next_8yr']))
        for thr in (200, 150):
            v = fc.loc[n, f'fc_decay_jd_{thr}']
            r[f'fc{thr}_err_d'] = None if math.isnan(v) else round(v - decay)
        rows.append(r)
    t = pd.DataFrame(rows)
    pd.set_option('display.width', 220)
    print(t.to_string(index=False))
    t.to_csv(f'{D}/midlife_summary.csv', index=False)
    print()
    for lab, sel in (("all", t), ("seen before (7)", t[t.group != "new (out-of-sample)"]), ("NEW out-of-sample", t[t.group == "new (out-of-sample)"])):
        if len(sel) == 0:
            continue
        line = f"{lab:<20} n={len(sel):2d}  OREM gave a date: {int(sel['orem_err_d'].notna().sum())}/{len(sel)}"
        oe = sel['orem_err_d'].dropna()
        if len(oe):
            line += f" (median |err| {oe.abs().median():.0f} d)"
        for thr in (200, 150):
            e = sel[f'fc{thr}_err_d'].dropna()
            line += f" | forecast {thr} km: median |err| {e.abs().median():.0f} d, mean {e.abs().mean():.0f} d, within 250 d {int((e.abs() <= 250).sum())}/{len(e)}, >1 yr off {int((e.abs() > 365).sum())}"
        print(line)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep': prep()
    elif cmd == 'orem': run_orem(int(sys.argv[2]))
    elif cmd == 'forecast': forecast()
    elif cmd == 'summary': summary()
