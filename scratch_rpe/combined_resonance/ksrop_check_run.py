"""Run KSROP's driver_KS (KS-regularized, full force) on every case of ksrop_check_cases.csv and record the decay epoch it reports.

Force model: geopotential degree NGEO (EGM2008 head), Sun and Moon (Legendre degrees 2 and 3 as in KSROP's GMAT-validated setting), atmospheric drag with
OREM's static ATM.DAT (SW-All.csv / ATM2D.DAT deliberately NOT supplied, so the legacy static-table path runs: same atmosphere family as the secular model),
no SRP.  Co-rotation factor FR (1.0 = on; the secular model has none).  Decay = driver's own stop at altitude < 80 km.
Start state: the case's Keplerian elements (a, e, i, RAAN, argp, M) -> Cartesian, EME2000, epoch = the case's jd0.

Resumable (skips cases already in the results CSV).  <= 3 worker processes (GitHub\\CLAUDE.md limit is 4; RAM is tight); each run is preceded by
E:\\GitHub\\scripts\\wait_for_cores.ps1.  Usage: python ksrop_check_run.py [--arm NAME] [--fr 1.0] [--ngeo 20] [--only substr] [--limit N]"""
import argparse
import concurrent.futures
import datetime as dt
import math
import os
import re
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd

D = 'E:/GitHub/OREM/scratch_rpe/combined_resonance'
KS = 'E:/GitHub/KSROP'
EXE = 'E:/claude/ksrop_check/bin/driver_KS.exe'
WORK = 'E:/claude/ksrop_check/work'
WAIT = 'E:/GitHub/scripts/wait_for_cores.ps1'
MU = 3.986004415e5
MAX_WORKERS = 3                      # <= 4 cores (GitHub\CLAUDE.md #1); 3 because the machine has ~1.7 GB free RAM
ISTEP = 360


def jd_to_dt(jd):
    return dt.datetime(1970, 1, 1) + dt.timedelta(days=jd - 2440587.5)


def epoch_to_jd(s):
    t = dt.datetime.strptime(s[:23], '%Y-%m-%dT%H:%M:%S.%f')
    return (t - dt.datetime(1970, 1, 1)).total_seconds() / 86400.0 + 2440587.5


def kep2cart(a, e, inc, raan, argp, M):
    i, O, w, M = (math.radians(x) for x in (inc, raan, argp, M))
    E = M
    for _ in range(60):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    xp, yp = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
    r = a * (1 - e * math.cos(E))
    n = math.sqrt(MU / a ** 3)
    vxp, vyp = -a * n * math.sin(E) / (1 - e * math.cos(E)), a * n * math.sqrt(1 - e * e) * math.cos(E) / (1 - e * math.cos(E))
    cO, sO, ci, si, cw, sw = math.cos(O), math.sin(O), math.cos(i), math.sin(i), math.cos(w), math.sin(w)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ np.array([xp, yp, 0.0]), R @ np.array([vxp, vyp, 0.0]), r


def write_inputs(wd, c, ngeo, fr, nrev):
    os.makedirs(f'{wd}/input', exist_ok=True); os.makedirs(f'{wd}/output', exist_ok=True)
    for f in ('ATM.DAT', 'EGM2008_to2190_TideFree'):
        if not os.path.exists(f'{wd}/input/{f}'):
            if f == 'ATM.DAT':
                shutil.copy(f'{KS}/input/ATM.DAT', f'{wd}/input/ATM.DAT')
            else:                                    # the driver streams only the first n(n+1)/2 rows: a 100-degree head is plenty
                with open(f'{KS}/input/{f}') as src, open(f'{wd}/input/{f}', 'w') as dst:
                    for _ in range(5200):
                        dst.write(src.readline())
    with open(f'{wd}/input/const_new.dat', 'w') as f:
        f.write('3.986004415D5 6378.1363D0 1.495978707d08 1.32712440018d11 4.902801076d3\n')
        f.write(f'{ngeo} 2 3\n4.56d-6\n')
    with open(f'{wd}/input/input.dat', 'w') as f:
        f.write(f'{nrev} {ISTEP} 1d-15\n1 1 1\n{c["BN_kg_m2"]:.10g} 1 7.2921150d-5 3.35281066d-3 {fr}\n1.2 0.01 0 1\n')
    x, v, _ = kep2cart(c['a'], c['e'], c['inc'], c['raan'], c['argp'], c['M'])
    ep = jd_to_dt(c['jd0']).strftime('%Y-%m-%dT%H:%M:%S.%f')[:23]
    with open(f'{wd}/input/input.opm', 'w') as f:
        f.write(f'CCSDS_OPM_VERS = 2.0\nCREATION_DATE  = {ep}\nORIGINATOR     = KSROP\n\nMETA_START\nOBJECT_NAME    = SATELLITE\nOBJECT_ID      = UNKNOWN\n'
                f'CENTER_NAME    = EARTH\nREF_FRAME      = EME2000\nTIME_SYSTEM    = UTC\nMETA_STOP\n\nSTATE_VECTOR\nEPOCH          = {ep}\n')
        for lab, val, u in (('X', x[0], 'km'), ('Y', x[1], 'km'), ('Z', x[2], 'km'), ('X_DOT', v[0], 'km/s'), ('Y_DOT', v[1], 'km/s'), ('Z_DOT', v[2], 'km/s')):
            f.write(f'{lab:<14} = {val:16.9f} [{u}]\n')


def run_case(args):
    c, ngeo, fr, years_margin = args
    wd = f'{WORK}/{os.getpid()}'
    while subprocess.run(['powershell', '-NoProfile', '-File', WAIT], capture_output=True).returncode != 0:
        time.sleep(60)
    cap_years = min(float(c['years_cap']), c['real_horizon_d'] / 365.25 + years_margin)
    P_days = 2 * math.pi * math.sqrt(c['a'] ** 3 / MU) / 86400.0
    nrev = int(cap_years * 365.25 / P_days) + 2
    write_inputs(wd, c, ngeo, fr, nrev)
    for pat in os.listdir(f'{wd}/output'):
        os.remove(f'{wd}/output/{pat}')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    t0 = time.time()
    try:
        r = subprocess.run([EXE], cwd=wd, capture_output=True, text=True, timeout=6 * 3600, env=env)
        out, rc = r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired:
        out, rc = 'TIMEOUT', -9
    wall = time.time() - t0
    status, end_jd = 'no_decay', float('nan')
    m = re.search(r'Re-entry occurred: altitude =\s*([-\d.]+) km \(< 80 km\) at epoch (\S+)', out)
    d = re.search(r'last valid epoch was (\S+)', out)
    if m:
        status, end_jd = 'decay', epoch_to_jd(m.group(2))
    elif d:
        status, end_jd = 'diverged', epoch_to_jd(d.group(1))
    elif rc != 0:
        status = f'error_rc{rc}'
    for pat in os.listdir(f'{wd}/output'):
        os.remove(f'{wd}/output/{pat}')
    return dict(case=c['case'], status=status, ksrop_end_jd=end_jd, nrev=nrev, wall_s=round(wall, 1), ngeo=ngeo, fr=fr, tail=out.strip()[-160:].replace('\n', ' | '))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', default='full'); ap.add_argument('--fr', type=float, default=1.0); ap.add_argument('--ngeo', type=int, default=20)
    ap.add_argument('--only', default=None); ap.add_argument('--limit', type=int, default=None); ap.add_argument('--margin', type=float, default=4.0)
    a = ap.parse_args()
    cases = pd.read_csv(f'{D}/ksrop_check_cases.csv')
    res_path = f'{D}/ksrop_check_results_{a.arm}.csv'
    done = set(pd.read_csv(res_path)['case']) if os.path.exists(res_path) else set()
    todo = [r for _, r in cases.iterrows() if r['case'] not in done and (a.only is None or a.only in r['case'])]
    if a.limit:
        todo = todo[:a.limit]
    print(f'{len(todo)} cases to run ({len(done)} done), arm={a.arm} ngeo={a.ngeo} fr={a.fr}', flush=True)
    os.makedirs(WORK, exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futs = [pool.submit(run_case, (r.to_dict(), a.ngeo, a.fr, a.margin)) for r in todo]
        for fu in concurrent.futures.as_completed(futs):
            row = fu.result()
            pd.DataFrame([row]).to_csv(res_path, mode='a', header=not os.path.exists(res_path), index=False)
            print(row['case'], row['status'], f"{row['wall_s']} s", row['tail'][-80:], flush=True)


if __name__ == '__main__':
    main()
