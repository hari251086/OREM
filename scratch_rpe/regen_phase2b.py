"""Phase 2b of the post-loader-fix regeneration (THEORY section 20): coverage pool, KSROP case table + full re-run of the
KSROP check (old results moved aside so every case is recomputed; the runner keeps <= 3 workers with its own CPU checks),
then the section 16-19 diagnostics. Logs in E:/claude/runs/regen/."""
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

CR = Path('E:/GitHub/OREM/scratch_rpe/combined_resonance')
LOG = Path('E:/claude/runs/regen')
WAIT = 'E:/GitHub/scripts/wait_for_cores.ps1'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
           MPLBACKEND='Agg', PYTHONIOENCODING='utf-8')


def wait_cores():
    subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', WAIT], capture_output=True)


def run(args, name):
    wait_cores(); t = time.time()
    with open(LOG / f'{name}.log', 'w', encoding='utf-8') as f:
        rc = subprocess.run([sys.executable, '-u'] + args, cwd=CR, env=ENV, stdout=f, stderr=subprocess.STDOUT).returncode
    print(f'{name}: rc={rc} {time.time() - t:.0f}s', flush=True)
    if rc != 0:
        print('STOPPING', flush=True); sys.exit(rc)


FROM_KSROP = '--from-ksrop' in sys.argv            # resume: pool and case table already regenerated
if not FROM_KSROP:
    run(['window_coverage_pool.py'], 'window_coverage_pool')
    run(['ksrop_check_cases.py'], 'ksrop_check_cases')
old = CR / 'ksrop_check_results_full.csv'
if old.exists() and not FROM_KSROP:
    shutil.move(str(old), str(LOG / 'ksrop_check_results_full.before.csv'))
run(['ksrop_check_run.py', '--arm', 'full'], 'ksrop_check_run_full')
run(['ksrop_check_analyze.py', 'full'], 'ksrop_check_analyze_full')

# section 16 diagnostics: 4 shards in parallel, then merge
procs = []
for k in range(4):
    wait_cores()
    f = open(LOG / f'diag_nodate_part{k}.log', 'w', encoding='utf-8')
    procs.append(subprocess.Popen([sys.executable, '-u', 'diag_nodate.py', str(k), '4'], cwd=CR, env=ENV,
                                  stdout=f, stderr=subprocess.STDOUT))
    time.sleep(5)
for p in procs:
    p.wait()
print('diag_nodate shards rc', [p.returncode for p in procs], flush=True)
import pandas as pd  # noqa: E402
parts = [CR / f'diag_nodate_part{k}.csv' for k in range(4)]
pd.concat([pd.read_csv(p) for p in parts]).sort_values(['norad', 'lead_yr']).to_csv(CR / 'diag_nodate.csv', index=False)
for p in parts:
    p.unlink()
run(['diag_nodate_bscan.py'], 'diag_nodate_bscan')
run(['diag_nodate_orembn.py'], 'diag_nodate_orembn')
run(['diag_mean_osc.py'], 'diag_mean_osc')
print('PHASE 2B DONE', flush=True)
