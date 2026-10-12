"""Phase 2a of the post-loader-fix regeneration (THEORY section 20): the 17 mid-life OREM runs at the moved cutoffs
(<= 4 concurrent, CPU check before each launch), then the dependent Python steps in order. Logs in E:/claude/runs/regen/."""
import os
import subprocess
import sys
import time
from pathlib import Path

CR = Path('E:/GitHub/OREM/scratch_rpe/combined_resonance')
LOG = Path('E:/claude/runs/regen'); LOG.mkdir(parents=True, exist_ok=True)
WAIT = 'E:/GitHub/scripts/wait_for_cores.ps1'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
           MPLBACKEND='Agg', PYTHONIOENCODING='utf-8')
sys.path.insert(0, str(CR))


def wait_cores():
    subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', WAIT], capture_output=True)


def run(args, name):
    wait_cores(); t = time.time()
    with open(LOG / f'{name}.log', 'w', encoding='utf-8') as f:
        rc = subprocess.run([sys.executable, '-u'] + args, cwd=CR, env=ENV, stdout=f, stderr=subprocess.STDOUT).returncode
    print(f'{name}: rc={rc} {time.time() - t:.0f}s', flush=True)
    if rc != 0:
        print('STOPPING', flush=True); sys.exit(rc)


# 1. mid-life OREM runs, <= 4 at a time
import midlife_test  # noqa: E402
procs = []
for n in midlife_test.ALL:
    while sum(p.poll() is None for p in procs) >= 4:
        time.sleep(10)
    wait_cores()
    f = open(LOG / f'midlife_orem_{n}.log', 'w', encoding='utf-8')
    procs.append(subprocess.Popen([sys.executable, '-u', 'midlife_test.py', 'orem', str(n)], cwd=CR, env=ENV,
                                  stdout=f, stderr=subprocess.STDOUT))
    print('launched OREM', n, flush=True)
    time.sleep(5)
for p in procs:
    p.wait()
print('mid-life OREM runs done; return codes', [p.returncode for p in procs], flush=True)

# 2. dependent Python steps
for args in (['midlife_test.py', 'forecast'], ['midlife_test.py', 'summary'], ['midlife_extra.py'], ['scan_second_hump.py'],
             ['drag_test.py'], ['make_fig5.py'], ['window_coverage_test.py'], ['perigee_descent_rates.py']):
    run(args, '_'.join(a.replace('.py', '') for a in args))
print('PHASE 2A DONE', flush=True)
