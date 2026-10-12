"""Re-run the issue #56 Python analysis pipeline in dependency order after the 3-line TLE loader fix (THEORY section 20).
Each step runs after the CPU-availability check (GitHub/CLAUDE.md section 14) with BLAS threads pinned to 1, one at a time.
Fortran-backed steps (OREM in midlife_test, the KSROP check, the restricted OREM campaigns) are NOT run here; the driver
stops before them so their inputs can be compared first. Log: E:/claude/runs/regen/<step>.log."""
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path('E:/GitHub/OREM/scratch_rpe')
CR = ROOT / 'combined_resonance'
LOG = Path('E:/claude/runs/regen'); LOG.mkdir(parents=True, exist_ok=True)
WAIT = 'E:/GitHub/scripts/wait_for_cores.ps1'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
           MPLBACKEND='Agg', PYTHONIOENCODING='utf-8')

STEPS = [
    (ROOT, ['resonance_campaign_sweep.py']),
    (ROOT, ['resonance_sweep_new_candidates.py']),
    (ROOT, ['make_postresonance_tles.py']),
    (ROOT / 'resonance_plots', ['prepare_plot_data.py']),
    (ROOT / 'resonance_plots', ['prepare_plot_data_fullpool.py']),
    (ROOT, ['validate_resonance_detector.py']),
    (CR, ['molniya_cycle.py']),
    (CR, ['validate_drive.py']),
    (CR, ['robust_stats.py']),
    (CR, ['validate_strength.py']),
    (CR, ['identify_drivers.py']),
    (CR, ['check_apogee_perigee.py']),
    (CR, ['trend_all28.py']),
    (CR, ['trend_phase.py']),
    (CR, ['forecast_test.py']),
    (CR, ['forecast_pool.py']),
    (CR, ['make_figures.py']),
    (CR, ['make_fig3.py']),
    (CR, ['cross_check_orem.py', '300']),
    (CR, ['cross_check_orem.py', '200']),
    (CR, ['cross_check_orem.py', '150']),
    (CR, ['make_fig4.py']),
    (CR, ['midlife_test.py', 'prep']),          # cutoffs only; the OREM step is not run here
]
only = sys.argv[1:]
for cwd, args in STEPS:
    name = '_'.join(a.replace('.py', '') for a in args)
    if only and name not in only:
        continue
    subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', WAIT], capture_output=True)
    t = time.time()
    with open(LOG / f'{name}.log', 'w', encoding='utf-8') as f:
        rc = subprocess.run([sys.executable, '-u'] + args, cwd=cwd, env=ENV, stdout=f, stderr=subprocess.STDOUT).returncode
    print(f'{name}: rc={rc} {time.time() - t:.0f}s', flush=True)
    if rc != 0:
        print('STOPPING: step failed, see', LOG / f'{name}.log', flush=True)
        sys.exit(rc)
print('PHASE 1 DONE', flush=True)
