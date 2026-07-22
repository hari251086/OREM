# OREM — Coding Standards

Rules for code contributed to this repository, independent of any single
issue or version. Architectural rationale lives in `ARCHITECTURE.md`;
build/test mechanics live in `README.md`. This file is for standing rules.

## 1. Resource limits

- **Max 4 cores/threads.** Any parallel or multi-threaded code path — OpenMP
  directives (`/Qopenmp`, `!$OMP`), MPI, spawned subprocesses, or a Python
  `multiprocessing`/`concurrent.futures` pool (e.g. in OREM-Watchlist tooling
  that calls into this repo) — must cap its own concurrency at 4, regardless
  of how many cores the host machine has available. Hardcode the cap or read
  it from a single named constant/config value; do not default to
  `os.cpu_count()` / `OMP_NUM_THREADS` unset / all available cores.
