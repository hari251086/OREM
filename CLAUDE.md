# CLAUDE.md

Cross-repo rules live in `GitHub\CLAUDE.md` — not repeated here, only
what's specific to OREM. Fortran-specific gotchas are in the shared,
tree-walked `.claude/rules/fortran-ks-gotchas.md`.

## Project

Re-entry prediction tool: zone-based orbit propagation (KSROP, consumed
as a git+tag fpm dependency, not vendored) + Response Surface Method +
Genetic Algorithm to estimate ballistic number and re-entry epoch from
TLE history. Feeds OREM-Watchlist (a sibling repo) as its prediction
engine, invoked as a subprocess.

## Build / Test

```bash
fpm build --compiler ifx --flag "-heap-arrays"
fpm test  --compiler ifx --flag "-heap-arrays"
fpm run orem --compiler ifx --flag "-heap-arrays" -- input/orem_42928.cfg
```
`-heap-arrays` (16 MB stack) is **required** for every target linking
`rsm.F`/`ga.F` — the `surfaces(5000,3,3)` arrays overflow the default
stack without it. On Linux: `-heap-arrays` + `ulimit -s unlimited`
(stack size is a shell/OS setting there, not a linker flag).
CI (`.github/workflows/ci.yml`) matrix-tests `ifx` and `gfortran`.

## Key code

- `src/` — RSM/GA zone-refinement logic (`rsm.F`, `ga.F`)
- KSROP propagator (`Subrouts.F`, `Legendre.F`, `TLEread.F`) is a git
  dependency (see `fpm.toml`), not vendored locally
- `input/*.cfg` — run config (TLE path, NORAD ID, zone/GA/force-model
  params — see README §4 for the exact field-by-field format)
- `output/OREM_<NORAD>_<DATE>.txt` — per-zone RPE report

## Always / never

- Never drop the `-heap-arrays` flag when building/testing/running any
  target that touches `rsm.F`/`ga.F` — it will overflow the stack.
- GA population size is fixed at `pop=20` in the shipped config format
  (v1.15) — don't change it casually without checking why that value was
  pinned (see README §4's config-file comments).
