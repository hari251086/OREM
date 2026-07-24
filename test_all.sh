#!/bin/bash
# test_all.sh -- issue #22: compile every OREM test suite with Intel
# ifx (Linux) and run them all. Used both locally and by
# .github/workflows/ci.yml. Exits nonzero on the first compile or
# test failure (set -e).
#
# ifx, not gfortran: ksrop/propagate_ks.F relies on several ifx-
# specific tolerances (implicit-real array dimensions, a function/
# array name collision on `R`, non-standard function-result syntax)
# that gfortran 13 rejects outright. Making that file portable is
# real, separate surgery on the most validated file in the repo --
# tracked as its own issue (#28), not bundled into CI setup. CI tests
# what the project actually builds and validates with.
#
# input/example_multi.tle.txt (94597-entry catalog, 13 MB) is
# deliberately gitignored -- test_tle_evolution.F's T42-T50 skip
# gracefully when it's absent (e.g. in CI); run this script with the
# full catalog present locally for that additional coverage.
set -e

ulimit -s unlimited 2>/dev/null || true

echo "=== Compiling ==="
ifx -heap-arrays test/test_propagate_ks.F src/propagate_ks.F src/Subrouts.F src/Legendre.F -o test_propagate_ks.exe
ifx -heap-arrays test/test_tle_evolution.F src/tle_evolution.F src/TLEread.F src/Subrouts.F src/Legendre.F -o test_tle_evolution.exe
ifx -heap-arrays test/test_zone_select.F src/zone_select.F src/tle_evolution.F src/TLEread.F src/Subrouts.F src/Legendre.F -o test_zone_select.exe
ifx -heap-arrays test/test_ga.F src/ga.F src/rsm.F src/tle_evolution.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F -o test_ga.exe
ifx -heap-arrays test/test_rsm.F src/rsm.F src/tle_evolution.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F src/ga.F -o test_rsm.exe
ifx -heap-arrays test/test_orem.F src/orem.F src/report.F src/rsm.F src/ga.F src/tle_evolution.F src/tle_filter.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F -o test_orem.exe
ifx -heap-arrays test/test_reentry.F src/orem.F src/rsm.F src/ga.F src/tle_evolution.F src/tle_filter.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F -o test_reentry.exe
ifx -heap-arrays test/test_e2e.F src/orem.F src/rsm.F src/ga.F src/tle_evolution.F src/tle_filter.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F -o test_e2e.exe
ifx -heap-arrays test/test_gmat.F src/rsm.F src/tle_evolution.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/Legendre.F src/TLEread.F src/ga.F -o test_gmat.exe
ifx -heap-arrays test/test_sw.F src/swx.F src/orem.F src/report.F src/rsm.F src/ga.F src/tle_evolution.F src/tle_filter.F src/zone_select.F src/propagate_ks.F src/Subrouts.F src/TLEread.F src/Legendre.F -o test_sw.exe
ifx -heap-arrays test/test_tle_filter.F src/tle_filter.F src/zone_select.F src/tle_evolution.F src/TLEread.F src/Subrouts.F src/Legendre.F -o test_tle_filter.exe

echo "=== Running ==="
fail=0
for exe in test_propagate_ks test_tle_evolution test_zone_select test_ga \
           test_rsm test_orem test_reentry test_e2e test_gmat test_sw \
           test_tle_filter; do
   echo "--- $exe ---"
   if ! "./${exe}.exe"; then
      echo "*** $exe FAILED ***"
      fail=1
   fi
done

if [ "$fail" -ne 0 ]; then
   echo "=== SOME SUITES FAILED ==="
   exit 1
fi

echo "=== ALL SUITES PASSED ==="
