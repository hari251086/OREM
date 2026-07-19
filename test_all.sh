#!/bin/bash
# test_all.sh -- issue #22: compile every OREM test suite with gfortran
# and run them all. Used both locally (Unix/gfortran) and by
# .github/workflows/ci.yml. Exits nonzero on the first compile or
# test failure (set -e).
#
# input/example_multi.tle.txt (94597-entry catalog, 13 MB) is
# deliberately gitignored -- test_tle_evolution.F's T42-T50 skip
# gracefully when it's absent (e.g. in CI); run this script with the
# full catalog present locally for that additional coverage.
set -e

ulimit -s unlimited 2>/dev/null || true

echo "=== Compiling ==="
gfortran test_propagate_ks.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F -o test_propagate_ks.exe
gfortran test_tle_evolution.F tle_evolution.F ksrop/TLEread.F ksrop/Subrouts.F ksrop/Legendre.F -o test_tle_evolution.exe
gfortran test_zone_select.F zone_select.F tle_evolution.F ksrop/TLEread.F ksrop/Subrouts.F ksrop/Legendre.F -o test_zone_select.exe
gfortran test_ga.F ga.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F -o test_ga.exe
gfortran test_rsm.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F ga.F -o test_rsm.exe
gfortran test_orem.F orem.F report.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F -o test_orem.exe
gfortran test_reentry.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F -o test_reentry.exe
gfortran test_e2e.F orem.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F -o test_e2e.exe
gfortran test_gmat.F rsm.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/Legendre.F ksrop/TLEread.F ga.F -o test_gmat.exe
gfortran test_sw.F swx.F orem.F report.F rsm.F ga.F tle_evolution.F zone_select.F ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F -o test_sw.exe

echo "=== Running ==="
fail=0
for exe in test_propagate_ks test_tle_evolution test_zone_select test_ga \
           test_rsm test_orem test_reentry test_e2e test_gmat test_sw; do
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
