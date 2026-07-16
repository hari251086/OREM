"""Parse gmat_xval_42928z0.script's ReportFile output into the Fortran
DATA blocks used by test_gmat.F (ref_apo0, ref_drop).

Usage: python parse_gmat_xval.py [path-to-gmat_xval_42928z0.txt]
Default path matches GMAT's default output directory for this install.
"""
import sys

path = sys.argv[1] if len(sys.argv) > 1 else \
    r"E:\Softwares\gmat-win-R2026a\output\gmat_xval_42928z0.txt"

rows = {}
with open(path) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 5 or not parts[0].isdigit():
            continue
        run_id, elapsed, radapo = int(parts[0]), float(parts[1]), float(parts[2])
        rows.setdefault(run_id, {})[round(elapsed / 604800)] = radapo

e_idx = [1, 2, 3]
bn_idx = [1, 2, 3]

print("      double precision ref_drop(3,3)")
print("      data ref_drop /")
for bi in bn_idx:
    vals = []
    for ei in e_idx:
        run = rows[ei * 10 + bi]
        drop = run[1] - run[0]
        vals.append(drop)
    tail = "," if bi < 3 else " /"
    label = {1: "BN=80", 2: "BN=120", 3: "BN=160"}[bi]
    print("     &  {:.4f}d0, {:.4f}d0, {:.4f}d0{}   ! e1,e2,e3: {}"
          .format(vals[0], vals[1], vals[2], tail, label))

print()
print("      double precision ref_apo0(3)")
apo0 = [rows[ei * 10 + 1][0] for ei in e_idx]
print("      data ref_apo0 / {:.4f}d0, {:.4f}d0, {:.4f}d0 /"
      .format(*apo0))
