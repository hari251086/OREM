"""Compare regenerated CSV outputs with their committed (HEAD) versions (issue #56, THEORY section 20).
Numeric columns: changed if |new - old| > atol + rtol*|old| (NaN == NaN). Other columns: exact string comparison.
Usage: regen_diff.py <path relative to repo root> [key column] ..."""
import io
import subprocess
import sys

import numpy as np
import pandas as pd

REPO = 'E:/GitHub/OREM'


def old_version(rel):
    out = subprocess.run(['git', 'show', f'HEAD:{rel}'], cwd=REPO, capture_output=True, text=True, encoding='utf-8')
    return pd.read_csv(io.StringIO(out.stdout)) if out.returncode == 0 else None


def compare(rel, key=None, rtol=1e-6, atol=1e-6, show=6):
    new = pd.read_csv(f'{REPO}/{rel}')
    old = old_version(rel)
    if old is None:
        print(f'{rel}: no committed version'); return
    print(f'== {rel}: rows old {len(old)} new {len(new)}')
    if key is None or key not in new.columns:
        if len(old) != len(new):
            print('   row count differs and no key column; cannot align'); return
        m = old.add_suffix('_old').join(new.add_suffix('_new')); kcols = []
    else:
        keys = [key] + ([c for c in ('lead_yr', 'case') if c in new.columns and c != key])
        m = old.merge(new, on=keys, how='outer', suffixes=('_old', '_new'), indicator=True)
        only = m['_merge'].value_counts().to_dict()
        if only.get('left_only') or only.get('right_only'):
            print(f"   rows only in old: {only.get('left_only', 0)}, only in new: {only.get('right_only', 0)}")
        m = m[m['_merge'] == 'both']; kcols = keys
    for c in [c for c in new.columns if c not in kcols and f'{c}_old' in m and f'{c}_new' in m]:
        a, b = m[f'{c}_old'], m[f'{c}_new']
        an, bn = pd.to_numeric(a, errors='coerce'), pd.to_numeric(b, errors='coerce')
        if an.notna().sum() + bn.notna().sum() > 0 and (an.notna() == a.notna()).all() and (bn.notna() == b.notna()).all():
            both_nan = an.isna() & bn.isna()
            diff = ~both_nan & ((an.isna() != bn.isna()) | ((an - bn).abs() > atol + rtol * an.abs()))
        else:
            diff = (a.astype(str) != b.astype(str)) & ~(a.isna() & b.isna())
        n = int(diff.sum())
        if n:
            ex = m.loc[diff, kcols + [f'{c}_old', f'{c}_new']].head(show)
            print(f'   {c}: {n} changed  e.g. ' + '; '.join(
                f"{' '.join(str(r[k]) for k in kcols)}: {r[f'{c}_old']} -> {r[f'{c}_new']}" for _, r in ex.iterrows()))


if __name__ == '__main__':
    args = sys.argv[1:]
    for i in range(0, len(args), 2):
        compare(args[i], args[i + 1] if i + 1 < len(args) else None)
