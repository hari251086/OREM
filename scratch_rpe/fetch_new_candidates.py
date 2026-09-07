"""One-time historical backfill of TLE data for the new GTO/HEO candidates
found by screen_new_candidates.py, beyond OREM's existing 97-object
campaign. Two sources, per the Space-Track usage policy's "store and
reuse local data" rule (never re-fetch what's already on disk):

1. 58 candidates already have a full history cached by OREM-Watchlist
   (`objects/<norad>/<norad>.tle.txt`, 3LE format) -- converted in
   place to OREM's plain 2-line TLE format, zero new Space-Track
   requests.
2. 122 candidates need a first-ever fetch. Uses the same shared,
   rate-limited client as every other script in this lineage
   (`heowatch.spacetrack_client.get_client()`, never a second client
   instance) and batches many NORAD IDs into ONE `gp_history` query
   each (the `spacetrack` library joins a Python list into a
   comma-separated predicate automatically) per
   SPACE_TRACK_USAGE_POLICY.md �4's "batch multiple objects into one
   comma-delimited query rather than looping per-object requests" --
   NOT the fetch_replacement_tles.py precedent's per-object loop,
   which predates that explicit requirement. `gp_history` is each
   object's full lifetime history requested once (never re-requested),
   matching the "1/lifetime per range" data class rule.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"C:\Users\hari2\OneDrive\Documents\GitHub\OREM-Watchlist\src")
from heowatch.spacetrack_client import get_client  # noqa: E402

DATA_ROOT = Path("E:/Research/1. R&D/Re-entry/2026/Data")
OBJECTS_DIR = DATA_ROOT / "objects"
INPUT_DIR = Path("E:/GitHub/OREM/input")
CHUNK_SIZE = 40


def to_plain_2line(text):
    """Strip any '0 <name>' line, keep only '1 ...'/'2 ...' pairs."""
    lines = [l for l in text.splitlines() if l.startswith('1 ') or l.startswith('2 ')]
    return '\n'.join(lines) + ('\n' if lines else '')


def convert_cached(candidates):
    written = []
    for c in candidates:
        norad = c['norad']
        src = OBJECTS_DIR / str(norad) / f"{norad}.tle.txt"
        text = src.read_text(encoding='utf-8', errors='replace')
        plain = to_plain_2line(text)
        n_tle = plain.count('\n1 ') + (1 if plain.startswith('1 ') else 0)
        dst = INPUT_DIR / f"example_{norad}.tle.txt"
        dst.write_text(plain, encoding='utf-8')
        written.append((norad, c['name'], n_tle, 'cache'))
        print(f"  {norad:>7} {c['name'][:30]:<30} {n_tle:5d} TLEs (from local cache)")
    return written


def fetch_batch(candidates):
    st = get_client()
    written = []
    norads = [c['norad'] for c in candidates]
    by_norad = {c['norad']: c for c in candidates}

    for i in range(0, len(norads), CHUNK_SIZE):
        chunk = norads[i:i + CHUNK_SIZE]
        print(f"\nfetching gp_history for {len(chunk)} objects "
              f"({i+1}-{i+len(chunk)} of {len(norads)})...")
        tle_text = st.gp_history(
            norad_cat_id=chunk, orderby="norad_cat_id,epoch", format="tle"
        )
        lines = tle_text.splitlines()
        # Group consecutive '1 '/'2 ' pairs by the NORAD ID embedded in
        # each line 1 (cols 3-7), since the response interleaves objects
        # in norad_cat_id order but doesn't repeat a header per object.
        per_object = {}
        k = 0
        while k < len(lines) - 1:
            l1, l2 = lines[k], lines[k + 1]
            if l1.startswith('1 ') and l2.startswith('2 '):
                norad = int(l1[2:7])
                per_object.setdefault(norad, []).append((l1, l2))
                k += 2
            else:
                k += 1

        for norad in chunk:
            pairs = per_object.get(norad, [])
            c = by_norad[norad]
            dst = INPUT_DIR / f"example_{norad}.tle.txt"
            if pairs:
                with open(dst, 'w', encoding='utf-8') as f:
                    for l1, l2 in pairs:
                        f.write(l1 + '\n')
                        f.write(l2 + '\n')
                written.append((norad, c['name'], len(pairs), 'fetched'))
                print(f"  {norad:>7} {c['name'][:30]:<30} {len(pairs):5d} TLEs (fetched)")
            else:
                written.append((norad, c['name'], 0, 'EMPTY'))
                print(f"  {norad:>7} {c['name'][:30]:<30}     0 TLEs -- no gp_history data returned")

        time.sleep(1)  # extra margin beyond the library's own 30/min throttle

    return written


def main():
    with open('scratch_rpe/new_candidates.json') as f:
        data = json.load(f)

    print(f"=== Converting {len(data['have_cache'])} cached objects ===")
    w1 = convert_cached(data['have_cache'])

    print(f"\n=== Fetching {len(data['need_fetch'])} new objects "
          f"in batches of {CHUNK_SIZE} ===")
    w2 = fetch_batch(data['need_fetch'])

    all_written = w1 + w2
    import csv
    with open('scratch_rpe/new_candidates_fetch_log.csv', 'w', newline='') as f:
        cw = csv.writer(f)
        cw.writerow(['norad', 'name', 'n_tle', 'source'])
        cw.writerows(all_written)

    n_ok = sum(1 for r in all_written if r[2] >= 3)
    n_empty = sum(1 for r in all_written if r[2] == 0)
    print(f"\n{len(all_written)} total: {n_ok} with >=3 TLEs (usable), "
          f"{n_empty} empty (no data returned -- likely too-old catalog "
          f"objects not in Space-Track's history for this NORAD range)")
    print("Written: scratch_rpe/new_candidates_fetch_log.csv")


if __name__ == '__main__':
    main()
