"""Merge just the 27 detector-flagged objects' pages (from the existing
pages_fullpool/obj_*.pdf, no re-render) into one small PDF for a second
manual visual review pass -- issue #56 follow-up after the full
253-object review. Ordered by the review-category verdict already
posted to the issue, not manifest/NORAD order, so re-scanning groups
like with like. A companion text index (page number -> category/object)
is written alongside since the pages themselves carry no divider (no
reportlab/Scilab re-render available in this pass)."""
import csv
from pypdf import PdfWriter

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots'
PAGES = f'{BASE}/pages_fullpool'

_CATEGORIES = [
    ("Genuine, strong U-turn", [40943, 59347, 35497, 10724]),
    ("Weaker but plausible", [37151, 27526, 18954]),
    ("Uncertain", [19881]),
    ("False positive (detector artifact)", [
        18983, 13112, 5980, 7096, 60328, 4570, 28140, 44911,
        6294, 6231, 29462, 5941, 26463, 68538, 37819,
    ]),
    ("Excluded independently (policy, not physics)", [
        39501, 44187, 42985, 13913,
    ]),
]

with open(f'{BASE}/manifest_fullpool.csv') as f:
    manifest = {row["norad"]: row for row in csv.DictReader(f)}

n_total = sum(len(ids) for _, ids in _CATEGORIES)
assert n_total == 27, n_total

writer = PdfWriter()
index_lines = []
page_no = 1
for label, norad_ids in _CATEGORIES:
    for norad in norad_ids:
        row = manifest[str(norad)]
        writer.append(f'{PAGES}/obj_{norad}.pdf')
        n_pages = len(writer.pages) - (page_no - 1)
        index_lines.append(
            f'p.{page_no:>3} [{label}]  {norad}  {row["name"]}'
        )
        page_no = len(writer.pages) + 1

out_path = f'{BASE}/OREM_resonance_flagged27_review.pdf'
with open(out_path, 'wb') as f:
    writer.write(f)

index_path = f'{BASE}/OREM_resonance_flagged27_review_index.txt'
with open(index_path, 'w') as f:
    f.write('\n'.join(index_lines) + '\n')

print(f'Merged {len(writer.pages)} pages ({n_total} objects) -> {out_path}')
print(f'Index -> {index_path}')
