"""Merge the Scilab-generated title/index/object pages for the FULL
253-object GTO/HEO pool into one PDF."""
import csv
from pypdf import PdfWriter

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots'
PAGES = f'{BASE}/pages_fullpool'

with open(f'{BASE}/manifest_fullpool.csv') as f:
    manifest = list(csv.DictReader(f))

npages_index = 11  # from plot_frontmatter_fullpool.sce's own run log

writer = PdfWriter()
writer.append(f'{PAGES}/_00_title.pdf')
for p in range(1, npages_index + 1):
    writer.append(f'{PAGES}/_01_index_p{p}.pdf')
for row in manifest:
    writer.append(f'{PAGES}/obj_{row["norad"]}.pdf')

out_path = f'{BASE}/OREM_solar_apsidal_resonance_plots_FULLPOOL.pdf'
with open(out_path, 'wb') as f:
    writer.write(f)

print(f'Merged {1 + npages_index + len(manifest)} pages -> {out_path}')
