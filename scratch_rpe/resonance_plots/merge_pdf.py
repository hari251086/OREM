"""Merge the Scilab-generated title/index/object pages into one PDF."""
import csv
from pypdf import PdfWriter

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots'
PAGES = f'{BASE}/pages'

with open(f'{BASE}/manifest.csv') as f:
    manifest = list(csv.DictReader(f))

writer = PdfWriter()
writer.append(f'{PAGES}/_00_title.pdf')
writer.append(f'{PAGES}/_01_index.pdf')
for row in manifest:
    writer.append(f'{PAGES}/obj_{row["norad"]}.pdf')

out_path = f'{BASE}/OREM_solar_apsidal_resonance_plots.pdf'
with open(out_path, 'wb') as f:
    writer.write(f)

print(f'Merged {2 + len(manifest)} pages -> {out_path}')
