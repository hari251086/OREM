"""Split the full-pool PDF (48.7MB, over the 30MB upload limit) into
several smaller parts, each comfortably under the limit."""
from pypdf import PdfReader, PdfWriter

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots'
SRC = f'{BASE}/OREM_solar_apsidal_resonance_plots_FULLPOOL.pdf'
MAX_BYTES = 25 * 1024 * 1024  # comfortable margin under the 30MB limit

reader = PdfReader(SRC)
n = len(reader.pages)

part = 1
writer = PdfWriter()
start_page = 0
parts_written = []

for i in range(n):
    writer.add_page(reader.pages[i])
    # Periodically check size; writing to a buffer each time is slow for
    # 265 pages, so only check every 10 pages (and always on the last).
    if (i - start_page + 1) % 5 == 0 or i == n - 1:
        import io
        buf = io.BytesIO()
        writer.write(buf)
        size = buf.tell()
        if size > MAX_BYTES or i == n - 1:
            out_path = f'{BASE}/OREM_solar_apsidal_resonance_plots_FULLPOOL_part{part}.pdf'
            with open(out_path, 'wb') as f:
                f.write(buf.getvalue())
            print(f'part{part}: pages {start_page+1}-{i+1} ({i-start_page+1} pages, '
                  f'{size/1024/1024:.1f} MB) -> {out_path}')
            parts_written.append(out_path)
            part += 1
            writer = PdfWriter()
            start_page = i + 1

print(f'\n{len(parts_written)} parts written, {n} pages total')
