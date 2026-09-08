// plot_frontmatter.sce -- title page + index/summary page for the
// resonance-plots PDF compile (OREM issue #56 follow-on). Each line is
// its own xstring() call at an explicit y-position -- Scilab's
// multi-row string-matrix xstring() spaces lines inconsistently with
// the 0-1 data_bounds used here, which caused overlapping text on the
// first attempt.

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots';

mlines = mgetl(BASE + '/manifest.csv');
nobj = size(mlines, 1) - 1;

function draw_lines(lines_list, y0, dy, fsize)
    for k = 1:size(lines_list, 1)
        xstring(0.05, y0 - (k-1)*dy, lines_list(k));
        t = gce();
        t.font_size = fsize;
    end
endfunction

// ---- Title page ----
f = scf();
f.figure_size = [1500, 950];
a = gca();
a.axes_visible = ['off', 'off', 'off'];
a.data_bounds = [0, 0; 1, 1];

xstring(0.05, 0.92, 'Solar Apsidal U-Turn Resonance in Real GTO/HEO Objects');
t = gce();
t.font_size = 4;

body = [ ...
  'A TLE-based reproduction of Wang & Gurfil (2017), Advances in Space Research'; ...
  '59, 2101-2116, applied to 23 real, independently tracked objects.'; ...
  ' '; ...
  'Resonance condition (paper Eq. 12):'; ...
  '   a (1-e^2)^(4/7)  =  lambda_crit(i)'; ...
  '   lambda_crit(i) = [ 3 J2 sqrt(mu) RE^2 (5cos^2 i - 2cos i - 1) / (4 dM_solar/dt) ]^(2/7)'; ...
  ' '; ...
  'Resonance angle plotted (proxy for the papers solar-azimuth U-turn):'; ...
  '   psi(t) = mod( RAAN(t) + argp(t) - L_sun(t), 360 deg )'; ...
  ' '; ...
  'a(t), e(t), i(t), RAAN(t), argp(t) computed directly from each objects real'; ...
  'Space-Track TLE mean elements -- no propagation.'; ...
  ' '; ...
  string(nobj) + ' objects shown: 16 from OREMs original 97-object RPE campaign'; ...
  '+ 7 found after expanding the candidate pool.'; ...
  ' '; ...
  'Source: OREM (github.com/hari251086/OREM), issue #56.' ...
];
draw_lines(body, 0.80, 0.045, 2);

xs2pdf(f.figure_id, BASE + '/pages/_00_title.pdf');
close(f);
disp('wrote title page');

// ---- Index page ----
f = scf();
f.figure_size = [1500, 950];
a = gca();
a.axes_visible = ['off', 'off', 'off'];
a.data_bounds = [0, 0; 1, 1];

xstring(0.05, 0.95, 'Index: 23 objects, in the order plotted');
t = gce();
t.font_size = 3;

xstring(0.05, 0.88, 'NORAD    Name                             i(deg)   Resonance yr   Observed decay yr');
t = gce();
t.font_size = 2;
t.font_style = 8;

rows_out = [];
for k = 1:nobj
    toks = strsplit(mlines(k + 1), ',');
    norad = toks(1);
    name = toks(2);
    i_mean = strtod(toks(4));
    resonance_year = strtod(toks(7));
    decay_year_s = toks(8);
    row = norad;
    while length(row) < 9
        row = row + ' ';
    end
    nm = name;
    while length(nm) < 33
        nm = nm + ' ';
    end
    ideg = string(round(i_mean*100)/100);
    while length(ideg) < 8
        ideg = ideg + ' ';
    end
    ryr = string(round(resonance_year*100)/100);
    while length(ryr) < 14
        ryr = ryr + ' ';
    end
    row = row + nm + ideg + ryr + decay_year_s;
    rows_out = [rows_out; row];
end
draw_lines(rows_out, 0.83, 0.033, 2);

xs2pdf(f.figure_id, BASE + '/pages/_01_index.pdf');
close(f);
disp('wrote index page');
quit
