// plot_frontmatter_fullpool.sce -- title page + paginated index for the
// full 253-object GTO/HEO pool PDF compile.

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots';

mlines = mgetl(BASE + '/manifest_fullpool.csv');
nobj = size(mlines, 1) - 1;

function draw_lines(lines_list, y0, dy, fsize)
    for k = 1:size(lines_list, 1)
        xstring(0.05, y0 - (k-1)*dy, lines_list(k));
        t = gce();
        t.font_size = fsize;
    end
endfunction

n_resonant = 0;
for k = 1:nobj
    toks = strsplit(mlines(k + 1), ',');
    if strtod(toks(8)) == 1 then
        n_resonant = n_resonant + 1;
    end
end

// ---- Title page ----
f = scf();
f.figure_size = [1500, 950];
a = gca();
a.axes_visible = ['off', 'off', 'off'];
a.data_bounds = [0, 0; 1, 1];

xstring(0.05, 0.92, 'Solar Apsidal Resonance Screen: Full GTO/HEO Candidate Pool');
t = gce(); t.font_size = 4;

body = [ ...
  'A TLE-based reproduction of Wang & Gurfil (2017), Advances in Space Research'; ...
  '59, 2101-2116, applied to every object in OREMs GTO/HEO candidate pool --'; ...
  'not only the ones an automated Eq.12 crossing check flagged as resonant.'; ...
  ' '; ...
  'This PDF is for MANUAL review: each page shows the same 6 panels regardless'; ...
  'of the automated DETECTOR verdict (shown in the info panel), so a genuine'; ...
  'U-turn the detector missed -- or a flagged crossing that is not visually a'; ...
  'real U-turn -- can be judged directly from the plots.'; ...
  ' '; ...
  'Resonance condition (paper Eq. 12):  a(1-e^2)^(4/7) = lambda_crit(i)'; ...
  'Resonance angle plotted:  psi(t) = mod( RAAN(t) + argp(t) - L_sun(t), 360 deg )'; ...
  ' '; ...
  string(nobj) + ' objects plotted (' + string(n_resonant) + ' auto-flagged RESONANT, ' + ...
     string(nobj - n_resonant) + ' not) out of 277 screened -- 24 excluded for having'; ...
  'fewer than 3 tracked TLEs (a trajectory cannot be drawn from 0-2 points).'; ...
  ' '; ...
  'a(t), e(t), i(t), RAAN(t), argp(t) computed directly from each objects real'; ...
  'Space-Track TLE mean elements -- no propagation. Long records (>800 TLEs)'; ...
  'are uniformly decimated to <=800 plotted points for file size/speed.'; ...
  ' '; ...
  'Source: OREM (github.com/hari251086/OREM), issue #56.' ...
];
draw_lines(body, 0.80, 0.040, 2);

xs2pdf(f.figure_id, BASE + '/pages_fullpool/_00_title.pdf');
close(f);
disp('wrote title page');

// ---- Paginated index ----
rows_per_page = 24;
npages = ceil(nobj / rows_per_page);

for p = 1:npages
    f = scf();
    f.figure_size = [1500, 950];
    a = gca();
    a.axes_visible = ['off', 'off', 'off'];
    a.data_bounds = [0, 0; 1, 1];

    xstring(0.05, 0.95, 'Index (' + string(p) + '/' + string(npages) + '): NORAD, name, pool, i, detector verdict, resonance yr, decay yr');
    t = gce(); t.font_size = 2;

    xstring(0.05, 0.895, 'NORAD    Name                        Pool          i(deg)  Verdict    Res.yr   Decay yr');
    t = gce(); t.font_size = 1; t.font_style = 8;

    i0 = (p - 1) * rows_per_page + 1;
    i1 = min(p * rows_per_page, nobj);
    rows_out = [];
    for k = i0:i1
        toks = strsplit(mlines(k + 1), ',');
        norad = toks(1);
        name = toks(2);
        pool = toks(3);
        i_mean = strtod(toks(6));
        is_res = strtod(toks(8));
        res_yr_s = toks(9);
        decay_yr = strtod(toks(11));

        row = norad;
        while length(row) < 9
            row = row + ' ';
        end
        nm = name;
        while length(nm) < 29
            nm = nm + ' ';
        end
        pl = pool;
        while length(pl) < 14
            pl = pl + ' ';
        end
        ideg = string(round(i_mean*100)/100);
        while length(ideg) < 8
            ideg = ideg + ' ';
        end
        verdict = 'RESONANT';
        if is_res <> 1 then
            verdict = '-';
        end
        while length(verdict) < 11
            verdict = verdict + ' ';
        end
        ryr = '-';
        if is_res == 1 then
            ryr = string(round(strtod(res_yr_s)*100)/100);
        end
        while length(ryr) < 9
            ryr = ryr + ' ';
        end
        row = row + nm + pl + ideg + verdict + ryr + string(round(decay_yr*100)/100);
        rows_out = [rows_out; row];
    end
    draw_lines(rows_out, 0.85, 0.033, 1);

    xs2pdf(f.figure_id, BASE + '/pages_fullpool/_01_index_p' + string(p) + '.pdf');
    close(f);
end
disp('wrote ' + string(npages) + ' index pages');
quit
