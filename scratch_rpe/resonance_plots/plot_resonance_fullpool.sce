// plot_resonance_fullpool.sce -- same 6-panel page as plot_resonance.sce,
// but for ALL 253 GTO/HEO objects in the candidate pool (not just the
// 27 my own Eq.12 detector flags as resonant) -- user wants to eyeball
// every object and pick the real U-turn cases manually. Resonance
// marking (red dashed line, info-panel resonance line) is drawn only
// when the manifest says is_resonant=1; every object still gets the
// resonance-condition/decay-condition curves in the (a,e) panel so a
// borderline or missed case is still visually checkable.

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots';
RE = 6378.1363;

function M = read_numeric_csv(path)
    lines = mgetl(path);
    n = size(lines, 1) - 1;
    M = [];
    for k = 2:n+1
        toks = strsplit(lines(k), ',');
        row = zeros(1, size(toks, 1));
        for j = 1:size(toks, 1)
            row(j) = strtod(toks(j));
        end
        M = [M; row];
    end
endfunction

mlines = mgetl(BASE + '/manifest_fullpool.csv');
nobj = size(mlines, 1) - 1;

mkdir(BASE + '/pages_fullpool');

for oi = 1:nobj
    toks = strsplit(mlines(oi + 1), ',');
    norad = toks(1);
    name = toks(2);
    pool = toks(3);
    i_mean = strtod(toks(6));
    lambda_crit = strtod(toks(7));
    is_resonant = strtod(toks(8));
    resonance_year_s = toks(9);
    decay_year = strtod(toks(11));
    first_year = strtod(toks(12));
    last_year = strtod(toks(13));
    note = toks(14);
    csv_rel = toks(15);

    data = read_numeric_csv(BASE + '/' + csv_rel);
    t_year = data(:, 2);
    a_km   = data(:, 3);
    e_val  = data(:, 4);
    hp_km  = data(:, 8);
    azim   = data(:, 9);

    f = scf();
    f.figure_size = [1500, 950];
    f.background = -2;

    // (1,1) resonance angle
    subplot(2, 3, 1);
    plot(t_year, azim, 'b-');
    if is_resonant == 1 then
        resonance_year = strtod(resonance_year_s);
        plot([resonance_year, resonance_year], [0, 360], 'r--');
    end
    xtitle('Resonance angle psi = (RAAN+AOP-L_sun) mod 360', ...
           'Year', 'psi [deg]');
    xgrid(1);

    // (1,2) perigee height
    subplot(2, 3, 2);
    plot(t_year, hp_km, 'b-');
    hpr = [min(hp_km) - 10, max(hp_km) + 10];
    if is_resonant == 1 then
        plot([resonance_year, resonance_year], hpr, 'r--');
    end
    plot([first_year, last_year], [100, 100], 'k:');
    xtitle('Perigee height h_p(t)  (red=resonance if detected, dotted=100km decay)', ...
           'Year', 'h_p [km]');
    xgrid(1);

    // (1,3) (a,e) plane
    subplot(2, 3, 3);
    e_grid = linspace(0.001, max(e_val) * 1.05 + 0.05, 200);
    a_res = lambda_crit ./ (1 - e_grid.^2).^(4/7);
    a_decay = (RE + 100) ./ (1 - e_grid);
    plot(e_val, a_km, 'b.-');
    plot(e_grid, a_res, 'r--');
    plot(e_grid, a_decay, 'k:');
    xtitle('(a,e) plane: orbit (blue), resonance cond. (red), decay cond. (dotted)', ...
           'eccentricity e', 'a [km]');
    xgrid(1);
    legend(['orbit trajectory'; 'resonance condition (Eq.12)'; 'decay condition (h_p=100km)'], 4);

    // (2,1) semi-major axis
    subplot(2, 3, 4);
    plot(t_year, a_km, 'b-');
    xtitle('Semi-major axis a(t)', 'Year', 'a [km]');
    xgrid(1);

    // (2,2) eccentricity
    subplot(2, 3, 5);
    plot(t_year, e_val, 'b-');
    xtitle('Eccentricity e(t)', 'Year', 'e');
    xgrid(1);

    // (2,3) info panel
    subplot(2, 3, 6);
    a2 = gca();
    a2.axes_visible = ['off', 'off', 'off'];
    a2.data_bounds = [0, 0; 1, 1];

    xstring(0.05, 0.92, 'NORAD ' + norad + '  ' + name);
    t = gce(); t.font_size = 3;
    xstring(0.05, 0.80, 'Pool: ' + pool);
    xstring(0.05, 0.72, 'Mean inclination: ' + string(round(i_mean*100)/100) + ' deg');
    xstring(0.05, 0.64, 'lambda_crit (Eq.12): ' + string(round(lambda_crit*100)/100));
    if is_resonant == 1 then
        xstring(0.05, 0.56, 'DETECTOR: RESONANT, epoch year ' + resonance_year_s);
    else
        xstring(0.05, 0.56, 'DETECTOR: not resonant (no Eq.12 crossing)');
    end
    xstring(0.05, 0.48, 'TLE record: ' + string(round(first_year*100)/100) + ' to ' + ...
            string(round(last_year*100)/100));
    xstring(0.05, 0.40, 'Observed decay: year ' + string(round(decay_year*100)/100));
    if note <> '' then
        xstring(0.05, 0.28, 'NOTE: ' + note);
        tn = gce(); tn.font_size = 1;
    end

    outpdf = BASE + '/pages_fullpool/obj_' + norad + '.pdf';
    xs2pdf(f.figure_id, outpdf);
    close(f);
    if modulo(oi, 20) == 0 then
        disp(string(oi) + ' / ' + string(nobj) + ' done');
    end
end

disp('ALL DONE');
quit
