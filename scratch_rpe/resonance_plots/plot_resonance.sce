// plot_resonance.sce -- Wang & Gurfil (2017)-style resonance plots for
// all 23 real objects found to show a genuine U-turn resonance crossing
// (OREM issue #56). Reads manifest.csv + one obj_<norad>.csv per object
// (both written by prepare_plot_data.py, which does all the orbital-
// mechanics work -- this script only plots).
//
// Per object, a 2x3-panel page:
//   (1,1) resonance angle psi = mod(RAAN+AOP-L_sun,360) vs time --
//         the U-turn itself (paper Fig 4/8/12 analog)
//   (1,2) perigee height h_p(t), resonance epoch + 100km decay
//         threshold marked (paper Fig 5/9/13 analog)
//   (1,3) trajectory in the (a,e) plane with the resonance-condition
//         curve a=lambda_crit/(1-e^2)^(4/7) and the 100km decay-
//         condition curve a=(RE+100)/(1-e) overlaid (paper Fig 7/11/15)
//   (2,1) semi-major axis a(t)
//   (2,2) eccentricity e(t)
//   (2,3) text info panel (NORAD/name/inclination/resonance date)
//
// One PDF page per object written to pages/obj_<norad>.pdf.

BASE = 'E:/GitHub/OREM/scratch_rpe/resonance_plots';
RE = 6378.1363;

function M = read_numeric_csv(path)
    lines = mgetl(path);
    n = size(lines, 1) - 1; // skip header
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

function s = safe_name(txt)
    s = txt;
endfunction

// ---- manifest: norad,name,n_tle,i_mean_deg,lambda_crit,resonance_jd,
//                resonance_year,decay_year,first_year,last_year,csv_path
mlines = mgetl(BASE + '/manifest.csv');
nobj = size(mlines, 1) - 1;

mkdir(BASE + '/pages');

for oi = 1:nobj
    toks = strsplit(mlines(oi + 1), ',');
    norad = toks(1);
    name = toks(2);
    i_mean = strtod(toks(4));
    lambda_crit = strtod(toks(5));
    resonance_year = strtod(toks(7));
    decay_year_s = toks(8);
    first_year = strtod(toks(9));
    last_year = strtod(toks(10));
    csv_rel = toks(11);

    data = read_numeric_csv(BASE + '/' + csv_rel);
    // columns: jd,year,a_km,e,i_deg,raan_deg,argp_deg,hp_km,azimuth_deg,lambda
    t_year = data(:, 2);
    a_km   = data(:, 3);
    e_val  = data(:, 4);
    hp_km  = data(:, 8);
    azim   = data(:, 9);
    lam    = data(:, 10);

    f = scf();
    f.figure_size = [1500, 950];
    f.background = -2;

    // (1,1) resonance angle U-turn
    subplot(2, 3, 1);
    plot(t_year, azim, 'b-');
    yr = [0, 360];
    plot([resonance_year, resonance_year], yr, 'r--');
    xtitle('Resonance angle psi = (RAAN+AOP-L_sun) mod 360', ...
           'Year', 'psi [deg]');
    xgrid(1);

    // (1,2) perigee height, resonance + decay threshold marked
    subplot(2, 3, 2);
    plot(t_year, hp_km, 'b-');
    hpr = [min(hp_km) - 10, max(hp_km) + 10];
    plot([resonance_year, resonance_year], hpr, 'r--');
    plot([first_year, last_year], [100, 100], 'k:');
    xtitle('Perigee height h_p(t)  (red=resonance, dotted=100km decay)', ...
           'Year', 'h_p [km]');
    xgrid(1);

    // (1,3) (a,e) plane trajectory + resonance/decay condition curves
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
    info = ['NORAD ' + norad + '  ' + name; ...
            'Mean inclination: ' + string(round(i_mean * 100) / 100) + ' deg'; ...
            'lambda_crit (Eq.12): ' + string(round(lambda_crit * 100) / 100); ...
            'Resonance epoch: year ' + string(round(resonance_year * 100) / 100); ...
            'TLE record: ' + string(round(first_year * 100) / 100) + ' to ' + ...
               string(round(last_year * 100) / 100); ...
            'Observed decay: year ' + decay_year_s];
    xstring(0.05, 0.85, info);

    outpdf = BASE + '/pages/obj_' + norad + '.pdf';
    xs2pdf(f.figure_id, outpdf);
    close(f);
    disp('wrote ' + outpdf);
end

disp('ALL DONE');
quit
