clear
close all
clc
format compact

yr_range = 1979:2009;

v_cutin = 3; % [m/s]
v_cutout = 25; % [m/s]
v_rated = 10.5;

% wind spd [m/s]; wind pwr [normalized]
pwr_curve = [
3	0.005993333
4	0.041946667
5	0.107866667
6	0.203746667
7	0.329586667
8	0.503373333
9	0.713106667
10	0.92884
10.61	1
];


%% wind speed barchat
for c = 5%1:5 % five locations
    switch c
        case 1
            subregion_name = 'Xilingol';
            lat = 43.9;
            lon = 116;
            color_code = [0 114 189]/255; % blue
            wst_yr = 8;
            bst_yr = 9;
        case 2
            subregion_name = 'Jiuquan';
            lat = 39.73;
            lon = 98.5;
            color_code = [217 83 25]/255; % red
            wst_yr = 22;
            bst_yr = 8;
        case 3
            subregion_name = 'Hami';
            lat = 42.9;
            lon = 93.44;
            color_code = [237 177 32]/255; % yellow
            wst_yr = 27;
            bst_yr = 7;
        case 4
            subregion_name = 'Zhundong';
            lat = 44.15;
            lon = 87.89;
            color_code = [126 47 142]/255; % purple
            wst_yr = 7;
            bst_yr = 19;
        case 5
            subregion_name = 'Ili';
            lat = 44;
            lon = 81.3;
            color_code = [119 172 48]/255; % green
            wst_yr = 29;
            bst_yr = 12;
    end
    
    i = 1;
    for iy = [1:10:length(yr_range), wst_yr, bst_yr] % 1979, 1989, 1999, 2009, worst yr, best yr
        yr = yr_range(iy);
        
        file_name = [subregion_name, '_', num2str(yr)];
        load(file_name);

        if iy == wst_yr
            extra_string = ' (worst year)';
            output_name = ['wsd_barchart_', subregion_name, '_worst'];
        elseif iy == bst_yr
            extra_string = ' (best year)';
            output_name = ['wsd_barchart_', subregion_name, '_best'];
        else
            extra_string = '';
            output_name = ['wsd_barchart_', subregion_name, '_', num2str(yr)];
        end

        % ==============================
        ctr = 0:1:25;
        cnt = hist(v, ctr);
        mn = mean(v);
        ma = median(v);
        
        % ==========
        figure(iy); clf; hold on; box on;
        bar(ctr, cnt/sum(cnt), 1, 'facec', color_code, 'edgecolor', 'w');
        xlim([-1 26]);
        ylim([0 0.3]);
        set(gca, 'fontsize', 8);
        set(gca, 'layer', 'bottom');
        set(gca, 'tickdir', 'out');
        set(gca, 'ygrid', 'on');
        set(gca, 'GridLineStyle', '-');
        xlabel('Wind Speed (m/s)');
        ylabel('Frequency (-)');
        title([subregion_name, '; ', num2str(yr), extra_string]);
        set(gcf, 'unit', 'inch', 'pos', [0.25+4.15*(i-1)    6.165    4.0000    3.0000]);
        
        hold on;
        plot([1 1]*mn, get(gca, 'ylim'), 'k', 'linewidth', 1);
        text(double(mn), 0.15, ['mean = ', num2str(mn, '%2.1f')], 'rotation', 90, 'verticalalignment', 'top', 'horizontalalignment', 'center', 'fontsize', 8);
        
        plot([1 1]*ma, get(gca, 'ylim'), 'r', 'linewidth', 1);
        text(double(ma), 0.15, ['median = ', num2str(ma, '%2.1f')], 'color', 'r', 'rotation', 90, 'verticalalignment', 'bottom', 'horizontalalignment', 'center', 'fontsize', 8);

        plot([1 1]*v_cutin, [0 0.3], '--', 'color', [0.7 0.7 0.7]);
        text(v_cutin, 0.15, ['Cut-In Speed (', num2str(v_cutin), ' m/s)'], 'rotation', 90, 'fontsize', 6, 'horizontalalignment', 'center');
        plot([1 1]*v_cutout, [0 0.3], '--', 'color', [0.7 0.7 0.7]);
        text(v_cutout, 0.15, ['Cut-Out Speed (', num2str(v_cutout), ' m/s)'], 'rotation', 90, 'fontsize', 6, 'horizontalalignment', 'center');
        plot([1 1]*v_rated, [0 0.3], '--', 'color', [0.7 0.7 0.7]);
        text(v_rated, 0.15, ['Rated Speed (', num2str(v_rated), ' m/s)'], 'rotation', 90, 'fontsize', 6, 'horizontalalignment', 'center');
        
        % export_fig(output_name, '-r300', '-painters');
        
        % ==============================
        ctr = 0:0.1:1;
        cnt = hist(p, ctr);
        mn = mean(p);
        ma = median(p);
        
        % ==========
        figure(100+i); clf;
        bar(ctr, cnt/sum(cnt), 1, 'facec', color_code, 'edgecolor', 'w');
%         axis equal;
        ylim([0 0.8]);
        xlim([-0.06 1.06]);
        set(gca, 'fontsize', 8);
        set(gca, 'layer', 'bottom');
        set(gca, 'tickdir', 'out');
        set(gca, 'ygrid', 'on');
        set(gca, 'GridLineStyle', '-');
        xlabel('Wind Power, \it{P} (normalized)');
        ylabel('Frequency (normalized)');
        title([subregion_name, '; ', num2str(yr)]);
        text(1.04, 0.95, '(bin size = 0.1)', 'fontsize', 7, 'horizontalalignment', 'right', 'color', [0.6 0.6 0.6]);
        set(gcf, 'unit', 'inch', 'pos', [0.25+4.15*(i-1)    2.130    4.0000    3.0000]);

        hold on;
        plot([1 1]*mn, get(gca, 'ylim'), 'k', 'linewidth', 1);
        text(double(mn), mean(get(gca, 'ylim')), ['mean = ', num2str(mn, '%2.2f')], 'rotation', 90, 'verticalalignment', 'top', 'horizontalalignment', 'center', 'fontsize', 8);
        
        plot([1 1]*ma, get(gca, 'ylim'), 'r', 'linewidth', 1);
        text(double(ma), mean(get(gca, 'ylim')), ['median = ', num2str(ma, '%2.2f')], 'color', 'r', 'rotation', 90, 'verticalalignment', 'bottom', 'horizontalalignment', 'center', 'fontsize', 8);
        
        % export_fig(output_name, '-r300', '-painters');
        
        i = i+1;
    end
end



