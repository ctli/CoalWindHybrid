% 2.5 GW wind at Xilingol in best year (1987) & 9240 MW of coal (660 MW x 14)
clear
close all
clc


%% supercritical air cooling
coal_nameplate = 660; % [MW]
coal_fuelrate_x = 0.4:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]


%% Fuel rate [g/kWh]
figure(1); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*392, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuelrate_y, 'linewidth', 1, 'color', [217 83 25]/255);
xlim([-15 710]);
ylim([296 392]);
% my_gridline([1 1 1]*0.85, 'front');
set(gca, 'layer', 'top');
set(gca, 'fontsize', 9);
xlabel('Output Power (MW)');
ylabel('Coal Consumption Rate (g/kWh)', 'fontweight', 'bold', 'fontsize', 12);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);
set(gca, 'pos', [0.1300    0.115    0.73    0.77]);
legend('Working Range', 'Fuel Rate');
set(legend, 'pos', [0.465    0.75    0.2890    0.0950]);
text(0.4*coal_nameplate, mean(get(gca, 'ylim')), '40% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);
text(1.0*coal_nameplate, mean(get(gca, 'ylim')), '100% of Nameplate Capacity    ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);

ax_pos = get(gca, 'pos');
x_lim = get(gca, 'xlim');
y_lim = get(gca, 'ylim');

ax2 = axes;
set(ax2, 'pos', ax_pos);
set(ax2, 'color', 'none');
set(ax2, 'fontsize', 9);
x_ticklabel = 0:0.15:1.5;
set(ax2, 'xaxislocation', 'top', 'xlim', x_lim, 'xtick', x_ticklabel*coal_nameplate, 'xticklabel', x_ticklabel);
set(ax2, 'yaxislocation', 'right', 'ylim', y_lim, 'ytick', (1000:500:11000)/21.8, 'yticklabel', 1000:500:11000);
xlabel('Output Power (Normalized)', 'fontsize', 9);
ylabel('Heat Rate (Btu/kWh)', 'fontweight', 'bold', 'fontsize', 12);


%% Fuel consumption [ton/h]
figure(2); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*202, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuel_y, 'linewidth', 1, 'color', [217 83 25]/255);
plot([coal_fuel_x(1),coal_fuel_x(end)],[coal_fuel_y(1),coal_fuel_y(end)], '--', 'color', [0.5 0.5 0.5]);
ylim([98 202]);
xlim([-15 710]);

xa = 507/(266*2)*660; % min fuel rate @ 0.953 (i.e. 419MW)
ya = interp1(coal_fuel_x, coal_fuel_y, xa);
plot(xa, ya, 'ko', 'markerf', 'k', 'markersize', 3);
text(xa+0.02, ya+7, 'Min. Fuel Rate      ', 'horizontalalignment', 'right');
text(xa+0.02, ya+2, ['@ ', num2str(xa, '%3.0f'), ' ( (y/x)^{\prime}=0 ) '], 'horizontalalignment', 'right');

xb = 507*2/(266*3*2)*660; % inflection point @ 0.6353 (i.e. 629MW)
yb = interp1(coal_fuel_x, coal_fuel_y, xb);
plot(xb, yb, 'ko', 'markerf', 'k', 'markersize', 3);
text(xb+0.02, yb-4, 'Inflection Point');
text(xb+0.02, yb-8, ['@ ', num2str(xb, '%3.0f'), ' ( y^{\prime\prime}=0 )']);

text(0.4*coal_nameplate, mean(get(gca, 'ylim')), '40% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);
text(1.0*coal_nameplate, mean(get(gca, 'ylim')), '100% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);

my_gridline([1 1 1]*0.85, 'front');
set(gca, 'fontsize', 9);
xlabel('Generation (MW)');
ylabel('Coal Consumption (ton/h)', 'fontweight', 'bold', 'fontsize', 12);
set(gcf, 'unit', 'inch', 'pos', [1.0    1.09    5.0000    3.7500]);
set(gca, 'pos', [0.1300    0.115    0.7750    0.77]);

ax_pos = get(gca, 'pos');
x_lim = get(gca, 'xlim');
y_lim = get(gca, 'ylim');
y_tick = get(gca, 'ytick');

ax2 = axes;
set(ax2, 'pos', ax_pos);
set(ax2, 'color', 'none');
set(ax2, 'fontsize', 9);
x_ticklabel = 0:0.15:1.5;
set(ax2, 'xaxislocation', 'top', 'xlim', x_lim, 'xtick', x_ticklabel*coal_nameplate, 'xticklabel', x_ticklabel);
set(ax2, 'yaxislocation', 'right', 'ylim', y_lim, 'ytick', y_tick, 'yticklabel', []);
xlabel('Output Power (Normalized)', 'fontsize', 9);


