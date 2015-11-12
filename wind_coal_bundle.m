% 2.5 GW wind at Xilingol in best year (1987) & 9240 MW of coal (660 MW x 14)
clear
close all
clc

p_HVDC = 9000; % [MW] line capacity of HVDC line: 9GW

yr_range = 1979:2009;

subregion_name = 'Xilingol';
lat = 43.9;
lon = 116;
color_code = [0 114 189]/255; % blue
wst_yr = 8;
bst_yr = 9;

v_cutin = 3; % [m/s]
v_cutout = 25; % [m/s]
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

%% supercritical air cooling
coal_nameplate = 660; % [MW]
coal_ramping = 0.01*coal_nameplate*60; % 396 MW/hr; (ramping: 1 precent/min)
coal_num = 2;

coal_fuelrate_x = 0.4:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]

% ==============================
% fuel rate
figure(1); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*392, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuelrate_y, 'linewidth', 1);
xlim([-15 710]);
ylim([298 392]);
my_gridline([1 1 1]*0.85, 'front');
set(gca, 'layer', 'top');
set(gca, 'fontsize', 9);
xlabel('Output Power (MW)');
ylabel('Coal Consumption Rate (g/kWh)', 'fontweight', 'bold', 'fontsize', 12);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);
set(gca, 'pos', [0.1300    0.115    0.7750    0.77]);
legend('Working Range', 'Fuel Rate');
set(legend, 'pos', [0.4924    0.75    0.2890    0.0950]);
text(0.4*coal_nameplate, mean(get(gca, 'ylim')), '40% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);
text(1.0*coal_nameplate, mean(get(gca, 'ylim')), '100% of Nameplate Capacity    ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);

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

% ==============================
% fuel consumption
figure(2); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*202, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuel_y, 'linewidth', 1);
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


%% two coal units
dx1 = 1000;
dx2 = 1001;
x_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
x_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
f_coal_unit1 = 266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542; %[g/kWh]
f_coal_unit2 = 266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542; %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

tic;
p_range = 450:50:1200;
opt_u1 = -1*ones(1,length(p_range));
opt_u2 = -1*ones(1,length(p_range));
opt_f = -1*ones(1,length(p_range));
for i = 1:length(p_range)
    p_level = p_range(i);
    u1 = linspace(min(x_coal_unit1),max(x_coal_unit1),1000);
    u1_out = u1 - 0.08*coal_nameplate;
    u2_out = p_level - u1_out;
    u2 = u2_out + 0.08*coal_nameplate;
    f = interp2(x_coal_unit1, x_coal_unit2, f_coal, u1, u2);
    if any(~isnan(f))
        [value,id] = min(f);
        opt_u1(i) = u1(id);
        opt_u2(i) = u2(id);
        opt_f(i) = value;
    end
end
exclud_id = find(opt_f==-1);
p_range(exclud_id) = [];
opt_u1(exclud_id) = [];
opt_u2(exclud_id) = [];
opt_f(exclud_id) = [];
toc;

%% ========================
dx1 = 20; % very course grid for plotting
dx2 = 21;
x_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
x_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
p_coal_unit1 = (linspace(0.4,1,dx1)-0.08) * coal_nameplate;
p_coal_unit2 = (linspace(0.4,1,dx2)-0.08) * coal_nameplate;
[c1,c2] = meshgrid(p_coal_unit1,p_coal_unit2);
p_coal = c1 + c2;
f_coal_unit1 = 266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542; %[g/kWh]
f_coal_unit2 = 266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542; %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

figure(3); clf;
mesh(x_coal_unit1, x_coal_unit2, p_coal); hold on;
plot3(opt_u2, opt_u1, p_range, 'o-', 'linewidth', 1);
xlabel('u2 (MW)');
ylabel('u1 (MW)');
zlabel('Combined Power Output (MW)');
view(315, 35);
set(gcf, 'unit', 'inch', 'pos', [6.2    5.85    5.0000    3.7500]);

figure(31); clf;
mesh(x_coal_unit1, x_coal_unit2, p_coal); hold on;
plot3(opt_u2, opt_u1, p_range, 'o-', 'linewidth', 1);
xlabel('u2 (MW)');
ylabel('u1 (MW)');
zlabel('Combined Power Output (MW)');
axis equal
xlim([250 680])
ylim([250 680])
view(0, 90);
set(gcf, 'units', 'inch', 'pos', [6.25  1.0900   5   3.75]);

figure(4); clf;
mesh(x_coal_unit1, x_coal_unit2, f_coal); hold on;
surface([0.4 1]*coal_nameplate, [0.4 1]*coal_nameplate, ones(2,2)*680, 'facec', [1 0.5 0.5], 'edgecolor', 'none');
alpha(0.7);
xlabel('u2 (MW)');
ylabel('u1 (MW)');
zlabel('Combined Fuel Consumption (ton/h)');
view(315+180, 35);
set(gcf, 'unit', 'inch', 'pos', [11.4    5.85    5.0000    3.7500]);

figure(41); clf;
mesh(x_coal_unit1, x_coal_unit2, f_coal); hold on;
plot3(opt_u1, opt_u2, opt_f, 'linewidth', 2);
xlabel('u2 (MW)');
ylabel('u1 (MW)');
zlabel('Combined Fuel Consumption (ton/h)');
view(115, 25);
set(gcf, 'unit', 'inch', 'pos', [11.4    1.0900    5.0000    3.7500]);

% ==================
figure(5); clf; hold on; box on;
plot([200 800], [200 800], '-', 'color', [0.7 0.7 0.7]);
plot(opt_u1, opt_u2, 'o-', 'markersize', 3);
xlabel('Optimal u1');
ylabel('Optimal u2');
set(gcf, 'unit', 'inch', 'pos', [8.6979    3.5625    5.0000    3.7500]);
axis equal;
axis([200 800 200 800]);
my_gridline;

