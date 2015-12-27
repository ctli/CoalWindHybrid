clear
close all
clc
format compact

x = [0.5, 0.75, 1];

y1=[ % subcritical air cooling (least efficient)
 3.65132743362832e+002
 3.25132743362832e+002
 3.08141592920354e+002];

y2 =[ % supercritical air cooling (use this one in modeling!)
 3.54867256637168e+002
 3.11327433628319e+002
 3.01061946902655e+002];

y3 =[ % subcritical water cooling
 3.49911504424779e+002
 3.05663716814159e+002
 2.95398230088496e+002];

y4 =[ % ultra supercritical air cooling
 3.44601769911504e+002
 3.02123893805310e+002
 2.84778761061947e+002];

y5 =[ % supercritical water cooling
 3.29734513274336e+002
 2.96460176991150e+002
 2.79823008849558e+002];

y6 =[ % ultra supper water cooling (most efficient)
 3.24778761061947e+002
 2.92920353982301e+002
 2.74513274336283e+002];


%%
figure(1); clf; hold on; box on;
plot(x, y1, '+-', 'linewidth', 1);
plot(x, y2, '*-', 'linewidth', 1);
plot(x, y3, 'o-', 'linewidth', 1);
plot(x, y4, '^-', 'linewidth', 1);
plot(x, y5, 'x-', 'linewidth', 1);
plot(x, y6, 's-', 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);
xlim([0.25 1.25]);
ylim([260 380]);
xlabel('Output Power (Normalized)');
ylabel('Coal Consumption Rate (g/kWh)');
title('Typical Coal Plants ');
set(gca, 'pos', [0.1300    0.1115    0.72    0.8135]);
box off;
% set(gca, 'ygrid', 'on');

legend('subcritical air cooling', ...
       'supercritical air cooling', ...
       'subcritical water cooling', ...
       'ultra supercritical air cooling', ...
       'supercritical water cooling', ...
       'ultra supper water cooling');

axis_pos = get(gca, 'pos');
x_lim = get(gca, 'xlim');
y_lim = get(gca, 'ylim');
y_ticklabel = 260:20:380;

     
% quadratic curve fit
XX = 0.45:0.05:1.05;

color_code = [0      0.5647 0.7412 %blue
              0.8510 0.3255 0.0980 %red
              0.9297 0.6941 0.1255 %yellow
              0.5000 0.1843 0.1255 %purple
              0.4667 0.6745 0.1882 %green
              0.3202 0.7451 0.7176 %light blue
              ];

% p1 = polyfit(x,y1',2);
% YY1 = XX.^2*p1(1) + XX.*p1(2) + p1(3);
% plot(XX, YY1, '--', 'color', color_code(1,:));
% text(0.3, 270, ['Eq: y1 = ', num2str(p1(1), '%3.0f'), 'x^2 ' num2str(p1(2), '%3.0f'), 'x + ', num2str(p1(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12);
% 
p2 = polyfit(x,y2',2);
YY2 = XX.^2*p2(1) + XX.*p2(2) + p2(3);
plot(XX, YY2, '--', 'color', color_code(2,:));
text(0.3, 272, ['Quadratic Eq: y2 = ', num2str(p2(1), '%3.0f'), 'x^2 ' num2str(p2(2), '%3.0f'), 'x + ', num2str(p2(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12, 'color', color_code(2,:));

[m,id] = min(YY2);
plot(XX(id), YY2(id), 'x', 'color', color_code(2,:), 'markersize', 6, 'linewidth', 2);

% p3 = polyfit(x,y3',2);
% YY3 = XX.^2*p3(1) + XX.*p3(2) + p3(3);
% plot(XX, YY3, '--', 'color', color_code(3,:));
% text(0.3, 290, ['Eq: y3 = ', num2str(p3(1), '%3.0f'), 'x^2 ' num2str(p3(2), '%3.0f'), 'x + ', num2str(p3(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12);
% 
% p4 = polyfit(x,y4',2); % too strong
% YY4 = XX.^2*p4(1) + XX.*p4(2) + p4(3);
% plot(XX, YY4, '--', 'color', color_code(4,:));
% text(0.3, 300, ['Eq: y4 = ', num2str(p4(1), '%3.0f'), 'x^2 ' num2str(p4(2), '%3.0f'), 'x + ', num2str(p4(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12);
% 
% p5 = polyfit(x,y5',2); % too strong
% YY5 = XX.^2*p5(1) + XX.*p5(2) + p5(3);
% plot(XX, YY5, '--', 'color', color_code(5,:));
% text(0.3, 310, ['Eq: y5 = ', num2str(p5(1), '%3.0f'), 'x^2 ' num2str(p5(2), '%3.0f'), 'x + ', num2str(p5(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12);
% 
% p6 = polyfit(x,y6',2);
% YY6 = XX.^2*p6(1) + XX.*p6(2) + p6(3);
% plot(XX, YY6, '--', 'color', color_code(6,:));
% text(0.3, 320, ['Eq: y6 = ', num2str(p6(1), '%3.0f'), 'x^2 ' num2str(p6(2), '%3.0f'), 'x + ', num2str(p6(3), '%3.0f')], 'FontName', 'Calibri', 'fontsize', 12);

ax2 = axes;
set(ax2, 'color', 'none', ...
         'pos', axis_pos, ...
         'yaxislocation', 'right', ...
         'ylim', y_lim, ...
         'ytick', (6000:500:8000)/21.8, ...
         'yticklabel', (6000:500:8000), ...
         'xaxislocation', 'top', ...
         'xlim', x_lim, ...
         'xticklabel', []);
ylabel(ax2, 'Heat Rate (Btu/kWh)');


%% fuel consumption
% supercritical air cooling
coal_fuelrate_x = 0.5:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]

coal_nameplate = 660; % [MW]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]

% ==============================
figure(2); clf; hold on; box on;
plot(coal_fuel_x, coal_fuel_y, 'color', color_code(2,:), 'linewidth', 1);
plot([coal_fuel_x(1),coal_fuel_x(end)],[coal_fuel_y(1),coal_fuel_y(end)], '--', 'color', [0.5 0.5 0.5]);

xa = 507/(266*2)*660 % min fuel rate @ 0.953 (i.e. 419MW)
ya = interp1(coal_fuel_x, coal_fuel_y, xa);
plot(xa, ya, 'ko', 'markerf', 'k', 'markersize', 3);
text(xa+0.02, ya-4, 'Min. Fuel Rate ');
text(xa+0.02, ya-8, ['@ ', num2str(xa, '%3.0f'), ' ( (y/x)^{\prime}=0 )']);

xb = 507*2/(266*3*2)*660 % inflection point @ 0.6353 (i.e. 629MW)
yb = interp1(coal_fuel_x, coal_fuel_y, xb);
plot(xb, yb, 'ko', 'markerf', 'k', 'markersize', 3);
text(xb+0.02, yb-4, 'Inflection Point');
text(xb+0.02, yb-8, ['@ ', num2str(xb, '%3.0f'), ' ( y^{\prime\prime}=0 )']);

grid on;
set(gcf, 'unit', 'inch', 'pos', [1.0    1.1    5.0000    3.7500]);
xlabel('Generation (MW)');
ylabel('Coal Consumption (ton/h)');
title('Coal Consumption of a Typical Coal Plants ');
legend('Supercritical air cooling');
set(legend, 'location', 'northwest');

% export_fig thermal_eff_2 -r600 -painters;


%% test convexity (extend x-axis)
coal_fuelrate_x = 0:0.01:1.2; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]

coal_nameplate = 660; % [MW]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]

% ==============================
figure(3); clf; hold on; %box on;
set(gcf, 'unit', 'inch', 'pos', [6.175    1.1    5.0000    3.7500]);
set(gca, 'pos', [0.115 0.11 0.735 0.78]);

plot(coal_fuel_x, coal_fuel_y, 'color', color_code(2,:), 'linewidth', 1);
xlim([0 800]);
ylim([0 250]);
set(gca, 'xtick', (0:0.2:1.2)*coal_nameplate);
set(gca, 'ytick', 0:50:250);
xlabel('Generation (MW)');
ylabel('Coal Consumption (ton/h)');
legend('Supercritical air cooling');
set(legend, 'location', 'northwest');
my_gridline;

plot([0.5 1]*coal_nameplate,interp1(coal_fuel_x, coal_fuel_y, [0.5 1]*coal_nameplate), '+-', 'color', [0.5 0.5 0.5]);
axx=[0.5 1]*coal_nameplate;
axy=[1 1]*205;
[arrowx,arrowy] = dsxy2figxy(gca, axx, axy);
annotation('doublearrow',arrowx,arrowy,'headsize',5,'headsize',5, 'color', [0.6 0.6 0.6]);
text(mean([0.5 1]*coal_nameplate), 212, 'Operation Range', 'color', [0.6 0.6 0.6], 'horizontalalignment', 'center', 'fontsize', 8);

xa = 507/(266*2)*coal_nameplate; % min fuel rate @ 0.953 (i.e. 419MW)
ya = interp1(coal_fuel_x, coal_fuel_y, xa);
plot(xa, ya, 'ko', 'markerf', 'k', 'markersize', 3);
text(xa-10, ya-11, 'Min. Fuel Rate ');
text(xa-10, ya-23, ['@ ', num2str(xa/coal_nameplate, '%0.3f'), ' ( (y/x)^{\prime}=0 )']);

xb = 507*2/(266*3*2)*coal_nameplate; % inflection point @ 0.6353 (i.e. 629MW)
yb = interp1(coal_fuel_x, coal_fuel_y, xb);
plot(xb, yb, 'ko', 'markerf', 'k', 'markersize', 3);
text(xb-10, yb-11, 'Inflection Point');
text(xb-10, yb-23, ['@ ', num2str(xb/coal_nameplate, '%0.3f'), ' ( y^{\prime\prime}=0 )']);

% text(775, 16, 'Extended x-axis', 'horizontalalignment', 'right', 'fontsize', 12);

ax1 = gca;

y_fuel = 0:50:250;
x_MW = interp1(coal_fuel_y, coal_fuel_x, y_fuel);
heat_rate = y_fuel.*(21.8)./(x_MW*1e3)*1e6; % [ton/MWh]] -> [Btu/kWh]; 1 tce = 27.778mmBT

ax2 = axes('pos', get(ax1, 'pos'), ...
           'color', 'none', ...
           'xaxislocation', 'top', 'xlim', get(ax1, 'xlim')/coal_nameplate, 'xtick', 0:0.2:1.2, ...
           'yaxislocation', 'right', 'ylim', get(ax1, 'ylim'), 'ytick', get(ax1, 'ytick'), 'yticklabel', num2str(heat_rate', '%7.0f'));
xlabel('Generation (Normalized)');
ylabel('Heat Rate (Btu/kWh)');

% export_fig thermal_eff_3 -r600 -painters;

