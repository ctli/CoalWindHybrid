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

%%
% supercritical air cooling
coal_nameplate = 660; % [MW]
coal_ramping = 0.01*coal_nameplate*60; % 396 MW/hr; (ramping: 1 precent/min)
coal_num = 2;

coal_fuelrate_x = 0.4:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]

% fuel rate
figure(1); clf; hold on; box on;
plot(coal_fuel_x, coal_fuelrate_y, 'linewidth', 1);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);
xlabel('Output Power (Normalized)');
ylabel('Coal Consumption Rate (g/kWh)');
title('Fuel Rate of Typical Coal Plants ');
grid on;

% fuel consumption
figure(2); clf; hold on; box on;
plot(coal_fuel_x, coal_fuel_y, 'linewidth', 1);
plot([coal_fuel_x(1),coal_fuel_x(end)],[coal_fuel_y(1),coal_fuel_y(end)], '--', 'color', [0.5 0.5 0.5]);

xa = 507/(266*2)*660; % min fuel rate @ 0.953 (i.e. 419MW)
ya = interp1(coal_fuel_x, coal_fuel_y, xa);
plot(xa, ya, 'ko', 'markerf', 'k', 'markersize', 3);
text(xa+0.02, ya-4, 'Min. Fuel Rate ');
text(xa+0.02, ya-8, ['@ ', num2str(xa, '%3.0f'), ' ( (y/x)^{\prime}=0 )']);

xb = 507*2/(266*3*2)*660; % inflection point @ 0.6353 (i.e. 629MW)
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


%% two coal units
x_coal_unit1 = linspace(0.4,1,10) * coal_nameplate;
x_coal_unit2 = linspace(0.4,1,11) * coal_nameplate;
p_coal_unit1 = (linspace(0.4,1,10)-0.08) * coal_nameplate;
p_coal_unit2 = (linspace(0.4,1,11)-0.08) * coal_nameplate;
[c1,c2] = meshgrid(p_coal_unit1,p_coal_unit2);
p_coal = c1 + c2;

f_coal_unit1 = 266*linspace(0.4,1,10).^2 -507*linspace(0.4,1,10) + 542; %[g/kWh]
f_coal_unit2 = 266*linspace(0.4,1,11).^2 -507*linspace(0.4,1,11) + 542; %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

figure(3); clf;
mesh(x_coal_unit1, x_coal_unit2, p_coal); hold on;
xlabel('u2 (11 grids)');
ylabel('u1 (10 grids)');

figure(4); clf;
mesh(x_coal_unit1, x_coal_unit2, f_coal); hold on;
surface(x_coal_unit1, x_coal_unit2, ones(size(f_coal))*680, 'facec', [0 1 0], 'edgecolor', [0 0.8 0]);
xlabel('u2 (11 grids)');
ylabel('u1 (10 grids)');

p_range = [600, 800, 1000];
opt_u1 = -1*ones(1,length(p_range));
opt_u2 = -1*ones(1,length(p_range));
for i = 1:length(p_range)
p_level = p_range(i);
isoquant1p = linspace(min(x_coal_unit1),max(x_coal_unit1),50);
isoquant2p = (p_level - isoquant1p*0.92) + 0.08*coal_nameplate;

out_of_bnd1 = isoquant2p>max(x_coal_unit2);
out_of_bnd2 = isoquant2p<min(x_coal_unit2);
out_of_bnd = out_of_bnd1 | out_of_bnd2;
if any(~out_of_bnd)
    isoquant1p(out_of_bnd) = [];
    isoquant2p(out_of_bnd) = [];
    
    isoquant3p = p_level * ones(size(isoquant1p));
    
    isoquant1f = isoquant1p;
    isoquant2f = isoquant2p;
    isoquant3f = interp2(x_coal_unit1, x_coal_unit2, f_coal, isoquant1f, isoquant2f);
    [value,id] = min(isoquant3f);
    opt_u1(i) = isoquant1p(id);
    opt_u2(i) = isoquant2p(id);
    
    figure(3); hold on;
    plot3(isoquant1p, isoquant2p, isoquant3p, 'linewidth', 2);
    figure(4); hold on;
    plot3(isoquant1f, isoquant2f, isoquant3f, 'linewidth', 2);
end

end


%%
% for t = 1%1:length(yr_range)
%     yr = yr_range(t);
%     file_name = [subregion_name, '_', num2str(yr)];
%     load(file_name);
%     
%     for coal_num = 2%4:14
%         for wind_capacity = 500%0:500:5000 % [MW]
%             p_wind = p*wind_capacity; % [MW]
%             
%             flag = ones(1,8760); % 1: plant-constrained; 2: line-constrained
%             
%             p_diff = p_HVDC - p_wind;
%             p_coal_unit = p_diff/coal_num + coal_inhouse; % how much each coal unit need to generate
%             flag(p_coal_unit<coal_nameplate) = 2; % line-constrained
%             p_coal_unit(p_coal_unit>coal_nameplate) = coal_nameplate;
%             
%             f_coal = interp1(coal_fuel_x, coal_fuel_y, p_coal_unit);
%             
%             % coal ramping
%             dp_coal = diff(p_coal_unit);
%             if any(abs(dp_coal) > coal_ramping)
%                 disp('ramping constraint is violated!');
%             end
%             
%             % bundled generation
%             p_coal_total = (p_coal_unit - coal_inhouse)*coal_num;
%             p_bundle = p_coal_total + p_wind;
%             dp_bundle = diff(p_bundle);
%             
%             p_coal_14 = (p_coal_unit - coal_inhouse)*ones(1,coal_num);
%         end
%     end
%     toc;
% end

