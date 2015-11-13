clear
close all
clc
format compact

coal_nameplate = 660; % [MW]

coal_fuelrate_x = 0.4:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]


p0 = 0;
q0 = 0;
f0 = 0;

p1 = u_coal_unit1;
q1 = p1 - 0.08*coal_nameplate;
f1 = f_coal_unit1;

p2 = uu_range;
q2 = vv_range;
f2 = opt_f;

% ========
q = 0:1:vv_range(end);
f = zeros(1,length(q));
u1_log = zeros(1, length(q));
u2_log = zeros(1, length(q));
f_u1 = zeros(1, length(q));
f_u2 = zeros(1, length(q));

id1 = (q>0 & q<=q1(end));
f(id1) = interp1(q1, f1, q(id1));
u1_log(id1) = interp1(q1, p1, q(id1));
f_u1(id1) = interp1(p1, f1, u1_log(id1));

id2 = (q>q1(end));
f(id2) = interp1(q2, f2, q(id2));
u1_log(id2) = interp1(q2, opt_u1, q(id2));
u2_log(id2) = interp1(q2, opt_u2, q(id2));
f_u1(id2) = interp1(p1, f1, u1_log(id2));
f_u2(id2) = interp1(p1, f1, u2_log(id2));

v1_log = u1_log - 0.08*coal_nameplate;
v2_log = zeros(1,length(q));
v2_log(u2_log>0) = u2_log(u2_log>0) - 0.08*coal_nameplate;

% ================================
figure(6); clf; hold on; box on;
plot(q0, f0, 'x', q1, f1, '-', q2, opt_f, '-', 'linewidth', 1);

xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Coal Consumption (ton/hr)');
ylim([0 400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend('No Unit is Commited', 'One Unit is Commited', 'Two Units are Commited');
set(legend, 'location', 'northwest');

% ================================
figure(7); clf; hold on; box on;
plot(q0, f0, 'x', q1, f1, '-', q2, opt_f, '-', 'linewidth', 1);
% area(q,f, 'facec', [0.8 0.95 1], 'edgecolor', 'none');
ha = area(q, [f_u1;f_u2]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');

xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Coal Consumption (ton/hr)');
ylim([0 400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend(ha, '1st Unit', '2nd Unit');
set(legend, 'location', 'northwest');

% ================================
figure(8); clf; hold on; box on;
ha = area(q, [u1_log;u2_log]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
plot(q, v1_log + v2_log, 'linewidth', 1);
xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Power Output (MW)');
ylim([0 1400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend('1st Unit', '2nd Unit', 'Useful Output');
set(legend, 'location', 'northwest');

