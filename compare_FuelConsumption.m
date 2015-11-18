clear
close all
clc
format compact

coal_nameplate = 660; % [MW]

% coal_fuelrate_x = 0.4:0.01:1; % [normalized]
% coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
% coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
% coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]


%% ========================================================================
p0 = 0;
q0 = 0;
f0 = 0;

% One unit
dx1 = 100;
p1 = linspace(0.4,1,dx1) * coal_nameplate;
q1 = p1 - 0.08*coal_nameplate;
f1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]

% Two units
load('TwoUnits');
p2 = opt_u1 + opt_u2;
q2 = opt_v1 + opt_v2;
f2 = opt_f;

% Three units
load('ThreeUnits');
p3 = opt_u1 + opt_u2 + opt_u3;
q3 = opt_v1 + opt_v2 + opt_v3;
f3 = opt_f;

% Four units
load('FourUnits');
p4 = opt_u1 + opt_u2  + opt_u3 + opt_u4;
q4 = opt_v1 + opt_v2 + opt_v3 + opt_v4;
f4 = opt_f;

figure(1); clf; hold on; box on;
plot(q0, f0, 'x');
plot(q1, f1, 'linewidth', 1);
plot(q2, f2, 'x-', 'linewidth', 1);
plot(q3, f3, 'linewidth', 1);
plot(q4, f4, 'linewidth', 1);
legend('No Unit is Commited', 'One Unit is Commited', 'Two Units is Commited', 'Three Units is Commited', 'Four Units is Commited');
set(legend, 'location', 'northwest');
xlabel('Power Output, MW (exclude in-house use)');
ylabel('Coal Consumption (ton/h)');

dx1 = 200; % 4.33GB RAM
dx2 = 200;
u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;

v_coal_unit1 = u_coal_unit1-0.08*coal_nameplate;
v_coal_unit2 = u_coal_unit2-0.08*coal_nameplate;

f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %[g/kWh]

% ==============================
dx = 200; % 4.33GB RAM
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; %[g/kWh]

% x = 0;
% y = 0;
% plot(x, y, 'x', 'linewidth', 1);

% x = v_coal_unit;
% y = f_coal_unit;
% plot(x, y, 'linewidth', 1);

x = v_coal_unit*2;
y = f_coal_unit*2;
plot(x, y, 'linewidth', 1);
% plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])

x = v_coal_unit*3;
y = f_coal_unit*3;
plot(x, y, 'linewidth', 1);
% plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])

x = v_coal_unit*4;
y = f_coal_unit*4;
plot(x, y, 'linewidth', 1);
% plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])


%%
% load('TwoUnits');
% figure(2); clf; hold on;
% 
% subplot(3,1,1);
% plot(v_unique, opt_u1);
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(3,1,2);
% plot(v_unique, opt_u2);
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(3,1,3);
% ha = area(v_unique, [opt_u1; opt_u2]');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% legend('u1', 'u2');
% set(legend, 'location', 'northwest');


%%
% load('ThreeUnits');
% figure(3); clf; hold on;
% 
% subplot(4,1,1);
% plot(v_unique, opt_u1);
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(4,1,2);
% plot(v_unique, opt_u2);
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(4,1,3);
% plot(v_unique, opt_u3);
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(4,1,4);
% ha = area(v_unique, [opt_u1; opt_u2; opt_u3]');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3');
% set(legend, 'location', 'northwest');


%%
% load('FourUnits');
% figure(4); clf; hold on;
% 
% subplot(5,1,1);
% plot(v_unique, opt_u1);
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(5,1,2);
% plot(v_unique, opt_u2);
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(5,1,3);
% plot(v_unique, opt_u3);
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(5,1,4);
% plot(v_unique, opt_u4);
% ylim([0 700]);
% ylabel('u4');
% 
% subplot(5,1,5);
% ha = area(v_unique, [opt_u1; opt_u2; opt_u3; opt_u4]', 'edgecolor', 'none');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3', 'u4');
% set(legend, 'location', 'northwest');


% % ========
% q = 0:1:v_unique(end);
% f = zeros(1,length(q));
% u1_log = zeros(1, length(q));
% u2_log = zeros(1, length(q));
% f_u1 = zeros(1, length(q));
% f_u2 = zeros(1, length(q));
% 
% id1 = (q>0 & q<=q1(end));
% f(id1) = interp1(q1, f1, q(id1));
% u1_log(id1) = interp1(q1, p1, q(id1));
% f_u1(id1) = interp1(p1, f1, u1_log(id1));
% 
% id2 = (q>q1(end));
% f(id2) = interp1(q2, f2, q(id2));
% u1_log(id2) = interp1(q2, opt_u1, q(id2));
% u2_log(id2) = interp1(q2, opt_u2, q(id2));
% f_u1(id2) = interp1(p1, f1, u1_log(id2));
% f_u2(id2) = interp1(p1, f1, u2_log(id2));
% 
% v1_log = u1_log - 0.08*coal_nameplate;
% v2_log = zeros(1,length(q));
% v2_log(u2_log>0) = u2_log(u2_log>0) - 0.08*coal_nameplate;

% 
% % ================================
% figure(9); clf; hold on; box on;
% plot(q0, f0, 'x', q1, f1, '-', q2, opt_f, '-', 'linewidth', 1);
% % area(q,f, 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% ha = area(q, [f_u1;f_u2]');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% 
% xlabel('Power Output, MW (exclude in-house generation)')
% ylabel('Coal Consumption (ton/hr)');
% ylim([0 400]);
% xlim([0 1250]);
% grid on;
% set(gca, 'layer', 'top');
% legend(ha, '1st Unit', '2nd Unit');
% set(legend, 'location', 'northwest');
% 
% % ================================
% figure(10); clf; hold on; box on;
% ha = area(q, [u1_log;u2_log]');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% plot(q, v1_log + v2_log, 'linewidth', 1);
% xlabel('Power Output, MW (exclude in-house generation)')
% ylabel('Power Output (MW)');
% ylim([0 1400]);
% xlim([0 1250]);
% grid on;
% set(gca, 'layer', 'top');
% legend('1st Unit', '2nd Unit', 'Useful Output');
% set(legend, 'location', 'northwest');
