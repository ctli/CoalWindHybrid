clear
close all
clc

coal_nameplate = 660; % [MW]


%% No units
% figure(1); clf; hold on; box on;
% x = 0;
% y = 0;
% plot(x, y, 'x', 'linewidth', 1);


%% One unit
dx = 200;
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; %[g/kWh]

% ==============================
v = v_coal_unit;
y = f_coal_unit;

% ==============================
% figure(1); clf; hold on; box on;
% plot(v, y, 'linewidth', 1);
% % plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');


%% Two units
dx = 1000;
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; %[g/kWh]
u_min = min(u_coal_unit)*ones(1,dx);
v_min = min(v_coal_unit)*ones(1,dx);
f_min = min(f_coal_unit)*ones(1,dx);

% ==============================
v1a = v_min;
f1a = f_min;
v2a = v_coal_unit;
f2a = f_coal_unit;
va = v1a + v2a;
fa = f1a + f2a;

v1b = v_coal_unit;
f1b = f_coal_unit;
v2b = v_coal_unit;
f2b = f_coal_unit;
vb = v1b + v2b;
fb = f1b + f2b;

% figure(2); clf; hold on; box on;
% plot(va, fa, 'linewidth', 1);
% plot(vb, fb, 'linewidth', 1);
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% legend('u1@minumum; u2 increases', 'u1 & u2 increase equally');
% set(legend, 'location', 'northwest');

% ==============================
v1a_extended = interp1(va, v1a, vb);
v2a_extended = interp1(va, v2a, vb);
fa_extended = interp1(va, fa, vb);
f = [fa_extended; fb];
[value, id_row] = min(f);

id_col = 1:length(id_row);
id = sub2ind([2,dx], id_row, id_col);

v1 = [v1a_extended; v1b];
v2 = [v2a_extended; v2b];
opt_f = value;
opt_v1 = v1(id);
opt_v2 = v2(id);

% =====
% figure(21); clf;
% area(vb, [opt_v1; opt_v2]');
% 
% figure(22); clf;
% subplot(2,1,1);
% plot(vb, opt_v1); hold on;
% ylim([0 700]);
% 
% subplot(2,1,2);
% plot(vb, opt_v2); hold on;
% ylim([0 700]);


%% Three units;
dx = 2000;
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; %[g/kWh]
u_min = min(u_coal_unit)*ones(1,dx);
v_min = min(v_coal_unit)*ones(1,dx);
f_min = min(f_coal_unit)*ones(1,dx);

% ==============================
x1 = v_min;
f1 = f_min;
x2 = v_min;
f2 = f_min;
x3 = v_coal_unit;
f3 = f_coal_unit;
xa = x1 + x2 + x3;
fa = f1 + f2 + f3;

x1 = v_min;
f1 = f_min;
x2 = v_coal_unit;
f2 = f_coal_unit;
x3 = v_coal_unit;
f3 = f_coal_unit;
xb = x1 + x2 + x3;
fb = f1 + f2 + f3;

x1 = v_coal_unit;
f1 = f_coal_unit;
x2 = v_coal_unit;
f2 = f_coal_unit;
x3 = v_coal_unit;
f3 = f_coal_unit;
xc = x1 + x2 + x3;
fc = f1 + f2 + f3;

figure(2); clf; hold on; box on;
plot(xa, fa, 'linewidth', 1);
plot(xb, fb, 'linewidth', 1);
plot(xc, fc, 'linewidth', 1);

legend('two plants @ minumum; one plants increases', 'one plant @ minumum; two plants increase equally', 'all three plants increase equally');
set(legend, 'location', 'northwest');

xlabel('Output power (in-house use excluded)');
ylabel('Coal Consumption (ton/h)');


%% Four units
% x = v_coal_unit1*3;
% y = f_coal_unit1*3;
% plot(x, y, 'linewidth', 1);
% % plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])
% 
% x = v_coal_unit1*4;
% y = f_coal_unit1*4;
% plot(x, y, 'linewidth', 1);
% % plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])
% 
% grid on;

