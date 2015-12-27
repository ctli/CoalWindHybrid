clear
close all
clc
format compact

coal_nameplate = 660; % [MW]


%% One unit
% dx = 200;
% u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
% v_coal_unit = u_coal_unit-0.08*coal_nameplate;
% f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; % [ton/h]
% 
% v = v_coal_unit;
% f = f_coal_unit;
% 
% figure(1); clf; hold on; box on;
% plot(v, f, 'linewidth', 1);
% % plot(x([1,end]), y([1,end]), '--', 'color', [0.7 0.7 0.7])
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');


%% Two units 
dx = 1000; % 0.115275 sec
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; % [ton/h]
u_min = min(u_coal_unit)*ones(1,dx);
v_min = min(v_coal_unit)*ones(1,dx);
f_min = min(f_coal_unit)*ones(1,dx);

% ==============================
u1a = u_min;
v1a = v_min;
f1a = f_min;
u2a = u_coal_unit;
v2a = v_coal_unit;
f2a = f_coal_unit;
ua = u1a + u2a;
va = v1a + v2a;
fa = f1a + f2a;

u1b = u_coal_unit;
v1b = v_coal_unit;
f1b = f_coal_unit;
u2b = u_coal_unit;
v2b = v_coal_unit;
f2b = f_coal_unit;
ub = u1b + u2b;
vb = v1b + v2b;
fb = f1b + f2b;

figure(2); clf; hold on; box on;
plot(va, fa, 'linewidth', 1);
plot(vb, fb, 'linewidth', 1);
xlabel('Output power (in-house use excluded)');
ylabel('Coal Consumption (ton/h)');
legend('u1@minumum; u2 increases', 'u1 & u2 increase equally');
set(legend, 'location', 'northwest');

% ==============================
u1a_extended = interp1(va, u1a, vb);
u2a_extended = interp1(va, u2a, vb);
v1a_extended = interp1(va, v1a, vb);
v2a_extended = interp1(va, v2a, vb);
fa_extended = interp1(va, fa, vb);
f = [fa_extended; fb];
[value, id_row] = min(f);

id_col = 1:length(id_row);
id = sub2ind([2,dx], id_row, id_col);

u1 = [u1a_extended; u1b];
u2 = [u2a_extended; u2b];
v1 = [v1a_extended; v1b];
v2 = [v2a_extended; v2b];
opt_f = value;
opt_u1 = u1(id);
opt_u2 = u2(id);
opt_v1 = v1(id);
opt_v2 = v2(id);
v_unique = vb;

% save('TwoUnitsClean', ...
%      'v_unique', 'opt_f', ...
%      'opt_u1', 'opt_u2', ...
%      'opt_v1', 'opt_v2', ...
%      'dx');

% ==========
figure(21); clf;
subplot(3,1,1);
plot(vb, opt_u1); hold on;
ylim([0 700]);
ylabel('u1');

subplot(3,1,2);
plot(vb, opt_u2); hold on;
ylim([0 700]);
ylabel('u2');

subplot(3,1,3);
ha = area(vb, [opt_u1; opt_u2]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
legend('u1', 'u2');
set(legend, 'location', 'northwest');


%% Three units;
% dx = 2000; % 0.115979 sec
% u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
% v_coal_unit = u_coal_unit-0.08*coal_nameplate;
% f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; % [ton/h]
% u_min = min(u_coal_unit)*ones(1,dx);
% v_min = min(v_coal_unit)*ones(1,dx);
% f_min = min(f_coal_unit)*ones(1,dx);
% 
% % ==============================
% u1a = u_min;
% v1a = v_min;
% f1a = f_min;
% u2a = u_min;
% v2a = v_min;
% f2a = f_min;
% u3a = u_coal_unit;
% v3a = v_coal_unit;
% f3a = f_coal_unit;
% ua = u1a + u2a + u3a;
% va = v1a + v2a + v3a;
% fa = f1a + f2a + f3a;
% 
% u1b = u_min;
% v1b = v_min;
% f1b = f_min;
% u2b = u_coal_unit;
% v2b = v_coal_unit;
% f2b = f_coal_unit;
% u3b = u_coal_unit;
% v3b = v_coal_unit;
% f3b = f_coal_unit;
% ub = u1b + u2b + u3b;
% vb = v1b + v2b + v3b;
% fb = f1b + f2b + f3b;
% 
% u1c = u_coal_unit;
% v1c = v_coal_unit;
% f1c = f_coal_unit;
% u2c = u_coal_unit;
% v2c = v_coal_unit;
% f2c = f_coal_unit;
% u3c = u_coal_unit;
% v3c = v_coal_unit;
% f3c = f_coal_unit;
% uc = u1c + u2c + u3c;
% vc = v1c + v2c + v3c;
% fc = f1c + f2c + f3c;
% 
% figure(3); clf; hold on; box on;
% plot(va, fa, 'linewidth', 1);
% plot(vb, fb, 'linewidth', 1);
% plot(vc, fc, 'linewidth', 1);
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% legend('2 plant2 @ minumum; 1 plants increase equally', '1 plant @ minumum; 2 plants increase equally', 'all three plants increase equally');
% set(legend, 'location', 'northwest');
% 
% % ==============================
% u1a_extended = interp1(va, u1a, vc);
% u2a_extended = interp1(va, u2a, vc);
% u3a_extended = interp1(va, u3a, vc);
% v1a_extended = interp1(va, v1a, vc);
% v2a_extended = interp1(va, v2a, vc);
% v3a_extended = interp1(va, v3a, vc);
% fa_extended = interp1(va, fa, vc);
% 
% u1b_extended = interp1(vb, u1b, vc);
% u2b_extended = interp1(vb, u2b, vc);
% u3b_extended = interp1(vb, u3b, vc);
% v1b_extended = interp1(vb, v1b, vc);
% v2b_extended = interp1(vb, v2b, vc);
% v3b_extended = interp1(vb, v3b, vc);
% fb_extended = interp1(vb, fb, vc);
% 
% f = [fa_extended; fb_extended; fc];
% [value, id_row] = min(f);
% 
% id_col = 1:length(id_row);
% id = sub2ind([3,dx], id_row, id_col);
% 
% u1 = [u1a_extended; u1b_extended; u1c];
% u2 = [u2a_extended; u2b_extended; u2c];
% u3 = [u3a_extended; u3b_extended; u3c];
% v1 = [v1a_extended; v1b_extended; v1c];
% v2 = [v2a_extended; v2b_extended; v2c];
% v3 = [v3a_extended; v3b_extended; v3c];
% opt_f = value;
% opt_u1 = u1(id);
% opt_u2 = u2(id);
% opt_u3 = u3(id);
% opt_v1 = v1(id);
% opt_v2 = v2(id);
% opt_v3 = v3(id);
% 
% v_unique = vc;
% 
% % save('ThreeUnitsClean', ...
% %      'v_unique', 'opt_f', ...
% %      'opt_u1', 'opt_u2', 'opt_u3', ...
% %      'opt_v1', 'opt_v2', 'opt_v3', ...
% %      'dx');
% 
% % ==========
% figure(31); clf;
% subplot(4,1,1);
% plot(vc, opt_u1, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(4,1,2);
% plot(vc, opt_u2, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(4,1,3);
% plot(vc, opt_u3, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(4,1,4);
% ha = area(vc, [opt_u1; opt_u2; opt_u3]', 'edgecolor', 'none');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3');
% set(legend, 'location', 'northwest');


%% Four units
% dx = 4000; % 0.132656 sec
% u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
% v_coal_unit = u_coal_unit-0.08*coal_nameplate;
% f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; % [ton/h]
% u_min = min(u_coal_unit)*ones(1,dx);
% v_min = min(v_coal_unit)*ones(1,dx);
% f_min = min(f_coal_unit)*ones(1,dx);
% 
% % ==============================
% u1a = u_min;
% v1a = v_min;
% f1a = f_min;
% u2a = u_min;
% v2a = v_min;
% f2a = f_min;
% u3a = u_min;
% v3a = v_min;
% f3a = f_min;
% u4a = u_coal_unit;
% v4a = v_coal_unit;
% f4a = f_coal_unit;
% ua = u1a + u2a + u3a + u4a;
% va = v1a + v2a + v3a + v4a;
% fa = f1a + f2a + f3a + f4a;
% 
% u1b = u_min;
% v1b = v_min;
% f1b = f_min;
% u2b = u_min;
% v2b = v_min;
% f2b = f_min;
% u3b = u_coal_unit;
% v3b = v_coal_unit;
% f3b = f_coal_unit;
% u4b = u_coal_unit;
% v4b = v_coal_unit;
% f4b = f_coal_unit;
% ub = u1b + u2b + u3b + u4b;
% vb = v1b + v2b + v3b + v4b;
% fb = f1b + f2b + f3b + f4b;
% 
% u1c = u_min;
% v1c = v_min;
% f1c = f_min;
% u2c = u_coal_unit;
% v2c = v_coal_unit;
% f2c = f_coal_unit;
% u3c = u_coal_unit;
% v3c = v_coal_unit;
% f3c = f_coal_unit;
% u4c = u_coal_unit;
% v4c = v_coal_unit;
% f4c = f_coal_unit;
% uc = u1c + u2c + u3c + u4c;
% vc = v1c + v2c + v3c + v4c;
% fc = f1c + f2c + f3c + f4c;
% 
% u1d = u_coal_unit;
% v1d = v_coal_unit;
% f1d = f_coal_unit;
% u2d = u_coal_unit;
% v2d = v_coal_unit;
% f2d = f_coal_unit;
% u3d = u_coal_unit;
% v3d = v_coal_unit;
% f3d = f_coal_unit;
% u4d = u_coal_unit;
% v4d = v_coal_unit;
% f4d = f_coal_unit;
% ud = u1d + u2d + u3d + u4d;
% vd = v1d + v2d + v3d + v4d;
% fd = f1d + f2d + f3d + f4d;
% 
% figure(4); clf; hold on; box on;
% plot(va, fa, 'linewidth', 1);
% plot(vb, fb, 'linewidth', 1);
% plot(vc, fc, 'linewidth', 1);
% plot(vd, fd, 'linewidth', 1);
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% legend('3 plants @ minumum; 1 plant increases', '2 plant2 @ minumum; 2 plants increase equally', '1 plant @ minumum; 3 plants increase equally', 'all three plants increase equally');
% set(legend, 'location', 'northwest');
% 
% % ==============================
% u1a_extended = interp1(va, u1a, vd);
% u2a_extended = interp1(va, u2a, vd);
% u3a_extended = interp1(va, u3a, vd);
% u4a_extended = interp1(va, u4a, vd);
% v1a_extended = interp1(va, v1a, vd);
% v2a_extended = interp1(va, v2a, vd);
% v3a_extended = interp1(va, v3a, vd);
% v4a_extended = interp1(va, v4a, vd);
% fa_extended = interp1(va, fa, vd);
% 
% u1b_extended = interp1(vb, u1b, vd);
% u2b_extended = interp1(vb, u2b, vd);
% u3b_extended = interp1(vb, u3b, vd);
% u4b_extended = interp1(vb, u4b, vd);
% v1b_extended = interp1(vb, v1b, vd);
% v2b_extended = interp1(vb, v2b, vd);
% v3b_extended = interp1(vb, v3b, vd);
% v4b_extended = interp1(vb, v4b, vd);
% fb_extended = interp1(vb, fb, vd);
% 
% u1c_extended = interp1(vc, u1c, vd);
% u2c_extended = interp1(vc, u2c, vd);
% u3c_extended = interp1(vc, u3c, vd);
% u4c_extended = interp1(vc, u4c, vd);
% v1c_extended = interp1(vc, v1c, vd);
% v2c_extended = interp1(vc, v2c, vd);
% v3c_extended = interp1(vc, v3c, vd);
% v4c_extended = interp1(vc, v4c, vd);
% fc_extended = interp1(vc, fc, vd);
% 
% f = [fa_extended; fb_extended; fc_extended; fd];
% [value, id_row] = min(f);
% 
% id_col = 1:length(id_row);
% id = sub2ind([4,dx], id_row, id_col);
% 
% u1 = [u1a_extended; u1b_extended; u1c_extended; u1d];
% u2 = [u2a_extended; u2b_extended; u2c_extended; u2d];
% u3 = [u3a_extended; u3b_extended; u3c_extended; u3d];
% u4 = [u4a_extended; u4b_extended; u4c_extended; u4d];
% v1 = [v1a_extended; v1b_extended; v1c_extended; v1d];
% v2 = [v2a_extended; v2b_extended; v2c_extended; v2d];
% v3 = [v3a_extended; v3b_extended; v3c_extended; v3d];
% v4 = [v4a_extended; v4b_extended; v4c_extended; v4d];
% opt_f = value;
% opt_u1 = u1(id);
% opt_u2 = u2(id);
% opt_u3 = u3(id);
% opt_u4 = u4(id);
% opt_v1 = v1(id);
% opt_v2 = v2(id);
% opt_v3 = v3(id);
% opt_v4 = v4(id);
% 
% v_unique = vd;
% 
% % save('FourUnitsClean', ...
% %      'v_unique', 'opt_f', ...
% %      'opt_u1', 'opt_u2', 'opt_u3', 'opt_u4', ...
% %      'opt_v1', 'opt_v2', 'opt_v3', 'opt_v4', ...
% %      'dx');
% 
% % ==========
% figure(41); clf;
% subplot(5,1,1);
% plot(vd, opt_u1, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(5,1,2);
% plot(vd, opt_u2, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(5,1,3);
% plot(vd, opt_u3, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(5,1,4);
% plot(vd, opt_u4, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u4');
% 
% subplot(5,1,5);
% ha = area(vd, [opt_u1; opt_u2; opt_u3; opt_u4]', 'edgecolor', 'none');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3', 'u4');
% set(legend, 'location', 'northwest');


%% Five units
% dx = 6000; % 0.138804 sec
% u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
% v_coal_unit = u_coal_unit-0.08*coal_nameplate;
% f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; % [ton/h]
% u_min = min(u_coal_unit)*ones(1,dx);
% v_min = min(v_coal_unit)*ones(1,dx);
% f_min = min(f_coal_unit)*ones(1,dx);
% 
% % ==============================
% u1a = u_min;
% v1a = v_min;
% f1a = f_min;
% u2a = u_min;
% v2a = v_min;
% f2a = f_min;
% u3a = u_min;
% v3a = v_min;
% f3a = f_min;
% u4a = u_min;
% v4a = v_min;
% f4a = f_min;
% u5a = u_coal_unit;
% v5a = v_coal_unit;
% f5a = f_coal_unit;
% ua = u1a + u2a + u3a + u4a + u5a;
% va = v1a + v2a + v3a + v4a + v5a;
% fa = f1a + f2a + f3a + f4a + f5a;
% 
% u1b = u_min;
% v1b = v_min;
% f1b = f_min;
% u2b = u_min;
% v2b = v_min;
% f2b = f_min;
% u3b = u_min;
% v3b = v_min;
% f3b = f_min;
% u4b = u_coal_unit;
% v4b = v_coal_unit;
% f4b = f_coal_unit;
% u5b = u_coal_unit;
% v5b = v_coal_unit;
% f5b = f_coal_unit;
% ub = u1b + u2b + u3b + u4b + u5b;
% vb = v1b + v2b + v3b + v4b + v5b;
% fb = f1b + f2b + f3b + f4b + f5b;
% 
% u1c = u_min;
% v1c = v_min;
% f1c = f_min;
% u2c = u_min;
% v2c = v_min;
% f2c = f_min;
% u3c = u_coal_unit;
% v3c = v_coal_unit;
% f3c = f_coal_unit;
% u4c = u_coal_unit;
% v4c = v_coal_unit;
% f4c = f_coal_unit;
% u5c = u_coal_unit;
% v5c = v_coal_unit;
% f5c = f_coal_unit;
% uc = u1c + u2c + u3c + u4c + u5c;
% vc = v1c + v2c + v3c + v4c + v5c;
% fc = f1c + f2c + f3c + f4c + f5c;
% 
% u1d = u_min;
% v1d = v_min;
% f1d = f_min;
% u2d = u_coal_unit;
% v2d = v_coal_unit;
% f2d = f_coal_unit;
% u3d = u_coal_unit;
% v3d = v_coal_unit;
% f3d = f_coal_unit;
% u4d = u_coal_unit;
% v4d = v_coal_unit;
% f4d = f_coal_unit;
% u5d = u_coal_unit;
% v5d = v_coal_unit;
% f5d = f_coal_unit;
% ud = u1d + u2d + u3d + u4d + u5d;
% vd = v1d + v2d + v3d + v4d + v5d;
% fd = f1d + f2d + f3d + f4d + f5d;
% 
% u1e = u_coal_unit;
% v1e = v_coal_unit;
% f1e = f_coal_unit;
% u2e = u_coal_unit;
% v2e = v_coal_unit;
% f2e = f_coal_unit;
% u3e = u_coal_unit;
% v3e = v_coal_unit;
% f3e = f_coal_unit;
% u4e = u_coal_unit;
% v4e = v_coal_unit;
% f4e = f_coal_unit;
% u5e = u_coal_unit;
% v5e = v_coal_unit;
% f5e = f_coal_unit;
% ue = u1e + u2e + u3e + u4e + u5e;
% ve = v1e + v2e + v3e + v4e + v5e;
% fe = f1e + f2e + f3e + f4e + f5e;
% 
% figure(5); clf; hold on; box on;
% plot(va, fa, 'linewidth', 1);
% plot(vb, fb, 'linewidth', 1);
% plot(vc, fc, 'linewidth', 1);
% plot(vd, fd, 'linewidth', 1);
% plot(ve, fe, 'linewidth', 1);
% xlabel('Output power (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% legend('4 plants @ minumum; 1 plant increases', '3 plants @ minumum; 1 plant increases', '2 plant2 @ minumum; 2 plants increase equally', '1 plant @ minumum; 3 plants increase equally', 'all three plants increase equally');
% set(legend, 'location', 'southeast');
% 
% % ==============================
% u1a_extended = interp1(va, u1a, ve);
% u2a_extended = interp1(va, u2a, ve);
% u3a_extended = interp1(va, u3a, ve);
% u4a_extended = interp1(va, u4a, ve);
% u5a_extended = interp1(va, u5a, ve);
% v1a_extended = interp1(va, v1a, ve);
% v2a_extended = interp1(va, v2a, ve);
% v3a_extended = interp1(va, v3a, ve);
% v4a_extended = interp1(va, v4a, ve);
% v5a_extended = interp1(va, v5a, ve);
% fa_extended = interp1(va, fa, ve);
% 
% u1b_extended = interp1(vb, u1b, ve);
% u2b_extended = interp1(vb, u2b, ve);
% u3b_extended = interp1(vb, u3b, ve);
% u4b_extended = interp1(vb, u4b, ve);
% u5b_extended = interp1(vb, u5b, ve);
% v1b_extended = interp1(vb, v1b, ve);
% v2b_extended = interp1(vb, v2b, ve);
% v3b_extended = interp1(vb, v3b, ve);
% v4b_extended = interp1(vb, v4b, ve);
% v5b_extended = interp1(vb, v5b, ve);
% fb_extended = interp1(vb, fb, ve);
% 
% u1c_extended = interp1(vc, u1c, ve);
% u2c_extended = interp1(vc, u2c, ve);
% u3c_extended = interp1(vc, u3c, ve);
% u4c_extended = interp1(vc, u4c, ve);
% u5c_extended = interp1(vc, u5c, ve);
% v1c_extended = interp1(vc, v1c, ve);
% v2c_extended = interp1(vc, v2c, ve);
% v3c_extended = interp1(vc, v3c, ve);
% v4c_extended = interp1(vc, v4c, ve);
% v5c_extended = interp1(vc, v5c, ve);
% fc_extended = interp1(vc, fc, ve);
% 
% u1d_extended = interp1(vd, u1d, ve);
% u2d_extended = interp1(vd, u2d, ve);
% u3d_extended = interp1(vd, u3d, ve);
% u4d_extended = interp1(vd, u4d, ve);
% u5d_extended = interp1(vd, u5d, ve);
% v1d_extended = interp1(vd, v1d, ve);
% v2d_extended = interp1(vd, v2d, ve);
% v3d_extended = interp1(vd, v3d, ve);
% v4d_extended = interp1(vd, v4d, ve);
% v5d_extended = interp1(vd, v5d, ve);
% fd_extended = interp1(vd, fd, ve);
% 
% f = [fa_extended; fb_extended; fc_extended; fd_extended; fe];
% [value, id_row] = min(f);
% 
% id_col = 1:length(id_row);
% id = sub2ind([5,dx], id_row, id_col);
% 
% u1 = [u1a_extended; u1b_extended; u1c_extended; u1d_extended; u1e];
% u2 = [u2a_extended; u2b_extended; u2c_extended; u2d_extended; u2e];
% u3 = [u3a_extended; u3b_extended; u3c_extended; u3d_extended; u3e];
% u4 = [u4a_extended; u4b_extended; u4c_extended; u4d_extended; u4e];
% u5 = [u5a_extended; u5b_extended; u5c_extended; u5d_extended; u5e];
% v1 = [v1a_extended; v1b_extended; v1c_extended; v1d_extended; v1e];
% v2 = [v2a_extended; v2b_extended; v2c_extended; v2d_extended; v2e];
% v3 = [v3a_extended; v3b_extended; v3c_extended; v3d_extended; v3e];
% v4 = [v4a_extended; v4b_extended; v4c_extended; v4d_extended; v4e];
% v5 = [v5a_extended; v5b_extended; v5c_extended; v5d_extended; v5e];
% opt_f = value;
% opt_u1 = u1(id);
% opt_u2 = u2(id);
% opt_u3 = u3(id);
% opt_u4 = u4(id);
% opt_u5 = u5(id);
% opt_v1 = v1(id);
% opt_v2 = v2(id);
% opt_v3 = v3(id);
% opt_v4 = v4(id);
% opt_v5 = v5(id);
% 
% v_unique = ve;
% 
% % save('FiveUnitsClean', ...
% %      'v_unique', 'opt_f', ...
% %      'opt_u1', 'opt_u2', 'opt_u3', 'opt_u4', 'opt_u5', ...
% %      'opt_v1', 'opt_v2', 'opt_v3', 'opt_v4', 'opt_v5', ...
% %      'dx');
% 
% % ==========
% figure(51); clf;
% subplot(6,1,1);
% plot(ve, opt_u1, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(6,1,2);
% plot(ve, opt_u2, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(6,1,3);
% plot(ve, opt_u3, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(6,1,4);
% plot(ve, opt_u4, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u4');
% 
% subplot(6,1,5);
% plot(ve, opt_u5, 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u5');
% 
% subplot(6,1,6);
% ha = area(ve, [opt_u1; opt_u2; opt_u3; opt_u4; opt_u5]', 'edgecolor', 'none');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% set(ha(4), 'facec', [0.95 0.95 0], 'edgecolor', 'none');
% set(ha(5), 'facec', [0.5 1 0.5], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3', 'u4', 'u5');
% set(legend, 'location', 'east');

