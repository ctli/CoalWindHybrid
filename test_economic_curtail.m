clear
close all
clc

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
coal_min = coal_nameplate*0.4; % 264 [MW] 
coal_num = 14;

% ====================
% Costs of coal power plants
% China:
coal_price = 57; % i.e. 380 [RMB/tonne] = 57 [$/tonne] (1RMB = 0.15USD)
coal_heatcontent = 21.8; % 5500 [kcal/kg] = 21.8 [mmBtu/tonne] (1kcal = 3.96567 Btu)

% ==========
% US:
startup = 38; % hot startup cost [$/MW]
startup_oth = 5.81; % other hot startup cost [$/MW]
coal_startup = (startup + startup_oth)*coal_nameplate; % 2.8915e+04 [$]

startup_fuel = 10.1; % hot start fuel [mmBtu/MW]
coal_startup_fuel = (startup_fuel*coal_nameplate)/coal_heatcontent*coal_price; % 1.7429e+04 [$]

coal_startup_cost = coal_startup + coal_startup_fuel; % 4.6344e+04 [$]

% Variable cost
coal_loadfollow = 1.72; % load following cost [$/MW]
coal_baseload = 3.22; % Baseload variable cost [$/MWh]


%% Without economic curtailment
load('DP_Xilingol_2009_debug');
cost_startup = sum(opt_cost_startup)
cost_base_vom = sum(opt_cost_base_vom)
cost_fuel = sum(opt_cost_fuel)
cost_ramp = sum(opt_cost_ramp)

disp('====================');
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp

disp('~~~~~~~~~~~~~~~');
% Check ramping
d_coal_pwr = [zeros(1,coal_num);
              diff(u_dispatch')]'; % [14]x[8760]
d_coal_pctg = d_coal_pwr/coal_nameplate;

% Check changes in commitment
d_cmt = [0, diff(cmt_dispatch)];
cmt_c = zeros(1, length(wind_pwr)); % Commition
cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
cmt_d = zeros(1, length(wind_pwr)); % Decommition
cmt_d(d_cmt<0) = d_cmt(d_cmt<0);

opt_cost_startup2 = cmt_c * coal_startup_cost;
opt_cost_base_vom2 = coal_dispatch * coal_baseload;
opt_cost_fuel2 = f_dispatch * coal_price;

x = abs(d_coal_pctg);
y = (x-0.3)*6.5/0.7+1.5;
ramp_scale = ones(size(x)); % Exceeding 30 precent nameplate ramping is scaled by 1.5-8 times
ramp_scale(x>0.3) = y(x>0.3);
opt_cost_ramp2 = abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow; % [$]

cost_startup = sum(opt_cost_startup2)
cost_base_vom = sum(opt_cost_base_vom2)
cost_fuel = sum(opt_cost_fuel2)
cost_ramp = sum(opt_cost_ramp2)

disp('====================');
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp
J_star(1)

% ========================================
figure(1); clf;
ax1 = subplot(2,2,1);
ha = area([coal_dispatch;wind_dispatch;wind_curtail]', 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.7 0.7]);
set(ha(2), 'facec', [0.6 1 0]);
set(ha(3), 'facec', [0 0.7 0]);
% ====================
ax2 = subplot(2,2,2);
plot(cmt_dispatch);

figure(2); clf;
plot(J_star, 'x-');


%% With economic curtailment
disp('========================================')
disp('========================================')
load('DP_Xilingol_2009_debug_new');
cost_startup = sum(opt_cost_startup)
cost_base_vom = sum(opt_cost_base_vom)
cost_fuel = sum(opt_cost_fuel)
cost_ramp = sum(opt_cost_ramp)

disp('====================');
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp


% ========================================
figure(1);
ax3 = subplot(2,2,3);
ha = area([coal_dispatch;wind_dispatch;wind_curtail]', 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.7 0.7]);
set(ha(2), 'facec', [0.6 1 0]);
set(ha(3), 'facec', [0 0.7 0]);
% ====================
ax4 = subplot(2,2,4);
plot(cmt_dispatch);

linkaxes([ax1, ax2, ax3, ax4], 'x');
xlim([8590 8760]);


figure(2); hold on;
plot(J_star, 'o-');

