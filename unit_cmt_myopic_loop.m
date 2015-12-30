clear
close all
clc
format compact

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
coal_min = coal_nameplate*0.4; % 264 [MW] 
coal_num = 14;

% ==============================
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

% ==============================
% Wind power
wind_file = 'Xilingol_2009';
load(wind_file);

wind_pwr_range = 0:500:5000; % [1x11]

% wind_ratio = linspace(1,1,3); save_name = ['Myopic_', wind_file, '_loop']; % Not allow economic wind curtialment
wind_ratio = linspace(1,0,101); save_name = ['Myopic_', wind_file, '_loop_new']; % Allow economic wind curtialment

target_pwr_range = 2500:500:8500; % [1x13]


%% Myopic unit commitment
% load  MyopicDispatch;
% id_range = 1:length(v_range);
% 
% table_coal_pwr_v   = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_coal_pwr_u   = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_wind_pwr     = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_wind_curtail = zeros(length(wind_pwr_range), length(target_pwr_range));
% 
% table_cost_total    = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_fuel     = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_startup  = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_base_vom = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_ramp     = zeros(length(wind_pwr_range), length(target_pwr_range));
% 
% tic;
% for w = 1:length(wind_pwr_range)
% wind_pwr = round(wind_pwr_range(w)*p)';
% 
% for tg = 1:length(target_pwr_range)
% disp(['w=', num2str(w), '; tg=', num2str(tg)]);
% target_pwr = target_pwr_range(tg);
% 
% % ====================
% % Myopic dispatch
% id_ratio = zeros(1, length(wind_pwr));
% id_dispatch = zeros(1, length(wind_pwr));
% wind_dispatch = zeros(1, length(wind_pwr));
% coal_dispatch = zeros(1, length(wind_pwr));
% for t = 1:length(wind_pwr)
%     wind_pwr_tmp = wind_pwr(t)*wind_ratio;
%     coal_pwr_tmp = target_pwr - wind_pwr_tmp;
%     coal_pwr_tmp(coal_pwr_tmp<coal_min) = coal_min;
% 
%     id_tmp = interp1(v_range, id_range, coal_pwr_tmp);
%     id_tmp = ceil(id_tmp);
%     f_tmp = interp1(id_range, f_myopic, id_tmp);
% 
%     cost_base_vom = coal_pwr_tmp*coal_baseload;
%     cost_fuel = f_tmp*coal_price;
%     [value, id_opt] = min(cost_base_vom + cost_fuel);
%     id_ratio(t) = id_opt;
%     id_dispatch(t) = id_tmp(id_opt);
%     coal_dispatch(t) = coal_pwr_tmp(id_opt);
%     wind_dispatch(t) = target_pwr - coal_dispatch(t);
% end
% cmt_dispatch = cmt_myopic(id_dispatch);
% f_dispatch = f_myopic(id_dispatch); % [ton/h]
% v_dispatch = v_myopic(:,id_dispatch);
% u_dispatch = u_myopic(:,id_dispatch);
% wind_curtail = wind_pwr - wind_dispatch;
% 
% % ====================
% % Check ramping
% d_coal_pwr = [zeros(1,coal_num); diff(u_dispatch')]'; % [14]x[8760]
% d_coal_pctg = d_coal_pwr/coal_nameplate;
% 
% % Check changes in commitment
% d_cmt = [0, diff(cmt_dispatch)];
% cmt_c = zeros(1, length(wind_pwr)); % Commition
% cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
% cmt_d = zeros(1, length(wind_pwr)); % Decommition
% cmt_d(d_cmt<0) = d_cmt(d_cmt<0);
% 
% % ====================
% % Costs
% opt_cost_startup = cmt_c * coal_startup_cost;
% opt_cost_base_vom = coal_dispatch * coal_baseload;
% opt_cost_fuel = f_dispatch * coal_price;
% 
% x = abs(d_coal_pctg);
% y = (x-0.3)*6.5/0.7+1.5;
% ramp_scale = ones(size(x)); % Exceeding 30 precent nameplate ramping is scaled by 1.5-8 times
% ramp_scale(x>0.3) = y(x>0.3);
% opt_cost_ramp = abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow;
% 
% cost_startup = sum(opt_cost_startup);
% cost_base_vom = sum(opt_cost_base_vom);
% cost_fuel = sum(opt_cost_fuel);
% cost_ramp = sum(opt_cost_ramp);
% cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp;
% 
% % ====================
% table_coal_pwr_v(w,tg) = sum(v_dispatch(:));
% table_coal_pwr_u(w,tg) = sum(u_dispatch(:));
% 
% table_wind_pwr(w,tg)     = sum(wind_pwr(:));
% table_wind_curtail(w,tg) = sum(wind_curtail(:));
% 
% table_cost_startup(w,tg)  = cost_startup;
% table_cost_base_vom(w,tg) = cost_base_vom;
% table_cost_fuel(w,tg)     = cost_fuel;
% table_cost_ramp(w,tg)     = cost_ramp;
% table_cost_total(w,tg)    = cost_total;
% 
% toc;
% end
% end
% 
% % save(save_name, ...
% %      'table_coal_pwr_v', 'table_coal_pwr_u', ...
% %      'table_wind_pwr', 'table_wind_curtail', ...
% %      'table_cost_startup', 'table_cost_base_vom', 'table_cost_fuel', 'table_cost_ramp', 'table_cost_total');

% load(['Myopic_', wind_file, '_loop']);
load(['Myopic_', wind_file, '_loop_new']);

pctg_startup = table_cost_startup./table_cost_total;
pctg_base_vom = table_cost_base_vom./table_cost_total;
pctg_fuel = table_cost_fuel./table_cost_total;
pctg_ramp = table_cost_ramp./table_cost_total;
pctg_wind_curtail = table_wind_curtail./(table_wind_pwr + table_wind_curtail);


%% Total generation cost
figure(1); clf;
surf(wind_pwr_range, target_pwr_range, table_cost_total'/1e9);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Total Generation Cost (Billion USD)');
ylim([2000 10000]);
zlim([0 1.8]);
set(gca, 'clim', [0.3 1.65]);
view(35, 10);

% ====================
figure(11); clf;
[C,h] = contourf(wind_pwr_range, target_pwr_range, table_cost_total'/1e9);
clabel(C, h);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
title('Countour: Total Generation Cost (Billion USD)');
set(gca, 'xtick', 0:1000:5000);
set(gca, 'clim', [0.3 1.65]);
grid on;

% ====================
figure(12); clf;
[C,h] = contourf(target_pwr_range, wind_pwr_range, table_cost_total/1e9);
clabel(C, h);
ylabel('Wind Capacity (MW)');
xlabel('Target Output (MW)');
title('Countour: Total Generation Cost (Billion USD)');
set(gca, 'ytick', 0:1000:5000);
set(gca, 'clim', [0.3 1.65]);
grid on;


%% Share of fuel cost
figure(2); clf;
surf(wind_pwr_range, target_pwr_range, pctg_fuel');
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Share of Fuel Cost (-)');
ylim([2000 10000]);
zlim([0 1]);
set(gca, 'clim', [0.66 0.86]);
view(30, 15);

% ====================
figure(21); clf;
[C,h] = contourf(wind_pwr_range, target_pwr_range, pctg_fuel');
clabel(C, h);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
title('Countour: Share of Fuel Cost (-)');
set(gca, 'xtick', 0:1000:5000);
set(gca, 'clim', [0.66 0.86]);
grid on;

% ====================
figure(22); clf;
[C,h] = contourf(target_pwr_range, wind_pwr_range, pctg_fuel);
clabel(C, h);
ylabel('Wind Capacity (MW)');
xlabel('Target Output (MW)');
title('Countour: Share of Fuel Cost (-)');
set(gca, 'ytick', 0:1000:5000);
set(gca, 'clim', [0.66 0.86]);
grid on;


%% Wind curtailment
figure(3); clf;
surf(wind_pwr_range, target_pwr_range, pctg_wind_curtail');
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Wind Curtailment (-)');
ylim([2000 10000]);
zlim([0 0.3]);
set(gca, 'clim', [0 0.3]);
view(30, 15);

% ====================
figure(31); clf;
[C,h] = contourf(wind_pwr_range, target_pwr_range, pctg_wind_curtail');
clabel(C, h);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
title('Countour: Wind Curtailment (-)');
set(gca, 'xtick', 0:1000:5000);
set(gca, 'clim', [0 0.3]);
grid on;

% ====================
figure(32); clf;
[C,h] = contourf(target_pwr_range, wind_pwr_range, pctg_wind_curtail);
clabel(C, h);
ylabel('Wind Capacity (MW)');
xlabel('Target Output (MW)');
title('Countour: Wind Curtailment (-)');
set(gca, 'ytick', 0:1000:5000);
set(gca, 'clim', [0 0.3]);
grid on;

