clear
close all
clc
format compact

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

% ====================
% Wind power
wind_file = 'Xilingol_2009';
load(wind_file);
wind_pwr_range = 0:500:5000; % [1x11]

% wind_ratio = linspace(1,1,3); save_name = ['DP_', wind_file, '_loop']; % Not allow economic wind curtialment
wind_ratio = linspace(1,0,11); save_name = ['DP_', wind_file, '_loop_new']; % Allow economic wind curtialment

target_pwr_range = 2500:500:8500; % [1x13]


%% Horizon-based unit commitment
% load FourteenUnits; % 'v_range', 'f_table', 'u_table', 'v_table', 'id_st', 'id_ed', 'v_st', 'v_ed'
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
% 
% target_pwr = target_pwr_range(tg);
% 
% % ====================
% % DP calculation
% N = length(wind_pwr);
% 
% % One state: number of coal plants commited
% id_ratio = zeros(1, length(wind_pwr));
% id_dispatch = zeros(1, length(wind_pwr));
% cmt_dispatch  = -1*ones(1, length(wind_pwr)); % Optimal commitment
% f_dispatch    = -1*ones(1, length(wind_pwr));
% u_dispatch    = -1*ones(coal_num, length(wind_pwr));
% v_dispatch    = -1*ones(coal_num, length(wind_pwr));
% wind_dispatch = -1*ones(1, length(wind_pwr));
% coal_dispatch = -1*ones(1, length(wind_pwr));
% J_star        = -1*ones(1, length(wind_pwr)); % Cummulative cost
% 
% opt_cost_startup  = -1*ones(1, length(wind_pwr));
% opt_cost_base_vom = -1*ones(1, length(wind_pwr));
% opt_cost_fuel     = -1*ones(1, length(wind_pwr));
% opt_cost_ramp     = -1*ones(1, length(wind_pwr));
% 
% t = N;
% wind_pwr_tmp = wind_pwr(t)*wind_ratio;
% coal_pwr_tmp = target_pwr - wind_pwr_tmp;
% coal_pwr_tmp(coal_pwr_tmp<coal_min) = coal_min;
% 
% id_tmp = interp1(v_range, id_range, coal_pwr_tmp);
% id_tmp = ceil(id_tmp);
% f_tmp = f_table(id_tmp,:); % [pwr range]x[cmt]
% 
% cost_base_vom_tmp = coal_pwr_tmp*coal_baseload;
% cost_fuel_tmp = f_tmp*coal_price;
% J_tmp = cost_fuel_tmp + repmat(cost_base_vom_tmp', 1, coal_num);
% [value, id_opt] = min(J_tmp(:));
% [id_row, id_col] = ind2sub(size(cost_fuel_tmp), id_opt);
% J_star(t) = value;
% id_ratio(t) = wind_ratio(id_row);
% id_dispatch(t) = id_tmp(id_row);
% cmt_dispatch(t) = id_col;
% f_dispatch(t) = f_tmp(id_opt);
% u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_dispatch(t));
% v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_dispatch(t));
% wind_dispatch(t) = wind_pwr_tmp(id_row);
% coal_dispatch(t) = coal_pwr_tmp(id_row);
% opt_cost_startup(t) = 0;
% opt_cost_base_vom(t) = cost_base_vom_tmp(id_row);
% opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
% opt_cost_ramp(t) = 0;
% 
% for t = N-1:-1:1
%     wind_pwr_tmp = wind_pwr(t)*wind_ratio;
%     coal_pwr_tmp = target_pwr - wind_pwr_tmp;
%     coal_pwr_tmp(coal_pwr_tmp<coal_min) = coal_min;
% 
%     id_tmp = interp1(v_range, id_range, coal_pwr_tmp);
%     id_tmp = ceil(id_tmp);
%     f_tmp = f_table(id_tmp,:); % [pwr range]x[cmt]
% 
%     cost_base_vom_tmp = coal_pwr_tmp * coal_baseload;
%     cost_fuel_tmp = f_tmp * coal_price;
%     
%     cost_startup_tmp = zeros(length(wind_ratio), coal_num);
%     cost_ramp_tmp = zeros(length(wind_ratio), coal_num);
%     for cmt = 1:coal_num
%         u_tmp = u_table(:,id_tmp,cmt);
%         v_tmp = v_table(:,id_tmp,cmt);
%         
%         c = cmt_dispatch(t+1) - cmt;
%         if c>0
%             cost_startup_tmp(:,cmt) = c*coal_startup_cost;
%         else
%             cost_startup_tmp(:,cmt) = 0;
%         end
%         
%         d_coal_pwr = u_tmp - repmat(u_dispatch(:,t+1), 1, length(wind_ratio)); % [14x10] = [14units]x[wind ratio]
%         d_coal_pctg = d_coal_pwr/coal_nameplate;
%         x = abs(d_coal_pctg);
%         y = (x-0.3)*6.5/0.7+1.5;
%         ramp_scale = ones(size(x));
%         ramp_scale(x>0.3) = y(x>0.3);
%         cost_ramp_tmp(:,cmt) = sum(abs(d_coal_pwr).*ramp_scale * coal_loadfollow);
%     end
%     J_tmp = J_star(t+1) + cost_fuel_tmp + repmat(cost_base_vom_tmp', 1, coal_num) + cost_startup_tmp + cost_ramp_tmp;
%     
%     [value, id_opt] = min(J_tmp(:));
%     [id_row, id_col] = ind2sub(size(cost_fuel_tmp), id_opt);
%     J_star(t) = value;
% 	id_ratio(t) = wind_ratio(id_row);
%     id_dispatch(t) = id_tmp(id_row);
%     cmt_dispatch(t) = id_col;
%     f_dispatch(t) = f_tmp(id_opt);
%     u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_dispatch(t));
%     v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_dispatch(t));
%     coal_dispatch(t) = coal_pwr_tmp(id_row);
%     wind_dispatch(t) = target_pwr - coal_dispatch(t);
%     opt_cost_startup(t) = cost_startup_tmp(id_opt);
%     opt_cost_base_vom(t) = cost_base_vom_tmp(id_row);
%     opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
%     opt_cost_ramp(t) = cost_ramp_tmp(id_opt);
% end
% wind_curtail = wind_pwr - wind_dispatch;
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

% save(save_name, ...
%      'table_coal_pwr_v', 'table_coal_pwr_u', ...
%      'table_wind_pwr', 'table_wind_curtail', ...
%      'table_cost_startup', 'table_cost_base_vom', 'table_cost_fuel', 'table_cost_ramp', 'table_cost_total');

% load(['DP_', wind_file, '_loop']);
load(['DP_', wind_file, '_loop_new']);

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

