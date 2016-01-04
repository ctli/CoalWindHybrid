% With parallelization; with economic wind curtailment

clear
close all
clc
format compact

load FourteenUnits; % 'v_range', 'f_table', 'u_table', 'v_table', 'id_st', 'id_ed', 'v_st', 'v_ed'
id_range = 1:length(v_range);

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
v_min = min(v_range); % 220 [MW] 
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

wind_ratio = linspace(1,1,3); save_name = ['DP_', wind_file, '_loop']; % Not allow economic wind curtialment
% wind_ratio = linspace(1,0,101); save_name = ['DP_', wind_file, '_loop_new'];  % Allow economic wind curtialment
wind_num = length(wind_ratio);

target_pwr_range = 2500:500:8500; % [1x13]


%% Horizon-based unit commitment
% table_coal_v        = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_coal_u        = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_coal_dispatch = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_wind_dispatch = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_wind_curtail  = zeros(length(wind_pwr_range), length(target_pwr_range));
% 
% table_cost_startup  = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_base_vom = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_fuel     = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_ramp     = zeros(length(wind_pwr_range), length(target_pwr_range));
% table_cost_total    = zeros(length(wind_pwr_range), length(target_pwr_range));
% 
% tic;
% for w = 1:length(wind_pwr_range)
% wind_pwr = round(wind_pwr_range(w)*p)';
% N = length(wind_pwr);
% 
% for tg = 1:length(target_pwr_range)
% disp(['w=', num2str(w), '; tg=', num2str(tg)]);
% target_pwr = target_pwr_range(tg);
% 
% %% One state: number of coal plants commited
% J_star         = nan*ones(coal_num, N); % [x]x[t]; Cummulative cost
% id_ratio       = zeros(coal_num, N);
% id_dispatch    = zeros(coal_num, N);
% cmt_dispatch   = nan*ones(coal_num, N); % Optimal commitment
% coal_dispatch  = -1*ones(coal_num, N);
% f_sum_dispatch = -1*ones(coal_num, N);
% f_dispatch     = -1*ones(coal_num, coal_num, N); % [14 units]x[x]x[t]
% u_dispatch     = nan*ones(coal_num, coal_num, N);
% v_dispatch     = nan*ones(coal_num, coal_num, N);
% 
% opt_cost_startup  = -1*ones(coal_num, N);
% opt_cost_base_vom = -1*ones(coal_num, N);
% opt_cost_fuel     = -1*ones(coal_num, N);
% opt_cost_ramp     = -1*ones(coal_num, N);
% 
% t = N;
% wind_pwr_tmp = wind_pwr(t)*wind_ratio; % [1x5]
% coal_pwr_tmp = target_pwr - wind_pwr_tmp;
% coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
% 
% id_tmp = interp1(v_range, id_range, coal_pwr_tmp); % [1x5]
% id_tmp = ceil(id_tmp);
% 
% f_sum_tmp = f_table_sum(:,id_tmp); % [x]x[5] = [14x5]
% 
% cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x5]
% cost_fuel_tmp = f_sum_tmp*coal_price; % [14x5]
% J_tmp = cost_fuel_tmp + repmat(cost_base_vom_tmp, coal_num, 1); % [x]x[5] = [14x5]
% [J_N, id_opt] = min(J_tmp, [], 2); % [14x1]
% J_star(:,t) = J_N; % [14x1]
% id_ratio(:,t) = id_opt; % [14x1]
% 
% id_r = 1:coal_num;
% id_c = id_opt';
% id_ex2d = sub2ind(size(J_tmp), id_r, id_c);
% id_dispatch(:,t) = id_tmp(id_c);
% cmt_dispatch(:,t) = id_r;
% coal_dispatch(:,t) = coal_pwr_tmp(id_c);
% f_sum_dispatch(:,t) = f_sum_tmp(id_ex2d);
% 
% id_rw = 1:coal_num;
% id_cl = 1:coal_num;
% [cc,rr] = meshgrid(id_rw, id_cl);
% id_pg = repmat(id_tmp(id_c), coal_num, 1);
% id_ex3d = sub2ind(size(f_table), rr, cc, id_pg);
% f_dispatch(:,:,t) = f_table(id_ex3d);
% u_dispatch(:,:,t) = u_table(id_ex3d);
% v_dispatch(:,:,t) = v_table(id_ex3d);
% 
% opt_cost_startup(:,t) = 0;
% opt_cost_base_vom(:,t) = cost_base_vom_tmp(id_c);
% opt_cost_fuel(:,t) = cost_fuel_tmp(id_ex2d);
% opt_cost_ramp(:,t) = 0;
% 
% for t = N-1:-1:1
%     wind_pwr_tmp = wind_pwr(t)*wind_ratio; % [1x5]
%     coal_pwr_tmp = target_pwr - wind_pwr_tmp;
%     coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
%     
%     id_tmp = interp1(v_range, id_range, coal_pwr_tmp); % [1x5]
%     id_tmp = ceil(id_tmp);
%     
%     f_sum_tmp = f_table_sum(:,id_tmp); % [x]x[5] = [14x5]
%     f_tmp = f_table(:,:,id_tmp); % [14 units]x[x]x[5] = [14x14x5]
%     u_tmp = u_table(:,:,id_tmp); % [14 units]x[x]x[5] = [14x14x5]
%     v_tmp = v_table(:,:,id_tmp); % [14 units]x[x]x[5] = [14x14x5]
% 
%     cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x5]
%     cost_base_vom_tmp = repmat(cost_base_vom_tmp, coal_num, 1); % % [x]x[5] = [14x5]
%     cost_base_vom_tmp = reshape(cost_base_vom_tmp, coal_num, 1, wind_num); % [14x5]->[14x1x5]
%     cost_base_vom_tmp = repmat(cost_base_vom_tmp, 1, coal_num, 1); % [14x1x5]->[14x14x5]
%     
%     cost_fuel_tmp = f_sum_tmp*coal_price; % [x]x[5] = [14x5]
%     cost_fuel_tmp = reshape(cost_fuel_tmp, coal_num, 1, wind_num); % [14x5]->[14x1x5]
%     cost_fuel_tmp = repmat(cost_fuel_tmp, 1, coal_num, 1); % [14x1x5]->[14x14x5]
%     
%     u = 1:coal_num; % 14 different commitments
%     x = (1:coal_num)';
%     [uu,xx] = meshgrid(x,u);
%     c = uu - xx;
%     c(c<0) = 0;
%     cost_startup_tmp = c*coal_startup_cost; % [x]x[u] = [14x14]
%     cost_startup_tmp = repmat(cost_startup_tmp, 1, 1, wind_num); % [x]x[u]x[u2] = [14x14x5]
%     
%     u_rep = reshape(u_tmp,14,14,1,wind_num); % [14x14x5]->[14x14x1x5]
%     u_rep2 = repmat(u_rep,1,1,14,1); % [14x14x14x5] = [14 units]x[x]x[u1]x[u2]
% 
%     u_flip = reshape(u_dispatch(:,:,t+1), 14,1,14); % [14x14]->[14x1x14]
%     u_flip1 = repmat(u_flip, 1,14,1); % [14x1x14]->[14x14x14]
%     u_flip2 = repmat(u_flip1,1,1,1, wind_num); % [14x14x14]->[14x14x14x5]
%     
%     d_coal_pwr = u_rep2 - u_flip2;
%     d_coal_pctg = d_coal_pwr/coal_nameplate;
%     x = abs(d_coal_pctg);
%     y = (x-0.3)*6.5/0.7+1.5;
%     ramp_scale = ones(size(x));
%     ramp_scale(x>0.3) = y(x>0.3);
%     cost_ramp_tmp = squeeze(sum(abs(d_coal_pwr).*ramp_scale * coal_loadfollow)); % [state]x[u1]x[u2]=[14x14x5]
%     
%     JJ = repmat(J_star(:,t+1)',coal_num,1,wind_num);
%     
%     J_tmp = repmat(J_star(:,t+1)',coal_num,1,wind_num) + ... % [14x14x5]
%             cost_fuel_tmp + ...
%             cost_base_vom_tmp + ...
%             cost_startup_tmp + ...
%             cost_ramp_tmp;
%     [value, id_opt] = min(J_tmp(:,:),[],2); 
%     J_star(:,t) = value;
% 
%     [id_r, id_c] = ind2sub(size(f_sum_tmp), id_opt);
%     id_dispatch(:,t) = id_tmp(id_c);
%     cmt_dispatch(:,t) = id_r;
%     coal_dispatch(:,t) = coal_pwr_tmp(id_c);
%     f_sum_dispatch(:,t) = f_sum_tmp(id_opt);
%     
%     id_rw = 1:coal_num;
%     id_cl = 1:coal_num;
%     [cc,rr] = meshgrid(id_rw, id_cl);
%     id_pg = repmat(id_tmp(id_c), coal_num, 1);
%     id_ex3d = sub2ind(size(f_table), rr, cc, id_pg);
%     f_dispatch(:,:,t) = f_table(id_ex3d);
%     u_dispatch(:,:,t) = u_table(id_ex3d);
%     v_dispatch(:,:,t) = v_table(id_ex3d);
%     
%     id_row = (1:coal_num)';
%     id_col = id_r;
%     id_page = id_c;
%     id_extract3d = sub2ind(size(cost_startup_tmp), id_row, id_col, id_page);
%     opt_cost_startup(:,t) = cost_startup_tmp(id_extract3d);
%     opt_cost_base_vom(:,t) = cost_base_vom_tmp(id_extract3d);
%     opt_cost_fuel(:,t) = cost_fuel_tmp(id_extract3d);
%     opt_cost_ramp(:,t) = cost_ramp_tmp(id_extract3d);
% end
% wind_dispatch = target_pwr - coal_dispatch;
% wind_curtail = repmat(wind_pwr, coal_num, 1) - wind_dispatch;
% 
% %% Extract one trajectory
% J_sim             = nan*ones(1,N);
% cmt_sim           = nan*ones(1,N);
% coal_dispatch_sim = nan*ones(1,N);
% f_sum_sim         = nan*ones(1,N);
% f_sim             = nan*ones(coal_num,N);
% u_sim             = nan*ones(coal_num,N);
% v_sim             = nan*ones(coal_num,N);
% cost_startup_sim  = nan*ones(1,N);
% cost_base_vom_sim = nan*ones(1,N);
% cost_fuel_sim     = nan*ones(1,N);
% cost_ramp_sim     = nan*ones(1,N);
% for t = 1:N
%     if t==1
%         [value, id] = min(J_star(:,t));
%         cmt_sim(t) = id;
%     else
%         cmt_sim(t) = cmt_dispatch(cmt_sim(t-1),t-1);
%     end
%     J_sim(t) = J_star(cmt_sim(t),t);
%     coal_dispatch_sim(t) = coal_dispatch(cmt_sim(t),t);
%     f_sum_sim(t) = f_sum_dispatch(cmt_sim(t),t);
%     f_sim(:,t) = f_dispatch(:,cmt_sim(t),t);
%     u_sim(:,t) = u_dispatch(:,cmt_sim(t),t);
%     v_sim(:,t) = v_dispatch(:,cmt_sim(t),t);
%     cost_startup_sim(t) = opt_cost_startup(cmt_sim(t),t);
%     cost_base_vom_sim(t) = opt_cost_base_vom(cmt_sim(t),t);
%     cost_fuel_sim(t) = opt_cost_fuel(cmt_sim(t),t);
%     cost_ramp_sim(t) = opt_cost_ramp(cmt_sim(t),t);
% end
% wind_dispatch_sim = target_pwr - coal_dispatch_sim;
% wind_curtail_sim = wind_pwr - wind_dispatch_sim;
% 
% %% Check ramping
% d_coal_pwr = [zeros(1,coal_num);
%               diff(u_sim')]'; % [14]x[8760]
% d_coal_pctg = d_coal_pwr/coal_nameplate;
% 
% % Check changes in commitment
% d_cmt = [0, diff(cmt_sim)]; % [1]x[t]
% cmt_c = zeros(1,N); % Commition
% cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
% cmt_d = zeros(1,N); % Decommition
% cmt_d(d_cmt<0) = d_cmt(d_cmt<0);
% 
% % Commit or decommit 2 or more coal plants
% id_jump = find(abs(d_cmt)>1) - 1;
% 
% %% Cost
% cost_startup = sum(cost_startup_sim);
% cost_base_vom = sum(cost_base_vom_sim);
% cost_fuel = sum(cost_fuel_sim);
% cost_ramp = sum(cost_ramp_sim);
% cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp;
% 
% %% ====================
% table_coal_v(w,tg) = sum(v_sim(:));
% table_coal_u(w,tg) = sum(u_sim(:));
% table_coal_dispatch(w,tg)= sum(coal_dispatch_sim);
% table_wind_dispatch(w,tg)= sum(wind_dispatch_sim);
% table_wind_curtail(w,tg) = sum(wind_curtail_sim);
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
%      'table_coal_v', 'table_coal_u', ...
%      'table_coal_dispatch', 'table_wind_dispatch', 'table_wind_curtail', ...
%      'table_cost_startup', 'table_cost_base_vom', 'table_cost_fuel', 'table_cost_ramp', 'table_cost_total');

load(save_name);

pctg_startup = table_cost_startup./table_cost_total;
pctg_base_vom = table_cost_base_vom./table_cost_total;
pctg_fuel = table_cost_fuel./table_cost_total;
pctg_ramp = table_cost_ramp./table_cost_total;
pctg_wind_curtail = table_wind_curtail./(table_wind_dispatch + table_wind_curtail);


%%
plot_switch = 'on'; % on/off
switch plot_switch
    case 'on'
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

%% Wind curtailment
figure(3); clf;
surf(wind_pwr_range, target_pwr_range, pctg_wind_curtail');
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Wind Curtailment (-)');
ylim([2000 10000]);
zlim([0 0.4]);
set(gca, 'clim', [0 0.4]);
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
end
