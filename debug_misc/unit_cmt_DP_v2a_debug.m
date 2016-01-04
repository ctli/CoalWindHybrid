% No parallelization; with economic wind curtailment

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
wind_pwr = round(p*5000)'; % [1x8760]
wind_pwr = wind_pwr(1:6);

% wind_ratio = linspace(1,1,3); save_name = ['DP_', wind_file, '_debug']; % Not allow economic wind curtialment
wind_ratio = linspace(1,0,5); save_name = ['DP_', wind_file, '_debug_new'];  % Allow economic wind curtialment
wind_num = length(wind_ratio);

target_pwr = 2500;


%% Horizon-based unit commitment
N = length(wind_pwr);

% One state: number of coal plants commited
J_star         = nan*ones(coal_num, length(wind_pwr)); % Cummulative cost
id_ratio       = zeros(coal_num, length(wind_pwr));
id_dispatch    = zeros(coal_num, length(wind_pwr));
cmt_dispatch   = nan*ones(coal_num, length(wind_pwr)); % Optimal commitment
coal_dispatch  = -1*ones(coal_num, length(wind_pwr));
f_sum_dispatch = -1*ones(coal_num, length(wind_pwr));
f_dispatch     = -1*ones(coal_num, coal_num, length(wind_pwr));
u_dispatch     = nan*ones(coal_num, coal_num, length(wind_pwr));
v_dispatch     = nan*ones(coal_num, coal_num, length(wind_pwr));

opt_cost_startup  = -1*ones(coal_num, length(wind_pwr));
opt_cost_base_vom = -1*ones(coal_num, length(wind_pwr));
opt_cost_fuel     = -1*ones(coal_num, length(wind_pwr));
opt_cost_ramp     = -1*ones(coal_num, length(wind_pwr));

t = N;
wind_pwr_tmp = wind_pwr(t)*wind_ratio; % [1x5]
coal_pwr_tmp = target_pwr - wind_pwr_tmp;
coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;

id_tmp = interp1(v_range, id_range, coal_pwr_tmp); % [1x5]
id_tmp = ceil(id_tmp);

f_sum_tmp = f_table_sum(:,id_tmp); % [x]x[5] = [14x5]

cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x5]
cost_fuel_tmp = f_sum_tmp*coal_price; % [14x5]
J_tmp = cost_fuel_tmp + repmat(cost_base_vom_tmp, coal_num, 1); % [x]x[5] = [14x5]
[J_N, id_opt] = min(J_tmp, [], 2); % [14x1]
J_star(:,t) = J_N; % [14x1]
id_ratio(:,t) = id_opt; % [14x1]

id_row = 1:coal_num;
id_col = id_opt';
id_extract = sub2ind(size(J_tmp), id_row, id_col);
id_dispatch(:,t) = id_tmp(id_col);
cmt_dispatch(:,t) = id_row;
coal_dispatch(:,t) = coal_pwr_tmp(id_col);
f_sum_dispatch(:,t) = f_sum_tmp(id_extract);

id_rw = 1:coal_num;
id_cl = 1:coal_num;
[cc,rr] = meshgrid(id_rw, id_cl);
id_pg = repmat(id_tmp(id_col), coal_num, 1);
id_extract3d = sub2ind(size(f_table), rr, cc, id_pg);
f_dispatch(:,:,t) = f_table(id_extract3d);
u_dispatch(:,:,t) = u_table(id_extract3d);
v_dispatch(:,:,t) = v_table(id_extract3d);

opt_cost_startup(:,t) = 0;
opt_cost_base_vom(:,t) = cost_base_vom_tmp(id_col);
opt_cost_fuel(:,t) = cost_fuel_tmp(id_extract);
opt_cost_ramp(:,t) = 0;

tic;
for t = N-1:-1:N-3%1
    J_star_w         = nan*ones(coal_num, wind_num); % [14x5]
    cmt_dispatch_w   = nan*ones(coal_num, wind_num); % [14x5]
    coal_pwr_tmp_w   = nan*ones(coal_num, wind_num); % [14x5]
    f_sum_dispatch_w = nan*ones(coal_num, wind_num); % [14x5]
    f_dispatch_w = nan*ones(coal_num, coal_num, wind_num); % [14x14x5]
    u_dispatch_w = nan*ones(coal_num, coal_num, wind_num); % [14x14x5]
    v_dispatch_w = nan*ones(coal_num, coal_num, wind_num); % [14x14x5]
    
    opt_cost_startup_w  = nan*ones(coal_num, wind_num); % [14x5]
    opt_cost_base_vom_w = nan*ones(coal_num, wind_num); % [14x5]
    opt_cost_fuel_w     = nan*ones(coal_num, wind_num); % [14x5]
    opt_cost_ramp_w     = nan*ones(coal_num, wind_num); % [14x5]
    
    cost_base_vom_tmp_wx = nan*ones(coal_num,coal_num,wind_num); % [14x14x5]=[x]x[u1]x[u2]
    cost_fuel_tmp_wx     = nan*ones(coal_num,coal_num,wind_num); % [14x14x5]
    cost_startup_tmp_wx  = nan*ones(coal_num,coal_num,wind_num); % [14x14x5]
    cost_ramp_tmp_wx     = nan*ones(coal_num,coal_num,wind_num); % [14x14x5]
    J_tmp_wx             = nan*ones(coal_num,coal_num,wind_num); % [14x14x5]
    u_rep2_wx  = nan*ones(14,14,14,5);
    u_flip2_wx = nan*ones(14,14,14,5);

    for w = 1:wind_num;
    wind_pwr_tmp = wind_pwr(t)*wind_ratio(w); % [1x1]
    coal_pwr_tmp = target_pwr - wind_pwr_tmp;
    coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
    coal_pwr_tmp_w(:,w) = coal_pwr_tmp;
    
    id_tmp = interp1(v_range, id_range, coal_pwr_tmp); % [1x1]
    id_tmp = ceil(id_tmp);

    f_sum_tmp = f_table_sum(:,id_tmp); % [x]x[1] = [14x1]

    cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x1]
    cost_fuel_tmp = f_sum_tmp*coal_price; % [x]x[1] = [14x1]
    
    for xx = 1:coal_num
%         if ~isnan(f_sum_tmp(xx))
            f_tmp = f_table(:,xx,id_tmp); % [14 units]x[1]
            u_tmp = u_table(:,xx,id_tmp); % [14 units]x[1]
            v_tmp = v_table(:,xx,id_tmp); % [14 units]x[1]
            
            u = 1:coal_num; % 14 different commitments
            c = u - xx;
            c(c<0) = 0;
            cost_startup_tmp = c*coal_startup_cost; % [1]x[u] = [1x14]
            cost_startup_tmp_wx(xx,:,w) = cost_startup_tmp;
            
            u_rep2_wx(:,xx,:,w) = repmat(u_tmp, 1, coal_num);
            u_flip2_wx(:,xx,:,w) = u_dispatch(:,:,t+1);

            d_coal_pwr = repmat(u_tmp, 1, coal_num) - u_dispatch(:,:,t+1); % [14x14]=[14 units]x[u]
            d_coal_pctg = d_coal_pwr/coal_nameplate;
            x = abs(d_coal_pctg);
            y = (x-0.3)*6.5/0.7+1.5;
            ramp_scale = ones(size(x));
            ramp_scale(x>0.3) = y(x>0.3);
            cost_ramp_tmp = sum(abs(d_coal_pwr).*ramp_scale * coal_loadfollow); % [1]x[u] = [1x14]
            cost_ramp_tmp_wx(xx,:,w) = cost_ramp_tmp;
            
            cost_fuel_tmp_wx(xx,:,w) = cost_fuel_tmp(xx);
            cost_base_vom_tmp_wx(xx,:,w) = cost_base_vom_tmp;
            J_tmp_wx(xx,:,w) = J_star(:,t+1);
            
            J_tmp = J_star(:,t+1)' + cost_fuel_tmp(xx) + cost_base_vom_tmp + cost_startup_tmp + cost_ramp_tmp; % [1]x[u] = [1x14]
            if ~all(isnan(J_tmp)) % at least some controls are legit
                [value, id_opt] = min(J_tmp);
                J_star_w(xx,w) = value;
                cmt_dispatch_w(xx,w) = id_opt;
                f_sum_dispatch_w(xx,w) = f_sum_tmp(xx);
                f_dispatch_w(:,xx,w) = f_tmp;
                u_dispatch_w(:,xx,w) = u_tmp;
                v_dispatch_w(:,xx,w) = v_tmp;
                
                opt_cost_startup_w(xx,w) = cost_startup_tmp(id_opt);
                opt_cost_base_vom_w(xx,w) = cost_base_vom_tmp;
                opt_cost_fuel_w(xx,w) = cost_fuel_tmp(xx);
                opt_cost_ramp_w(xx,w) = cost_ramp_tmp(id_opt);
            end
%         end
    end
    
    end
    [value, id] = min(J_star_w, [], 2);
    J_star(:,t) = value;
    id_ratio(:,t) = id;
    
    id_row = 1:coal_num;
    id_col = id';
    id_extract = sub2ind(size(f_sum_dispatch_w), id_row, id_col);
    cmt_dispatch(:,t) = cmt_dispatch_w(id_extract);
    coal_dispatch(:,t) = coal_pwr_tmp_w(id_extract);
    f_sum_dispatch(:,t) = f_sum_dispatch_w(id_extract);

    id_rw = 1:coal_num;
    id_cl = 1:coal_num;
    [cc,rr] = meshgrid(id_rw, id_cl);
    id_pg = repmat(id', coal_num, 1);
    id_extract3d = sub2ind(size(f_dispatch_w), rr, cc, id_pg);
    f_dispatch(:,:,t) = f_dispatch_w(id_extract3d);
    u_dispatch(:,:,t) = u_dispatch_w(id_extract3d);
    v_dispatch(:,:,t) = v_dispatch_w(id_extract3d);

    opt_cost_startup(:,t) = opt_cost_startup_w(id_extract);
    opt_cost_base_vom(:,t) = opt_cost_base_vom_w(id_extract);
    opt_cost_fuel(:,t) = opt_cost_fuel_w(id_extract);
    opt_cost_ramp(:,t) = opt_cost_ramp_w(id_extract);
end
toc;
wind_dispatch = target_pwr - coal_dispatch;
wind_curtail = repmat(wind_pwr, coal_num, 1) - wind_dispatch;

u_dispatch_wx = u_dispatch;
save('DP_v3a_debug', ...
     'J_tmp_wx', 'cost_base_vom_tmp_wx', 'cost_fuel_tmp_wx', 'cost_startup_tmp_wx', 'cost_ramp_tmp_wx', ...
     'u_rep2_wx', 'u_flip2_wx', ...
     'u_dispatch_wx');


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
% % load DP_v3a;
% 
% 
% %% ========================================================================
% % Check ramping
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
% % Double check costs (not identical, but extremely close)
% cost_startup = sum(cost_startup_sim)
% cost_base_vom = sum(cost_base_vom_sim)
% cost_fuel = sum(cost_fuel_sim)
% cost_ramp = sum(cost_ramp_sim)
% cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp
% 
% dJ_pctg = (J_sim(1) - cost_total)/cost_total
% 



