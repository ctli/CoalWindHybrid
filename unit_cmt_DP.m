clear
close all
clc
format compact

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
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
wind_pwr = p*2500;

target_pwr = 8500;
coal_pwr = target_pwr - wind_pwr; % Use coal to make up deficit

wind_curtail = (coal_pwr + wind_pwr) - target_pwr;


%% Horizon-based unit commitment
load FourteenUnits; % 'v_range', 'f_table', 'u_table', 'v_table', 'id_st', 'id_ed', 'v_st', 'v_ed'

N = length(coal_pwr);

% One state: number of coal plants commited
J_star = -1*ones(1, length(coal_pwr)); % Cummulative cost
cmt_star = -1*ones(1, length(coal_pwr)); % Optimal commitment
u_dispatch = -1*ones(coal_num, length(coal_pwr));
v_dispatch = -1*ones(coal_num, length(coal_pwr));

opt_cost_startup = -1*ones(1, length(coal_pwr));
opt_cost_base_vom = -1*ones(1, length(coal_pwr));
opt_cost_fuel = -1*ones(1, length(coal_pwr));
opt_cost_ramp = -1*ones(1, length(coal_pwr));

t = N;
id_dispatch(t) = find(v_range>=coal_pwr(t), 1, 'first');
f_extract = f_table(id_dispatch(t),:);
cost_fuel_tmp = f_extract * coal_price;
cost_base_vom_tmp = coal_pwr(t) * coal_baseload;
J_tmp = cost_fuel_tmp + cost_base_vom_tmp;
[value, id_opt] = min(J_tmp);
J_star(t) = value;
cmt_star(t) = id_opt;
u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_star(t));
v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_star(t));
opt_cost_startup(t) = 0;
opt_cost_base_vom(t) = cost_base_vom_tmp;
opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
opt_cost_ramp(t) = 0;

tic;
for t = N-1:-1:1
    id_dispatch(t) = find(v_range>=coal_pwr(t), 1, 'first');
    f_extract = f_table(id_dispatch(t),:);
    cost_fuel_tmp = f_extract * coal_price;
    cost_base_vom_tmp = coal_pwr(t) * coal_baseload;
    
    cost_startup_tmp = zeros(1, coal_num);
    cost_ramp_tmp = zeros(1, coal_num);
    for cmt = 1:coal_num
        u_tmp = u_table(:,id_dispatch(t),cmt);
        v_tmp = v_table(:,id_dispatch(t),cmt);
        
        c = cmt - cmt_star(t+1);
        if c>0, cost_startup_tmp(cmt) = c*coal_startup_cost;
        else cost_startup_tmp(cmt) = 0;
        end
        
        d_coal_pwr = u_tmp - u_dispatch(:,t+1);
        d_coal_pctg = d_coal_pwr/coal_nameplate;
        x = abs(d_coal_pctg);
        y = (x-0.3)*6.5/0.7+1.5;
        ramp_scale = zeros(size(x));
        ramp_scale(x>0.3) = y(x>0.3);
        cost_ramp_tmp(cmt) = sum(abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow);
    end
    J_tmp = J_star(t+1) + cost_fuel_tmp + cost_base_vom_tmp + cost_startup_tmp + cost_ramp_tmp;
    
    [value, id_opt] = min(J_tmp);
    J_star(t) = value;
    cmt_star(t) = id_opt;
    v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_star(t));
    u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_star(t));
    opt_cost_startup(t) = cost_startup_tmp(id_opt);
    opt_cost_base_vom(t) = cost_base_vom_tmp;
    opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
    opt_cost_ramp(t) = cost_ramp_tmp(id_opt);
end
toc;

cost_startup = sum(opt_cost_startup)
cost_base_vom = sum(opt_cost_base_vom)
cost_fuel = sum(opt_cost_fuel)
cost_ramp = sum(opt_cost_ramp)
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp

% WHY coal_total not equal to J_star(1)???
J_star(1)

% save(['DP_', wind_file, '_nominal'], ...
%      'id_dispatch', 'v_dispatch', 'u_dispatch', ...
%      'wind_pwr', 'wind_curtail', ...
%      'cost_startup', 'cost_base_vom', 'cost_fuel', 'cost_ramp', 'cost_total');



