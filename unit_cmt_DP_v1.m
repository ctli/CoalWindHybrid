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
wind_pwr = round(p*2500)'; % [1x8760]

% wind_ratio = linspace(1,1,3); save_name = ['DP_', wind_file, '_debug']; % Not allow economic wind curtialment
% wind_ratio = linspace(1,0,5); save_name = ['DP_', wind_file, '_debug_new'];  % Allow economic wind curtialment

target_pwr = 8500;


%% Horizon-based unit commitment
N = length(wind_pwr);

% One state: number of coal plants commited
J_star         = nan*ones(coal_num, length(wind_pwr)); % Cummulative cost
id_ratio       = zeros(1, length(wind_pwr));
id_dispatch    = zeros(1, length(wind_pwr));
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
wind_pwr_tmp = wind_pwr(t)*1; % [1x1]
coal_pwr_tmp = target_pwr - wind_pwr_tmp;
coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
coal_dispatch(:,t) = coal_pwr_tmp;

id_tmp = interp1(v_range, id_range, coal_pwr_tmp); % [1x1]
id_tmp = ceil(id_tmp);
id_dispatch(t) = id_tmp;

f_tmp = f_table_sum(:,id_tmp); % [x]x[1] = [14x1]

cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x1]
cost_fuel_tmp = f_tmp*coal_price; % [14x1]
J_N = cost_fuel_tmp + cost_base_vom_tmp; % [14x1]
J_star(:,t) = J_N; % [1x14]

cmt_dispatch(:,t) = 1:coal_num;
f_sum_dispatch(:,t) = f_tmp;

f_dispatch(:,:,t) = f_table(:,:,id_dispatch(t));
u_dispatch(:,:,t) = u_table(:,:,id_dispatch(t));
v_dispatch(:,:,t) = v_table(:,:,id_dispatch(t));

opt_cost_startup(:,t) = 0;
opt_cost_base_vom(:,t) = cost_base_vom_tmp;
opt_cost_fuel(:,t) = cost_fuel_tmp;
opt_cost_ramp(:,t) = 0;

tic;
for t = N-1:-1:1
    wind_pwr_tmp = wind_pwr(t)*1;
    coal_pwr_tmp = target_pwr - wind_pwr_tmp;
    coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
    coal_dispatch(:,t) = coal_pwr_tmp;

    id_tmp = interp1(v_range, id_range, coal_pwr_tmp);
    id_tmp = ceil(id_tmp);
    id_dispatch(t) = id_tmp;

    f_sum_tmp = f_table_sum(:,id_tmp); % [x]x[1] = [14x1]

    cost_base_vom_tmp = coal_pwr_tmp * coal_baseload; % [1x1]
    cost_fuel_tmp = f_sum_tmp*coal_price; % [14x1]
    
    cost_ramp_tmp = zeros(1, coal_num); % [1]x[x] = [1x14]
    for xx = 1:coal_num
        if ~isnan(f_sum_tmp(xx))
            f_tmp = f_table(:,xx,id_tmp);
            u_tmp = u_table(:,xx,id_tmp); % [14x1]
            v_tmp = v_table(:,xx,id_tmp);
            
            u = 1:coal_num; % 14 different commitments
            c = u - xx;
            c(c<0) = 0;
            cost_startup_tmp = c*coal_startup_cost; % [1]x[u] = [1x14]
            
            d_coal_pwr = repmat(u_tmp, 1, coal_num) - u_dispatch(:,:,t+1); % [14x14]=[14 units]x[u]
            d_coal_pctg = d_coal_pwr/coal_nameplate;
            x = abs(d_coal_pctg);
            y = (x-0.3)*6.5/0.7+1.5;
            ramp_scale = ones(size(x));
            ramp_scale(x>0.3) = y(x>0.3);
            cost_ramp_tmp = sum(abs(d_coal_pwr).*ramp_scale * coal_loadfollow); % [1]x[u] = [1x14]
            
            J_tmp = J_star(:,t+1)' + cost_fuel_tmp(xx) + cost_base_vom_tmp + cost_startup_tmp + cost_ramp_tmp;
            if ~all(isnan(J_tmp)) % at least some controls are legit
                [value, id_opt] = min(J_tmp);
                J_star(xx,t) = value;
                cmt_dispatch(xx,t) = id_opt;
                f_sum_dispatch(xx,t) = f_sum_tmp(xx);
                
                f_dispatch(:,xx,t) = f_tmp;
                u_dispatch(:,xx,t) = u_tmp;
                v_dispatch(:,xx,t) = v_tmp;
                
                opt_cost_startup(xx,t) = cost_startup_tmp(id_opt);
                opt_cost_base_vom(xx,t) = cost_base_vom_tmp;
                opt_cost_fuel(xx,t) = cost_fuel_tmp(xx);
                opt_cost_ramp(xx,t) = cost_ramp_tmp(id_opt);
            end
        end
    end
end
toc;
wind_dispatch = target_pwr - coal_dispatch;
wind_curtail = repmat(wind_pwr, coal_num, 1) - wind_dispatch;


%% Extract one trajectory
J_sim = zeros(1,N);
f_sum_sim = zeros(1,N);
cmt_sim = zeros(1,N);
f_sim = zeros(coal_num,N);
u_sim = zeros(coal_num,N);
v_sim = zeros(coal_num,N);
cost_startup_sim  = -1*ones(1,N);
cost_base_vom_sim = -1*ones(1,N);
cost_fuel_sim     = -1*ones(1,N);
cost_ramp_sim     = -1*ones(1,N);

coal_dispatch_sim = -1*ones(1,N);
wind_dispatch_sim = -1*ones(1,N);

for t = 1:N
    if t==1
        [value, id] = min(J_star(:,t));
        cmt_sim(t) = id;
    else
        cmt_sim(t) = cmt_dispatch(cmt_sim(t-1),t-1);
    end
    J_sim(t) = J_star(cmt_sim(t),t);
    f_sum_sim(t) = f_sum_dispatch(cmt_sim(t),t);
    f_sim(:,t) = f_dispatch(:,cmt_sim(t),t);
    u_sim(:,t) = u_dispatch(:,cmt_sim(t),t);
    v_sim(:,t) = v_dispatch(:,cmt_sim(t),t);
    cost_startup_sim(t) = opt_cost_startup(cmt_sim(t),t);
    cost_base_vom_sim(t) = opt_cost_base_vom(cmt_sim(t),t);
    cost_fuel_sim(t) = opt_cost_fuel(cmt_sim(t),t);
    cost_ramp_sim(t) = opt_cost_ramp(cmt_sim(t),t);
    
    coal_dispatch_sim(t) = coal_dispatch(cmt_sim(t),t);
    wind_dispatch_sim(t) = wind_dispatch(cmt_sim(t),t);
end
wind_curtail_sim = wind_pwr - wind_dispatch_sim;


%% ========================================================================
% Check ramping
d_coal_pwr = [zeros(1,coal_num);
              diff(u_sim')]'; % [14]x[8760]
d_coal_pctg = d_coal_pwr/coal_nameplate;

% Check changes in commitment
d_cmt = [0, diff(cmt_sim)]; % [1]x[t]
cmt_c = zeros(1,N); % Commition
cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
cmt_d = zeros(1,N); % Decommition
cmt_d(d_cmt<0) = d_cmt(d_cmt<0);

% Commit or decommit 2 or more coal plants
id_jump = find(abs(d_cmt)>1) - 1;

% Double check costs (not identical, but extremely close)
cost_startup = sum(cost_startup_sim)
cost_base_vom = sum(cost_base_vom_sim)
cost_fuel = sum(cost_fuel_sim)
cost_ramp = sum(cost_ramp_sim)
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp

dJ_pctg = (J_sim(1) - cost_total)/cost_total

save('DP_v1', ...
     'J_star', 'id_ratio', 'id_dispatch', ...
     'cmt_dispatch', 'f_sum_dispatch', 'f_dispatch', 'u_dispatch', 'v_dispatch', ...
     'coal_dispatch', 'wind_dispatch', 'wind_pwr', 'wind_curtail', ...
     'opt_cost_startup', 'opt_cost_base_vom', 'opt_cost_fuel', 'opt_cost_ramp', ...
     'J_sim', ...
     'cmt_sim', 'f_sim', 'u_sim', 'v_sim', ...
     'coal_dispatch_sim', 'wind_dispatch_sim', 'wind_curtail_sim', ...
     'cost_startup', 'cost_base_vom', 'cost_fuel', 'cost_ramp', 'cost_total');


%% ========================================================================
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
v_unique = zeros(1, length(wind_pwr));
u_unique = zeros(1, length(wind_pwr));
for t = 1:length(wind_pwr)
    vt = v_sim(:,t);
    if length(unique(vt))>1, 
        v_unique(t) = vt(1);
    else
        v_unique(t) = unique(vt);
    end

    ut = u_sim(:,t);
    if length(unique(ut))>1
        u_unique(t) = ut(1);
    else
        u_unique(t) = unique(ut);
    end
end
%% ==============================
figure(3); clf; hold on;
set(gcf, 'units', 'inch', 'pos', [2.9792    1.4583    5.8333    6.75]);

% ====================
ax1 = subplot(4,1,1:2); hold on; box on;
ha = area([coal_dispatch_sim; wind_dispatch_sim; wind_curtail_sim]', 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.7 0.7]);
set(ha(2), 'facec', [0.6 1 0]);
set(ha(3), 'facec', [0 0.7 0]);

area(v_sim', 'facec', 'none', 'edgecolor', [1 1 1]);
for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], coal_dispatch_sim(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end

ylim([0 9000]);
xlim([0 24*7*4]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 0:1500:9000);
set(gca, 'layer', 'top');
ylabel('Total Output Power (MW)');
title('Time Window: 4 Weeks');
% ylim([8000 9000]);
% set(gca, 'ytick', 8000:200:9000);
xlim([0 30]);

% ====================
ax2 = subplot(4,1,3); hold on; box on;
plot(cmt_sim, 'color', [1 1 1]*0.6, 'linewidth', 1);
if ~isempty(id_jump)
for i = 1:length(id_jump)
h = plot(id_jump(i)+[0 1], cmt_sim(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end
legend(h, 'Steep Change');
set(legend, 'location', 'southwest');
end

xlim([0 168*7]);
% ylim([9.5 14.5]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 1:1:14);
ylabel({'Unit Commited', '(count)'});
my_gridline;

% ====================
ax3 = subplot(4,1,4); hold on; box on;
plot(u_unique);
plot(v_unique);

ylabel({'Output Power of One', 'Coal Plant (MW)'});
ylim([0 690]);
set(gca, 'ytick', 0:330/2:660);
set(gca, 'xtick', 0:168:24*7*4);
xlabel('Time (hr)');
legend('Generation', 'Useful Output (Exclude In-House Use)');
set(legend, 'location', 'southwest');
my_gridline;
linkaxes([ax1, ax2, ax3], 'x');

% export_fig DP -painters
% export_fig DP_new -painters

end


