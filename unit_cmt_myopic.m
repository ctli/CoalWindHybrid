clear
close all
clc
format compact

load FourteenUnits; % 'v_range', 'f_table', 'u_table', 'v_table', 'id_st', 'id_ed', 'v_st', 'v_ed'
id_range = 1:length(v_range);

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
u_min = coal_nameplate*0.4; % 264 [MW]
v_min = min(v_range); % 220 [MW]; 0.33 of capacity
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
wind_pwr = round(p*2500)';

wind_ratio = linspace(1,1,3); save_name = ['Myopic_', wind_file, '_nominal']; % Not allow economic wind curtialment
% wind_ratio = linspace(1,0,101); save_name = ['Myopic_', wind_file, '_nominal_new']; % Allow economic wind curtialment

target_pwr = 8500;


%% Myopic unit commitment
% cmt_myopic = zeros(1,length(v_range));
% f_sum_myopic = zeros(1,length(v_range));
% f_myopic = zeros(coal_num,length(v_range));
% v_myopic = zeros(coal_num,length(v_range));
% u_myopic = zeros(coal_num,length(v_range));
% tic;
% for vv = 1:length(v_range)
%     f_column = f_table_sum(:,vv);
%     [f_sum_myopic(vv), cmt_myopic(vv)] = min(f_column);
%     f_myopic(:,vv) = f_table(:,cmt_myopic(vv),vv);
%     v_myopic(:,vv) = v_table(:,cmt_myopic(vv),vv);
%     u_myopic(:,vv) = u_table(:,cmt_myopic(vv),vv);
% end
% toc;
% % save('MyopicDispatch', 'v_range', 'cmt_myopic', 'f_sum_myopic', 'f_myopic', 'v_myopic', 'u_myopic');
% 
% % ====================
% figure(1); clf; hold on; box on; % Fuel consumption
% plot(0,0,'x');
% plot(v_range, f_table_sum, 'linewidth', 1);
% for n = 1:coal_num
%     not_nan = ~isnan(f_table_sum(n,:));
%     id_st = find(not_nan==1, 1, 'first');
%     id_ed = find(not_nan==1, 1, 'last');
%     v_st = v_range(id_st);
%     v_ed = v_range(id_ed);
%     if n==1
%         text(v_ed, max(f_table_sum(n,:)), [' ', num2str(n), ' Unit is commited'], 'fontsize', 7);
%     elseif n<11
%         text(v_ed, max(f_table_sum(n,:)), [' ', num2str(n), ' Units are commited'], 'fontsize', 7);
%     else
%         text(v_st, min(f_table_sum(n,:)), [' ', num2str(n), ' Units are commited '], 'fontsize', 7, 'horizontalalignment', 'right');
%     end
% end
% xlabel('Output Power, MW (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% my_gridline;
% h1 = plot(v_range, f_sum_myopic, 'color', [1 1 1]*0, 'linewidth', 0.35);
% legend(h1, 'Myopic Dispatch');
% set(legend, 'location', 'northwest');
% 
% % ====================
% figure(2); clf; % Unit commitment
% plot(v_range, cmt_myopic);
% ylim([0 15]);
% set(gca, 'ytick', 0:3:15);
% xlabel('Output Power, MW (in-house use excluded)');
% ylabel('Number of Plants Dispatched (Count)');
% my_gridline;

load  MyopicDispatch; % 'v_range', 'cmt_myopic', 'f_sum_myopic', 'f_myopic', 'v_myopic', 'u_myopic'


%% Myopic dispatch
% id_ratio = zeros(1, length(wind_pwr));
% id_dispatch = zeros(1, length(wind_pwr));
% coal_dispatch = zeros(1, length(wind_pwr));
% tic;
% for t = 1:length(wind_pwr)
%     wind_pwr_tmp = wind_pwr(t)*wind_ratio;
%     coal_pwr_tmp = target_pwr - wind_pwr_tmp;
%     coal_pwr_tmp(coal_pwr_tmp<v_min) = v_min;
%     
%     id_tmp = interp1(v_range, id_range, coal_pwr_tmp);
%     id_tmp = ceil(id_tmp);
%     f_tmp = interp1(id_range, f_sum_myopic, id_tmp);
%     
%     cost_base_vom = coal_pwr_tmp*coal_baseload;
%     cost_fuel = f_tmp*coal_price;
%     [value, id_opt] = min(cost_base_vom + cost_fuel);
%     id_ratio(t) = id_opt;
%     id_dispatch(t) = id_tmp(id_opt);
%     coal_dispatch(t) = coal_pwr_tmp(id_opt);
% end
% toc;
% cmt_dispatch = cmt_myopic(id_dispatch);
% f_sum_dispatch = f_sum_myopic(id_dispatch); % [ton/h]
% f_dispatch = f_myopic(:,id_dispatch);
% v_dispatch = v_myopic(:,id_dispatch);
% u_dispatch = u_myopic(:,id_dispatch);
% wind_dispatch = target_pwr - coal_dispatch;
% wind_curtail = wind_pwr - wind_dispatch;

load(['Myopic_', wind_file, '_nominal']);


%% ========================================================================
% Check ramping
d_coal_pwr = [zeros(1,coal_num); diff(u_dispatch')]'; % [8760]x[14]
d_coal_pctg = d_coal_pwr/coal_nameplate;

% Check changes in commitment
d_cmt = [0, diff(cmt_dispatch)];
cmt_c = zeros(1, length(wind_pwr)); % Commition
cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
cmt_d = zeros(1, length(wind_pwr)); % Decommition
cmt_d(d_cmt<0) = d_cmt(d_cmt<0);

% Commit or decommit 2 or more coal plants
id_jump = find(abs(d_cmt)>1) - 1;


%% Costs
opt_cost_startup = cmt_c * coal_startup_cost;
opt_cost_base_vom = coal_dispatch * coal_baseload;
opt_cost_fuel = f_sum_dispatch * coal_price;

% opt_cost_ramp = abs(d_coal_pwr(:)) * coal_loadfollow;
x = abs(d_coal_pctg);
y = (x-0.3)*6.5/0.7+1.5;
ramp_scale = ones(size(x)); % Exceeding 30 precent nameplate ramping is scaled by 1.5-8 times
ramp_scale(x>0.3) = y(x>0.3);
opt_cost_ramp = abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow; % [$]

cost_startup = sum(opt_cost_startup)
cost_base_vom = sum(opt_cost_base_vom)
cost_fuel = sum(opt_cost_fuel)
cost_ramp = sum(opt_cost_ramp)

disp('====================');
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp

disp('====================');
pctg_startup = cost_startup/cost_total
pctg_base_vom = cost_base_vom/cost_total
pctg_fuel = cost_fuel/cost_total
pctg_ramp = cost_ramp/cost_total

% save(save_name, ...
%      'id_ratio', 'id_dispatch', 'cmt_dispatch', 'f_sum_dispatch', ...
%      'f_dispatch', 'v_dispatch', 'u_dispatch', ...
%      'coal_dispatch', 'wind_dispatch', 'wind_curtail', 'wind_pwr', ...
%      'opt_cost_startup', 'opt_cost_base_vom', 'opt_cost_fuel', 'opt_cost_ramp', ...
%      'cost_startup', 'cost_base_vom', 'cost_fuel', 'cost_ramp', 'cost_total');


%% ========================================================================
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
v_unique = zeros(1, length(wind_pwr));
u_unique = zeros(1, length(wind_pwr));
for t = 1:length(wind_pwr)
    vt = v_dispatch(1:cmt_dispatch(t),t);
    if length(unique(vt))>1
        v_unique(t) = vt(1);
    else
        v_unique(t) = unique(vt);
    end

    ut = u_dispatch(1:cmt_dispatch(t),t);
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
ha = area([coal_dispatch; wind_dispatch; wind_curtail]', 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.7 0.7]);
set(ha(2), 'facec', [0.6 1 0]);
set(ha(3), 'facec', [0 0.7 0]);

area(v_dispatch', 'facec', 'none', 'edgecolor', [1 1 1]);
for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], coal_dispatch(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end

ylim([0 9000]);
xlim([0 24*7*4]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 0:1500:9000);
set(gca, 'layer', 'top');
ylabel('Total Output Power (MW)');
title('Time Window: 4 Weeks');
% legend('Coal', 'Wind', 'Wind Curtailment');
% set(legend, 'location', 'northwest');

% ====================
ax2 = subplot(4,1,3); hold on; box on;
plot(cmt_dispatch, 'color', [1 1 1]*0.6, 'linewidth', 1);
for i = 1:length(id_jump)
h = plot(id_jump(i)+[0 1], cmt_dispatch(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end
legend(h, 'Steep Change');
set(legend, 'location', 'southwest');

xlim([0 168*7]);
ylim([min(cmt_dispatch)-0.5, max(cmt_dispatch)+0.5]);
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

% export_fig myopic -painters

end

