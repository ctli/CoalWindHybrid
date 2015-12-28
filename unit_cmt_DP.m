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
wind_pwr = round(p*2500)'; % [1x8760]

target_pwr = 8500;
coal_pwr = target_pwr - wind_pwr; % Use coal to make up deficit

wind_curtail = (coal_pwr + wind_pwr) - target_pwr;


%% Horizon-based unit commitment
% load FourteenUnits; % 'v_range', 'f_table', 'u_table', 'v_table', 'id_st', 'id_ed', 'v_st', 'v_ed'
% 
% N = length(coal_pwr);
% 
% % One state: number of coal plants commited
% J_star       = -1*ones(1, length(coal_pwr)); % Cummulative cost
% cmt_dispatch = -1*ones(1, length(coal_pwr)); % Optimal commitment
% f_dispatch   = -1*ones(1, length(coal_pwr));
% u_dispatch   = -1*ones(coal_num, length(coal_pwr));
% v_dispatch   = -1*ones(coal_num, length(coal_pwr));
% 
% opt_cost_startup  = -1*ones(1, length(coal_pwr));
% opt_cost_base_vom = -1*ones(1, length(coal_pwr));
% opt_cost_fuel     = -1*ones(1, length(coal_pwr));
% opt_cost_ramp     = -1*ones(1, length(coal_pwr));
% 
% t = N;
% id_dispatch(t) = find(v_range>=coal_pwr(t), 1, 'first');
% f_extract = f_table(id_dispatch(t),:);
% cost_fuel_tmp = f_extract * coal_price;
% cost_base_vom_tmp = coal_pwr(t) * coal_baseload;
% J_tmp = cost_fuel_tmp + cost_base_vom_tmp;
% [value, id_opt] = min(J_tmp);
% J_star(t) = value;
% cmt_dispatch(t) = id_opt;
% f_dispatch(t) = f_extract(id_opt);
% u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_dispatch(t));
% v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_dispatch(t));
% opt_cost_startup(t) = 0;
% opt_cost_base_vom(t) = cost_base_vom_tmp;
% opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
% opt_cost_ramp(t) = 0;
% 
% tic;
% for t = N-1:-1:1
%     id_dispatch(t) = find(v_range>=coal_pwr(t), 1, 'first');
%     f_extract = f_table(id_dispatch(t),:);
%     cost_fuel_tmp = f_extract * coal_price;
%     cost_base_vom_tmp = coal_pwr(t) * coal_baseload;
%     
%     cost_startup_tmp = zeros(1, coal_num);
%     cost_ramp_tmp = zeros(1, coal_num);
%     for cmt = 1:coal_num
%         u_tmp = u_table(:,id_dispatch(t),cmt);
%         v_tmp = v_table(:,id_dispatch(t),cmt);
%         
%         c = cmt_dispatch(t+1) - cmt;
%         if c>0
%             cost_startup_tmp(cmt) = c*coal_startup_cost;
%         else
%             cost_startup_tmp(cmt) = 0;
%         end
%         
%         d_coal_pwr = u_tmp - u_dispatch(:,t+1);
%         d_coal_pctg = d_coal_pwr/coal_nameplate;
%         x = abs(d_coal_pctg);
%         y = (x-0.3)*6.5/0.7+1.5;
%         ramp_scale = ones(size(x));
%         ramp_scale(x>0.3) = y(x>0.3);
%         cost_ramp_tmp(cmt) = sum(abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow);
%     end
%     J_tmp = J_star(t+1) + cost_fuel_tmp + cost_base_vom_tmp + cost_startup_tmp + cost_ramp_tmp;
%     
%     [value, id_opt] = min(J_tmp);
%     J_star(t) = value;
%     cmt_dispatch(t) = id_opt;
%     f_dispatch(t) = f_table(id_dispatch(t), id_opt);
%     u_dispatch(:,t) = u_table(:,id_dispatch(t),cmt_dispatch(t));
%     v_dispatch(:,t) = v_table(:,id_dispatch(t),cmt_dispatch(t));
%     opt_cost_startup(t) = cost_startup_tmp(id_opt);
%     opt_cost_base_vom(t) = cost_base_vom_tmp;
%     opt_cost_fuel(t) = cost_fuel_tmp(id_opt);
%     opt_cost_ramp(t) = cost_ramp_tmp(id_opt);
% end
% toc;

load(['DP_', wind_file, '_nominal']);


%%
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

% save(['DP_', wind_file, '_nominal'], ...
%      'id_dispatch', 'cmt_dispatch', 'f_dispatch', 'v_dispatch', 'u_dispatch', ...
%      'wind_pwr', 'wind_curtail', ...
%      'opt_cost_startup', 'opt_cost_base_vom', 'opt_cost_fuel', 'opt_cost_ramp', ...
%      'cost_startup', 'cost_base_vom', 'cost_fuel', 'cost_ramp', 'cost_total', ...
%      'J_star');


%% ========================================================================
% Check ramping
d_coal_pwr = [zeros(1,coal_num);
              diff(u_dispatch')]'; % [14]x[8760]
d_coal_pctg = d_coal_pwr/coal_nameplate;

% Check changes in commitment
d_cmt = [0, diff(cmt_dispatch)];
cmt_c = zeros(1, length(coal_pwr)); % Commition
cmt_c(d_cmt>0) = d_cmt(d_cmt>0);
cmt_d = zeros(1, length(coal_pwr)); % Decommition
cmt_d(d_cmt<0) = d_cmt(d_cmt<0);

% Commit or decommit 2 or more coal plants
id_jump = find(abs(d_cmt)>1) - 1;


%% ========================================================================
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
v_unique = zeros(1, length(coal_pwr));
u_unique = zeros(1, length(coal_pwr));
for t = 1:length(coal_pwr)
    vt = v_dispatch(1:cmt_dispatch(t),t);
    if length(unique(vt))>1
        disp([num2str(t), ': output power not equally distributed']); % this will triger errors
        v_unique(t) = vt(1);
    else
        v_unique(t) = unique(vt);
    end

    ut = u_dispatch(1:cmt_dispatch(t),t);
    if length(unique(ut))>1
        disp([num2str(t), ': output power not equally distributed']); % this will triger errors
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
ha = area([coal_pwr;wind_pwr]', 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.7 0.7]);
set(ha(2), 'facec', [0.6 1 0]);

area(v_dispatch', 'facec', 'none', 'edgecolor', [1 1 1]);
for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], coal_pwr(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end

ylim([0 9000]);
xlim([0 24*7*4]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 0:1500:9000);
set(gca, 'layer', 'top');
ylabel('Total Output Power (MW)');
title('Time Window: 4 Weeks');

% ====================
ax2 = subplot(4,1,3); hold on; box on;
plot(cmt_dispatch, 'color', [1 1 1]*0.6, 'linewidth', 1);
if ~isempty(id_jump)
for i = 1:length(id_jump)
h = plot(id_jump(i)+[0 1], cmt_dispatch(id_jump(i)+[0 1]), 'ko-', 'markersize', 2, 'markerf', 'k');
end
legend(h, 'Steep Change');
set(legend, 'location', 'southwest');
end

xlim([0 168*7]);
ylim([9.5 14.5]);
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

end


%% Check consistency: not identical, but extremely close
% figure(4); clf; hold on; box on;
% z1 = downsample(J_star, 500);
% plot(z1, 'x-');
% 
% opt_cost_total = opt_cost_startup + opt_cost_base_vom + opt_cost_fuel + opt_cost_ramp;
% J_alternative = fliplr(cumsum(fliplr(opt_cost_total)));
% z2 = downsample(J_alternative, 500);
% plot(z2, 'o-');
% title('Check Consistency: J\_star & Total Cumulative Cost');
% legend('J\_star', 'Total Cumulative Cost');
% 
% dJ = J_star - J_alternative;
% dJ_pctg = dJ./J_star;
% 
% disp('========================================');
% disp('Double check costs (identical)');
% opt_cost_startup2 = cmt_c * coal_startup_cost;
% opt_cost_base_vom2 = coal_pwr * coal_baseload;
% opt_cost_fuel2 = f_dispatch * coal_price;
% 
% x = abs(d_coal_pctg);
% y = (x-0.3)*6.5/0.7+1.5;
% ramp_scale = ones(size(x)); % Exceeding 30 precent nameplate ramping is scaled by 1.5-8 times
% ramp_scale(x>0.3) = y(x>0.3);
% opt_cost_ramp2 = abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow; % [$]
% 
% cost_startup = sum(opt_cost_startup2)
% cost_base_vom = sum(opt_cost_base_vom2)
% cost_fuel = sum(opt_cost_fuel2)
% cost_ramp = sum(opt_cost_ramp2)
% 
% disp('====================');
% cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp
% 
% disp('====================');
% pctg_startup = cost_startup/cost_total
% pctg_base_vom = cost_base_vom/cost_total
% pctg_fuel = cost_fuel/cost_total
% pctg_ramp = cost_ramp/cost_total

