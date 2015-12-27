clear
close all
clc
format compact

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
coal_num = 14;


%% Myopic unit commitment
% load OptTable;
% 
% f_myopic = zeros(1,length(v_range));
% cmt_myopic = zeros(1,length(v_range));
% v_myopic = zeros(coal_num,length(v_range));
% u_myopic = zeros(coal_num,length(v_range));
% tic;
% for vv = 1:length(v_range)
%     f_column = f_table(vv,:);
%     [f_myopic(vv), cmt_myopic(vv)] = min(f_column);
%     v_myopic(:,vv) = v_table(:,vv,cmt_myopic(vv));
%     u_myopic(:,vv) = u_table(:,vv,cmt_myopic(vv));
% end
% toc;
% save('MyopicDispatch', 'v_range', 'f_myopic', 'cmt_myopic', 'v_myopic', 'u_myopic');
% 
% % ====================
% figure(1); clf; hold on; box on; % Fuel consumption
% plot(0,0,'x');
% plot(v_range, f_table, 'linewidth', 1);
% for n = 1:coal_num
%     if n==1
%         text(v_ed(n), max(f_table(:,n)), [' ', num2str(n), ' Unit is commited'], 'fontsize', 7);
%     elseif n<11
%         text(v_ed(n), max(f_table(:,n)), [' ', num2str(n), ' Units are commited'], 'fontsize', 7);
%     else
%         text(v_st(n), min(f_table(:,n)), [' ', num2str(n), ' Units are commited '], 'fontsize', 7, 'horizontalalignment', 'right');
%     end
% end
% xlabel('Output Power, MW (in-house use excluded)');
% ylabel('Coal Consumption (ton/h)');
% my_gridline;
% h1 = plot(v_range, f_myopic, 'color', [1 1 1]*0, 'linewidth', 0.35);
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

load  MyopicDispatch;


%% Wind power
load Xilingol_2009;
wind_pwr = p*2500;
d_wind_pwr = diff(p);

target_pwr = 8500;
coal_pwr = target_pwr - wind_pwr; % Use coal to make up deficit


%% Myopic dispatch
tic;
id_dispatch = zeros(1, length(coal_pwr));
for t = 1:length(coal_pwr)
    id_dispatch(t) = find(v_range>=coal_pwr(t), 1, 'first');
end
f_dispatch = f_myopic(id_dispatch); % [ton/h]
cmt_dispatch = cmt_myopic(id_dispatch);
v_dispatch = v_myopic(:,id_dispatch);
u_dispatch = u_myopic(:,id_dispatch);
toc;

% Check ramping
d_coal_pwr = diff(u_dispatch'); % [8759]x[14]
d_coal_pctg = d_coal_pwr/coal_nameplate;

% Check changes in commitment
d_cmt = diff(cmt_dispatch);
cmt_c = d_cmt(d_cmt>0); % Commition
cmt_d = d_cmt(d_cmt<0); % Decommition

% Commit or decommit 2 or more coal plants
id_jump = find(abs(d_cmt)>1);

plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
% ====================
figure(3); clf; hold on;
ax1 = subplot(3,1,1:2); hold on; box on;
ha = area([coal_pwr, wind_pwr], 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.6 0.6]);
set(ha(2), 'facec', [0.6 1 0]);

area(v_dispatch', 'facec', 'none', 'edgecolor', [0.5 0 0]);

ylim([0 9000]);
xlim([0 24*7*4]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 0:1500:9000);
set(gca, 'layer', 'top');
ylabel('Total Output Power (MW)');
title('Time Window: 4 Weeks');
box on;

% ==========
ax2 = subplot(3,1,3); hold on; box on;
v_unique = zeros(1, length(coal_pwr));
for t = 1:length(coal_pwr)
    vt = v_dispatch(1:cmt_dispatch(t),t);
    if length(unique(vt))>1
        disp('output power not equally distributed'); % this will triger errors
    end
    v_unique(t) = unique(vt);
end
plot(v_unique);
ylabel({'Output Power of', 'Each Coal Plant (MW)'});
xlabel('Time (hr)');
ylim([0 660]);
set(gca, 'xtick', 0:168:24*7*4);
linkaxes([ax1, ax2], 'x');

% ====================
figure(4); clf;
ax1 = subplot(3,1,1:2); hold on; box on;
ha = area([coal_pwr, wind_pwr], 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.6 0.6]);
set(ha(2), 'facec', [0.6 1 0]);
for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], coal_pwr(id_jump(i)+[0 1]), 'kx-', 'markersize',  3);
end

ylim([0 9000]);
xlim([0 24*7*4]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 0:1500:9000);
legend(fliplr(ha), 'Wind', 'Coal');
set(legend, 'location', 'southwest');
ylabel('Output Power (MW)');
title('Time Window: 4 Weeks');
my_gridline('front');

% ==========
ax2 = subplot(3,1,3); hold on; box on;
plot(cmt_dispatch);
for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], cmt_dispatch(id_jump(i)+[0 1]), 'kx-', 'markersize', 3);
end

xlim([0 168*7]);
ylim([min(cmt_dispatch)-0.5, max(cmt_dispatch)+0.5]);
set(gca, 'xtick', 0:168:24*7*4);
set(gca, 'ytick', 1:1:14);
ylabel({'Unit Commited', '(count)'});
xlabel('Time (hr)');
linkaxes([ax1, ax2], 'x');
my_gridline;
end


%% Costs of coal power plants
% China:
coal_price = 57; % i.e. 380 [RMB/tonne] = 57 [$/tonne] (1RMB = 0.15USD)
coal_heatcontent = 21.8; % 5500 [kcal/kg] = 21.8 [mmBtu/tonne] (1kcal = 3.96567 Btu)

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


%%
cost_startup = sum(cmt_c * coal_startup_cost) % [$]
cost_base_vom = sum(coal_pwr * coal_baseload)
cost_fuel = sum(f_dispatch * coal_price) % [$]

% cost_ramp = sum(abs(d_coal_pwr(:)) * coal_loadfollow)
x = abs(d_coal_pctg);
y = (x-0.3)*6.5/0.7+1.5;
ramp_scale = zeros(size(x)); % Exceeding 30 precent nameplate ramping is scaled by 1.5-8 times
ramp_scale(x>0.3) = y(x>0.3);
cost_ramp = sum(abs(d_coal_pwr(:)).*ramp_scale(:) * coal_loadfollow)

disp('====================');
cost_total = cost_startup + cost_base_vom + cost_fuel + cost_ramp

disp('====================');
pctg_startup = cost_startup/cost_total
pctg_base_vom = cost_base_vom/cost_total
pctg_fuel = cost_fuel/cost_total
pctg_ramp = cost_ramp/cost_total



