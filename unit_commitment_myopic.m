clear
close all
clc
format compact

coal_nameplate = 660; % [MW]
coal_useable = coal_nameplate*0.92; % 607.2 [MW]
coal_num = 14;


%% Costs of coal power plants
% China:
coal_price = 57; % i.e. 380 [RMB/tonne] = 57 [$/tonne] (1RMB = 0.15USD)
coal_heatcontent = 21.8; % 5500 [kcal/kg] = 21.8 [mmBtu/tonne] (1kcal = 3.96567 Btu)

% US:
coal_startup = 38; % hot startup cost [$/MW]
coal_startup_oth = 5.81; % other hot startup cost [$/MW]
coal_startup = coal_startup + coal_startup_oth;

coal_loadfollow = 1.72; % load following cost [$/MW]
coal_baseload = 3.22; % Baseload variable cost [$/MWh]
coal_startup_fuel = 10.1; % hot start fuel [mmBtu/MW]


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
% figure(1); clf; hold on; box on;
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
% % export_fig myopic -r300
% 
% figure(2); clf;
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
f_dispatch = f_myopic(id_dispatch);
cmt_dispatch = cmt_myopic(id_dispatch);
v_dispatch = v_myopic(:,id_dispatch);
u_dispatch = u_myopic(:,id_dispatch);
toc;

% Check ramping
d_coal_pwr = diff(coal_pwr);

% Check changes in commitment
d_cmt = diff(cmt_dispatch);

% Commit or decommit 2 or more coal plants
id_jump = find(abs(d_cmt)>1);


%% ========================================================================
figure(4); clf; hold on;
ax1 = subplot(3,1,1:2); hold on; box on;
ha = area([coal_pwr, wind_pwr], 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.6 0.6]);
set(ha(2), 'facec', [0.6 1 0]);

area(v_dispatch', 'facec', 'none', 'edgecolor', [0.5 0 0]);

ylim([0 9000]);
xlim([0 168*7]);
set(gca, 'xtick', 0:168:168*7);
set(gca, 'ytick', 0:1500:9000);
set(gca, 'layer', 'top');
ylabel('Total Output Power (MW)');
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
set(gca, 'xtick', 0:168:168*7);
linkaxes([ax1, ax2], 'x');

%% ====================
figure(3); clf;
ax1 = subplot(3,1,1:2); hold on; box on;
ha = area([coal_pwr, wind_pwr], 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.6 0.6]);
set(ha(2), 'facec', [0.6 1 0]);

for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], coal_pwr(id_jump(i)+[0 1]), 'kx-', 'markersize',  3);
end

ylim([0 9000]);
xlim([0 168*7]);
set(gca, 'xtick', 0:168:168*7);
set(gca, 'ytick', 0:1500:9000);
legend(fliplr(ha), 'Wind', 'Coal');
set(legend, 'location', 'southwest');
ylabel('Output Power (MW)');
my_gridline('front');

% ==========
ax2 = subplot(3,1,3); hold on; box on;
plot(cmt_dispatch);

for i = 1:length(id_jump)
plot(id_jump(i)+[0 1], cmt_dispatch(id_jump(i)+[0 1]), 'kx-', 'markersize', 3);
end

xlim([0 168*7]);
ylim([min(cmt_dispatch)-0.5, max(cmt_dispatch)+0.5]);
set(gca, 'xtick', 0:168:168*7);
set(gca, 'ytick', 1:1:14);
ylabel({'Unit Commited', '(count)'});
xlabel('Time (hr)');
linkaxes([ax1, ax2], 'x');
my_gridline;


%%


