clear
close all
clc

wind_file = 'Xilingol_2009';

wind_pwr_range = 0:500:5000; % [1x11]
target_pwr_range = 2500:500:8500; % [1x13]


%% Myopic disaptch; myopic dispatch does not do "economic" wind curtailment
% m1 = load(['Myopic_', wind_file, '_nominal']);
% m2 = load(['Myopic_', wind_file, '_nominal_new']);
% [a1,b1] = unique(m1.wind_curtail);
% [a2,b2] = unique(m2.wind_curtail);
% isequal(b1,b2) % Identical

% % ==============================
% myopic = load(['myopic_', wind_file, '_loop']);
% myopic_new = load(['myopic_', wind_file, '_loop_new']);
% [a1,b1] = unique(myopic.table_wind_curtail);
% [a2,b2] = unique(myopic_new.table_wind_curtail);
% isequal(b1,b2) % Identical
% 
% % Total generation cost
% figure(1); clf; hold on; grid on;
% mesh(wind_pwr_range, target_pwr_range, myopic_new.table_cost_total'/1e9);
% % mesh(wind_pwr_range, target_pwr_range, myopic.table_cost_total'/1e9);
% xlabel('Wind Capacity (MW)');
% ylabel('Target Output (MW)');
% zlabel('Total Generation Cost (Billion USD)');
% ylim([2000 10000]);
% zlim([0 1.8]);
% view(35, 10);
% alpha(0.5);


%% Horizon dispatch; this has economic wind curtailment but very tiny
% dp = load(['DP_', wind_file, '_loop']);
% dp_new = load(['DP_', wind_file, '_loop_new']);
% 
% % ==============================
% % With & without economic wind curtailment
% figure(2); clf;
% mesh(wind_pwr_range, target_pwr_range, dp.table_cost_total'/1e9); hold on;
% mesh(wind_pwr_range, target_pwr_range, dp_new.table_cost_total'/1e9);
% xlabel('Wind Capacity (MW)');
% ylabel('Target Output (MW)');
% zlabel('Total Generation Cost (Billion USD)');
% ylim([2000 10000]);
% zlim([0 1.8]);
% view(35, 10);
% alpha(0.5);
% text(wind_pwr_range(end), target_pwr_range(1), dp.table_cost_total(end,1)'/1e9, '\leftarroww/o Economic Wind Curtailment', 'color', 'r');
% text(wind_pwr_range(end), target_pwr_range(1), dp_new.table_cost_total(end,1)'/1e9, '\leftarroww/ Economic Wind Curtailment');
% 
% % Cost improvement
% dc = (dp_new.table_cost_total)./(dp.table_cost_total);
% d = 1-dc;
% figure(21); clf;
% surf(wind_pwr_range, target_pwr_range, dc');
% xlabel('Wind Capacity (MW)');
% ylabel('Target Output (MW)');
% zlabel('Cost Difference (-)');
% title ('With & Without Economic Wind Curtailment');
% ylim([2000 10000]);
% view(25, 15);


%% Myopic vs. DP
myopic = load(['myopic_', wind_file, '_loop']);
dp = load(['DP_', wind_file, '_loop_new']);

% ====================
figure(3); clf; hold on; grid on;
mesh(wind_pwr_range, target_pwr_range, myopic.table_cost_total'/1e9);
mesh(wind_pwr_range, target_pwr_range, dp.table_cost_total'/1e9);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Total Generation Cost (Billion USD)');
ylim([2000 10000]);
zlim([0 1.8]);
view(35, 10);
alpha(0.5);

text(wind_pwr_range(end), target_pwr_range(end), myopic.table_cost_total(end,end)'/1e9, '\leftarrowMyopic');
text(wind_pwr_range(end), target_pwr_range(end), dp.table_cost_total(end,end)'/1e9, '\leftarrowDP');

% ====================
% Cost improvement
dc = (dp.table_cost_total)./(myopic.table_cost_total);
% figure(31); clf;
% surf(wind_pwr_range, target_pwr_range, dc');
% xlabel('Wind Capacity (MW)');
% ylabel('Target Output (MW)');
% zlabel('Cost Difference (-)');
% title ('Myopic vs. Horizon-Based (DP) Dispatch');
% ylim([2000 10000]);
% view(40, 25);

figure(32); clf;
[C,h] = contourf(wind_pwr_range, target_pwr_range, dc');
clabel(C, h);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Cost Difference (-)');
title ('Cost Difference: Myopic vs. Horizon-Based (DP) Dispatch');
grid on;
set(gca, 'ytick', 2500:1000:8500);
set(gca, 'xtick', 0:1000:5000);


%% ====================
% Cost makeups
pctg_fuel1 = dp.table_cost_fuel./dp.table_cost_total*100;
pctg_base_vom1 = dp.table_cost_base_vom./dp.table_cost_total*100;
pctg_startup1 = dp.table_cost_startup./dp.table_cost_total*100;
pctg_ramp1 = dp.table_cost_ramp./dp.table_cost_total*100;

pctg_fuel2 = myopic.table_cost_fuel./myopic.table_cost_total*100;
pctg_base_vom2 = myopic.table_cost_base_vom./myopic.table_cost_total*100;
pctg_startup2 = myopic.table_cost_startup./myopic.table_cost_total*100;
pctg_ramp2 = myopic.table_cost_ramp./myopic.table_cost_total*100;

id_r_c = [1,1
          1,13
          11,13
          11,1];
id = sub2ind(size(pctg_fuel1), id_r_c(:,1), id_r_c(:,2));

dp_bar = [pctg_fuel1(id), pctg_base_vom1(id), pctg_startup1(id), pctg_ramp1(id)];
myopic_bar = [pctg_fuel2(id), pctg_base_vom2(id), pctg_startup2(id), pctg_ramp2(id)];

% figure(4); clf;
% i = 1;
% data = [myopic_bar(i,:); dp_bar(i,:)];
% bar(1:2, data, 'stacked');
% colormap summer
% title('2500MW wind; Target Output@2500MW');
% goldenratio;
% set(gca, 'xtick', 1:2, 'xticklabel', {'Myopic', 'DP'});
% 
% figure(41); clf;
% i = 2;
% data = [myopic_bar(i,:); dp_bar(i,:)];
% bar(1:2, data, 'stacked');
% colormap summer
% title('2500MW wind; Target Output@8500MW');
% goldenratio;
% set(gca, 'xtick', 1:2, 'xticklabel', {'Myopic', 'DP'});
% 
% figure(42); clf;
% i = 3;
% data = [myopic_bar(i,:); dp_bar(i,:)];
% bar(1:2, data, 'stacked');
% colormap summer
% title('5000MW wind; Target Output@2500MW');
% goldenratio;
% set(gca, 'xtick', 1:2, 'xticklabel', {'Myopic', 'DP'});

figure(43); clf;
i = 4;
data = [myopic_bar(i,:); dp_bar(i,:)];
hb = bar(1:2, data, 0.7,'stacked', 'edgecolor', 'w');
set(hb(4), 'facec', [0.95 0.95 0]);
colormap summer
title('5000MW wind; Target Output@8500MW');
set(gca, 'xtick', 1:2, 'xticklabel', {'Myopic', 'DP'});
ylabel('Cost Shares (%)');
ylim([0 102]);
box off;
set(gcf, 'units', 'inch', 'pos', [5.8438    2.3542    3.85    2.35]);

text(1, myopic_bar(i,1)/2, ['Fuel: ', num2str(myopic_bar(i,1), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(1, myopic_bar(i,1)+myopic_bar(i,2)/2, ['O&M: ', num2str(myopic_bar(i,2), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(1, myopic_bar(i,1)+myopic_bar(i,2)+myopic_bar(i,3)/2, ['Startup: ', num2str(myopic_bar(i,3), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(1, myopic_bar(i,1)+myopic_bar(i,2)+myopic_bar(i,3)+myopic_bar(i,4)/2, ['Ramping: ', num2str(myopic_bar(i,4), '%2.0f'), '%'], 'horizontalalignment', 'center');

text(2, dp_bar(i,1)/2, [num2str(dp_bar(i,1), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(2, dp_bar(i,1)+dp_bar(i,2)/2, [ num2str(dp_bar(i,2), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(2, dp_bar(i,1)+dp_bar(i,2)+dp_bar(i,3)/2, [num2str(dp_bar(i,3), '%2.0f'), '%'], 'horizontalalignment', 'center');
text(2, dp_bar(i,1)+dp_bar(i,2)+dp_bar(i,3)+dp_bar(i,4)/2, [num2str(dp_bar(i,4), '%2.0f'), '%'], 'horizontalalignment', 'center');

%% Wind Curtailment

