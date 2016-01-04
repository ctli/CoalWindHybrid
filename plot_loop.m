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

% Cost improvement
dc = (dp.table_cost_total)./(myopic.table_cost_total);
figure(31); clf;
surf(wind_pwr_range, target_pwr_range, dc');
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Cost Difference (-)');
title ('Myopic vs. Horizon-Based (DP) Dispatch');
ylim([2000 10000]);
view(40, 25);

figure(32); clf;
[C,h] = contourf(wind_pwr_range, target_pwr_range, dc');
clabel(C, h);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Cost Difference (-)');
title ('Myopic vs. Horizon-Based (DP) Dispatch');
grid on;

