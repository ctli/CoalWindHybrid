clear
close all
clc

wind_file = 'Xilingol_2009';

wind_pwr_range = 0:500:5000; % [1x11]
target_pwr_range = 2500:500:8500; % [1x13]


%% Myopic disaptch; myopic dispatch does not do "economic" wind curtailment
% myopic = load(['myopic_', wind_file, '_loop']);
% myopic_new = load(['myopic_', wind_file, '_loop_new']);
% [a1,b1] = unique(myopic.table_wind_curtail);
% [a2,b2] = unique(myopic_new.table_wind_curtail);
% isequal(b1,b2) % Identical
% 
% % ==============================
% % Total generation cost
% figure(1); clf; hold on; grid on;
% mesh(wind_pwr_range, target_pwr_range, myopic_new.table_cost_total'/1e9);
% xlabel('Wind Capacity (MW)');
% ylabel('Target Output (MW)');
% zlabel('Total Generation Cost (Billion USD)');
% ylim([2000 10000]);
% zlim([0 1.8]);
% view(35, 10);
% alpha(0.5);


%% Horizon dispatch; this has economic wind curtailment
dp1 = load(['DP_', wind_file, '_nominal']);
dp2 = load(['DP_', wind_file, '_nominal_new']);
isequal(dp1.J_star, dp2.J_star)  % not same

% dp = load(['DP_', wind_file, '_loop']);
% dp_new = load(['DP_', wind_file, '_loop_new']);

