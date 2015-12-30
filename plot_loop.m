clear
close all
clc

wind_file = 'Xilingol_2009';

wind_pwr_range = 0:500:5000; % [1x11]
target_pwr_range = 2500:500:8500; % [1x13]

myopic = load(['myopic_', wind_file, '_loop']);
myopic_new = load(['myopic_', wind_file, '_loop_new']);

% dp = load(['DP_', wind_file, '_loop']);
% dp_new = load(['DP_', wind_file, '_loop_new']);


%% Total generation cost
figure(1); clf; hold on; grid on;
mesh(wind_pwr_range, target_pwr_range, myopic.table_cost_total'/1e9);
mesh(wind_pwr_range, target_pwr_range, myopic_new.table_cost_total'/1e9);
xlabel('Wind Capacity (MW)');
ylabel('Target Output (MW)');
zlabel('Total Generation Cost (Billion USD)');
ylim([2000 10000]);
zlim([0 1.8]);
view(35, 10);
alpha(0.5);

text(wind_pwr_range(end), target_pwr_range(end), myopic.table_cost_total(end, end)/1e9, '\leftarrowMyopic', 'fontsize', 8);
text(wind_pwr_range(end), target_pwr_range(end), myopic_new.table_cost_total(end, end)/1e9, '\leftarrowMyopic (wind curtail allowed)', 'fontsize', 8);


isequal(myopic.table_cost_total, myopic_new.table_cost_total)
z = myopic.table_cost_total - myopic_new.table_cost_total;