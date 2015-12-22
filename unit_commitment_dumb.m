clear
close all
clc
format compact

load OptTable;

load Xilingol_2009;
wind_pwr = p*2500;

target_pwr = 8500;
coal_pwr = target_pwr - wind_pwr; % Use coal to make up deficit

%% Myopic unit commitment
for t = 1:10%length(coal_pwr)
    id = find(v_st>coal_pwr(t), 1, 'first');
end