clear
close all
clc


yr_range = 1979:2009;
file_prefix = 'Xilingol_';

p_log = nan*ones(8760, length(yr_range));
v_log = nan*ones(8760, length(yr_range));
for iy = 1:length(yr_range)
    yr = yr_range(iy);
    file_name = [file_prefix, num2str(yr)];
    load(file_name);
    if length(p)>=8760;
        p_log(:,iy) = p(1:8760);
    end
    if length(v)>=8760;
        v_log(:,iy) = v(1:8760);
    end
end
p = p_log(:);
mu = mean(p)
rho = 0.16%std(p)
% Haessig (2013): autocorrelation paper: sigma_P = 0.16 of capacity

phi1 = 0.8871;
theta1 = -0.1541;

rng(1); % Fix the seed of random generator
rho_esp = rho*sqrt(1-phi1^2)
esp = normrnd(0,rho_esp, [length(p), 1]);
fcst = phi1*p(1:end-1) + theta1*esp(1:end-1) + esp(2:end); % Fail->std in the white noise is not right

figure(1); clf;
plot(p, 'x-'); hold on;
plot(fcst);
xlim([0 168/4]);
