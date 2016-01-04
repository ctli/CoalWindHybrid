clear
close all
clc
format compact

%% =====================================================================
% data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
% obs = data.data(:,3);  %obs = reshape(obs, 24, length(obs)/24);
% fcst = data.data(:,4); %fcst = reshape(fcst, 24, length(fcst)/24);
% 
% scale = 800/max(fcst);
% fcst = fcst*scale;
% obs = obs*scale;
% 
% % ====================
% figure(2); clf;
% hist(obs);
% title('Wind Power Distribution');


%% =====================================================================
% wind_file = 'Xilingol_2009';
% load(wind_file);
% wind_pwr = round(p*2500)';

% % ====================
% figure(1); clf;
% set(gcf, 'units', 'inch', 'pos', [2.9792    1.4583    5.8333    6.75]);
% 
% subplot(3,1,1);
% bin = 0:2.5:25;
% n = hist(v, bin);
% bar(bin, n/sum(n), 1, 'facec', [0.3 0.65 1]);
% xlabel('Wind Speed (m/s)')
% ylabel('Probability (-)');
% xlim([-1.5 26.5]);
% 
% subplot(3,1,2);
% bin = 0:0.1:1;
% n = hist(p, bin);
% bar(bin, n/sum(n), 1, 'facec', [0.3 0.65 1]);
% xlabel('Wind Power (-)')
% ylabel('Probability (-)');
% xlim([-0.06 1.06]);
% 
% subplot(3,1,3);
% lag_number = 48;
% autocorr(p,lag_number);
% set(gca, 'xtick', 0:4:48);
% xlim([-1 49]);
% ylim([-0.1 1.1]);
% set(gca, 'ytick', -0:0.25:1);


%% =====================================================================
% wind_file = 'Xilingol_2009';
% load(wind_file);
% 
% n = 5;
% x1 = p(1:end-4); % [Nx1]; column vector
% x2 = p(2:end-3);
% x3 = p(3:end-2);
% x4 = p(4:end-1);
% x5 = p(5:end);
% x = [x1, x2, x3, x4, x5];
% mu = mean(x);
% 
% N = length(x1);
% sig = zeros(n,n);
% for r = 1:n % Row
%     for c = 1:n % Column
%         xi = x(:,r);
%         xj = x(:,c);
%         mi = mean(xi);
%         mj = mean(xj);
%         cov_ij = 1/N*sum((xi-mi).*(xj-mj));
%         sig(r,c) = cov_ij;
%     end
% end
% 
% figure(2); clf; hold on; box on;
% plot(x1, 'linewidth', 2);
% plot(x2);
% plot(x3);
% plot(x4);
% plot(x5);
% xlim([0 24]);
% ylim([-0.02 1.02]);
% set(gca, 'xtick', 0:4:8760);
% grid on;
% 
% % ====================
% rng(1); % Fix the seed of random generator
% xr = mvnrnd(mu,sig,200);
% 
% figure(3); clf; hold on; box on;
% plot(xr');
% plot(xr(1,:), 'k', 'linewidth', 2);
% xlim([0 24]);


%% =====================================================================
wind_file = 'Xilingol_2009';
load(wind_file);

n = 24;
x = zeros(length(p)-n+1, n);
for nn = 1:n
    st = nn;
    ed = n-st;
    xn = p(st:end-ed);
    x(:,nn) = xn;
end
mu = mean(x);

% ====================
N = length(x);
sig = zeros(n,n);
for r = 1:n % Row
    for c = 1:n % Column
        xi = x(:,r);
        xj = x(:,c);
        mi = mean(xi);
        mj = mean(xj);
        cov_ij = 1/N*sum((xi-mi).*(xj-mj));
        sig(r,c) = cov_ij;
    end
end

% ====================
rng(1); % Fix the seed of random generator
xr = mvnrnd(mu,sig,50);

figure(3); clf; hold on; box on;
plot(xr');
plot(xr(1,:), 'k', 'linewidth', 2);
xlim([0 25]);
title('50 Realizations of 24-Hr Long Random Variables');

