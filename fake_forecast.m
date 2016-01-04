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
wind_file = 'Xilingol_2009';
load(wind_file);
wind_pwr = round(p*2500)';

% ====================
figure(1); clf;
set(gcf, 'units', 'inch', 'pos', [2.9792    1.4583    5.8333    6.75]);

subplot(3,1,1);
bin = 0:2.5:25;
n = hist(v, bin);
bar(bin, n/sum(n), 1, 'facec', [0.3 0.65 1]);
xlabel('Wind Speed (m/s)')
ylabel('Probability (-)');
xlim([-1.5 26.5]);

subplot(3,1,2);
bin = 0:0.1:1;
n = hist(p, bin);
bar(bin, n/sum(n), 1, 'facec', [0.3 0.65 1]);
xlabel('Wind Power (-)')
ylabel('Probability (-)');
xlim([-0.06 1.06]);

subplot(3,1,3);
lag_number = 48;
autocorr(p,lag_number);
set(gca, 'xtick', 0:4:48);
xlim([-1 49]);
ylim([-0.1 1.1]);
set(gca, 'ytick', -0:0.25:1);


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
% X = mvnrnd(zeros(1,n),sig,200);
% 
% figure(3); clf; hold on; box on;
% plot(X');
% plot(X(1,:), 'k', 'linewidth', 2);
% xlim([0 24]);


%% =====================================================================
% wind_file = 'Xilingol_2009';
% load(wind_file);

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
p = v_log(:);

n = 24;
x = zeros(length(p)-n+1, n);
for nn = 1:n
    st = nn;
    ed = n-st;
    xn = p(st:end-ed);
    x(:,nn) = xn;
end

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
sig_log = sig;

figure(100); clf;
surf(sig_log);


%% ====================
rng(1); % Fix the seed of random generator
sig = sig_log;
X = mvnrnd(zeros(1,n),sig,50);

figure(3); clf; hold on; box on;
plot(X');
plot(X(1,:), 'k', 'linewidth', 2);
xlim([0 25]);
title('X (Normal Distribution: -inf~+inf)');

%% ===========
% Y = sqrt(2)*erfinv(2*X-1); % Uniform distribution->normal distribution (http://www.mathworks.com/matlabcentral/answers/35281-transforming-uniform-variables-to-normal-variables)
% figure(4); clf; hold on; box on;
% plot(Y');
% plot(Y(1,:), 'k', 'linewidth', 2);
% xlim([0 25]);
% title('50 Realizations of 24-Hr Long Random Variables');

% Normal distribution->uniform distribution (http://math.stackexchange.com/questions/153793/how-to-transform-normally-distributed-random-sequence-n0-1-to-uniformly-distri)
Y = 0.5*(erf(X/sqrt(2))+1);
figure(5); clf; hold on; box on;
plot(Y');
plot(Y(1,:), 'k', 'linewidth', 2);
xlim([0 25]);
title('Y (Uniform Distribution: 0~1)');

%% 
load Xilingol_probability;
p_realization = interp1(Fv,xv,Y);

figure(6); clf; hold on; box on;
plot(p_realization');
plot(p_realization(1,:), 'k', 'linewidth', 2);
title('p (Realization)');
