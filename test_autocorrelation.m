clear
close all
clc

yr_range = 1979:2009;
file_prefix = 'Xilingol_';

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
