clear
close all
clc
format compact


%% Random walk with continuous step sizes
sample_cnt = 50;
esp = rand(24,sample_cnt); % 0~1
s = esp-0.5; % -0.5~0.5
w = cumsum(s); % Walk

% ====================
figure(1); clf; hold on; box on;
plot(w);
ylabel('Position (-)');
xlabel('Horizon (Step Counts)');
title([num2str(sample_cnt), ' Radom Walk Trails']);
ylim([-5 5]);
xlim([0 25]);
grid on;

%% Scaled random walk
r0 = 0.16/0.5;
rn = 0.5/5;
r = repmat(linspace(r0,rn,24)',1,sample_cnt);
ww = w.*0.1;

figure(2); clf;
plot(ww);
ylabel('Position (-)');
xlabel('Horizon (Step Counts)');
title([num2str(sample_cnt), ' Scaled Radom Walk Trails']);
ylim([-0.5 0.5]);
xlim([0 25]);
grid on;


%%
wind_file = 'Xilingol_2009';
load(wind_file);

% ====================
i = 240;
f24 = p((1:24)+i); % Forecast
p24 = repmat(f24,1,sample_cnt) + ww; % Fake realizeations
p24(p24<0) = 0;
p24(p24>1) = 1;

figure(3); clf; hold on; box on;
plot(p24);
plot(f24, 'k', 'linewidth', 2);
ylim([0 1]);
ylim([-0.02 1.02]);
grid on;
xlabel('Time (hr)');
ylabel('Wind Power (normalized)');
title(['Forecast and 50 Random Realization in Day ', num2str(i)]);

% ====================
i = 625;
f24 = p((1:24)+i); % Forecast
p24 = repmat(f24,1,sample_cnt) + ww; % Fake realizeations
p24(p24<0) = 0;
p24(p24>1) = 1;

figure(31); clf; hold on; box on;
plot(p24);
plot(f24, 'k', 'linewidth', 2);
ylim([0 1]);
ylim([-0.02 1.02]);
grid on;
xlabel('Time (hr)');
ylabel('Wind Power (normalized)');
title(['Forecast and 50 Random Realization in Day ', num2str(i)]);

% ====================
i = 1033;
f24 = p((1:24)+i); % Forecast
p24 = repmat(f24,1,sample_cnt) + ww; % Fake realizeations
p24(p24<0) = 0;
p24(p24>1) = 1;

figure(32); clf; hold on; box on;
plot(p24);
plot(f24, 'k', 'linewidth', 2);
ylim([0 1]);
ylim([-0.02 1.02]);
grid on;
xlabel('Time (hr)');
ylabel('Wind Power (normalized)');
title(['Forecast and 50 Random Realization in Day ', num2str(i)]);

