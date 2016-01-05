clear
clc
close all

%% Unity step size
% sample_cnt = 10000;
% esp = rand(24,sample_cnt); % 0~1
% resp = round(esp);  % 0 or 1
% s = -1 + 2*resp;    %-1 or 1
% w = cumsum(s);
% mu = mean(w);
% rho = std(w);
% w0 = w(1,:); % First position
% wn = w(end,:); % Last position
% 
% % ==============================
% figure(1); clf; hold on; box on;
% plot(w);
% ylabel('Position (-)');
% xlabel('Horizon (Step Counts)');
% title([num2str(sample_cnt), ' Radom Walk Trails']);
% 
% % ====================
% figure(11); clf;
% subplot(3,1,1);
% plot(mu);
% ylabel('Mean (-)');
% xlabel('Samples');
% 
% % ====================
% subplot(3,1,2);
% plot(rho);
% ylabel('STD (-)');
% xlabel('Samples');
% 
% % ====================
% subplot(3,1,3);
% bin = -18:2:18;
% cnt = hist(wn, bin);
% bar(bin,cnt,1,'facec', [0.35 0.65 1]);
% xlabel('Distribution of Final Position');


%% Continous step size
sample_cnt = 500;
esp = rand(24,sample_cnt); % 0~1
s = esp-0.5; % -0.5~0.5
w = cumsum(s); % Walk
mu = mean(w);
rho = std(w);
w0 = w(1,:); % First position
wn = w(end,:); % Last position

% ==============================
figure(2); clf; hold on; box on;
plot(w);
ylabel('Position (-)');
xlabel('Horizon (Step Counts)');
title([num2str(sample_cnt), ' Radom Walk Trails']);
ylim([-5 5]);
xlim([0 25]);

% ====================
figure(21); clf;
subplot(3,1,1);
plot(mu);
ylabel('Mean (-)');
xlabel('Samples');

% ====================
subplot(3,1,2);
plot(rho);
ylabel('STD (-)');
xlabel('Samples');

% ====================
subplot(3,1,3);
bin = -5:1:5;
cnt = hist(wn, bin);
bar(bin,cnt,1,'facec', [0.35 0.65 1]);
xlabel('Distribution of Final Position');
ylabel('Counts');

