% Recursive estimation of mean & covariance
clear
close all
clc
format compact

rng(1); % Fix the seed
x = randn(1,100);
m1 = mean(x(1))
sig1 = cov(x(1))

log1 = zeros(1,100);
log2 = zeros(1,100);

sig_old = sig1;
m_old = m1;
for t = 2:100
disp('==============================');
x2 = x(1:t);

% ====================
m_now = mean(x2);
m_now2 = m_old + 1/t*(x2(end)-m_old); % recursive mean

% ====================
sig_now = cov(x2)
log1(t) = sig_now;

% sig_now2 = (t-2)/(t-1)*sig_old + 1/(t-1)*(x2(end)-m_old)*(x2(end)-m_old)' % Not as good as the following estimation
sig_now2 = (t-1)/(t)*sig_old + 1/(t)*(x2(end)-m_old)*(x2(end)-m_old)' % Only approximation, not idential to actual covariance
log2(t) = sig_now2;

% ====================
m_old = m_now;
sig_old = sig_now;
end

figure(1); clf; hold on; box on;
plot(log1);
plot(log2);
