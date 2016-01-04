clear
clc

rng(1);
esp = rand(10000,1); % 0~1
resp = round(esp);  % 0 or 1
s = -1 + 2*resp;    %-1 or 1
walk = cumsum(s);

figure(1); clf;
plot(walk);

mu = mean(walk)
rho = std(walk)
