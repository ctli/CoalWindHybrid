% Test white noise

clear
close all
clc

y = wgn(1,10000,0); % white Gaussian noise

figure(1); clf;
hist(y);
