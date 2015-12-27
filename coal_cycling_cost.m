%% Cycling costs of coal power plants (super critial)
clear;
close all;
clc;
format compact;


%% Units NOT designed for cycling
% Start cost [$/MW]
% coal_startup = 39;  %hot start; 25th centile
% coal_startup = 54;  %hot start; median
% coal_startup = 63;  %hot start; 75th centile
% coal_startup = 54;  %warm start; 25th centile
% coal_startup = 64;  %warm start; median
% coal_startup = 89;  %warm start; 75th centile
% coal_startup = 73;  %cold start; 25th centile
% coal_startup = 104; %cold start; median
% coal_startup = 120; %cold start; 75th centile

% Load following [$/MW]
% coal_loadfollow = 1.96; %25th centile
% coal_loadfollow = 1.52; %median
% coal_loadfollow = 2.38; %75th centile

% Baseload variable cost [$/MWh]
% coal_baseload = 2.48; %25th centile
% coal_baseload = 2.96; %median
% coal_baseload = 3.40; %75th centile

%% Units designed for cycling
% Start cost [$/MW]
coal_startup = 38; %median; hot start
% coal_startup = 56; %median; warm start
% coal_startup = 99; %median; cold start

% Load following [$/MW]
coal_loadfollow = 1.72; %median

% Baseload variable cost [$/MWh]
coal_baseload = 3.22; %median

%% Other start up costs [$/MW]
coal_startup_oth = 5.81; %hot start
% coal_startup_oth = 8.62; %warm start
% coal_startup_oth = 11.58; %cold start

coal_startup = coal_startup + coal_startup_oth;


%% Starup fuel [mmBtu/MW]
% Hot start
coal_startup_fuel = 10.1; %hot start
% coal_startup_fuel = 17.1; %warm start
% coal_startup_fuel = 20.1; % cold start

%% Coal price (the two number blow are consistent)
% http://www.smh.com.au/business/energy/coal-prices-fall-to-12year-lows-as-china-india-join-demand-slowdown-20150819-gj2jk6.html
% coal_price = 52.85; %[$/metric ton]

% EIA short-term energy outlook: US annual coal price to electricity averaged at $2.37/MMBtu in 2014
coal_price = 2.37; %[$/MMBtu]

% Heat content of US coal
coal_heatcontent = 22.22; % [mmBtu/metric ton]

