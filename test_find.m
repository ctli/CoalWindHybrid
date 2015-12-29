clear
close all
clc
format compact


%% Coal dispatch
load  MyopicDispatch;


%% Wind power
wind_file = 'Xilingol_2009';
load(wind_file);
wind_pwr = round(p*2500)';


%% Myopic dispatch
target_pwr = 8500;
coal_pwr = target_pwr - wind_pwr;
coal_pwr(coal_pwr<0) = 0; % [1x8760];

% cc = downsample(coal_pwr, 4000);
% [AA,BB] = meshgrid(coal_pwr, v_range);

rr = repmat(v_range, 6,1);


%% http://www.mathworks.com/matlabcentral/answers/108413-how-to-find-first-nonzero-element-first-1-per-row-and-set-other-elements-to-zero-without-loops-in
% input
  A = round(rand(5,6,2)) 
% engine
  isOne = A == 1 ;
  B = isOne & cumsum(isOne,2) == 1