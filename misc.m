p_HVDC = 9000; % [MW] line capacity of HVDC line: 9GW

yr_range = 1979:2009;

subregion_name = 'Xilingol';
lat = 43.9;
lon = 116;
color_code = [0 114 189]/255; % blue
wst_yr = 8;
bst_yr = 9;

v_cutin = 3; % [m/s]
v_cutout = 25; % [m/s]
pwr_curve = [
3	0.005993333
4	0.041946667
5	0.107866667
6	0.203746667
7	0.329586667
8	0.503373333
9	0.713106667
10	0.92884
10.61	1
];
