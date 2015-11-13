clear
close all
clc
format compact

coal_nameplate = 660; % [MW]

plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
coal_fuelrate_x = 0.4:0.01:1; % [normalized]
coal_fuelrate_y = 266*coal_fuelrate_x.^2 -507*coal_fuelrate_x + 542; %[g/kWh]
coal_fuel_x  = coal_fuelrate_x*coal_nameplate; % [MW]
coal_fuel_y = (coal_fuelrate_y*1000).*(coal_fuelrate_x*coal_nameplate)/1e6; % [g/h] -> [ton/h]

% ==============================
% fuel rate
figure(1); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*392, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuelrate_y, 'linewidth', 1);
xlim([-15 710]);
ylim([298 392]);
my_gridline([1 1 1]*0.85, 'front');
set(gca, 'layer', 'top');
set(gca, 'fontsize', 9);
xlabel('Output Power (MW)');
ylabel('Coal Consumption Rate (g/kWh)', 'fontweight', 'bold', 'fontsize', 12);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);
set(gca, 'pos', [0.1300    0.115    0.7750    0.77]);
legend('Working Range', 'Fuel Rate');
set(legend, 'pos', [0.4924    0.75    0.2890    0.0950]);
text(0.4*coal_nameplate, mean(get(gca, 'ylim')), '40% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);
text(1.0*coal_nameplate, mean(get(gca, 'ylim')), '100% of Nameplate Capacity    ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);

ax_pos = get(gca, 'pos');
x_lim = get(gca, 'xlim');
y_lim = get(gca, 'ylim');
y_tick = get(gca, 'ytick');

ax2 = axes;
set(ax2, 'pos', ax_pos);
set(ax2, 'color', 'none');
set(ax2, 'fontsize', 9);
x_ticklabel = 0:0.15:1.5;
set(ax2, 'xaxislocation', 'top', 'xlim', x_lim, 'xtick', x_ticklabel*coal_nameplate, 'xticklabel', x_ticklabel);
set(ax2, 'yaxislocation', 'right', 'ylim', y_lim, 'ytick', y_tick, 'yticklabel', []);
xlabel('Output Power (Normalized)', 'fontsize', 9);

% ==============================
% fuel consumption
figure(2); clf; hold on;
area([0.4 1]*coal_nameplate, [1 1]*202, 'facec', [0.85 0.98 1], 'edgecolor', 'none');
plot(coal_fuel_x, coal_fuel_y, 'linewidth', 1);
plot([coal_fuel_x(1),coal_fuel_x(end)],[coal_fuel_y(1),coal_fuel_y(end)], '--', 'color', [0.5 0.5 0.5]);
ylim([98 202]);
xlim([-15 710]);

xa = 507/(266*2)*660; % min fuel rate @ 0.953 (i.e. 419MW)
ya = interp1(coal_fuel_x, coal_fuel_y, xa);
plot(xa, ya, 'ko', 'markerf', 'k', 'markersize', 3);
text(xa+0.02, ya+7, 'Min. Fuel Rate      ', 'horizontalalignment', 'right');
text(xa+0.02, ya+2, ['@ ', num2str(xa, '%3.0f'), ' ( (y/x)^{\prime}=0 ) '], 'horizontalalignment', 'right');

xb = 507*2/(266*3*2)*660; % inflection point @ 0.6353 (i.e. 629MW)
yb = interp1(coal_fuel_x, coal_fuel_y, xb);
plot(xb, yb, 'ko', 'markerf', 'k', 'markersize', 3);
text(xb+0.02, yb-4, 'Inflection Point');
text(xb+0.02, yb-8, ['@ ', num2str(xb, '%3.0f'), ' ( y^{\prime\prime}=0 )']);

text(0.4*coal_nameplate, mean(get(gca, 'ylim')), '40% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);
text(1.0*coal_nameplate, mean(get(gca, 'ylim')), '100% of Nameplate Capacity     ', 'rotation', 90, 'horizontalalignment', 'center', 'color', [0.5 0.5 0.5], 'fontsize', 9);

my_gridline([1 1 1]*0.85, 'front');
set(gca, 'fontsize', 9);
xlabel('Generation (MW)');
ylabel('Coal Consumption (ton/h)', 'fontweight', 'bold', 'fontsize', 12);
set(gcf, 'unit', 'inch', 'pos', [1.0    1.09    5.0000    3.7500]);
set(gca, 'pos', [0.1300    0.115    0.7750    0.77]);

ax_pos = get(gca, 'pos');
x_lim = get(gca, 'xlim');
y_lim = get(gca, 'ylim');
y_tick = get(gca, 'ytick');

ax2 = axes;
set(ax2, 'pos', ax_pos);
set(ax2, 'color', 'none');
set(ax2, 'fontsize', 9);
x_ticklabel = 0:0.15:1.5;
set(ax2, 'xaxislocation', 'top', 'xlim', x_lim, 'xtick', x_ticklabel*coal_nameplate, 'xticklabel', x_ticklabel);
set(ax2, 'yaxislocation', 'right', 'ylim', y_lim, 'ytick', y_tick, 'yticklabel', []);
xlabel('Output Power (Normalized)', 'fontsize', 9);
    case 'off'
        ;
end


%% two coal units
% coal_num = 2;
% dx1 = 1000;
% dx2 = 1001;
% u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
% u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
% f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
% f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %[g/kWh]
% [f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
% f_coal = f1 + f2;
% 
% tic;
% uu_range = (0.4*coal_nameplate*coal_num):1:(coal_nameplate*coal_num);
% vv_range = uu_range - 0.08*coal_nameplate*coal_num;
% opt_u1 = -1*ones(1,length(uu_range));
% opt_u2 = -1*ones(1,length(uu_range));
% opt_f = -1*ones(1,length(uu_range));
% flag = -1*ones(1,length(uu_range)); % 1: pass; -1: fail
% for i = 1:length(uu_range)
%     q = vv_range(i);
%     u1 = linspace(min(u_coal_unit1),max(u_coal_unit1),2000);
%     v1 = u1 - 0.08*coal_nameplate;
%     v2 = q - v1;
%     u2 = v2 + 0.08*coal_nameplate;
%     f = interp2(u_coal_unit1, u_coal_unit2, f_coal, u1, u2);
%     if any(~isnan(f))
%         [value,id] = min(f);
%         opt_u1(i) = u1(id);
%         opt_u2(i) = u2(id);
%         opt_f(i) = value;
%         flag(i) = 1;
%     end
% end
% exclud_id = find(flag==-1);
% uu_range(exclud_id) = [];
% vv_range(exclud_id) = [];
% opt_u1(exclud_id) = [];
% opt_u2(exclud_id) = [];
% opt_f(exclud_id) = [];
% toc;

% ========================================================================
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on' % don't plot
dx1 = 10; % very course grid for plotting
dx2 = 15;
u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
[p1,p2] = meshgrid(u_coal_unit1,u_coal_unit2);
p_coal = p1 + p2;
f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*(linspace(0.4,1,dx1)*coal_nameplate/1e3); %[g/kWh]
f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*(linspace(0.4,1,dx2)*coal_nameplate/1e3); %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

% ==============================
figure(3); clf;
mesh(u_coal_unit1, u_coal_unit2, p_coal); hold on;
plot3(opt_u1, opt_u2, uu_range, 'o-', 'linewidth', 1, 'markersize', 4);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Power Output (MW)');
view(315, 35);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);

figure(31); clf;
mesh(u_coal_unit1, u_coal_unit2, p_coal); hold on;
plot3(opt_u1, opt_u2, uu_range, 'o-', 'linewidth', 1, 'markersize', 4);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Power Output (MW)');
axis equal
xlim([250 680])
ylim([250 680])
view(0, 90);
set(gcf, 'unit', 'inch', 'pos', [1.0    1.09    5.0000    3.7500]);

% ==============================
figure(4); clf;
mesh(u_coal_unit1, u_coal_unit2, f_coal); hold on;
surface([0.4 1]*coal_nameplate, [0.4 1]*coal_nameplate, ones(2,2)*299.4, 'facec', [1 0.5 0.5], 'edgecolor', 'none');
plot3([0.44 1]*coal_nameplate, [1 0.44]*coal_nameplate, ones(1,2)*299.4, 'color', 'k')
alpha(0.7);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Fuel Consumption (ton/h)');
% view(315, 45);
view(-110, 25);
set(gcf, 'unit', 'inch', 'pos', [6.2    5.85    5.0000    3.7500]);

figure(41); clf;
mesh(u_coal_unit1, u_coal_unit2, f_coal); hold on;
surface([0.4 1]*coal_nameplate, [0.4 1]*coal_nameplate, ones(2,2)*299.4, 'facec', [1 0.5 0.5], 'edgecolor', 'none');
plot3([0.44 1]*coal_nameplate, [1 0.44]*coal_nameplate, ones(1,2)*400, 'color', 'k', 'linewidth', 1)
alpha(0.7);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Fuel Consumption (ton/h)');
axis equal
xlim([250 680])
ylim([250 680])
view(0, 90);
set(gcf, 'units', 'inch', 'pos', [6.25  1.0900   5   3.75]);

figure(42); clf;
mesh(u_coal_unit1, u_coal_unit2, f_coal); hold on;
plot3(opt_u1, opt_u2, opt_f, 'o-', 'linewidth', 1, 'markersize', 4);
alpha(0.8);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Fuel Consumption (ton/h)');
view(315, 15);
set(gcf, 'unit', 'inch', 'pos', [11.4    5.85    5.0000    3.7500]);

% ==============================
figure(5); clf; hold on; box on;
plot([200 800], [200 800], '-', 'color', [0.7 0.7 0.7]);
plot(opt_u1, opt_u2, 'o-', 'markersize', 3);
xlabel('Optimal u1');
ylabel('Optimal u2');
set(gcf, 'unit', 'inch', 'pos', [11.4    1.0900    5.0000    3.7500]);
axis equal;
axis([200 800 200 800]);
my_gridline;
    otherwise % don't plot
        ;
end


% ========================================================================
% one unit vs. two units
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on'
p0 = 0;
q0 = 0;
f0 = 0;

p1 = u_coal_unit1;
q1 = p1 - 0.08*coal_nameplate;
f1 = f_coal_unit1;

p2 = uu_range;
q2 = vv_range;
f2 = opt_f;

% ========
q = 0:1:vv_range(end);
f = zeros(1,length(q));
u1_log = zeros(1, length(q));
u2_log = zeros(1, length(q));
f_u1 = zeros(1, length(q));
f_u2 = zeros(1, length(q));

id1 = (q>0 & q<=q1(end));
f(id1) = interp1(q1, f1, q(id1));
u1_log(id1) = interp1(q1, p1, q(id1));
f_u1(id1) = interp1(p1, f1, u1_log(id1));

id2 = (q>q1(end));
f(id2) = interp1(q2, f2, q(id2));
u1_log(id2) = interp1(q2, opt_u1, q(id2));
u2_log(id2) = interp1(q2, opt_u2, q(id2));
f_u1(id2) = interp1(p1, f1, u1_log(id2));
f_u2(id2) = interp1(p1, f1, u2_log(id2));

v1_log = u1_log - 0.08*coal_nameplate;
v2_log = zeros(1,length(q));
v2_log(u2_log>0) = u2_log(u2_log>0) - 0.08*coal_nameplate;

% ================================
figure(6); clf; hold on; box on;
plot(q0, f0, 'x', q1, f1, '-', q2, opt_f, '-', 'linewidth', 1);

xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Coal Consumption (ton/hr)');
ylim([0 400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend('No Unit is Commited', 'One Unit is Commited', 'Two Units are Commited');
set(legend, 'location', 'northwest');

% ================================
figure(7); clf; hold on; box on;
plot(q0, f0, 'x', q1, f1, '-', q2, opt_f, '-', 'linewidth', 1);
% area(q,f, 'facec', [0.8 0.95 1], 'edgecolor', 'none');
ha = area(q, [f_u1;f_u2]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');

xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Coal Consumption (ton/hr)');
ylim([0 400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend(ha, '1st Unit', '2nd Unit');
set(legend, 'location', 'northwest');

% ================================
figure(8); clf; hold on; box on;
ha = area(q, [u1_log;u2_log]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
plot(q, v1_log + v2_log, 'linewidth', 1);
xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Power Output (MW)');
ylim([0 1400]);
xlim([0 1250]);
grid on;
set(gca, 'layer', 'top');
legend('1st Unit', '2nd Unit', 'Useful Output');
set(legend, 'location', 'northwest');
    case 'off'
        ;
end


%% three coal units
coal_num = 3;
dx1 = 100;
dx2 = 100;
dx3 = 100;
u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
u_coal_unit3 = linspace(0.4,1,dx3) * coal_nameplate;
[u1_mesh, u2_mesh, u3_mesh] = meshgrid(u_coal_unit1, u_coal_unit2, u_coal_unit3);

v_coal_unit1 = u_coal_unit1-0.08*coal_nameplate;
v_coal_unit2 = u_coal_unit2-0.08*coal_nameplate;
v_coal_unit3 = u_coal_unit3-0.08*coal_nameplate;
[v1_mesh, v2_mesh, v3_mesh] = meshgrid(v_coal_unit1, v_coal_unit2, v_coal_unit3);
v = v1_mesh + v2_mesh + v3_mesh;
id2 = v2_mesh>v1_mesh;
id3 = v3_mesh>v2_mesh;

f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %create tiny difference to eliminate non-unique solutions
f_coal_unit3 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %create tiny differenct
[f1,f2,f3] = meshgrid(f_coal_unit1,f_coal_unit2,f_coal_unit3);
f_coal = f1 + f2 + f3;

v_unique = unique(v(:)); % 51x1
opt_u1 = -1*ones(1, length(v_unique));
opt_u2 = -1*ones(1, length(v_unique));
opt_u3 = -1*ones(1, length(v_unique));
opt_f = -1*ones(1, length(v_unique));

tic;
for i = 1:length(v_unique)
    id_extract = find(v(:)==v_unique(i));
    
    if length(id_extract) > 1
        f = f_coal(id_extract);
        u1_extract = u1_mesh(id_extract);
        u2_extract = u2_mesh(id_extract);
        u3_extract = u3_mesh(id_extract);
        
        id2 = u2_extract >= u1_extract;
        id3 = u3_extract >= u2_extract;
        check_valid = id2 & id3;
        id_valid = find(check_valid == 1);
        
        if isempty(id_valid)
            disp('error!!!');
        elseif length(id_valid) == 1
            opt_f(i) = f(id_valid);
            opt_u1(i) = u1_mesh(id_extract(id_valid));
            opt_u2(i) = u2_mesh(id_extract(id_valid));
            opt_u3(i) = u3_mesh(id_extract(id_valid));
        else % length(id_valid) > 1
            f_valid = f(id_valid);
            [value, id_opt] = min(f_valid);
            opt_f(i) = value;
            opt_u1(i) = u1_mesh(id_extract(id_valid(id_opt)));
            opt_u2(i) = u2_mesh(id_extract(id_valid(id_opt)));
            opt_u3(i) = u3_mesh(id_extract(id_valid(id_opt)));
        end
    else
        opt_f(i) = f_coal(id_extract);
        opt_u1(i) = u1_mesh(id_extract);
        opt_u2(i) = u2_mesh(id_extract);
        opt_u3(i) = u3_mesh(id_extract);
    end
end
tt = toc;

exclud_id = find(opt_f==-1);
v_unique(exclud_id) = [];
opt_u1(exclud_id) = [];
opt_u2(exclud_id) = [];
opt_u3(exclud_id) = [];
opt_f(exclud_id) = [];

%%
figure(9); clf;
ha = area(v_unique, [opt_u1; opt_u2; opt_u3]', 'edgecolor', 'none');
set(ha(1), 'facec', [0.5 0.85 1]);
set(ha(2), 'facec', [1 0.6 0.6]);
set(ha(3), 'facec', [1 0.8 0]);
xlim([600 1850]);
title([num2str(dx1), 'x' num2str(dx2), 'x' num2str(dx3), ' (', num2str(tt, '%10.1f'), ' sec)']);
xlabel('Power Output (MW) (exclude in-house use)');
ylabel('Power Production (MW)');
legend('u1', 'u2', 'u3');
set(legend, 'location', 'northwest');

