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
        % Nothing
end


%% two coal units
coal_num = 2;
dx1 = 200; % 1000x1000 = 4.33GB RAM
dx2 = 200;
u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
[u1_mesh, u2_mesh] = meshgrid(u_coal_unit1, u_coal_unit2);
u = u1_mesh + u2_mesh;

f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

u_unique = (0.4*coal_nameplate*coal_num):1:(coal_nameplate*coal_num);
v_unique = u_unique - 0.08*coal_nameplate*coal_num;

tic;
opt_f = -1*ones(1, length(v_unique));
opt_u1 = -1*ones(1, length(v_unique));
opt_u2 = -1*ones(1, length(v_unique));
opt_v1 = -1*ones(1, length(v_unique));
opt_v2 = -1*ones(1, length(v_unique));
for i = 1:length(v_unique)
    q = v_unique(i);
    u1 = linspace(min(u_coal_unit1),max(u_coal_unit1),2000);
    v1 = u1 - 0.08*coal_nameplate;
    v2 = q - v1;
    u2 = v2 + 0.08*coal_nameplate;
    f = interp2(u_coal_unit1, u_coal_unit2, f_coal, u1, u2);
    if any(~isnan(f))
        [value,id] = min(f);
        opt_u1(i) = u1(id);
        opt_u2(i) = u2(id);
        opt_v1(i) = v1(id);
        opt_v2(i) = v2(id);
        opt_f(i) = value;
    end
end
exclud_id = find(opt_f==-1);
u_unique(exclud_id) = [];
v_unique(exclud_id) = [];
opt_u1(exclud_id) = [];
opt_u2(exclud_id) = [];
opt_f(exclud_id) = [];
toc;

% save('TwoUnits', ...
%      'v_unique', 'opt_f', ...
%      'opt_u1', 'opt_u2', ...
%      'opt_v1', 'opt_v2', ...
%      'dx1', 'dx2');

figure(3); clf; hold on; box on;
ha = area(v_unique, [opt_u1; opt_u2]');
plot(v_unique, v_unique, 'linewidth', 1);
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
xlabel('Power Output, MW (exclude in-house generation)')
ylabel('Power Production (MW)');
xlim([400 1250]);
legend('u1', 'u2', 'Useful Output');
set(legend, 'location', 'northwest');

figure(31); clf; hold on; box on;
opt_f1 = interp1(u_coal_unit1, f_coal_unit1, opt_u1);
opt_f2 = interp1(u_coal_unit2, f_coal_unit2, opt_u2);
ha = area(v_unique, [opt_f1; opt_f2]');
set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
xlabel('Output, MWh (exclude in-house generation)')
ylabel('Hourly Coal Consumption (ton)');
xlim([400 1250]);
legend(ha, 'u1', 'u2');
set(legend, 'location', 'northwest');


%% ========================================================================
plot_switch = 'off'; % on/off
switch plot_switch
    case 'on' % don't plot
dx1 = 10; % very course grid for plotting
dx2 = 15;
u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
[u1_mesh, u2_mesh] = meshgrid(u_coal_unit1, u_coal_unit2);
u = u1_mesh + u2_mesh;
f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*(linspace(0.4,1,dx1)*coal_nameplate/1e3); %[g/kWh]
f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*(linspace(0.4,1,dx2)*coal_nameplate/1e3); %[g/kWh]
[f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
f_coal = f1 + f2;

% ==============================
figure(4); clf;
mesh(u_coal_unit1, u_coal_unit2, u); hold on;
plot3(opt_u1, opt_u2, u_unique, 'o-', 'linewidth', 1, 'markersize', 4);
xlabel('u1 (MWh)');
ylabel('u2 (MWh)');
zlabel('Combined Power Output (MWh)');
view(315, 35);
set(gcf, 'unit', 'inch', 'pos', [1.0    5.85    5.0000    3.7500]);

figure(41); clf;
mesh(u_coal_unit1, u_coal_unit2, u); hold on;
plot3(opt_u1, opt_u2, u_unique, 'o-', 'linewidth', 1, 'markersize', 4);
xlabel('u1 (MW)');
ylabel('u2 (MW)');
zlabel('Combined Power Output (MW)');
axis equal
xlim([250 680])
ylim([250 680])
view(0, 90);
set(gcf, 'unit', 'inch', 'pos', [1.0    1.09    5.0000    3.7500]);

% ==============================
figure(5); clf;
mesh(u_coal_unit1, u_coal_unit2, f_coal); hold on;
surface([0.4 1]*coal_nameplate, [0.4 1]*coal_nameplate, ones(2,2)*299.4, 'facec', [1 0.5 0.5], 'edgecolor', 'none');
plot3([0.44 1]*coal_nameplate, [1 0.44]*coal_nameplate, ones(1,2)*299.4, 'color', 'k')
alpha(0.7);
xlabel('u1 (MWh)');
ylabel('u2 (MWh)');
zlabel('Combined Fuel Consumption (ton)');
% view(315, 45);
view(-110, 25);
set(gcf, 'unit', 'inch', 'pos', [6.2    5.85    5.0000    3.7500]);

figure(51); clf;
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

figure(52); clf;
mesh(u_coal_unit1, u_coal_unit2, f_coal); hold on;
plot3(opt_u1, opt_u2, opt_f, 'o-', 'linewidth', 1, 'markersize', 4);
alpha(0.8);
xlabel('u1 (MWh)');
ylabel('u2 (MWh)');
zlabel('Combined Fuel Consumption (ton)');
view(315, 15);
set(gcf, 'unit', 'inch', 'pos', [11.4    5.85    5.0000    3.7500]);

% ==============================
figure(6); clf; hold on; box on;
plot([200 800], [200 800], '--', 'color', [1 1 1]*0.8);
plot(opt_u1, opt_u2, 'o-', 'markersize', 3, 'color', [0 0.5647 0.7412]);
xlabel('Optimal u1');
ylabel('Optimal u2');
set(gcf, 'unit', 'inch', 'pos', [11.4    1.0900    5.0000    3.7500]);
axis equal;
axis([200 800 200 800]);
my_gridline;
    otherwise % don't plot
        % Nothing
end


%% Another way to optimize the two-units case (but with dx1 = 1000, solution still zigzags because the cost difference is too small)
% dx1 = 1000; % 4.33GB RAM
% dx2 = 1000;
% u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
% u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
% [u1_mesh, u2_mesh] = meshgrid(u_coal_unit1, u_coal_unit2);
% u = u1_mesh + u2_mesh;
% 
% v_coal_unit1 = u_coal_unit1-0.08*coal_nameplate;
% v_coal_unit2 = u_coal_unit2-0.08*coal_nameplate;
% [v1_mesh, v2_mesh] = meshgrid(v_coal_unit1, v_coal_unit2);
% v = v1_mesh + v2_mesh;
% v_unique = unique(v);
% 
% f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
% f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3; %[g/kWh]
% [f1,f2] = meshgrid(f_coal_unit1,f_coal_unit2);
% f_coal = f1 + f2;
% 
% tic;
% for i = 1:length(v_unique) 
%     id_extract = find(v(:)==v_unique(i));
%     
%     if length(id_extract) > 1
%         f = f_coal(id_extract);
%         u1_extract = u1_mesh(id_extract);
%         u2_extract = u2_mesh(id_extract);
%         v1_extract = v1_mesh(id_extract);
%         v2_extract = v2_mesh(id_extract);
%         
%         id2 = u2_extract >= u1_extract;
%         check_valid = id2;
%         id_valid = find(check_valid == 1);
%         
%         if isempty(id_valid)
%             disp('error!!!');
%         elseif length(id_valid) == 1
%             opt_f(i) = f(id_valid);
%             opt_u1(i) = u1_mesh(id_extract(id_valid));
%             opt_u2(i) = u2_mesh(id_extract(id_valid));
%             opt_v1(i) = v1_mesh(id_extract(id_valid));
%             opt_v2(i) = v2_mesh(id_extract(id_valid));
%         else % length(id_valid) > 1
%             f_valid = f(id_valid);
%             [value, id_opt] = min(f_valid);
%             opt_f(i) = value;
%             opt_u1(i) = u1_mesh(id_extract(id_valid(id_opt)));
%             opt_u2(i) = u2_mesh(id_extract(id_valid(id_opt)));
%             opt_v1(i) = v1_mesh(id_extract(id_valid(id_opt)));
%             opt_v2(i) = v2_mesh(id_extract(id_valid(id_opt)));
%         end
%     else
%         opt_f(i) = f_coal(id_extract);
%         opt_u1(i) = u1_mesh(id_extract);
%         opt_u2(i) = u2_mesh(id_extract);
%         opt_v1(i) = v1_mesh(id_extract);
%         opt_v2(i) = v2_mesh(id_extract);
%     end    
% end
% tt = toc
% 
% figure(7); clf; hold on; box on;
% ha = area(v_unique, [opt_u1; opt_u2]');
% plot(v_unique, v_unique, 'linewidth', 1);
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% xlabel('Power Output, MW (exclude in-house generation)')
% ylabel('Power Production (MW)');
% xlim([400 1250]);
% legend('u1', 'u2', 'Useful Output');
% set(legend, 'location', 'northwest');


%% three coal units
% dx1 = 100;
% dx2 = 100;
% dx3 = 100;
% u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
% u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
% u_coal_unit3 = linspace(0.4,1,dx3) * coal_nameplate;
% u_coal_unit1 = roundn(u_coal_unit1, -1);
% u_coal_unit2 = roundn(u_coal_unit2, -1);
% u_coal_unit3 = roundn(u_coal_unit3, -1);
% [u1_mesh, u2_mesh, u3_mesh] = meshgrid(u_coal_unit1, u_coal_unit2, u_coal_unit3);
% 
% v_coal_unit1 = u_coal_unit1-0.08*coal_nameplate;
% v_coal_unit2 = u_coal_unit2-0.08*coal_nameplate;
% v_coal_unit3 = u_coal_unit3-0.08*coal_nameplate;
% v_coal_unit1 = roundn(v_coal_unit1, -1);
% v_coal_unit2 = roundn(v_coal_unit2, -1);
% v_coal_unit3 = roundn(v_coal_unit3, -1);
% [v1_mesh, v2_mesh, v3_mesh] = meshgrid(v_coal_unit1, v_coal_unit2, v_coal_unit3);
% v = v1_mesh + v2_mesh + v3_mesh;
% 
% f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
% f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3;
% f_coal_unit3 = (266*linspace(0.4,1,dx3).^2 -507*linspace(0.4,1,dx3) + 542).*linspace(0.4,1,dx3)*coal_nameplate/1e3;
% [f1, f2, f3] = meshgrid(f_coal_unit1, f_coal_unit2, f_coal_unit3);
% f_coal = f1 + f2 + f3;
% 
% v_unique = unique(v(:));
% opt_f = -1*ones(1, length(v_unique));
% opt_u1 = -1*ones(1, length(v_unique));
% opt_u2 = -1*ones(1, length(v_unique));
% opt_u3 = -1*ones(1, length(v_unique));
% opt_v1 = -1*ones(1, length(v_unique));
% opt_v2 = -1*ones(1, length(v_unique));
% opt_v3 = -1*ones(1, length(v_unique));
% 
% tic;
% for i = 1:length(v_unique)
%     id_extract = find(v(:)==v_unique(i));
%     
%     if length(id_extract) > 1
%         f = f_coal(id_extract);
%         u1_extract = u1_mesh(id_extract);
%         u2_extract = u2_mesh(id_extract);
%         u3_extract = u3_mesh(id_extract);
%         v1_extract = v1_mesh(id_extract);
%         v2_extract = v2_mesh(id_extract);
%         v3_extract = v3_mesh(id_extract);
%         
%         id2 = u2_extract >= u1_extract;
%         id3 = u3_extract >= u2_extract;
%         check_valid = id2 & id3;
%         id_valid = find(check_valid == 1);
%         
%         if isempty(id_valid)
%             disp('error!!!');
%         elseif length(id_valid) == 1
%             opt_f(i) = f(id_valid);
%             opt_u1(i) = u1_mesh(id_extract(id_valid));
%             opt_u2(i) = u2_mesh(id_extract(id_valid));
%             opt_u3(i) = u3_mesh(id_extract(id_valid));
%             opt_v1(i) = v1_mesh(id_extract(id_valid));
%             opt_v2(i) = v2_mesh(id_extract(id_valid));
%             opt_v3(i) = v3_mesh(id_extract(id_valid));
%         else % length(id_valid) > 1
%             f_valid = f(id_valid);
%             [value, id_opt] = min(f_valid);
%             opt_f(i) = value;
%             opt_u1(i) = u1_mesh(id_extract(id_valid(id_opt)));
%             opt_u2(i) = u2_mesh(id_extract(id_valid(id_opt)));
%             opt_u3(i) = u3_mesh(id_extract(id_valid(id_opt)));
%             opt_v1(i) = v1_mesh(id_extract(id_valid(id_opt)));
%             opt_v2(i) = v2_mesh(id_extract(id_valid(id_opt)));
%             opt_v3(i) = v3_mesh(id_extract(id_valid(id_opt)));
%         end
%     else
%         opt_f(i) = f_coal(id_extract);
%         opt_u1(i) = u1_mesh(id_extract);
%         opt_u2(i) = u2_mesh(id_extract);
%         opt_u3(i) = u3_mesh(id_extract);
%         opt_v1(i) = v1_mesh(id_extract);
%         opt_v2(i) = v2_mesh(id_extract);
%         opt_v3(i) = v3_mesh(id_extract);
%     end
% end
% tt = toc;
% 
% exclud_id = find(opt_f==-1);
% v_unique(exclud_id) = [];
% opt_u1(exclud_id) = [];
% opt_u2(exclud_id) = [];
% opt_u3(exclud_id) = [];
% opt_v1(exclud_id) = [];
% opt_v2(exclud_id) = [];
% opt_v3(exclud_id) = [];
% opt_f(exclud_id) = [];
% 
% if (dx1 == 100)
% save('ThreeUnits', ...
%      'v_unique', 'opt_f', ...
%      'opt_u1', 'opt_u2', 'opt_u3', ...
%      'opt_v1', 'opt_v2', 'opt_v3', ...
%      'dx1', 'dx2', 'dx3');
% end
% 
% % =================================
% figure(8); clf;  hold on; box on;
% ha = area(v_unique, [opt_u1; opt_u2; opt_u3]', 'edgecolor', 'none');
% plot(v_unique, opt_v1+opt_v2+opt_v3);
% set(ha(1), 'facec', [0.5 0.85 1]);
% set(ha(2), 'facec', [1 0.6 0.6]);
% set(ha(3), 'facec', [1 0.8 0]);
% xlim([600 1850]);
% title([num2str(dx1), 'x' num2str(dx2), 'x' num2str(dx3), ' (', num2str(tt, '%10.1f'), ' sec)']);
% xlabel('Power Output (MW) (exclude in-house use)');
% ylabel('Power Production (MW)');
% legend('u1', 'u2', 'u3', 'Useful Output');
% set(legend, 'location', 'northwest');


%% Four units
% tic;
% dx1 = 80; % use 8.6 GB RAM
% dx2 = 80;
% dx3 = 80;
% dx4 = 80;
% u_coal_unit1 = linspace(0.4,1,dx1) * coal_nameplate;
% u_coal_unit2 = linspace(0.4,1,dx2) * coal_nameplate;
% u_coal_unit3 = linspace(0.4,1,dx3) * coal_nameplate;
% u_coal_unit4 = linspace(0.4,1,dx4) * coal_nameplate;
% u_coal_unit1 = roundn(u_coal_unit1, -1);
% u_coal_unit2 = roundn(u_coal_unit2, -1);
% u_coal_unit3 = roundn(u_coal_unit3, -1);
% u_coal_unit4 = roundn(u_coal_unit4, -1);
% [u1_mesh, u2_mesh, u3_mesh, u4_mesh] = ndgrid(u_coal_unit1, u_coal_unit2, u_coal_unit3, u_coal_unit4);
% 
% v_coal_unit1 = u_coal_unit1-0.08*coal_nameplate;
% v_coal_unit2 = u_coal_unit2-0.08*coal_nameplate;
% v_coal_unit3 = u_coal_unit3-0.08*coal_nameplate;
% v_coal_unit4 = u_coal_unit4-0.08*coal_nameplate;
% v_coal_unit1 = roundn(v_coal_unit1, -1);
% v_coal_unit2 = roundn(v_coal_unit2, -1);
% v_coal_unit3 = roundn(v_coal_unit3, -1);
% v_coal_unit4 = roundn(v_coal_unit4, -1);
% [v1_mesh, v2_mesh, v3_mesh, v4_mesh] = ndgrid(v_coal_unit1, v_coal_unit2, v_coal_unit3, v_coal_unit4);
% v = v1_mesh + v2_mesh + v3_mesh + v4_mesh;
% 
% f_coal_unit1 = (266*linspace(0.4,1,dx1).^2 -507*linspace(0.4,1,dx1) + 542).*linspace(0.4,1,dx1)*coal_nameplate/1e3; %[g/kWh]
% f_coal_unit2 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx2) + 542).*linspace(0.4,1,dx2)*coal_nameplate/1e3;
% f_coal_unit3 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx3) + 542).*linspace(0.4,1,dx3)*coal_nameplate/1e3;
% f_coal_unit4 = (266*linspace(0.4,1,dx2).^2 -507*linspace(0.4,1,dx4) + 542).*linspace(0.4,1,dx4)*coal_nameplate/1e3;
% [f1, f2, f3, f4] = ndgrid(f_coal_unit1, f_coal_unit2, f_coal_unit3, f_coal_unit4);
% f_coal = f1 + f2 + f3 + f4;
% toc;
% 
% v_unique = unique(v(:));
% opt_u1 = -1*ones(1, length(v_unique));
% opt_u2 = -1*ones(1, length(v_unique));
% opt_u3 = -1*ones(1, length(v_unique));
% opt_u4 = -1*ones(1, length(v_unique));
% opt_v1 = -1*ones(1, length(v_unique));
% opt_v2 = -1*ones(1, length(v_unique));
% opt_v3 = -1*ones(1, length(v_unique));
% opt_v4 = -1*ones(1, length(v_unique));
% opt_f = -1*ones(1, length(v_unique));
% 
% tic;
% for i = 1:length(v_unique)
%     id_extract = find(v(:)==v_unique(i));
%     
%     if length(id_extract) > 1
%         f = f_coal(id_extract);
%         u1_extract = u1_mesh(id_extract);
%         u2_extract = u2_mesh(id_extract);
%         u3_extract = u3_mesh(id_extract);
%         u4_extract = u4_mesh(id_extract);
%         v1_extract = v1_mesh(id_extract);
%         v2_extract = v2_mesh(id_extract);
%         v3_extract = v3_mesh(id_extract);
%         v4_extract = v4_mesh(id_extract);
%         
%         id2 = u2_extract >= u1_extract;
%         id3 = u3_extract >= u2_extract;
%         id4 = u4_extract >= u3_extract;
%         check_valid = id2 & id3 & id4;
%         id_valid = find(check_valid == 1);
%         
%         if isempty(id_valid)
%             disp('error!!!');
%         elseif length(id_valid) == 1
%             opt_f(i) = f(id_valid);
%             opt_u1(i) = u1_mesh(id_extract(id_valid));
%             opt_u2(i) = u2_mesh(id_extract(id_valid));
%             opt_u3(i) = u3_mesh(id_extract(id_valid));
%             opt_u4(i) = u4_mesh(id_extract(id_valid));
%             opt_v1(i) = v1_mesh(id_extract(id_valid));
%             opt_v2(i) = v2_mesh(id_extract(id_valid));
%             opt_v3(i) = v3_mesh(id_extract(id_valid));
%             opt_v4(i) = v4_mesh(id_extract(id_valid));
%         else % length(id_valid) > 1
%             f_valid = f(id_valid);
%             [value, id_opt] = min(f_valid);
%             opt_f(i) = value;
%             opt_u1(i) = u1_mesh(id_extract(id_valid(id_opt)));
%             opt_u2(i) = u2_mesh(id_extract(id_valid(id_opt)));
%             opt_u3(i) = u3_mesh(id_extract(id_valid(id_opt)));
%             opt_u4(i) = u4_mesh(id_extract(id_valid(id_opt)));
%             opt_v1(i) = v1_mesh(id_extract(id_valid(id_opt)));
%             opt_v2(i) = v2_mesh(id_extract(id_valid(id_opt)));
%             opt_v3(i) = v3_mesh(id_extract(id_valid(id_opt)));
%             opt_v4(i) = v4_mesh(id_extract(id_valid(id_opt)));
%         end
%     else % length(id_extract) == 1
%         opt_f(i) = f_coal(id_extract);
%         opt_u1(i) = u1_mesh(id_extract);
%         opt_u2(i) = u2_mesh(id_extract);
%         opt_u3(i) = u3_mesh(id_extract);
%         opt_u4(i) = u4_mesh(id_extract);
%         opt_v1(i) = v1_mesh(id_extract);
%         opt_v2(i) = v2_mesh(id_extract);
%         opt_v3(i) = v3_mesh(id_extract);
%         opt_v4(i) = v4_mesh(id_extract);
%     end
% end
% tt = toc
% 
% exclud_id = find(opt_f==-1);
% v_unique(exclud_id) = [];
% opt_u1(exclud_id) = [];
% opt_u2(exclud_id) = [];
% opt_u3(exclud_id) = [];
% opt_u4(exclud_id) = [];
% opt_v1(exclud_id) = [];
% opt_v2(exclud_id) = [];
% opt_v3(exclud_id) = [];
% opt_v4(exclud_id) = [];
% opt_f(exclud_id) = [];
% 
% if (dx1==80)
% save('FourUnits', ...
%      'v_unique', 'opt_f', ...
%      'opt_u1', 'opt_u2', 'opt_u3', 'opt_u4', ...
%      'opt_v1', 'opt_v2', 'opt_v3', 'opt_v4', ...
%      'dx1', 'dx2', 'dx3', 'dx4');
% end
% 
% % =================================
% figure(9); clf; hold on; box on;
% ha = area(v_unique, [opt_u1; opt_u2; opt_u3; opt_u4]', 'edgecolor', 'none');
% plot(v_unique, opt_v1+opt_v2+opt_v3+opt_v4);
% set(ha(1), 'facec', [0.5 0.85 1]);
% set(ha(2), 'facec', [1 0.6 0.6]);
% set(ha(3), 'facec', [1 0.8 0]);
% xlim([800 2500]);
% title([num2str(dx1), 'x' num2str(dx2), 'x' num2str(dx3), ' (', num2str(tt, '%10.1f'), ' sec)']);
% xlabel('Power Output (MW) (exclude in-house use)');
% ylabel('Power Production (MW)');
% legend('u1', 'u2', 'u3', 'u4', 'Useful Output');
% set(legend, 'location', 'northwest');


