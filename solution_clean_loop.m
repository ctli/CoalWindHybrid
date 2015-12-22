clear
close all
clc
format compact

coal_nameplate = 660; % [MW]

coal_num = 14; % 3 units

dx = 6000;
u_coal_unit = linspace(0.4,1,dx) * coal_nameplate;
v_coal_unit = u_coal_unit-0.08*coal_nameplate;
f_coal_unit = (266*linspace(0.4,1,dx).^2 -507*linspace(0.4,1,dx) + 542).*linspace(0.4,1,dx)*coal_nameplate/1e3; %[g/kWh]
u_min = min(u_coal_unit);
v_min = min(v_coal_unit);
f_min = min(f_coal_unit);


%% ========================================================================
% 211.2~8500.8MW -> round to 220-8500MW
v_range = 220:0.1:8500; % [1x415]
u_table = zeros(coal_num, length(v_range), coal_num);
v_table = zeros(coal_num, length(v_range), coal_num);
f_table = zeros(length(v_range), coal_num);
id_st = zeros(coal_num, 1); % Starting entry that is not nan
id_ed = zeros(coal_num, 1); % Ending entry that is not nan

figure(3); clf; hold on; box on;
% No unit
plot(0, 0, 'x');

% ====================
% One unit
n = 1;
u1_extended = interp1(v_coal_unit, u_coal_unit, v_range);
u_table(n,:,n) =  u1_extended;
v1_extended = interp1(v_coal_unit, v_coal_unit, v_range);
v_table(n,:,n) =  v1_extended;
f1_extended = interp1(v_coal_unit, f_coal_unit, v_range);
f_table(:,n) =  f1_extended;
id_st(n) = find(~isnan(f1_extended)==1, 1, 'first');
id_ed(n) = find(~isnan(f1_extended)==1, 1, 'last');

plot(v_coal_unit, f_coal_unit);
text(v_coal_unit(end), f_coal_unit(end), ' 1 Unit is commited', 'fontsize', 7);

% ====================
% Multiple units
% 2 units; 0.008600 sec
% 3 units; 0.013865 sec
% 4 units; 0.021214 sec
% 5 units; 0.030380 sec
% 10 units; 0.116095 sec
% 14 units; 0.241221 sec
tic;
for n = 2:coal_num
    v_long = v_coal_unit*n;

    u = ones(n, dx) * u_min;
    v = ones(n, dx) * v_min;
    f = ones(n, dx) * f_min;
    
    u_extended = zeros(n, dx, n);
    v_extended = zeros(n, dx, n);
    f_extended = zeros(n, dx, n);
    u_sum = zeros(n, dx);
    v_sum = zeros(n, dx);
    f_sum = zeros(n, dx);
    for cmt = 1:n
        u(cmt,:) = u_coal_unit;
        v(cmt,:) = v_coal_unit;
        f(cmt,:) = f_coal_unit;
        
        u_sum(cmt,:) = sum(u);
        v_sum(cmt,:) = sum(v);
        f_sum(cmt,:) = sum(f);
        
        for nn = 1:n
            u_extended(nn,:,cmt) = interp1(v_sum(cmt,:),u(nn,:),v_long);
            v_extended(nn,:,cmt) = interp1(v_sum(cmt,:),v(nn,:),v_long);
            f_extended(nn,:,cmt) = interp1(v_sum(cmt,:),f(nn,:),v_long);
        end
    end
    f_combined = squeeze(sum(f_extended))';
    [value, id_row] = min(f_combined);
    
    id_col = 1:length(id_row);
    id = sub2ind([n,dx], id_row, id_col);
    opt_f = value; % [1]x[dx]
    opt_u = zeros(n,dx);
    opt_v = zeros(n,dx);
    for nn = 1:n
        u1_extended = squeeze(u_extended(nn,:,:))';
        opt_u(nn,:) = u1_extended(id);
        v1_extended = squeeze(u_extended(nn,:,:))';
        opt_v(nn,:) = v1_extended(id);
        
        % Save opt solutions to tables
        u1_extended = interp1(v_long, opt_u(nn,:), v_range);
        u_table(nn,:,n) = u1_extended;
        v1_extended = interp1(v_long, opt_v(nn,:), v_range);
        v_table(nn,:,n) = v1_extended;
    end
    f_extended = interp1(v_long, opt_f, v_range);
    f_table(:,n) = f_extended;
    id_st(n) = find(~isnan(f_extended)==1, 1, 'first');
    id_ed(n) = find(~isnan(f_extended)==1, 1, 'last');
    
    id_bad = isnan(opt_f);
    v_long(id_bad) = [];
    opt_f(id_bad) = [];
    opt_u(:,id_bad) = [];
    opt_v(:,id_bad) = [];
    plot(v_long, opt_f);
    if n<11
    text(v_long(end), opt_f(end), [' ', num2str(n), ' Units are commited'], 'fontsize', 7);
    else
    text(v_long(1), opt_f(1), [' ', num2str(n), ' Units are commited '], 'fontsize', 7, 'horizontalalignment', 'right');
    end
    toc;
end
v_st = v_range(id_st);
v_ed = v_range(id_ed);
xlabel('Output Power, MW (in-house use excluded)');
ylabel('Coal Consumption (ton/h)');
my_gridline;



% ====================
% flag_badrow = boolean(zeros(1,length(v_range)));
% for vv = 1:length(v_range)
%     f_row = f_table(vv,:);
%     if all(isnan(f_row))
%         flag_badrow(vv) = 1;
%     end
% end
% f_table(flag_badrow,:) = [];
% u_table(:,flag_badrow,:) = [];
% v_table(:,flag_badrow,:) = [];
% v_range(flag_badrow) = [];
% save('OptTable', ...
%      'v_range', ...
%      'f_table', 'u_table', 'v_table', ...
%      'id_st', 'id_ed', 'v_st', 'v_ed');

% ====================
% figure(1); clf;
% opt_u = flipud(opt_u);
% for nn = 1:n
%     subplot(n+1,1,nn);
%     plot(v_long, opt_u(nn,:));
%     ylim([0 700]);
%     set(gca, 'ytick', 0:350:700);
%     ylabel(['u', num2str(nn)]);
% end
% subplot(n+1,1,n+1);
% area(v_long, opt_u', 'edgecolor', 'none');
% ylim([0 4000]);


%% Double check 3-units result to be identical to solution_clean.m
% figure(2); clf;
% subplot(4,1,1); hold on; box on;
% plot(v_long, opt_u(1,:), 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u1');
% 
% subplot(4,1,2); hold on; box on;
% plot(v_long, opt_u(2,:), 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u2');
% 
% subplot(4,1,3); hold on; box on;
% plot(v_long, opt_u(3,:), 'linewidth', 1); hold on;
% ylim([0 700]);
% ylabel('u3');
% 
% subplot(4,1,4); hold on; box on;
% ha = area(v_long, opt_u', 'edgecolor', 'none');
% set(ha(1), 'facec', [0.8 0.95 1], 'edgecolor', 'none');
% set(ha(2), 'facec', [1 0.8 0.8], 'edgecolor', 'none');
% set(ha(3), 'facec', [1 0.8 0], 'edgecolor', 'none');
% legend('u1', 'u2', 'u3');
% set(legend, 'location', 'northwest');
% 
% load ThreeUnitsClean;
% subplot(4,1,3);
% plot(v_unique, opt_u1, '--'); hold on;
% subplot(4,1,2);
% plot(v_unique, opt_u2, '--'); hold on;
% subplot(4,1,1);
% plot(v_unique, opt_u3, '--'); hold on;


