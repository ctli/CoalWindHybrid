clear
close all
clc
format compact


%% Coal power plants
coal_num = 14;
load OptTable;

% Myopic unit commitment
f_myopic = zeros(1,length(v_range));
cmt_myopic = zeros(1,length(v_range));
tic;
for vv = 1:length(v_range)
    f_column = f_table(vv,:);
    [f_myopic(vv), cmt_myopic(vv)] = min(f_column);
end
toc;

figure(1); clf; hold on; box on;
plot(0,0,'x');
plot(v_range, f_table, 'linewidth', 1);
for n = 1:coal_num
    if n==1
        text(v_ed(n), max(f_table(:,n)), [' ', num2str(n), ' Unit is commited'], 'fontsize', 7);
    elseif n<11
        text(v_ed(n), max(f_table(:,n)), [' ', num2str(n), ' Units are commited'], 'fontsize', 7);
    else
        text(v_st(n), min(f_table(:,n)), [' ', num2str(n), ' Units are commited '], 'fontsize', 7, 'horizontalalignment', 'right');
    end
end
xlabel('Output Power, MW (in-house use excluded)');
ylabel('Coal Consumption (ton/h)');
my_gridline;
% h1 = plot(v_range, f_myopic, 'color', [1 1 1]*0, 'linewidth', 0.35);
% legend(h1, 'Myopic Dispatch');
% set(legend, 'location', 'northwest');
% 
% % export_fig myopic -r300
% 
% figure(2); clf;
% plot(v_range, cmt_num);
% ylim([0 15]);
% xlabel('Output Power, MW (in-house use excluded)');
% ylabel('Number of Plants Dispatched (Count)');
% my_gridline;


%%
load Xilingol_2009;
wind_pwr = p*2500;

target_pwr = 8499;
coal_pwr = target_pwr - wind_pwr; % Use coal to make up deficit


%% Myopic dispatch
f_min = interp1(v_range, f_myopic, coal_pwr);
cmt_min = zeros(1,length(coal_pwr));
v_min = zeros(1,length(coal_pwr));
for t = 1:length(coal_pwr)
    id = find(v_range>=coal_pwr(t), 1, 'first');
    f_min(t) = f_myopic(id);
    cmt_min(t) = cmt_myopic(id);
    v_min(t) = v_range(id);
end
plot(coal_pwr, f_min, 'x');


%%
figure(2); clf;
ha = area([coal_pwr, wind_pwr], 'edgecolor', 'none');
set(ha(1), 'facec', [1 0.4 0.4]);
set(ha(2), 'facec', [0.6 1 0]);
xlim([ 0 168]);


