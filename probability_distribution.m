% Derive probability distribution
clear
close all
clc

yr_range = 1979:2009;
file_prefix = 'Xilingol_';

p_log = nan*ones(8760, length(yr_range));
v_log = nan*ones(8760, length(yr_range));
for iy = 1:length(yr_range)
    yr = yr_range(iy);
    file_name = [file_prefix, num2str(yr)];
    load(file_name);
    if length(p)>=8760;
        p_log(:,iy) = p(1:8760);
    end
    if length(v)>=8760;
        v_log(:,iy) = v(1:8760);
    end
end


%% Wind speed probability distribution
v_bin = linspace(0,25,21);

figure(1); clf;
set(gcf, 'units', 'inch', 'pos', [0.45    1.4583    5.8333    6.75]);

% ====================
% PDF
subplot(3,1,1);
v_cnt = hist(v_log(:), v_bin);
bar(v_bin, v_cnt/sum(v_cnt), 1, 'facec', [0.35 0.65 1]);
xlim([-1 26]);
ylabel('Probability (-)');
set(gca, 'layer', 'top');
set(gca, 'ygrid', 'on');

% ====================
% CDF
subplot(3,1,2);
xv = sort(v_log(:));
Fv = (1:length(v_log(:)))/length(v_log(:));
plot(xv,Fv, 'linewidth', 1);
xlim([-1 26]);
ylim([-0.05 1.05]);
set(gca, 'ytick', 0:0.25:1);
grid on;
set(gca, 'tickdir', 'out');
ylabel('CDF (-)');

% ====================
% Realization
subplot(3,1,3); hold on; box on;
Fv_r = linspace(0,1,10000);
xv_r = interp1(Fv,xv,Fv_r);
xv_r = xv_r(~isnan(xv_r));
vr_cnt = hist(xv_r, v_bin);
bar(v_bin, vr_cnt/sum(vr_cnt), 1, 'facec', [1 0.4 0.4]);

Fv_r = linspace(0,1,50);
xv_r = interp1(Fv,xv,Fv_r);
xv_r = xv_r(~isnan(xv_r));
vr_cnt = hist(xv_r, v_bin);
bar(v_bin, vr_cnt/sum(vr_cnt), 0.5, 'facec', [0 0.8 0]);

xlim([-1 26]);
ylim([0 0.15]);
xlabel('Wind Speed (m/s)');
ylabel('Probability (-)');
set(gca, 'layer', 'top');
set(gca, 'ygrid', 'on');
legend('Reconstruction with 10000 Samples', 'Reconstruction with 50 Samples');


%% Wind power probability distribution
p_bin = linspace(0,1,21);

figure(2); clf;
set(gcf, 'units', 'inch', 'pos', [6.5    1.4583    5.8333    6.75]);

% ====================
% PDF
subplot(3,1,1);
p_cnt = hist(p_log(:), p_bin);
bar(p_bin, p_cnt/sum(p_cnt), 1, 'facec', [0.35 0.65 1]);
xlim([-0.06 1.06]);
set(gca, 'xtick', 0:0.25:1);
ylim([0 0.25]);
set(gca, 'ytick', 0:0.05:0.25);
grid on;
ylabel('Probability (-)');
set(gca, 'layer', 'top');
set(gca, 'ygrid', 'on');

% ====================
% CDF
subplot(3,1,2);
xp = sort(p_log(:));
Fp = (1:length(p_log(:)))/length(p_log(:));
plot(xp,Fp);
xlim([-0.06 1.06]);
ylim([-0.05 1.05]);
set(gca, 'xtick', 0:0.25:1);
set(gca, 'ytick', 0:0.25:1);
grid on;
set(gca, 'tickdir', 'out');
ylabel('CDF (-)');

% ====================
% Realization
subplot(3,1,3); hold on; box on;
Fp_r = linspace(0,1,10000);
xp_r = interp1(Fp,xp,Fp_r);
xp_r = xp_r(~isnan(xp_r));
pr_cnt = hist(xp_r, p_bin);
bar(p_bin, pr_cnt/sum(pr_cnt), 1, 'facec', [1 0.4 0.4]);

Fp_r = linspace(0,1,50);
xp_r = interp1(Fp,xp,Fp_r);
xp_r = xp_r(~isnan(xp_r));
pr_cnt = hist(xp_r, p_bin);
bar(p_bin, pr_cnt/sum(pr_cnt), 0.5, 'facec', [0 0.8 0]);

xlim([-0.06 1.06]);
ylim([0 0.25]);
set(gca, 'ytick', 0:0.05:0.25);
xlabel('Wind Power (normalized)');
ylabel('Probability (-)');
set(gca, 'layer', 'top');
set(gca, 'ygrid', 'on');
legend('Reconstruction with 10000 Samples', 'Reconstruction with 50 Samples');

save_name = [file_prefix, 'probability'];
save(save_name, 'xp', 'Fp', 'xv', 'Fv');
