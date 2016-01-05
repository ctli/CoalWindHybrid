clear
close all
clc

sample_cnt_range = [50 100 500 1000 5000 10000 20000];

median_log = zeros(1, length(sample_cnt_range));
Q1_log     = zeros(1, length(sample_cnt_range));
Q3_log     = zeros(1, length(sample_cnt_range));
top_log    = zeros(1, length(sample_cnt_range));
btm_log    = zeros(1, length(sample_cnt_range));
figure(1); clf; hold on; box on;
for sc = 1:length(sample_cnt_range)
    save_name = ['stochastic_Xilingol_', num2str(sc)];
    load(save_name);
    tmp = stochastic_cost_total(~isnan(stochastic_cost_total))/1e6;
    
    median_log(sc) = median(tmp);
    Q1_log(sc) = quantile(tmp, 0.25);
    Q3_log(sc) = quantile(tmp, 0.75);
    top_log(sc) = quantile(tmp, 0.975);
    btm_log(sc) = quantile(tmp, 0.025);
    
    dx1 = 0.12;
    dx2 = 0.15;
    plot([-dx1 dx1]+sc, [1 1]*top_log(sc), 'k');
    plot([-dx1 dx1]+sc, [1 1]*btm_log(sc), 'k');
    
    plot([1 1]*sc, [btm_log(sc), Q1_log(sc)], 'k-');
    plot([1 1]*sc, [Q3_log(sc), top_log(sc)], 'k-');
    
    plot([-dx2 dx2]+sc, [1 1]*Q1_log(sc), 'b');
    plot([-dx2 dx2]+sc, [1 1]*Q3_log(sc), 'b');
    plot([-dx2 -dx2]+sc, [Q1_log(sc) Q3_log(sc)], 'b');
    plot([+dx2 +dx2]+sc, [Q1_log(sc) Q3_log(sc)], 'b');
    
    plot([-dx2 dx2]+sc, [1 1]*median_log(sc), 'r');
end
% ylim([0 4]); title('Day 1033');
ylim([2.9 3.4]); title('Zoom In');
set(gca, 'xtick', 1:length(sample_cnt_range), 'xticklabel', sample_cnt_range);
set(gca, 'ygrid', 'on');
xlabel('Number of Samples');
ylabel('Distribution of Total Cost (bn USD)');


