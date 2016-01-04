clear
close all
clc


wind_file = 'Xilingol_2009';
load(wind_file);

% yr_range = 1979:2009;
% file_prefix = 'Xilingol_';
% p_log = nan*ones(8760, length(yr_range));
% v_log = nan*ones(8760, length(yr_range));
% for iy = 1:length(yr_range)
%     yr = yr_range(iy);
%     file_name = [file_prefix, num2str(yr)];
%     load(file_name);
%     if length(p)>=8760;
%         p_log(:,iy) = p(1:8760);
%     end
%     if length(v)>=8760;
%         v_log(:,iy) = v(1:8760);
%     end
% end
% p = v_log(:);

n = 24;
x = zeros(length(p)-n+1, n);
for nn = 1:n
    st = nn;
    ed = n-st;
    xn = p(st:end-ed);
    x(:,nn) = xn;
end

% ====================
N = length(x);
sig = zeros(n,n);
for r = 1:n % Row
    for c = 1:n % Column
        xi = x(:,r);
        xj = x(:,c);
        mi = mean(xi);
        mj = mean(xj);
        cov_ij = 1/N*sum((xi-mi).*(xj-mj));
        sig(r,c) = cov_ij;
    end
end
sig_log = sig;

figure(100); clf;
surf(sig_log);
