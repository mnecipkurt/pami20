clear variables; clc; close all;

load PowerGrid_Data.mat;
N = size(X,1); % size of the training data
d = size(X,2); % dimensionality
load IEEE57bus.mat;
L = m;
sigma = 0.1;

% Space Partitioning
K = 16;
W = 256;
indices = randi(d,K+1,1);   % dimensions to seperate the observation space
indices = unique(indices);  % make sure: length = K
%zheta = (randn(K,1)>=0);
zheta = zeros(K,1);
thresholds = zeros(K,1);    % space partitioning thresholds determined by the available nominal data
X_rem = X;
for i = 1:K
    N_i = N/K;
    ind = indices(i);
    tmp_vec = X_rem(:,ind);
    [B,I] = sort(tmp_vec,'ascend');
    if (zheta(i) == 0)
        thresholds(i) = B(N_i); % upper threshold for subregion i
        ind_range = I(N_i+1:end);
        X_rem = X_rem(ind_range,:);
    elseif (zheta(i) == 1)
        thresholds(i) = B(N - i*N_i + 1); % lower threshold for subregion i
        ind_range = I(1:N-i*N_i);
        X_rem = X_rem(ind_range,:);
    end
end

% Online Attack Detection -- FAP calculations
h = (0.2:0.1:38)';
no_trials = 1000; % number of trials
sum_fap = zeros(length(h),1);
for nn = 1:no_trials
    nn
    alarm_flag = zeros(length(h),1);
    t = 1;
    % Initial Window
    sliding_window_subregions = zeros(W,1); % shows which observations belong to which subregion among K subregions
    for i = 1:W
        nominal_point = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
        if (nominal_point(indices(1)) <= thresholds(1))
            sliding_window_subregions(i) = 1;
        elseif (nominal_point(indices(2)) <= thresholds(2))
            sliding_window_subregions(i) = 2;
        elseif (nominal_point(indices(3)) <= thresholds(3))
            sliding_window_subregions(i) = 3;
        elseif (nominal_point(indices(4)) <= thresholds(4))
            sliding_window_subregions(i) = 4;
        elseif (nominal_point(indices(5)) <= thresholds(5))
            sliding_window_subregions(i) = 5;
        elseif (nominal_point(indices(6)) <= thresholds(6))
            sliding_window_subregions(i) = 6;
        elseif (nominal_point(indices(7)) <= thresholds(7))
            sliding_window_subregions(i) = 7;
        elseif (nominal_point(indices(8)) <= thresholds(8))
            sliding_window_subregions(i) = 8;
        elseif (nominal_point(indices(9)) <= thresholds(9))
            sliding_window_subregions(i) = 9;
        elseif (nominal_point(indices(10)) <= thresholds(10))
            sliding_window_subregions(i) = 10;
        elseif (nominal_point(indices(11)) <= thresholds(11))
            sliding_window_subregions(i) = 11;
        elseif (nominal_point(indices(12)) <= thresholds(12))
            sliding_window_subregions(i) = 12;
        elseif (nominal_point(indices(13)) <= thresholds(13))
            sliding_window_subregions(i) = 13;
        elseif (nominal_point(indices(14)) <= thresholds(14))
            sliding_window_subregions(i) = 14;
        elseif (nominal_point(indices(15)) <= thresholds(15))
            sliding_window_subregions(i) = 15;
        else
            sliding_window_subregions(i) = 16;
        end
    end
    while (alarm_flag(length(h)) == 0)  
        sliding_window_subregions(1:W-1) = sliding_window_subregions(2:W);
        nominal_point = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
        if (nominal_point(indices(1)) <= thresholds(1))
            sliding_window_subregions(W) = 1;
        elseif (nominal_point(indices(2)) <= thresholds(2))
            sliding_window_subregions(W) = 2;
        elseif (nominal_point(indices(3)) <= thresholds(3))
            sliding_window_subregions(W) = 3;
        elseif (nominal_point(indices(4)) <= thresholds(4))
            sliding_window_subregions(W) = 4;
        elseif (nominal_point(indices(5)) <= thresholds(5))
            sliding_window_subregions(W) = 5;
        elseif (nominal_point(indices(6)) <= thresholds(6))
            sliding_window_subregions(W) = 6;
        elseif (nominal_point(indices(7)) <= thresholds(7))
            sliding_window_subregions(W) = 7;
        elseif (nominal_point(indices(8)) <= thresholds(8))
            sliding_window_subregions(W) = 8;
        elseif (nominal_point(indices(9)) <= thresholds(9))
            sliding_window_subregions(W) = 9;
        elseif (nominal_point(indices(10)) <= thresholds(10))
            sliding_window_subregions(W) = 10;
        elseif (nominal_point(indices(11)) <= thresholds(11))
            sliding_window_subregions(W) = 11;
        elseif (nominal_point(indices(12)) <= thresholds(12))
            sliding_window_subregions(W) = 12;
        elseif (nominal_point(indices(13)) <= thresholds(13))
            sliding_window_subregions(W) = 13;
        elseif (nominal_point(indices(14)) <= thresholds(14))
            sliding_window_subregions(W) = 14;
        elseif (nominal_point(indices(15)) <= thresholds(15))
            sliding_window_subregions(W) = 15;
        else
            sliding_window_subregions(W) = 16;
        end
        nums = zeros(K,1);
        for i = 1:K
            nums(i) = sum(sliding_window_subregions == i);
        end
        prob = 1/K;
        g = sum((nums-W*prob).^2)/(W*prob); % test statistic
        sum_fap = sum_fap + t*(g >= h).*(alarm_flag == 0);
        alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
        % time update   
        t = t + 1;
    end
end
mean_fap = sum_fap/no_trials;
false_positive_rate = 1./mean_fap;  % false alarm rate = reciprocal of the mean false period

% Online Attack Detection -- ADD calculations
%tau = 1: attack launch time (corresponds to the worst case detection delay)
sum_add = zeros(length(h),1);
num_trials = 1e4;
max_tolerable_delay1 = 10;
cnt_detected1 = zeros(length(h),1);  % compute how many times detected within the delay bound, for each threshold level
max_tolerable_delay2 = 15;
cnt_detected2 = zeros(length(h),1);
max_tolerable_delay3 = 20;
cnt_detected3 = zeros(length(h),1);
max_tolerable_delay4 = 25;
cnt_detected4 = zeros(length(h),1);

for nn = 1:num_trials
    nn
    alarm_flag = zeros(length(h),1);
    t = 0;
    % Initial Window
    sliding_window_subregions = zeros(W,1); % shows which observations belong to which subregion among K subregions
    for i = 1:W
        nominal_point = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
        if (nominal_point(indices(1)) <= thresholds(1))
            sliding_window_subregions(i) = 1;
        elseif (nominal_point(indices(2)) <= thresholds(2))
            sliding_window_subregions(i) = 2;
        elseif (nominal_point(indices(3)) <= thresholds(3))
            sliding_window_subregions(i) = 3;
        elseif (nominal_point(indices(4)) <= thresholds(4))
            sliding_window_subregions(i) = 4;
        elseif (nominal_point(indices(5)) <= thresholds(5))
            sliding_window_subregions(i) = 5;
        elseif (nominal_point(indices(6)) <= thresholds(6))
            sliding_window_subregions(i) = 6;
        elseif (nominal_point(indices(7)) <= thresholds(7))
            sliding_window_subregions(i) = 7;
        elseif (nominal_point(indices(8)) <= thresholds(8))
            sliding_window_subregions(i) = 8;
        elseif (nominal_point(indices(9)) <= thresholds(9))
            sliding_window_subregions(i) = 9;
        elseif (nominal_point(indices(10)) <= thresholds(10))
            sliding_window_subregions(i) = 10;
        elseif (nominal_point(indices(11)) <= thresholds(11))
            sliding_window_subregions(i) = 11;
        elseif (nominal_point(indices(12)) <= thresholds(12))
            sliding_window_subregions(i) = 12;
        elseif (nominal_point(indices(13)) <= thresholds(13))
            sliding_window_subregions(i) = 13;
        elseif (nominal_point(indices(14)) <= thresholds(14))
            sliding_window_subregions(i) = 14;
        elseif (nominal_point(indices(15)) <= thresholds(15))
            sliding_window_subregions(i) = 15;
        else
            sliding_window_subregions(i) = 16;
        end
    end
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4)       
        sliding_window_subregions(1:W-1) = sliding_window_subregions(2:W);
        FD = -0.14 + 0.28*rand(1,L);   % false data: Uniform[-,-] for each entry
        anomalous_point = mvnrnd(H*phases,sigma^2*eye(L)) + FD;  % FDI Attack
        if (anomalous_point(indices(1)) <= thresholds(1))
            sliding_window_subregions(W) = 1;
        elseif (anomalous_point(indices(2)) <= thresholds(2))
            sliding_window_subregions(W) = 2;
        elseif (anomalous_point(indices(3)) <= thresholds(3))
            sliding_window_subregions(W) = 3;
        elseif (anomalous_point(indices(4)) <= thresholds(4))
            sliding_window_subregions(W) = 4;
        elseif (anomalous_point(indices(5)) <= thresholds(5))
            sliding_window_subregions(W) = 5;
        elseif (anomalous_point(indices(6)) <= thresholds(6))
            sliding_window_subregions(W) = 6;
        elseif (anomalous_point(indices(7)) <= thresholds(7))
            sliding_window_subregions(W) = 7;
        elseif (anomalous_point(indices(8)) <= thresholds(8))
            sliding_window_subregions(W) = 8;
        elseif (anomalous_point(indices(9)) <= thresholds(9))
            sliding_window_subregions(W) = 9;
        elseif (anomalous_point(indices(10)) <= thresholds(10))
            sliding_window_subregions(W) = 10;
        elseif (anomalous_point(indices(11)) <= thresholds(11))
            sliding_window_subregions(W) = 11;
        elseif (anomalous_point(indices(12)) <= thresholds(12))
            sliding_window_subregions(W) = 12;
        elseif (anomalous_point(indices(13)) <= thresholds(13))
            sliding_window_subregions(W) = 13;
        elseif (anomalous_point(indices(14)) <= thresholds(14))
            sliding_window_subregions(W) = 14;
        elseif (anomalous_point(indices(15)) <= thresholds(15))
            sliding_window_subregions(W) = 15;
        else
            sliding_window_subregions(W) = 16;
        end
        nums = zeros(K,1);
        for i = 1:K
            nums(i) = sum(sliding_window_subregions == i);
        end
        prob = 1/K;
        g = sum((nums-W*prob).^2)/(W*prob); % test statistic
        sum_add = sum_add + t*(g >= h).*(alarm_flag == 0);
        alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
       
        if (t == max_tolerable_delay1)
           cnt_detected1 = cnt_detected1 + alarm_flag; % increment detection cases if detected within the given delay bound
        end
       
        if (t == max_tolerable_delay2)
           cnt_detected2 = cnt_detected2 + alarm_flag; 
        end
       
        if (t == max_tolerable_delay3)
           cnt_detected3 = cnt_detected3 + alarm_flag; 
        end
       
        if (t == max_tolerable_delay4)
           cnt_detected4 = cnt_detected4 + alarm_flag; 
        end
       
        % time update
        t = t + 1;
    end
end
mean_add = sum_add/num_trials;
recall_10 = cnt_detected1/num_trials;   % true positive rate when the max. allowed detection delay is 10
recall_15 = cnt_detected2/num_trials;   % true positive rate when the max. allowed detection delay is 15
recall_20 = cnt_detected3/num_trials;
recall_25 = cnt_detected4/num_trials;

save('Smart_grid_QuantTree_results');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Online Attack Detection -- Sample Path
% 
% % Initial Window
% sliding_window_subregions = zeros(W,1);
% for i = 1:W
%     nominal_point = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
%     if (nominal_point(indices(1)) <= thresholds(1))
%         sliding_window_subregions(i) = 1;
%     elseif (nominal_point(indices(2)) <= thresholds(2))
%         sliding_window_subregions(i) = 2;
%     elseif (nominal_point(indices(3)) <= thresholds(3))
%         sliding_window_subregions(i) = 3;
%     elseif (nominal_point(indices(4)) <= thresholds(4))
%         sliding_window_subregions(i) = 4;
%     elseif (nominal_point(indices(5)) <= thresholds(5))
%         sliding_window_subregions(i) = 5;
%     elseif (nominal_point(indices(6)) <= thresholds(6))
%         sliding_window_subregions(i) = 6;
%     elseif (nominal_point(indices(7)) <= thresholds(7))
%         sliding_window_subregions(i) = 7;
%     elseif (nominal_point(indices(8)) <= thresholds(8))
%         sliding_window_subregions(i) = 8;
%     elseif (nominal_point(indices(9)) <= thresholds(9))
%         sliding_window_subregions(i) = 9;
%     elseif (nominal_point(indices(10)) <= thresholds(10))
%         sliding_window_subregions(i) = 10;
%     elseif (nominal_point(indices(11)) <= thresholds(11))
%         sliding_window_subregions(i) = 11;
%     elseif (nominal_point(indices(12)) <= thresholds(12))
%         sliding_window_subregions(i) = 12;
%     elseif (nominal_point(indices(13)) <= thresholds(13))
%         sliding_window_subregions(i) = 13;
%     elseif (nominal_point(indices(14)) <= thresholds(14))
%         sliding_window_subregions(i) = 14;
%     elseif (nominal_point(indices(15)) <= thresholds(15))
%         sliding_window_subregions(i) = 15;
%     else
%         sliding_window_subregions(i) = 16;
%     end
% end
% 
% tmax = 500; % max time
% tau = 200; % attack launch time
% decision_stat = zeros(tmax,1);
% t = 1;
% g = 0;
% while (t <= tmax)
%     if (t < tau)
%         sliding_window_subregions(1:W-1) = sliding_window_subregions(2:W);
%         nominal_point = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
%         if (nominal_point(indices(1)) <= thresholds(1))
%         sliding_window_subregions(W) = 1;
%         elseif (nominal_point(indices(2)) <= thresholds(2))
%         sliding_window_subregions(W) = 2;
%         elseif (nominal_point(indices(3)) <= thresholds(3))
%         sliding_window_subregions(W) = 3;
%         elseif (nominal_point(indices(4)) <= thresholds(4))
%         sliding_window_subregions(W) = 4;
%         elseif (nominal_point(indices(5)) <= thresholds(5))
%         sliding_window_subregions(W) = 5;
%         elseif (nominal_point(indices(6)) <= thresholds(6))
%         sliding_window_subregions(W) = 6;
%         elseif (nominal_point(indices(7)) <= thresholds(7))
%         sliding_window_subregions(W) = 7;
%         elseif (nominal_point(indices(8)) <= thresholds(8))
%         sliding_window_subregions(W) = 8;
%         elseif (nominal_point(indices(9)) <= thresholds(9))
%         sliding_window_subregions(W) = 9;
%         elseif (nominal_point(indices(10)) <= thresholds(10))
%         sliding_window_subregions(W) = 10;
%         elseif (nominal_point(indices(11)) <= thresholds(11))
%         sliding_window_subregions(W) = 11;
%         elseif (nominal_point(indices(12)) <= thresholds(12))
%         sliding_window_subregions(W) = 12;
%         elseif (nominal_point(indices(13)) <= thresholds(13))
%         sliding_window_subregions(W) = 13;
%         elseif (nominal_point(indices(14)) <= thresholds(14))
%         sliding_window_subregions(W) = 14;
%         elseif (nominal_point(indices(15)) <= thresholds(15))
%         sliding_window_subregions(W) = 15;
%         else
%         sliding_window_subregions(W) = 16;
%         end
%         nums = zeros(K,1);
%         for i = 1:K
%         nums(i) = sum(sliding_window_subregions == i);
%         end
%         prob = 1/K;
%         g = sum((nums-W*prob).^2)/(W*prob); % test statistic
%     elseif (t >= tau)
%         sliding_window_subregions(1:W-1) = sliding_window_subregions(2:W);
%         FD = -0.12 + 0.24*rand(1,L);   % false data: Uniform[-,-] for each entry
%         anomalous_point = mvnrnd(H*phases,sigma^2*eye(L)) + FD;  % FDI Attack
%         if (anomalous_point(indices(1)) <= thresholds(1))
%         sliding_window_subregions(W) = 1;
%         elseif (anomalous_point(indices(2)) <= thresholds(2))
%         sliding_window_subregions(W) = 2;
%         elseif (anomalous_point(indices(3)) <= thresholds(3))
%         sliding_window_subregions(W) = 3;
%         elseif (anomalous_point(indices(4)) <= thresholds(4))
%         sliding_window_subregions(W) = 4;
%         elseif (anomalous_point(indices(5)) <= thresholds(5))
%         sliding_window_subregions(W) = 5;
%         elseif (anomalous_point(indices(6)) <= thresholds(6))
%         sliding_window_subregions(W) = 6;
%         elseif (anomalous_point(indices(7)) <= thresholds(7))
%         sliding_window_subregions(W) = 7;
%         elseif (anomalous_point(indices(8)) <= thresholds(8))
%         sliding_window_subregions(W) = 8;
%         elseif (anomalous_point(indices(9)) <= thresholds(9))
%         sliding_window_subregions(W) = 9;
%         elseif (anomalous_point(indices(10)) <= thresholds(10))
%         sliding_window_subregions(W) = 10;
%         elseif (anomalous_point(indices(11)) <= thresholds(11))
%         sliding_window_subregions(W) = 11;
%         elseif (anomalous_point(indices(12)) <= thresholds(12))
%         sliding_window_subregions(W) = 12;
%         elseif (anomalous_point(indices(13)) <= thresholds(13))
%         sliding_window_subregions(W) = 13;
%         elseif (anomalous_point(indices(14)) <= thresholds(14))
%         sliding_window_subregions(W) = 14;
%         elseif (anomalous_point(indices(15)) <= thresholds(15))
%         sliding_window_subregions(W) = 15;
%         else
%         sliding_window_subregions(W) = 16;
%         end
%         nums = zeros(K,1);
%         for i = 1:K
%         nums(i) = sum(sliding_window_subregions == i);
%         end
%         prob = 1/K;
%         g = sum((nums-W*prob).^2)/(W*prob); % test statistic
%     end    
%     decision_stat(t) = g;
%     t = t + 1;
% end
% 
% plot(decision_stat, 'linewidth',2);
% hold on;
% xlabel('$t$','interpreter','latex','fontsize',14);
% ylabel('$g_t$','interpreter','latex','fontsize',14);
% 
% 
% 
