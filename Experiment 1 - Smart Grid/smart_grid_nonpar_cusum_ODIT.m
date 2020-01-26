close all; clear variables; clc;
 
load baseline;
nominal_mean = mean(baseline_distances);

alpha = 0.2;
thresh_index = round((length(baseline_distances))*(1-alpha));
sorted_arr = sort(baseline_distances);
odit_thresh = sorted_arr(thresh_index);

% Online Attack Detection -- FAP calculations
h = (0.01:0.05:15.96)';
h_odit = (0.001:0.01:1.021)';
no_trials = 100; % number of trials
sum_fap = zeros(length(h),1);
sum_fap_odit = zeros(length(h_odit),1);

for nn = 1:no_trials
    nn
    alarm_flag = zeros(length(h),1);
    alarm_flag_odit = zeros(length(h_odit),1);
    g = 0; % initial decision stat
    g_odit = 0;
    t = 1;
    while (alarm_flag(length(h)) == 0 || alarm_flag_odit(length(h_odit)) == 0)  
       datum_t = mvnrnd(H*phases,sigma^2*eye(L)); % nominal data point
       tmp_dist = zeros(num,1);
       for j = 1:num
           datum_j = C(j,:);
           tmp_dist(j) = norm(datum_t-datum_j,2);  % Euclidean distance
       end
       sort_dist = sort(tmp_dist,'ascend');
       sum_kNN = sum(sort_dist(1:k));
       
       g = max(0, g + sum_kNN - nominal_mean);
       sum_fap = sum_fap + t*(g >= h).*(alarm_flag == 0);
       alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
       
       g_odit = max(0, g_odit + sum_kNN - odit_thresh);
       sum_fap_odit = sum_fap_odit + t*(g_odit >= h_odit).*(alarm_flag_odit == 0);
       alarm_flag_odit = alarm_flag_odit + (g_odit >= h_odit).*(alarm_flag_odit == 0);
       
       % time update   
       t = t + 1;
    end
end
mean_fap = sum_fap/no_trials;
false_positive_rate = 1./mean_fap;  % false alarm rate = reciprocal of the mean false period

mean_fap_odit = sum_fap_odit/no_trials;
false_positive_rate_odit = 1./mean_fap_odit;

% Online Attack Detection -- ADD calculations
% tau = 1: % attack launch time (corresponds to the worst case scenerio)
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

sum_add_odit = zeros(length(h_odit),1);
cnt_detected1_odit = zeros(length(h_odit),1);  % compute how many times detected within the delay bound, for each threshold level
cnt_detected2_odit = zeros(length(h_odit),1);
cnt_detected3_odit = zeros(length(h_odit),1);
cnt_detected4_odit = zeros(length(h_odit),1);

for nn = 1:num_trials
    nn
    alarm_flag = zeros(length(h),1);
    alarm_flag_odit = zeros(length(h_odit),1);
    g = 0; % initial decision stat
    g_odit = 0;
    t = 0;
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4 || alarm_flag_odit(length(h_odit)) == 0)       
       FD = -0.14 + 0.28*rand(1,L);   % false data: Uniform[-,-] for each entry
       %JD = mvnrnd(zeros(1,L),0.5*sigma^2*eye(L));  % AWGN jamming data
       datum_t = mvnrnd(H*phases,sigma^2*eye(L)) + FD;  % FDI Attack
       %datum_t = mvnrnd(H*phases,sigma^2*eye(L)) + JD;  % Jamming Attack
       tmp_dist = zeros(num,1);
       for j = 1:num
           datum_j = C(j,:);
           tmp_dist(j) = norm(datum_t-datum_j,2);  % Euclidean distance
       end
       sort_dist = sort(tmp_dist,'ascend');
       sum_kNN = sum(sort_dist(1:k));
       
       g = max(0, g + sum_kNN - nominal_mean);
       sum_add = sum_add + t*(g >= h).*(alarm_flag == 0);
       alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
       
       g_odit = max(0, g_odit + sum_kNN - odit_thresh);
       sum_add_odit = sum_add_odit + t*(g_odit >= h_odit).*(alarm_flag_odit == 0);
       alarm_flag_odit = alarm_flag_odit + (g_odit >= h_odit).*(alarm_flag_odit == 0);
       
       if (t == max_tolerable_delay1)
           cnt_detected1 = cnt_detected1 + alarm_flag; % increment detection cases if detected within the given delay bound
           cnt_detected1_odit = cnt_detected1_odit + alarm_flag_odit;
       end
       
       if (t == max_tolerable_delay2)
           cnt_detected2 = cnt_detected2 + alarm_flag; 
           cnt_detected2_odit = cnt_detected2_odit + alarm_flag_odit; 
       end
       
       if (t == max_tolerable_delay3)
           cnt_detected3 = cnt_detected3 + alarm_flag; 
           cnt_detected3_odit = cnt_detected3_odit + alarm_flag_odit; 
       end
       
       if (t == max_tolerable_delay4)
           cnt_detected4 = cnt_detected4 + alarm_flag; 
           cnt_detected4_odit = cnt_detected4_odit + alarm_flag_odit; 
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

mean_add_odit = sum_add_odit/num_trials;
recall_10_odit = cnt_detected1_odit/num_trials;  
recall_15_odit = cnt_detected2_odit/num_trials;  
recall_20_odit = cnt_detected3_odit/num_trials;
recall_25_odit = cnt_detected4_odit/num_trials;

save('Smart_grid_nonpar_cusum_ODIT_results');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Online Attack Detection -- Sample Path
% tmax = 250; % max time
% tau = 200; % attack launch time
% decision_stat = zeros(tmax,1);
% t = 1;
% g = 0;
% while (t <= tmax)
%     if (t < tau)
%        datum_t = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
%        tmp_dist = zeros(num,1);
%        for j = 1:num
%            datum_j = C(j,:);
%            dist_tj = norm(datum_t-datum_j,2);  % Euclidean distance
%            tmp_dist(j) = dist_tj;
%        end
%        sort_dist = sort(tmp_dist,'ascend');
%        sum_kNN = sum(sort_dist(1:k));
%        %g = max(0, g + (sum_kNN - nominal_mean));
%        g = max(0, g + (sum_kNN - odit_thresh));
%     elseif (t >= tau)
%        FD = -0.1 + 0.2*rand(1,L);   % false data: Uniform[-,-] for each entry
%        datum_t = mvnrnd(H*phases,sigma^2*eye(L)) + FD;  % FDI Attack
%        tmp_dist = zeros(num,1);
%        for j = 1:num
%            datum_j = C(j,:);
%            dist_tj = norm(datum_t-datum_j,2);  % Euclidean distance
%            tmp_dist(j) = dist_tj;
%        end
%        sort_dist = sort(tmp_dist,'ascend');
%        sum_kNN = sum(sort_dist(1:k));
%        %g = max(0, g + (sum_kNN - nominal_mean));
%        g = max(0, g + (sum_kNN - odit_thresh));
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
