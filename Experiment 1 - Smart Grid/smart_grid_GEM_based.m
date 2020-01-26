close all; clear variables; clc;

load baseline;

% Online Attack Detection -- FAP calculations
alpha = 0.2;
h = (0.01:0.02:11.8)';
no_trials = 400; % number of trials
sum_fap = zeros(length(h),1);
for nn = 1:no_trials
    nn
    alarm_flag = zeros(length(h),1);
    g = 0; % initial decision stat
    t = 1;
    while (alarm_flag(length(h)) == 0)  
       datum_t = mvnrnd(H*phases,sigma^2*eye(L)); % nominal data point
       tmp_dist = zeros(num,1);
       for j = 1:num
           datum_j = C(j,:);
           tmp_dist(j) = norm(datum_t-datum_j,2);  % Euclidean distance
       end
       sort_dist = sort(tmp_dist,'ascend');
       sum_kNN = sum(sort_dist(1:k));
       tail_prob = sum(baseline_distances > sum_kNN)/(N-M);
       if (tail_prob == 0)   %  modification for small-size datasets
           tail_prob = 1/(N-M);
       end       
       g = g + log(alpha/tail_prob);
       if (g < 0)
           g = 0;
       end
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
    g = 0; % initial decision stat
    t = 0;
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4)       
       FD = -0.14 + 0.28*rand(1,L);    % false data: Uniform[-,-] for each entry
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
       tail_prob = sum(baseline_distances > sum_kNN)/(N-M);
       if (tail_prob == 0)   %  modification for small-size datasets
           tail_prob = 1/(N-M);
       end 
       g = g + log(alpha/tail_prob);
       if (g < 0)
           g = 0;
       end  
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

save('Smart_grid_GEM_based_results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online Attack Detection -- Sample Path
tmax = 250; % max time
tau = 200; % attack launch time
decision_stat = zeros(tmax,1);
g = 0; % initial decision stat
t = 1;
while (t <= tmax)
    if (t < tau)
       datum_t = mvnrnd(H*phases,sigma^2*eye(L));   % nominal data point
       tmp_dist = zeros(num,1);
       for j = 1:num
           datum_j = C(j,:);
           dist_tj = norm(datum_t-datum_j,2);  % Euclidean distance
           tmp_dist(j) = dist_tj;
       end
       sort_dist = sort(tmp_dist,'ascend');
       sum_kNN = sum(sort_dist(1:k));
       tail_prob = sum(baseline_distances > sum_kNN)/(N-M);
       if (tail_prob == 0)   %  modification
           tail_prob = 1/(N-M);
       end 
       g = g + log(alpha/tail_prob);
       if (g < 0)
           g = 0;
       end
    elseif (t >= tau)
       FD = -0.14 + 0.28*rand(1,L);   % false data: Uniform[-,-] for each entry
       datum_t = mvnrnd(H*phases,sigma^2*eye(L)) + FD;  % FDI Attack
       tmp_dist = zeros(num,1);
       for j = 1:num
           datum_j = C(j,:);
           dist_tj = norm(datum_t-datum_j,2);  % Euclidean distance
           tmp_dist(j) = dist_tj;
       end
       sort_dist = sort(tmp_dist,'ascend');
       sum_kNN = sum(sort_dist(1:k));
       tail_prob = sum(baseline_distances > sum_kNN)/(N-M);
       if (tail_prob == 0)   %  modification
           tail_prob = 1/(N-M);
       end 
       g = g + log(alpha/tail_prob);
       if (g < 0)
           g = 0;
       end 
    end    
    decision_stat(t) = g;
    t = t + 1;
end

plot(decision_stat,'linewidth',2);
hold on;
plot(12*ones(tmax,1),'r--');
plot(tau*ones(1.5*tmax,1),1:1.5*tmax,'k--');
xlabel('$t$','interpreter','latex','fontsize',14);
ylabel('$g_t$','interpreter','latex','fontsize',14);
ylim([0 350])


