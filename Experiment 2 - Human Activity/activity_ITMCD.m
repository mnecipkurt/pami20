clear variables; clc; close all;

load X_static.mat;
X = X_static';
N = size(X,1); % size of the training data
d = size(X,2); % dimensionality
p = 100;
f = 20;
k = 4;

load X_active.mat;
X_anomalous = X_active';
N_anomalous = size(X_anomalous,1); % size of the training data

% Online Detection -- FAP calculations
no_trials = 300; % number of trials
h = (0.1:0.2:79.7)';
sum_fap = zeros(length(h),1);
for mm = 1:no_trials
    mm
    % initial windows (free of anomaly, randomly chosen)
    indexes_rand1 = randi(N,p,1);
    past_window = X(indexes_rand1,:);
    indexes_rand2 = randi(N,f,1);
    future_window = X(indexes_rand2,:);
    %%%
    alarm_flag = zeros(length(h),1);
    t = 1;
    while (alarm_flag(length(h)) == 0)        
        past_window(1:p-1,:) = past_window(2:p,:);
        past_window(p,:) = future_window(1,:);
        ind_rand = randi(N,1);
        new_point = X(ind_rand,:);
        future_window(1:f-1,:) = future_window(2:f,:);
        future_window(f,:) = new_point;
        %%%
        dist2past = zeros(p,1);
        dist2future = zeros(p,1);
        for i = 1:p
            X_i = past_window(i,:);
            %%%
            tmp_dist1 = zeros(p,1);            
            for j = 1:p
                X_j = past_window(j,:);
                tmp_dist1(j) = norm(X_i-X_j,2);
            end
            tmp_dist1 = sort(tmp_dist1,'ascend');
            dist2past(i) = tmp_dist1(k+1);
            %%%
            tmp_dist2 = zeros(f,1);
            for j = 1:f
                X_j = future_window(j,:);
                tmp_dist2(j) = norm(X_i-X_j,2);
            end
            tmp_dist2 = sort(tmp_dist2,'ascend');
            dist2future(i) = tmp_dist2(k);            
        end
         
        dist3 = zeros(f,1);
        dist4 = zeros(f,1);
        for i = 1:f
            X_i = future_window(i,:);            
            %%%
            tmp_dist1 = zeros(f,1);
            for j = 1:f
                X_j = future_window(j,:);
                tmp_dist1(j) = norm(X_i-X_j,2);
            end
            tmp_dist1 = sort(tmp_dist1,'ascend');
            dist3(i) = tmp_dist1(k+1);
            %%%
            tmp_dist2 = zeros(p,1);
            for j = 1:p
                X_j = past_window(j,:);
                tmp_dist2(j) = norm(X_i-X_j,2);
            end
            tmp_dist2 = sort(tmp_dist2,'ascend');
            dist4(i) = tmp_dist2(k);            
        end
        
       g = (d/p)*(sum(log(dist2future./dist2past))) + log(f/(p-1)) ...
           + (d/f)*(sum(log(dist4./dist3))) + log(p/(f-1)); % decision stat
       sum_fap = sum_fap + t*(g >= h).*(alarm_flag == 0);
       alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
       % time update   
       t = t + 1;
    end
end
mean_fap = sum_fap/no_trials;
false_positive_rate = 1./mean_fap;  % false alarm rate: reciprocal of the mean false period

% Online Detection -- ADD calculations
tau = 1; % attack launch time
sum_add = zeros(length(h),1);
num_trials = 1e3;
max_tolerable_delay1 = 10;
cnt_detected1 = zeros(length(h),1);  % compute how many times detected within the delay bound, for each threshold level
max_tolerable_delay2 = 15;
cnt_detected2 = zeros(length(h),1);
max_tolerable_delay3 = 20;
cnt_detected3 = zeros(length(h),1);
max_tolerable_delay4 = 25;
cnt_detected4 = zeros(length(h),1);

for mm = 1:num_trials
    mm
    % initial windows (free of anomaly)
    indexes_rand1 = randi(N,p,1);
    past_window = X(indexes_rand1,:);
    indexes_rand2 = randi(N,f,1);
    future_window = X(indexes_rand2,:);
    %%%
    alarm_flag = zeros(length(h),1);
    t = 0;
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4)
        past_window(1:p-1,:) = past_window(2:p,:);
        past_window(p,:) = future_window(1,:);
        ind_rand = randi(N_anomalous,1);
        new_point = X_anomalous(ind_rand,:);
        future_window(1:f-1,:) = future_window(2:f,:);
        future_window(f,:) = new_point;
        
        dist2past = zeros(p,1);
        dist2future = zeros(p,1);
        for i = 1:p
            tmp_dist1 = zeros(p,1);
            X_i = past_window(i,:);
            for j = 1:p
                X_j = past_window(j,:);
                tmp_dist1(j) = norm(X_i-X_j,2);
            end
            tmp_dist1 = sort(tmp_dist1,'ascend');
            dist2past(i) = tmp_dist1(k+1);
            
            tmp_dist2 = zeros(f,1);
            for j = 1:f
                X_j = future_window(j,:);
                tmp_dist2(j) = norm(X_i-X_j,2);
            end
            tmp_dist2 = sort(tmp_dist2,'ascend');
            dist2future(i) = tmp_dist2(k);            
        end
         
        dist3 = zeros(f,1);
        dist4 = zeros(f,1);
        for i = 1:f
            tmp_dist1 = zeros(f,1);
            X_i = future_window(i,:);
            for j = 1:f
                X_j = future_window(j,:);
                tmp_dist1(j) = norm(X_i-X_j,2);
            end
            tmp_dist1 = sort(tmp_dist1,'ascend');
            dist3(i) = tmp_dist1(k+1);
            
            tmp_dist2 = zeros(p,1);
            for j = 1:p
                X_j = past_window(j,:);
                tmp_dist2(j) = norm(X_i-X_j,2);
            end
            tmp_dist2 = sort(tmp_dist2,'ascend');
            dist4(i) = tmp_dist2(k);            
        end
        
       g = (d/p)*(sum(log(dist2future./dist2past))) + log(f/(p-1)) ...
           + (d/f)*(sum(log(dist4./dist3))) + log(p/(f-1));
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
recall_10 = cnt_detected1/num_trials;  % true positive rate when the max. allowed detection delay is 10
recall_15 = cnt_detected2/num_trials;
recall_20 = cnt_detected3/num_trials;
recall_25 = cnt_detected4/num_trials;

save('activity_ITMCD_results');


