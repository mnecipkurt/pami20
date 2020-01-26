clear variables; clc; close all;

load nominal_thermostat.mat;
X_nominal = nominal_thermostat';
N = size(X_nominal,2); % size of the training data
num = 9e3;  
indexes_rand = randi(N,num,1);  
indexes_1 = unique(indexes_rand);  
Xstat_PCA = X_nominal(:,indexes_1);
% Xstat_PCA = X_static;

num = 9e3;
indexes_rand = randi(N,num,1);  
indexes_1 = unique(indexes_rand);  
M = length(indexes_1);
aa = 1:N;
aa(indexes_1) = 0;
aa = sort(aa,'ascend');
indexes_2 = aa(M+1:N);
Xstat_edf = X_nominal(:,indexes_1);
% Xstat_test = X_static(:,indexes_2);
Xstat_test = X_nominal;

L = size(X_nominal,1); % L-dimensional vector is observed at each time

% Determine the principal subspace using Xstat_PCA
X_1 = Xstat_PCA';
M = size(X_1,1);
mean_vec = mean(X_1);
new_X = X_1 - ones(M,1)*mean_vec;
S = (new_X'*new_X)/M;
[V,D] = eig(S);
[d,ind] = sort(diag(D),'descend');
Vs = V(:,ind);
dd = svd(S);
plot(dd,'linewidth',2); xlabel('Dimension'); ylabel('Eigenvalue');
gamma_hat = sum(dd(1:5))/sum(dd);
P = Vs(:,1:5);
Proj = P*P';
Res_Proj = eye(L) - Proj;

% Compute baseline distances (residual magnitudes) using X_2 and the principal subspace
% X_2 = Xstat_edf';
X_2 = X_nominal';
td = length(X_2);
new_data2 = X_2 - ones(td,1)*mean_vec;
tmp_mtrx = Res_Proj*new_data2'; % we need to compute the norm of each column
baseline_distances = sqrt(sum(tmp_mtrx.*tmp_mtrx)); 
figure; hist(baseline_distances,100);
nominal_mean = mean(baseline_distances);

% Online Detection -- FAP calculations
h = ([1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5 501000 1e6 1501000 2001000 2501000 3001000 3501000 4001000 4501000 5001000 6001000])';
no_trials = 1000; % number of trials
sum_fap = zeros(length(h),1);
ss = size(Xstat_test,2);
for j = 1:no_trials
    j
    alarm_flag = zeros(length(h),1);
    g = 0; % initial decision stat
    t = 1;
    while (alarm_flag(length(h)) == 0)       
       index_rand = randi(ss,1);
       datum_t = (Xstat_test(:,index_rand))'; % nominal data point
       new_datum = datum_t - mean_vec;
       res_datum = Res_Proj*new_datum';
       dist = norm(res_datum,2);
       
       g = max(0, g + dist - nominal_mean);
       sum_fap = sum_fap + t*(g >= h).*(alarm_flag == 0);
       alarm_flag = alarm_flag + (g >= h).*(alarm_flag == 0);
       
       % time update   
       t = t + 1;
    end
end
mean_fap = sum_fap/no_trials;
false_positive_rate = 1./mean_fap;  % false alarm rate = reciprocal of the mean false period

% Attack Data
load junk_thermostat.mat;
X_3 = junk_thermostat;
rr = length(X_3);

% Online Detection -- ADD calculations
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

for j = 1:num_trials
    j
    alarm_flag = zeros(length(h),1);
    g = 0; % initial decision stat
    t = 0;
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4)       
       index_rand = randi(rr,1);
       datum_t = X_3(index_rand,:); 
       new_datum = datum_t - mean_vec;
       res_datum = Res_Proj*new_datum';
       dist = norm(res_datum,2);
       
       g = max(0, g + dist - nominal_mean);
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

save('botnet_nonpar_cusum_results');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Online Attack Detection -- Sample Path
% tmax = 250; % max time
% tau = 200; % attack launch time
% decision_stat = zeros(tmax,1);
% t = 1;
% g = 0;
% while (t <= tmax)
%     if (t < tau)
%        index_rand = randi(ss,1);
%        datum_t = (Xstat_test(:,index_rand))'; % nominal data point
%        new_datum = datum_t - mean_vec;
%        res_datum = Res_Proj*new_datum';
%        dist = norm(res_datum,2);
%        g = max(0, g + (dist - nominal_mean));
%     elseif (t >= tau)
%        index_rand = randi(rr,1);
%        datum_t = X_3(index_rand,:); 
%        new_datum = datum_t - mean_vec;
%        res_datum = Res_Proj*new_datum';
%        dist = norm(res_datum,2);
%        g = max(0, g + (dist - nominal_mean));
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
