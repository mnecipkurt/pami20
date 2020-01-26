clear variables; clc; close all;

load PowerGrid_Data.mat;
N = size(X,1); % size of the training data
d = size(X,2); % dimensionality

load IEEE57bus.mat;
L = m;
sigma = 0.1;

K = 50; n0 = 10;
k = 10;

% Online Attack Detection -- FAP calculations
h = (-2:0.02:3.87)';
no_trials = 100; % number of trials
sum_fap = zeros(length(h),1);
for nn = 1:no_trials
    nn
    sliding_data_window = X(1:K,:);  % Initial Nominal Data Window
    alarm_flag = zeros(length(h),1);
    t = 1;
    while (alarm_flag(length(h)) == 0)  
        sliding_data_window(1:K-1,:) = sliding_data_window(2:K,:);
        sliding_data_window(K,:) = mvnrnd(H*phases,sigma^2*eye(L)); % nominal data point
        NN_mtrx = zeros(K,K);   % k nearest neighbor matrix: NNs take 1, others take 0 for each data point
        for i = 1:K
            for j = 1:K
                NN_mtrx(i,j) = norm((sliding_data_window(i,:)-sliding_data_window(j,:)),2);
            end
            NN_mtrx(i,i) = inf;
            tmpp = NN_mtrx;
            [B,I] = sort(tmpp(i,:),'ascend');
            NN_mtrx(i,I(1:k)) = 1;
            NN_mtrx(i,I(k+1:end)) = 0;
        end

        sum1 = 0;
        for i = 1:K
            for j = 1:K
                sum1 = sum1 + NN_mtrx(i,j) * NN_mtrx(j,i);
            end
        end

        sum2 = 0;
        for i = 1:K
            for j = 1:K
                for ell = 1:K
                    sum2 = sum2 + NN_mtrx(j,i) * NN_mtrx(ell,i);
                end
            end
        end

        g = -inf;
        for mm = n0:K-n0
        % data_range1 = 1:mm; data_range2 = mm+1:K;   
            R_mm = 0;
            for i = 1:K
                for j = 1:K
                    R_mm = R_mm + (NN_mtrx(i,j) + NN_mtrx(j,i))*((i<=mm && j>mm) || (i>mm && j<=mm));
                end
            end
            Exp_R_mm = 4*k*(K-mm)*mm/(K-1);

            hh = 4*(K-mm-1)*(mm-1)/((K-2)*(K-3));
            Var_R_mm = (4*(K-mm)*mm/(K-1))*(hh*(sum1/K + k - (2*(k^2))/(K-1)) + (1-hh)*(sum2/K - k^2));

            g_mm = (-R_mm + Exp_R_mm)/sqrt(Var_R_mm);

            if (g_mm > g)   % max operation over mm's
                g = g_mm;
            end
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
num_trials = 2e3;
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
    sliding_data_window = X(1:K,:);  % Initial Nominal Data Window
    alarm_flag = zeros(length(h),1);
    t = 0;
    while (alarm_flag(length(h)) == 0 || t <= max_tolerable_delay4)       
        FD = -0.14 + 0.28*rand(1,L);    % false data: Uniform[-,-] for each entry
        sliding_data_window(1:K-1,:) = sliding_data_window(2:K,:);
        sliding_data_window(K,:) = mvnrnd(H*phases,(sigma^2)*eye(L)) + FD;  % FDI Attack
        NN_mtrx = zeros(K,K);   % k nearest neighbor matrix: NNs take 1, others take 0 for each data point
        for i = 1:K
            for j = 1:K
                NN_mtrx(i,j) = norm((sliding_data_window(i,:)-sliding_data_window(j,:)),2);
            end
            NN_mtrx(i,i) = inf;
            tmpp = NN_mtrx;
            [B,I] = sort(tmpp(i,:),'ascend');
            NN_mtrx(i,I(1:k)) = 1;
            NN_mtrx(i,I(k+1:end)) = 0;
        end

        sum1 = 0;
        for i = 1:K
            for j = 1:K
                sum1 = sum1 + NN_mtrx(i,j) * NN_mtrx(j,i);
            end
        end

        sum2 = 0;
        for i = 1:K
            for j = 1:K
                for ell = 1:K
                    sum2 = sum2 + NN_mtrx(j,i) * NN_mtrx(ell,i);
                end
            end
        end

        g = -inf;
        for mm = n0:K-n0
        % data_range1 = 1:mm; data_range2 = mm+1:K;   
            R_mm = 0;
            for i = 1:K
                for j = 1:K
                    R_mm = R_mm + (NN_mtrx(i,j) + NN_mtrx(j,i))*((i<=mm && j>mm) || (i>mm && j<=mm));
                end
            end
            Exp_R_mm = 4*k*(K-mm)*mm/(K-1);

            hh = 4*(K-mm-1)*(mm-1)/((K-2)*(K-3));
            Var_R_mm = (4*(K-mm)*mm/(K-1))*(hh*(sum1/K + k - (2*(k^2))/(K-1)) + (1-hh)*(sum2/K - k^2));

            g_mm = (-R_mm + Exp_R_mm)/sqrt(Var_R_mm);

            if (g_mm > g)
                g = g_mm;
            end
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

save('Smart_grid_NN_heavy_results');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Online Attack Detection -- Sample Path
% tmax = 250; % max time
% tau = 200; % attack launch time
% decision_stat = zeros(tmax,1);
% t = 1;
% % Initial Nominal Data Window
% sliding_data_window = X(1:K,:);
% while (t <= tmax)
%     t
%     if (t < tau) 
%         sliding_data_window(1:K-1,:) = sliding_data_window(2:K,:);
%         sliding_data_window(K,:) = mvnrnd(H*phases,sigma^2*eye(L)); % nominal data point
%         NN_mtrx = zeros(K,K);   % k nearest neighbor matrix: NNs take 1, others take 0 for each data point
%         for i = 1:K
%             for j = 1:K
%                 NN_mtrx(i,j) = norm((sliding_data_window(i,:)-sliding_data_window(j,:)),2);
%             end
%             NN_mtrx(i,i) = inf;
%             tmpp = NN_mtrx;
%             [B,I] = sort(tmpp(i,:),'ascend');
%             NN_mtrx(i,I(1:k)) = 1;
%             NN_mtrx(i,I(k+1:end)) = 0;
%         end
% 
%         sum1 = 0;
%         for i = 1:K
%             for j = 1:K
%                 sum1 = sum1 + NN_mtrx(i,j) * NN_mtrx(j,i);
%             end
%         end
% 
%         sum2 = 0;
%         for i = 1:K
%             for j = 1:K
%                 for ell = 1:K
%                     sum2 = sum2 + NN_mtrx(j,i) * NN_mtrx(ell,i);
%                 end
%             end
%         end
% 
%         g = -inf;
%         for mm = n0:K-n0
%         % data_range1 = 1:mm; data_range2 = mm+1:K;   
%             R_mm = 0;
%             for i = 1:K
%                 for j = 1:K
%                     R_mm = R_mm + (NN_mtrx(i,j) + NN_mtrx(j,i))*((i<=mm && j>mm) || (i>mm && j<=mm));
%                 end
%             end
%             Exp_R_mm = 4*k*(K-mm)*mm/(K-1);
% 
%             hh = 4*(K-mm-1)*(mm-1)/((K-2)*(K-3));
%             Var_R_mm = (4*(K-mm)*mm/(K-1))*(hh*(sum1/K + k - (2*(k^2))/(K-1)) + (1-hh)*(sum2/K - k^2));
% 
%             g_mm = (-R_mm + Exp_R_mm)/sqrt(Var_R_mm);
% 
%             if (g_mm > g)
%                 g = g_mm;
%             end
%         end
%     elseif (t >= tau)       
%         sliding_data_window(1:K-1,:) = sliding_data_window(2:K,:);
%         FD = -0.3 + 0.2*rand(1,L);   % false data: Uniform[-,-] for each entry
%         sliding_data_window(K,:) = mvnrnd(H*phases,(sigma^2)*eye(L)) + FD;  % FDI Attack
%         NN_mtrx = zeros(K,K);   % k nearest neighbor matrix: NNs take 1, others take 0 for each data point
%         for i = 1:K
%             for j = 1:K
%                 NN_mtrx(i,j) = norm((sliding_data_window(i,:)-sliding_data_window(j,:)),2);
%             end
%             NN_mtrx(i,i) = inf;
%             tmpp = NN_mtrx;
%             [B,I] = sort(tmpp(i,:),'ascend');
%             NN_mtrx(i,I(1:k)) = 1;
%             NN_mtrx(i,I(k+1:end)) = 0;
%         end
% 
%         sum1 = 0;
%         for i = 1:K
%             for j = 1:K
%                 sum1 = sum1 + NN_mtrx(i,j) * NN_mtrx(j,i);
%             end
%         end
% 
%         sum2 = 0;
%         for i = 1:K
%             for j = 1:K
%                 for ell = 1:K
%                     sum2 = sum2 + NN_mtrx(j,i) * NN_mtrx(ell,i);
%                 end
%             end
%         end
% 
%         g = -inf;
%         for mm = n0:K-n0
%         % data_range1 = 1:mm; data_range2 = mm+1:K;   
%             R_mm = 0;
%             for i = 1:K
%                 for j = 1:K
%                     R_mm = R_mm + (NN_mtrx(i,j) + NN_mtrx(j,i))*((i<=mm && j>mm) || (i>mm && j<=mm));
%                 end
%             end
%             Exp_R_mm = 4*k*(K-mm)*mm/(K-1);
% 
%             hh = 4*(K-mm-1)*(mm-1)/((K-2)*(K-3));
%             Var_R_mm = (4*(K-mm)*mm/(K-1))*(hh*(sum1/K + k - (2*(k^2))/(K-1)) + (1-hh)*(sum2/K - k^2));
% 
%             g_mm = (-R_mm + Exp_R_mm)/sqrt(Var_R_mm);
% 
%             if (g_mm > g)
%                 g = g_mm;
%             end
%         end
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
