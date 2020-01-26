clear variables; clc; close all;

load X_static.mat;
X_static = X_static';
load V_activity.mat;
data_reduced_dim = X_static*P;

% L-dimensional vector is observed at each time
N = size(data_reduced_dim,1); % size of the training data
k = 4; % compute sum of distances to k nearest neighbors

num = 880;  
indexes_rand = randi(N,num,1);  
indexes_1 = unique(indexes_rand);  
M = length(indexes_1);  % randomly partition the training data of sizes M and N-M
L = size(data_reduced_dim,2);

aa = 1:N;
aa(indexes_1) = 0;
aa = sort(aa,'ascend');
indexes_2 = aa(M+1:N);

% Partition the training data into two parts
X = data_reduced_dim;
X_1 = X(indexes_1,:);
X_2 = X(indexes_2,:);

C = X_1;
num = M;

% Compute the sum of the kNN distances from data points in X_2 to data points in compact subset C of X_1 (Offline)
baseline_distances = zeros(N-M,1);
for i = 1:(N-M)
    i
    datum_i = X_2(i,:);
    tmp_dist = zeros(num,1);
    for j = 1:num
        datum_j = C(j,:);
        dist_ij = norm(datum_i-datum_j,2);  % Euclidean distance
        tmp_dist(j) = dist_ij;
    end
    sort_dist = sort(tmp_dist,'ascend');
    sum_kNN = sum(sort_dist(1:k));
    baseline_distances(i) = sum_kNN;
end

save('activity_baseline');
