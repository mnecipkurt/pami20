
% L-dimensional vector is observed at each time
N = 1e5; % size of the training data
k = 4; % compute sum of distances to k nearest neighbors

num = 2e3 + 24;  
indexes_rand = randi(N,num,1);  
indexes_1 = unique(indexes_rand);  
M = length(indexes_1);  % randomly partition the training data of sizes M and N-M

aa = 1:N;
aa(indexes_1) = 0;
aa = sort(aa,'ascend');
indexes_2 = aa(M+1:N);

% Training data
load IEEE57bus.mat;
L = m;
sigma = 0.1;
X = mvnrnd(H*phases,sigma^2*eye(L),N);  

% Partition the training data into two parts
X_1 = X(indexes_1,:);
X_2 = X(indexes_2,:);

% % Find a compact subset of X_1 (For each point in X_1, find the kNN in X_2 and compute the distance to the kNN)
% distances_abcd = zeros(M,1);
% for i = 1:M
%     i
%     datum_i = X_1(i,:);
%     tmp_dist = zeros(N-M,1);
%     for j = 1:(N-M)
%         datum_j = X_2(j,:);
%         dist_ij = norm(datum_i-datum_j,2);  % Euclidean distance
%         tmp_dist(j) = dist_ij;
%     end
%     sort_dist = sort(tmp_dist,'ascend');
%     sum_kNN = sum(sort_dist(1:k));
%     distances_abcd(i) = sum_kNN;
% end
% 
% [a,b] = sort(distances_abcd,'ascend');
% alpha = 0.1;
% num = floor((1-alpha)*M);
% C = X_1(b(1:num),:);  % most compact set in X_1 with the level of alpha
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

save('baseline');
