close all; clear variables;

% ADD vs. FAP

load activity_GEM_based_reduced_dim_results.mat
figure; plot(mean_fap([63 97 127 145 164 175 182 188]), mean_add([63 97 127 145 164 175 182 188]),'k-s','linewidth',1.8); hold on; grid on;

load activity_PCA_based_results.mat
plot(mean_fap([419 570 640 646 826 949 1019 1072]), mean_add([419 570 640 646 826 949 1019 1072]),'m-v','linewidth',1.8);

load activity_ITMCD_results.mat;
plot(mean_fap([2 73 172 260 308 331 348 362]), mean_add([2 73 172 260 308 331 348 362]),'b-d','linewidth',1.8);

load activity_nonpar_cusum_results.mat;
plot(mean_fap([5 13 19 27 33 38 42]), mean_add([5 13 19 27 33 38 42]),'r-+','linewidth',1.8);

load activity_ODIT_results.mat;
plot(mean_fap([92 141 162 177 196 208 216 224]), mean_add([92 141 162 177 196 208 216 224]),'g-o','linewidth',1.8);

load activity_QuantTree_results.mat;
plot(mean_fap([32 80 110 139 156 165 170 175]), mean_add([32 80 110 139 156 165 170 175]),'c-^','linewidth',1.8);

load activity_NN_heavy_results.mat;
plot(mean_fap([19 35 48 62 68 73 76 78]), mean_add([19 35 48 62 68 73 76 78]),'y-*','linewidth',1.8);

xlabel('$\rm{E}_\infty[\Gamma]$','interpreter','latex','fontsize',14);
ylabel('$\rm{E}_\tau \big[(\Gamma-\tau)^+\big]$','interpreter','latex','fontsize',14);
xlim([1 1e4]);
legend('Algorithm 1*','Algorithm 2','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','northwest');

axes('position',[0.42 0.34 0.46 0.26]); box on; 
load activity_GEM_based_reduced_dim_results.mat
plot(mean_fap([63 97 127 145 164 175 182 188]), mean_add([63 97 127 145 164 175 182 188]),'k-s','linewidth',1.8); hold on; grid on;
load activity_PCA_based_results.mat
plot(mean_fap([419 570 640 646 826 949 1019 1072]), mean_add([419 570 640 646 826 949 1019 1072]),'m-v','linewidth',1.8);
load activity_ODIT_results.mat;
plot(mean_fap([92 141 162 177 196 208 216 224]), mean_add([92 141 162 177 196 208 216 224]),'g-o','linewidth',1.8);
load activity_nonpar_cusum_results.mat;
plot(mean_fap([5 13 19 27 33 38 42]), mean_add([5 13 19 27 33 38 42]),'r-+','linewidth',1.8);
axis tight; grid on;

% Recall vs. FAR

load activity_GEM_based_reduced_dim_results.mat
figure; semilogx(1./mean_fap([1 30 63 97 127 164 188]), recall_10([1 30 63 97 127 164 188]),'k-s','linewidth',1.8); hold on; grid on;

load activity_PCA_based_results.mat
semilogx(1./mean_fap([1 200 400 600 800 940 1072]), recall_10([1 200 400 600 800 940 1072]),'m-v','linewidth',1.8);

load activity_ITMCD_results.mat;
semilogx(1./mean_fap([1 73 140 200 270 320 362]), recall_10([1 73 140 200 270 320 362]),'b-d','linewidth',1.8);

load activity_nonpar_cusum_results.mat;
semilogx(1./mean_fap([1 2 5 13 27 42]), recall_10([1 2 5 13 27 42]),'r-+','linewidth',1.8);

load activity_ODIT_results.mat;
semilogx(1./mean_fap([1 46 92 141 177 224]), recall_10([1 46 92 141 177 224]),'g-o','linewidth',1.8);

load activity_QuantTree_results.mat;
semilogx(1./mean_fap([1 40 55 70 87 110 139 162 175]), recall_10([1 40 55 70 87 110 139 162 175]),'c-^','linewidth',1.8);

load activity_NN_heavy_results.mat;
semilogx(1./mean_fap([1 30 48 62 70 78]), recall_10([1 30 48 62 70 78]),'y-*','linewidth',1.8);

xlabel('FAR');
ylabel('TPR');
xlim([1e-4 1]);
legend('Algorithm 1*','Algorithm 2','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','northwest');

% FAP, Approximation, Lower Bound
load activity_GEM_based_reduced_dim_results.mat
figure; semilogy(h,mean_fap,'k','linewidth',2); hold on;

load activity_PCA_based_results.mat
theta = lambertw(alpha*log(alpha))/log(alpha);
w_0 = theta - 1;
lower_bnd_fap = exp(-w_0*h);
apprx_fap = lower_bnd_fap*10.1;
Wald_apprx_fap = (h + (exp(-w_0*h)-1)/w_0)/(1+log(alpha));
semilogy(h,mean_fap,'m','linewidth',2); hold on;
semilogy(h,apprx_fap,'b-.','linewidth',2);
semilogy(h,Wald_apprx_fap,'c-.','linewidth',2);
semilogy(h,lower_bnd_fap,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period (Alg.1*)','False Alarm Period (Alg.2)','Our Approximation',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;
xlim([2 12])
