close all; clear variables;

% ADD vs. FAP

load botnet_GEM_based_reduced_dim_results.mat
figure; plot(mean_fap([21 41 100 144 162 178 199 209]), mean_add([21 41 100 144 162 178 199 209]),'k-s','linewidth',1.8); hold on; grid on;

load botnet_PCA_based_results.mat
plot(mean_fap([21 41 100 144 162 184 199 212]), mean_add([21 41 100 144 162 184 199 212]),'m-v','linewidth',1.8);

load botnet_ITMCD_results.mat;
plot(mean_fap([20 60 100 140 180 200 212 223]), mean_add([20 60 100 140 180 200 212 223]),'b-d','linewidth',1.8);

load botnet_nonpar_cusum_results_v3.mat;
plot(mean_fap([8 13 17 20]), mean_add([8 13 17 20]),'r-+','linewidth',1.8);

load botnet_ODIT_results.mat;
plot(mean_fap(28:35), mean_add(28:35),'g-o','linewidth',1.8);

load botnet_QuantTree_results.mat;
plot(mean_fap([72 103 120 140 152 158 162 166]), mean_add([72 103 120 140 152 158 162 166]),'c-^','linewidth',1.8);

load botnet_NN_heavy_results.mat;
plot(mean_fap([74 166 220 280 332 360 381 396]), mean_add([74 166 220 280 332 360 381 396]),'y-*','linewidth',1.8);

xlabel('$\rm{E}_\infty[\Gamma]$','interpreter','latex','fontsize',14);
ylabel('$\rm{E}_\tau \big[(\Gamma-\tau)^+\big]$','interpreter','latex','fontsize',14);
legend('Algorithm 1*','Algorithm 2','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','northeast');
xlim([1 1e4]);

axes('position',[0.4 0.3 0.46 0.26]); box on; 
load botnet_GEM_based_reduced_dim_results.mat
plot(mean_fap([21 41 100 144 162 178 199 209]), mean_add([21 41 100 144 162 178 199 209]),'k-s','linewidth',1.8); hold on; grid on;
load botnet_PCA_based_results.mat
plot(mean_fap([21 41 100 144 162 184 199 212]), mean_add([21 41 100 144 162 184 199 212]),'m-v','linewidth',1.8);
load botnet_ITMCD_results.mat;
plot(mean_fap([20 60 100 140 180 200 212 223]), mean_add([20 60 100 140 180 200 212 223]),'b-d','linewidth',1.8);
load botnet_NN_heavy_results.mat;
plot(mean_fap([74 166 220 280 332 360 381 396]), mean_add([74 166 220 280 332 360 381 396]),'y-*','linewidth',1.8);
axis tight; grid on;



% Recall vs. FAR

load botnet_GEM_based_reduced_dim_results.mat
figure; semilogx(1./mean_fap([1 29:20:209]), recall_10([1 29:20:209]),'k-s','linewidth',1.8); hold on; grid on;

load botnet_PCA_based_results.mat
semilogx(1./mean_fap([1 15 32:20:212]), recall_10([1 15 32:20:212]),'m-v','linewidth',1.8);

load botnet_ITMCD_results.mat;
semilogx(1./mean_fap([1 30 60 95 128 160 184 200 212 223]), recall_10([1 30 60 95 128 160 184 200 212 223]),'b-d','linewidth',1.8);

load botnet_nonpar_cusum_results_v3.mat;
semilogx(1./mean_fap([8 13 17 20]), recall_10([8 13 17 20]),'r-+','linewidth',1.8);

load botnet_ODIT_results.mat;
semilogx(1./mean_fap([1 14 20 24 25 26 30 40 50 64]), recall_10([1 14 20 24 25 26 30 40 50 64]),'g-o','linewidth',1.8);

load botnet_QuantTree_results.mat;
semilogx(1./mean_fap([1 50 72 78 86 94 103 110 120 130 140 152 158 162 166]), recall_10([1 50 72 78 86 94 103 110 120 130 140 152 158 162 166]),'c-^','linewidth',1.8);

load botnet_NN_heavy_results.mat;
semilogx(1./mean_fap([1 74 140 180 220 250 280 332 360 396]), recall_10([1 74 140 180 220 250 280 332 360 396]),'y-*','linewidth',1.8);

xlabel('FAR');
ylabel('TPR');
xlim([1e-4 1]);
legend('Algorithm 1*','Algorithm 2','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','southeast');


% FAP, Approximation, Lower Bound
load botnet_GEM_based_reduced_dim_results.mat;
figure; semilogy(h,mean_fap,'k','linewidth',2); hold on;

load botnet_PCA_based_results.mat;
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

