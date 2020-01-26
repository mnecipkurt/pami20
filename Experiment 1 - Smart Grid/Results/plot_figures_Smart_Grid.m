clear variables; close all;

% ADD vs. FAP

load Smart_grid_GEM_based_results.mat
figure; plot(mean_fap([40 272 385 456 495 513 526 539]), mean_add([40 272 385 456 495 513 526 539]),'k-s','linewidth',1.8); hold on; grid on;

load Smart_grid_ITMCD_results.mat;
plot(mean_fap([10 230 600 950 1042 1118 1138 1168]), mean_add_FDI([10 230 600 950 1042 1118 1138 1168]),'b-d','linewidth',1.8);

load Smart_grid_nonpar_cusum_ODIT_results.mat;
plot(mean_fap([18 77 124 167 220 251 280 320]), mean_add([18 77 124 167 220 251 280 320]),'r-+','linewidth',1.8);
plot(mean_fap_odit([20 52 67 77 84 90 94 98]), mean_add_odit([20 52 67 77 84 90 94 98]),'g-o','linewidth',1.8);

load Smart_grid_QuantTree_results.mat;
plot(mean_fap([70 158 211 262 290 315 330 343 351]), mean_add([70 158 211 262 290 315 330 343 351]),'c-^','linewidth',1.8);

load Smart_grid_NN_heavy_results.mat;
plot(mean_fap([169 255 275 281 286 289 293 294]), mean_add([169 255 275 281 286 289 293 294]),'y-*','linewidth',1.8);

xlabel('$\rm{E}_\infty[\Gamma]$','interpreter','latex','fontsize',14);
ylabel('$\rm{E}_\tau \big[(\Gamma-\tau)^+\big]$','interpreter','latex','fontsize',14);
xlim([1 1e4]);
legend('Algorithm 1','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','northwest');

axes('position',[0.42 0.34 0.46 0.26]); box on; 
load Smart_grid_GEM_based_results.mat
plot(mean_fap([40 272 385 456 495 513 526 539]), mean_add([40 272 385 456 495 513 526 539]),'k-s','linewidth',1.8); hold on; grid on;
load Smart_grid_nonpar_cusum_ODIT_results.mat;
plot(mean_fap_odit([20 52 67 77 84 90 94 98]), mean_add_odit([20 52 67 77 84 90 94 98]),'g-o','linewidth',1.8);
axis tight; grid on;

% Recall vs. FAR
load Smart_grid_GEM_based_results.mat
figure; semilogx(1./mean_fap([1 70 150 272 385 456 539]), recall_10([1 70 150 272 385 456 539]),'k-s','linewidth',1.8); hold on; grid on;

load Smart_grid_ITMCD_results.mat;
semilogx(1./mean_fap([1 110 230 400 600 780 950 1118 1168]), recall_10([1 110 230 400 600 780 950 1118 1168]),'b-d','linewidth',1.8);

load Smart_grid_nonpar_cusum_ODIT_results.mat;
semilogx(1./mean_fap([1 18 77 130 148 160 180 270 320]), recall_10([1 18 77 130 148 160 180 270 320]),'r-+','linewidth',1.8);
semilogx(1./mean_fap_odit([1 20 32 52 67 77 90 98]), recall_10_odit([1 20 32 52 67 77 90 98]),'g-o','linewidth',1.8);

load Smart_grid_QuantTree_results.mat;
semilogx(1./mean_fap([1 100 132 158 204 260 290 322 351]), recall_10([1 100 132 158 204 260 290 322 351]),'c-^','linewidth',1.8);

load Smart_grid_NN_heavy_results.mat;
semilogx(1./mean_fap([1 160 200 255 275 281 286 289 292 294]), recall_10([1 160 200 255 275 281 286 289 292 294]),'y-*','linewidth',1.8);

xlabel('FAR');
ylabel('TPR');
xlim([1e-4 1]);
legend('Algorithm 1','ITMCD','Nonpar. CUSUM','ODIT','QuantTree','NN-Based','location','southeast');

% FAP, Approximation, Lower Bound
load Smart_grid_GEM_based_results.mat;
theta = lambertw(alpha*log(alpha))/log(alpha);
w_0 = theta - 1;
lower_bnd_fap = exp(-w_0*h);
apprx_fap = lower_bnd_fap*10.1;
Wald_apprx_fap = (h + (exp(-w_0*h)-1)/w_0)/(1+log(alpha));
figure; semilogy(h,mean_fap,'k','linewidth',2); hold on;
semilogy(h,apprx_fap,'b-.','linewidth',2);
semilogy(h,Wald_apprx_fap,'c-.','linewidth',2);
semilogy(h,lower_bnd_fap,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period (Alg.1)','Our Approximation',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;
xlim([2 12])
