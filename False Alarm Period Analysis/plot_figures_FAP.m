clear variables; clc;

h_slopes = [1040 101 51 34.1 21.8 12.1 9.9 10.1 13 25.8 67 109 230 980];

% load results_alpha_0p001.mat;
% low_bnd_vec1 = theory_lower_bnd_fap;
% ARL_vec1 = simulation_avg_fap;
% app_vec1 = theory_apprx_fap;

load results_alpha_0p01.mat;
low_bnd_vec2 = theory_lower_bnd_fap;
ARL_vec2 = simulation_avg_fap;
app_vec2 = theory_apprx_fap;

% load results_alpha_0p02.mat;
% low_bnd_vec3 = theory_lower_bnd_fap;
% ARL_vec3 = simulation_avg_fap;
% app_vec3 = theory_apprx_fap;
% 
% load results_alpha_0p03.mat;
% low_bnd_vec4 = theory_lower_bnd_fap;
% ARL_vec4 = simulation_avg_fap;
% app_vec4 = theory_apprx_fap;

load results_alpha_0p05.mat;
low_bnd_vec5 = theory_lower_bnd_fap;
ARL_vec5 = simulation_avg_fap;
app_vec5 = theory_apprx_fap;

load ARL_alpha_0p1.mat;
low_bnd_vec6 = theory_lower_bnd_fap;
ARL_vec6 = simulation_avg_fap;
app_vec6 = theory_apprx_fap;
h_6 = h';

load ARL_alpha_0p15.mat;
low_bnd_vec7 = theory_lower_bnd_fap;
ARL_vec7 = simulation_avg_fap;
app_vec7 = theory_apprx_fap;

load ARL_alpha_0p2.mat;
low_bnd_vec8 = theory_lower_bnd_fap;
ARL_vec8 = simulation_avg_fap;
app_vec8 = theory_apprx_fap;
h_8 = h';

load ARL_alpha_0p25.mat;
low_bnd_vec9 = theory_lower_bnd_fap;
ARL_vec9 = simulation_avg_fap;
app_vec9 = theory_apprx_fap;

load ARL_alpha_0p3.mat;
low_bnd_vec10 = theory_lower_bnd_fap;
ARL_vec10 = simulation_avg_fap;
app_vec10 = theory_apprx_fap;
h_10 = h';

% load ARL_alpha_0p33.mat;
% low_bnd_vec11 = theory_lower_bnd_fap;
% ARL_vec11 = simulation_avg_fap;
% app_vec11 = theory_apprx_fap;
% 
% load ARL_alpha_0p34.mat;
% low_bnd_vec12 = theory_lower_bnd_fap;
% ARL_vec12 = simulation_avg_fap;
% app_vec12 = theory_apprx_fap;

load ARL_alpha_0p35.mat;
low_bnd_vec13 = theory_lower_bnd_fap;
ARL_vec13 = simulation_avg_fap;
app_vec13 = theory_apprx_fap;
h_13 = h';

% load ARL_alpha_0p36.mat;
% low_bnd_vec14 = theory_lower_bnd_fap;
% ARL_vec14 = simulation_avg_fap;
% app_vec14 = theory_apprx_fap;

figure; hold on;
% plot(low_bnd_vec1,ARL_vec1,'linewidth',3); 
% plot(low_bnd_vec2,ARL_vec2,'linewidth',3);
% plot(low_bnd_vec3,ARL_vec3,'linewidth',3);
% plot(low_bnd_vec4,ARL_vec4,'linewidth',3);
plot(low_bnd_vec5,ARL_vec5,'linewidth',3);
plot(low_bnd_vec6,ARL_vec6,'linewidth',3);
plot(low_bnd_vec7,ARL_vec7,'linewidth',3);
plot(low_bnd_vec8,ARL_vec8,'linewidth',3);
plot(low_bnd_vec9,ARL_vec9,'linewidth',3);
plot(low_bnd_vec10,ARL_vec10,'linewidth',3);
% plot(low_bnd_vec11,ARL_vec11,'linewidth',3);
% plot(low_bnd_vec12,ARL_vec12,'linewidth',3);
% plot(low_bnd_vec13,ARL_vec13,'linewidth',3);
% plot(low_bnd_vec14,ARL_vec14,'linewidth',3);

tmpp_vecc = 1:5:1e4;
% plot(tmpp_vecc,tmpp_vecc*h_slopes(1),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(2),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(3),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(4),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(5),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(6),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(7),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(8),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(9),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(10),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(11),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(12),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(13),'b-o','linewidth',1.1);
% plot(tmpp_vecc,tmpp_vecc*h_slopes(14),'b-o','linewidth',1.1);


xlabel('Lower Bound');
ylabel('False Alarm Period');
leg1 = legend('$\alpha = 0.05$', '$\alpha = 0.1$', '$\alpha = 0.15$','$\alpha = 0.2$','$\alpha = 0.25$','$\alpha = 0.3$','location','northwest');
set(leg1,'Interpreter','latex');
grid on; box on;

% leg1 = legend('$\alpha = 0.001$','$\alpha = 0.01$', '$\alpha = 0.02$', '$\alpha = 0.03$', '$\alpha = 0.05$', '$\alpha = 0.1$', ...
% '$\alpha = 0.15$','$\alpha = 0.2$','$\alpha = 0.25$','$\alpha = 0.3$','$\alpha = 0.33$','$\alpha = 0.34$','$\alpha = 0.35$', ...
% '$\alpha = 0.36$','location','northwestoutside');
% set(leg1,'Interpreter','latex');

xlim([0 1e4]);
% ylim([0 4e4])

% Wald's approximation Curves
figure; semilogy(h_6,ARL_vec6,'k','linewidth',2); hold on;
semilogy(h_6,app_vec6,'c-.','linewidth',2);
semilogy(h_6,low_bnd_vec6,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period','Wald Approximation','Lower Bound','location','northwest');
grid on; box on;
%%%%%%%%%%%%%%%%%%%%%%%%
figure; semilogy(h_8,ARL_vec8,'k','linewidth',2); hold on;
semilogy(h_8,app_vec8,'c-.','linewidth',2);
semilogy(h_8,low_bnd_vec8,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period','Wald Approximation','Lower Bound','location','northwest');
grid on; box on;

%%%%%%%%%%%%%%%%%%%%%%%%
figure; semilogy(h_10,ARL_vec10,'k','linewidth',2); hold on;
semilogy(h_10,app_vec10,'c-.','linewidth',2);
semilogy(h_10,low_bnd_vec10,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period','Wald Approximation','Lower Bound','location','northwest');
grid on; box on;

%%%%%%%%%%%%%%%%%%%%%%%%
figure; semilogy(h_13(91:end),ARL_vec13(91:end),'k','linewidth',2); hold on;
semilogy(h_13(91:end),app_vec13(91:end),'c-.','linewidth',2);
semilogy(h_13(91:end),low_bnd_vec13(91:end),'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
legend('False Alarm Period','Wald Approximation','Lower Bound','location','northwest');
grid on; box on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
subplot('position',[0.04 0.57 0.44 0.37]); 
semilogy(h_6,ARL_vec6,'k','linewidth',2); hold on;
semilogy(h_6,app_vec6,'c-.','linewidth',2);
semilogy(h_6,low_bnd_vec6,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
title('$\alpha = 0.1$','Interpreter','latex','fontsize',14);
legend('False Alarm Period',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;

subplot('position',[0.53 0.57 0.44 0.37]);   
semilogy(h_8,ARL_vec8,'k','linewidth',2); hold on;
semilogy(h_8,app_vec8,'c-.','linewidth',2);
semilogy(h_8,low_bnd_vec8,'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
title('$\alpha = 0.2$','Interpreter','latex','fontsize',14);
legend('False Alarm Period',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;

subplot('position',[0.04 0.08 0.44 0.37]); 
semilogy(h_10(301:end),ARL_vec10(301:end),'k','linewidth',2); hold on;
semilogy(h_10(301:end),app_vec10(301:end),'c-.','linewidth',2);
semilogy(h_10(301:end),low_bnd_vec10(301:end),'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
title('$\alpha = 0.3$','Interpreter','latex','fontsize',14);
legend('False Alarm Period',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;

subplot('position',[0.53 0.08 0.44 0.37]); 
semilogy(h_13(91:end),ARL_vec13(91:end),'k','linewidth',2); hold on;
semilogy(h_13(91:end),app_vec13(91:end),'c-.','linewidth',2);
semilogy(h_13(91:end),low_bnd_vec13(91:end),'r-.','linewidth',2);
xlabel('$h$','Interpreter','latex','fontsize',14);
title('$\alpha = 0.35$','Interpreter','latex','fontsize',14);
legend('False Alarm Period',"Wald's Approximation",'Lower Bound','location','northwest');
grid on; box on;


 