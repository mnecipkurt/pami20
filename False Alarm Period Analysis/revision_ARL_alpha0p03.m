% compare the ARL bound and approximation to the Monte Carlo simulation
% in the asymptotic regime where N_2 goes to infinity (tail prob. is uniform)

clear variables; clc;
alpha = 0.03;
h = 1:0.01:7;
theta = lambertw(alpha*log(alpha))/log(alpha);
w_0 = theta - 1;
theory_apprx_fap = (h + (exp(-w_0*h)-1)/w_0)/(1+log(alpha));
theory_lower_bnd_fap = exp(-w_0*h);

no_trials = 500;
sum_alarm_times = zeros(1,length(h));

for i=1:no_trials
    i
    g = 0;
    k = 1;
    false_alarm_flag = zeros(1,length(h));
    while (false_alarm_flag(length(h)) == 0)
        u = rand(1);
        s = log(alpha/u);
        g = max(0, g+s); 
        sum_alarm_times = sum_alarm_times + k*(g >= h).*(false_alarm_flag == 0);
        false_alarm_flag = false_alarm_flag + (g >= h).*(false_alarm_flag == 0);
        k = k+1;
    end      
end
simulation_avg_fap = sum_alarm_times/no_trials;

figure; plot(h,simulation_avg_fap,'b-o','linewidth',2); hold on;
plot(h,theory_lower_bnd_fap,'r-o','linewidth',2);
plot(h,theory_apprx_fap,'k-o','linewidth',2);
legend('Monte Carlo','Lower Bound','Approximation','location','northwest')

figure; plot(theory_lower_bnd_fap,simulation_avg_fap,'g-o','linewidth',2);
xlabel('Lower Bound'); ylabel('Monte Carlo')

save('results_alpha_0p03');

