%% HA3 Sundell & Ekman
close all
clear all  
clc
tau = load ('coal_mine_disasters.mat').T; 

histogram(tau,'BinWidth',5)
title('Histogram of events over the years (tau)')
% this^ plot is used in the report. 

%% Hybrid MCMC

% this entire section is the initial sketch for the algorithm. 
% For the majority of this script, the function sample.m was used, which 
% is largely the same implementation as presented in this section. 

rho = 0.03; % initially chosen quite randomly, this is refined in following tasks.
psi = 30; % same as with rho.   
d = 2; 
N = 10000; 
burn_in = N/5;
M = N + burn_in;
first = 1658;
last = 1980;
t = ones(d+1,M);
s = linspace(first,last,d+1)';
t = s.*t; 

theta = ones(M,1);
lambda = ones(d,M);

for j = 1:M-1
   theta(j+1) = gamrnd(2*d+2,1/(psi + sum(lambda(:,j)))); 
   for i = 2:d
      R = rho*(t(i+1,j)-t(i-1,j+1));
      t_cand = t(i,j)+(rand(1)-1/2)*2*R;
      if (t_cand <= t(i+1,j) && t_cand >= t(i-1,j+1))
          
          f_cand = f_t([lambda(1,j);lambda(2:i-1,j+1);lambda(i:end,j)],[t(1:i-1,j+1);t_cand;t(i+1:end,j)], tau);
          f_stay = f_t([lambda(1,j);lambda(2:i-1,j+1);lambda(i:end,j)],[t(1:i-1,j+1); t(i:end,j)], tau);
          alpha = min(1,f_cand/f_stay);

          if rand <= alpha
              t(i,j+1)=t_cand;
          else
              t(i,j+1) = t(i,j);
          end
      else
          t(i,j+1) = t(i,j);
      end
      lambda(i,j+1) = gamrnd(n_i(tau,t(i,j+1),t(i+1,j))+2,1/(theta(j+1)+t(i+1,j)-t(i,j+1)));      
   end
    lambda(1,j+1) = gamrnd(n_i(tau,t(1,j+1),t(2,j+1))+2,1/(theta(j+1)+t(2,j)-t(1,j+1)));
end
close all
plot(t(:,burn_in:end)')
legend() % this command labels the data, however no easy way of dynamically 
% naming all data series, hence "data1" etc. 
figure
plot(lambda(:,burn_in:end)')
legend()

%% 1c) behaviour for different d. 
close all
rhos = [0.015 0.02 0.026 0.03 0.05]; 
% these rhos were found to yield acceptance rate close to 30 % for each d. 
psi = 30; % this is the choice of psi that seemed suitable, see report. 

% plots for each d, see report. 
for d = 2:6
   [lambda,~,t,burn_in,acc_ratio] = sample(psi,rhos(d-1),tau,d); 
   figure
   plot(t(:,burn_in:end)')
   title(strcat('Breakpoints, including endpoints, d = ',string(d)))
   legend()
   figure
   for b = 1:d-1
       histogram(t(b+1,burn_in:end)','BinWidth',1)
       hold on
   end
   title(strcat('Histogram of breakpoint placement, d = ', string(d)))
   legend()
   hold off
   figure
   plot(lambda(:,burn_in:end)')
   title(strcat('Intensities, d = ', string(d)))
   legend()
%    d 
%    acc_ratio
end


%% Checking out best rho for each d
close all
psi = 30;
n = 20;
rhos = linspace(0.001,0.5,n);

% this is a very computationally demanding section to run. 
% however, proved very useful to find the most suitable rhos for each d.
for d = 2:6
   collect_accs = zeros(1,n);
   for i= 1:n
      [~,~,~,~,acc_ratio] = sample(psi,rhos(i),tau,d);
      collect_accs(i) = acc_ratio;
   end
   figure
   plot(rhos,collect_accs)
   title(strcat('Acceptence rate over Rho, d = ', string(d)))
end


%% 1d) depending on Psi
close all
d = 5;
rho = 0.03;
n = 100;
psis = linspace(0.1,40,n);

mean_lambdas = zeros(d,n);
var_lambdas = zeros(d,n);

mean_thetas = zeros(n,1);
var_thetas = zeros(n,1);

mean_ts = zeros(d+1,n);
var_ts = zeros(d+1,n);
acc_rates = zeros(1,n);
% pretty self-explanatory, collecting some interesting quantities for the 
% different values of psi. note that d = 5 was fixed. 
for i = 1:n 
   [lambda, theta, t,burn_in, acc_ratio] = sample(psis(i), rho, tau, d);
   mean_lambdas(:,i) = mean(lambda(:,burn_in:end)')';
   var_lambdas(:,i) = var(lambda(:,burn_in:end)')';
   
   mean_thetas(i) = mean(theta(burn_in:end));
   var_thetas(i) = var(theta(burn_in:end));
   
   mean_ts(:,i) = mean(t(:,burn_in:end)')';
   var_ts(:,i) = var(t(:,burn_in:end)')';
   acc_rates(i) = acc_ratio;
end



%% plotting for Psis
close all
plot(psis,mean_lambdas')
title('Mean of lambdas over psi')
figure 
plot(psis,var_lambdas')
title('Variance of lambdas over psi')
figure 
plot(psis,mean_thetas)
title('Mean of theta over psi')
figure
plot(psis,var_thetas)
title('Variance of theta over psi')
figure
plot(psis, mean_ts')
title('Mean of breakpoints over psi')
figure 
plot(psis, var_ts')
title('Variance of breakpoints over psi')
figure 
plot(psis,acc_rates)
title('Acceptance rate over psi')

%% 1 e) Dependence on Rho
% identical setup as in 1.d), but varying rho in stead of psi. 

close all
d = 5;
psi = 30;
n = 100;
rhos = linspace(0.001,2,n);

mean_lambdas = zeros(d,n);
var_lambdas = zeros(d,n);

mean_thetas = zeros(n,1);
var_thetas = zeros(n,1);

mean_ts = zeros(d+1,n);
var_ts = zeros(d+1,n);
acc_rates = zeros(1,n);

for i = 1:n 
   [lambda, theta, t,burn_in, acc_ratio] = sample(psi, rhos(i), tau, d);
   mean_lambdas(:,i) = mean(lambda(:,burn_in:end)')';
   var_lambdas(:,i) = var(lambda(:,burn_in:end)')';
   
   mean_thetas(i) = mean(theta(burn_in:end));
   var_thetas(i) = var(theta(burn_in:end));
   
   mean_ts(:,i) = mean(t(:,burn_in:end)')';
   var_ts(:,i) = var(t(:,burn_in:end)')';
   acc_rates(i) = acc_ratio;
end

%% plotting for Rhos:

close all
plot(rhos,mean_lambdas')
title('Mean of lambdas over rho')
figure 
plot(rhos,var_lambdas')
title('Variance of lambdas over rho')
figure 
plot(rhos,mean_thetas)
title('Mean of theta over rho')
figure
plot(rhos,var_thetas)
title('Variance of theta over rho')
figure
plot(rhos, mean_ts')
title('Mean of breakpoints over rho')
figure 
plot(rhos, var_ts')
title('Variance of breakpoints over rho')
figure 
plot(rhos,acc_rates)
title('Acceptance rate over rho')
% acc ratio as indicator means rho should be around 0.03

%% 2 - bootstrapping

% All of task 2 is given in detail in the report. 

load atlantic.txt
alpha = 0.05;
[beta_est,mu_est] = est_gumbel(atlantic);

n = length(atlantic);
B = 1000;

F_inv = @(u,beta,mu) mu - beta.*log(-log(u));

boot = zeros(2,B);

for b = 1:B
   u = rand(1,n);
   y_sim = F_inv(u,beta_est,mu_est);
   [beta,mu] = est_gumbel(y_sim);
   boot(1,b) = beta;
   boot(2,b) = mu;
end
delta_beta = sort(boot(1,:)-beta_est);
delta_mu = sort(boot(2,:)-mu_est);
bias_beta = mean(delta_beta);
bias_mu = mean(delta_mu);
L_beta = beta_est - bias_beta-delta_beta(ceil((1-alpha/2)*B));
U_beta = beta_est-bias_beta -delta_beta(ceil(alpha*B/2));
L_mu = mu_est-bias_mu - delta_mu(ceil((1-alpha/2)*B));
U_mu = mu_est-bias_mu - delta_mu(ceil(alpha*B/2));


%% 2c)
T = 3*14*100;

wave_est = F_inv(1-1/T,beta_est,mu_est);
wave_boot= zeros(1,B);

for i = 1:B
   wave_boot(i) = F_inv(1-1/T,boot(1,i),boot(2,i)); 
end

delta_wave = sort(wave_boot-wave_est);
bias_wave = mean(delta_wave);
U_wave = wave_est-bias_wave - delta_wave(ceil(alpha*B));

%% some plots - for the curious. 
% Mainly to see that the deltas don't seem very skewed. Bias is already
% handled. 
close all
histogram(delta_beta)
title(strcat('Histogram of error for beta, mean: ', string(mean(delta_beta))))
figure
histogram(delta_mu)
title(strcat('Histogram of error for mu, mean: ', string(mean(delta_mu))))
figure 
histogram(delta_wave)
title(strcat('Histogram of error for wave, mean: ', string(mean(delta_wave))))
