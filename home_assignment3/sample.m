function [lambda,theta, t, burn_in, acc_ratio] = sample(psi,rho, tau, d, N)

if nargin < 5 
    N = 10000; % default
end
burn_in = N/5;
M = N + burn_in; % Ensuring we always analyze N samples. note that this means 
first = 1658;% that we always simulate 20 % more samples. 
last = 1980;
t = ones(d+1,M);
s = linspace(first,last,d+1)';
t = s.*t; % initializing breakpoints as evenly spaced. 

accepted = 0; % for calculating mean acceptance rate across all the breakpoints. 

theta = zeros(M,1);
theta(1) = gamrnd(2,1/psi); % first theta initialized by its prior. 

lambda = zeros(d,M); % OBS osäkert med proportionaliteten. SNIS?
for k = 1:d
    lambda(k,1) = gamrnd(n_i(tau,t(k,1),t(k+1,1))+2,1/(theta(k)+t(k+1,1)-t(k,1)));
end
%^initializing with equal (very low) intensities.
for j = 1:M-1
   theta(j+1) = gamrnd(2*d+2,1/(psi + sum(lambda(:,j)))); 
   for i = 2:d
      R = rho*(t(i+1,j)-t(i-1,j+1));
      t_cand = t(i,j)+(rand(1)-1/2)*2*R; % candidate generated.
      if (t_cand <= t(i+1,j) && t_cand >= t(i-1,j+1)) % checking for violations of ordering of breakpoints. 
          %See documentation on f_t.m
          f_cand = f_t([lambda(1,j);lambda(2:i-1,j+1);lambda(i:end,j)],[t(1:i-1,j+1);t_cand;t(i+1:end,j)], tau);
          f_stay = f_t([lambda(1,j);lambda(2:i-1,j+1);lambda(i:end,j)],[t(1:i-1,j+1); t(i:end,j)], tau);
          alpha = min(1,f_cand/f_stay); % MH acceptance probability. 

          if rand <= alpha
              t(i,j+1)=t_cand;
              accepted =  accepted + 1;
          else
              t(i,j+1) = t(i,j);
          end
      else
          t(i,j+1) = t(i,j);
      end
      lambda(i,j+1) = gamrnd(n_i(tau,t(i,j+1),t(i+1,j))+2,1/(theta(j+1)+t(i+1,j)-t(i,j+1))); 
   end % gibbs sampling on lambdas. note lambda_1 is sampled separately as there is always one more lambda than breakpoints. 
    lambda(1,j+1) = gamrnd(n_i(tau,t(1,j+1),t(2,j+1))+2,1/(theta(j+1)+t(2,j)-t(1,j+1)));
end

acc_ratio = accepted/(M*(d-1)); 
% broad measurement of acceptance rate for the MH part, i.e. sampling 
% the breakpoints. 
% note not accounting for every breakpoint, but as an average across
% breakpoints. 
end

