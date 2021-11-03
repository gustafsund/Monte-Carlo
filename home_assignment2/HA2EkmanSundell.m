%% HA 2 Monte Carlo 

clear all
close all
clc
n_max = 40; % nbr steps
N = 1000; % nbr simulations
n_runs = 10; % nbr runs for the outer loop, looking at variance in A,mu and gamma.
n_max69 = 10;% nbr steps for questions 6 and 9. 
%% 3. 
% rng(1)

ratio_est = zeros(n_max,1); % Estimation of SAW_ratio (see report). 
T = zeros(2,n_max+1,N); % storing the trajectories in parallel (hence the ndarray)
for k = 1:n_max
successes = ones(n_max,N); % Counting successes   
for i = 1:N % generating trajectories
   u = rand(1); % Taking the random step.
   if u <0.25 %x+1
       T(:,k+1,i) = T(:,k,i) + [1;0];
   elseif u<0.5%x-1
       T(:,k+1,i) = T(:,k,i) + [-1;0];
   elseif u<0.75%y+1
       T(:,k+1,i) = T(:,k,i) + [0;1];
   else %y-1
       T(:,k+1,i) = T(:,k,i) + [0;-1];
   end
   if length(T(:,1:k+1,i)) ~= length(unique(T(:,1:k+1,i)','rows'))
       successes(k,i) = 0;
       break
   end
end
ratio_est(k) = mean(successes(k,:)); % here counting the success rate. 
end

max_pos = zeros(n_max,1); % This vector is filled with k_n(2), see report. 
for i = 1:n_max 
   max_pos(i) = 4.^i; 
end

cn_naive = ratio_est.*max_pos % estimating cn
%% 4. SIS
% rng(1)
cn_sis = zeros(n_max,1); % for storing estimates of cn
T = zeros(2,n_max+1,N);    
w = ones(n_max+1,N); 
s = ones(n_max,N);
for k = 1:n_max 
    for i = 1:N 
       possible = [T(:,k,i) + [1;0], T(:,k,i) - [1;0], T(:,k,i) + [0;1], T(:,k,i) - [0;1]];
       % all the neighbouring points.
       ind = ~ismember(possible',T(:,1:k,i)','rows');
       free = possible(:,ind);
       % The available points^
       s(k,i) = sum(ind); % s as in report. 
       if s(k,i) == 0 % handling if nowhere left to go
           T(:,k+1,i) = T(:,k,i);
           w(k+1,i) = 0; % weight set to zero. 
       else % if SA, we keep going
           l = size(free);
           j = randi(l(2)); %Randomize uniformly between the available points.
           T(:,k+1,i) = free(:,j); % Take step
           w(k+1,i) = w(k,i).*s(k,i); % update weights. 
       end

    end

end

for n = 1:n_max
    cn_sis(n) = mean(w(n+1,:));% estimate normalizing constant! 
end

%% 5. SISR
% rng(1)
% roughly the same setup, many steps are identical to the once taken in
% question 4, comments ommitted for those.
T = zeros(2,n_max+1,N);    
w = ones(n_max+1,N); 
s = ones(n_max,N);

for k = 1:n_max     
    for i = 1:N 
       possible = [T(:,k,i) + [1;0], T(:,k,i) - [1;0], T(:,k,i) + [0;1], T(:,k,i) - [0;1]];
       ind = ~ismember(possible',T(:,1:k,i)','rows');
       free = possible(:,ind);
       s(k,i) = sum(ind);
       if s(k,i) == 0
           T(:,k+1,i) = T(:,k,i);  
           w(k,i) = 0;
       else
           l = size(free);
           j = randi(l(2));
           T(:,k+1,i) = free(:,j);
           w(k,i) = s(k,i); % Implicitly setting previous weights to 1 and multiplying
       end
    end
    w_probs = w(k,:)./sum(w(k,:)); % probabilities of resampling 
    draws = randsample(N,N,true,w_probs); % resampling indexes
    T(:,:,:) = T(:,:,draws); % resampling trajectory.
end

cn_sisr = zeros(n_max,1); % for storing estimates
cn_sisr(1) = mean(w(1,:)); % recursively estimating cn
for r = 2:n_max
    cn_sisr(r) = cn_sisr(r-1)*mean(w(r,:));
end
%% 6. parameter estimation in d = 2 
% rng(1)
collect_est6 = zeros(3,n_runs); % storing n_runs estimates of 3 parameters. 

% this is the same algorithm as in 5, repeated n_runs times. 
for p = 1:n_runs
    T = zeros(2,n_max69+1,N);    
    w = ones(n_max69+1,N); 
    s = ones(n_max69,N);
    for k = 1:n_max69 
        for i = 1:N 
           possible = [T(:,k,i) + [1;0], T(:,k,i) - [1;0], T(:,k,i) + [0;1], T(:,k,i) - [0;1]];
           ind = ~ismember(possible',T(:,1:k,i)','rows');
           free = possible(:,ind);
           s(k,i) = sum(ind);
           if s(k,i) == 0
               T(:,k+1,i) = T(:,k,i);
               w(k,i) = 0;
           else
               l = size(free);
               j = randi(l(2));
               T(:,k+1,i) = free(:,j);
               w(k,i) = s(k,i); % Implicitly setting previous weights to 1 and multiplying
           end
        end
        w_probs = w(k+1,:)./sum(w(k+1,:)); 
        draws = randsample(N,N,true,w_probs);
        T(:,:,:) = T(:,:,draws);
    end

    cn_sisr6 = zeros(n_max69,1);
    cn_sisr6(1) = mean(w(1,:));
    for r = 2:n_max69
        cn_sisr6(r) = cn_sisr6(r-1)*mean(w(r,:));
    end
    
    % Setting up the least square estimate (linear regression)
    
    X = [ones(n_max69,1) (1:n_max69)' log(1:n_max69)'];

    Y = log(cn_sisr6);
    
    %Simultaneous estimation:
    beta = (X'*X)\X'*Y;

    % Save estimations
    collect_est6(1,p) = exp(beta(1)); % A_est, de-transforming from beta
    collect_est6(2,p) = exp(beta(2)); % mu_est, de-transforming from beta
    collect_est6(3,p) = beta(3)+1; % gamma_est, de-transforming from beta

end
% Present estimations with mean and variances, see table in report. 
collect_est6
var(collect_est6(1,:))
var(collect_est6(2,:))
var(collect_est6(3,:))
mean(collect_est6(1,:))
mean(collect_est6(2,:))
mean(collect_est6(3,:))
43/32 % for comparison with theoretical gamma_2 value. 

%% 9. More parameter estimation with higher dimensionality
% rng(1)
d = 5; % dynamic code, can try different dimensions of problem. 

% Once again re-using SISR algorithm from question 5. 
% Almost everything in this section is identical to 6., the only difference
% being that one may choose dimensionality. 
T = zeros(d,n_max69+1,N);    
w = ones(n_max69+1,N); 
s = ones(n_max69,N);

% I stores the steps possible to take. 
I = [eye(d), -eye(d)];
possible = zeros(d,2*d);
collect_est9 = zeros(3,n_runs);
com_mu_bounds = zeros(2,n_runs);
for p = 1:n_runs
    for k = 1:n_max69 
        for i = 1:N 
           for c = 1:2*d % more generalised way of looking at possible
               possible(:,c) = T(:,k,i) + I(:,c);
           end
           ind = ~ismember(possible',T(:,1:k,i)','rows');
           free = possible(:,ind);
           s(k,i) = sum(ind);
           if s(k,i) == 0
               T(:,k+1,i) = T(:,k,i);
               w(k,i) = 0;
           else
               l = size(free);
               j = randi(l(2));
               T(:,k+1,i) = free(:,j);
               w(k,i) = s(k,i); % Implicitly setting previous weights to 1 and multiplying
           end
        end
        w_probs = w(k,:)./sum(w(k,:)); 
        draws = randsample(N,N,true,w_probs);
        T(:,:,:) = T(:,:,draws);
    end

    cn_sisr9 = zeros(n_max69,1);
    cn_sisr9(1) = mean(w(1,:));
    for r = 2:n_max69
        cn_sisr9(r) = cn_sisr9(r-1)*mean(w(r,:));
    end

    X = [ones(n_max69,1) (1:n_max69)' log(1:n_max69)'];

    Y = log(cn_sisr9);
    beta = (X'*X)\X'*Y;

    collect_est9(1,p) = exp(beta(1)); %A_est
    collect_est9(2,p) = exp(beta(2)); % Mu_est  
    collect_est9(3,p)= beta(3)+1; % gamma_est


end 
mean(collect_est9(1,:))
mean(collect_est9(2,:))
mean(collect_est9(3,:))

% Finally, the approximation to compare to. Big O omitted of course. 
asym_bound_mu = 2*d-1 - 1/(2*d) - 3/((2*d).^2) - 16/((2*d).^3) 

