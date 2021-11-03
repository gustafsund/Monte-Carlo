%% Setup
close all
clear all
clc
load powercurve_V164

% Agnes Ekman (ag8720ek-s)
% Gustaf Sundell (gu0147su-s)


n = 10000;
% sample size for all simulations.
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6]';
k = [2 2 2 1.9 1.9 1.9 1.9 1.9 2 1.9 2 2]';
% Weibull parameters for the 12 months of the year. 

z = 1.96; % approximate Normal quantile of interest.
a = 3.5;
b = 25; % [a,b] makes up the power producing interval.

%% Setup och standard method 2.a)
close all
tauN_stand = zeros(12,3); 
% saving fit and CI in matrices like this one for each method across all
% months. 

% rng(1)
% for recreating the exact numbers in the report.

for m = 1:12
    wind = wblrnd(lambda(m),k(m),n,1);
    power = P(wind);
    mean_power = mean(power);
    sigma = std(power);
    tauN_stand(m,2) = mean_power;
    tauN_stand(m,3) = mean_power + z*sigma/sqrt(n);
    tauN_stand(m,1) = mean_power - z*sigma/sqrt(n);  
end

plot(tauN_stand(:,2),'b')
hold on 
plot(tauN_stand(:,1),'-.r')
plot(tauN_stand(:,3),'-.r')
hold off
width_stand = tauN_stand(:,3) - tauN_stand(:,1);
% width and inside for table in report. 

inside_stand = zeros(12,1); % standard is of course inside its own CI
%% Truncated version 2.a)
% rng(1);

tauN_trunc = zeros(12,3);

for m = 1:12
u = rand(n,1);
ind = u~=1 & u~= 0;% Generating x in [3.5,25] by means of inverse method. 
wind = wblinv((wblcdf(b, lambda(m), k(m)) - wblcdf(a,lambda(m), k(m)))*u(ind) + wblcdf(a,lambda(m),k(m)),lambda(m),k(m));
power = P(wind);
mean_power = mean(power).*(wblcdf(b,lambda(m),k(m)) - wblcdf(a,lambda(m),k(m)));
sigma = std(power);
tauN_trunc(m,2) = mean_power;
tauN_trunc(m,3) = mean_power + z*sigma/sqrt(n);%Bör ta length(u(ind)) ist för n?
tauN_trunc(m,1) = mean_power - z*sigma/sqrt(n);
end
figure
plot(tauN_trunc(:,2),'b')
hold on 
plot(tauN_trunc(:,1),'-.r')
plot(tauN_trunc(:,3), '-.r')
width_trunc= tauN_trunc(:,3) - tauN_trunc(:,1);

inside_trunc = zeros(12,1); % checking for overlapping CI's
for m = 1:12
if tauN_trunc(m,1) > tauN_stand(m,3)
inside_trunc(m) = 1;
       
elseif tauN_trunc(m,3) < tauN_stand(m,1)
inside_trunc(m) = - 1;
    
end
   
end
inside_trunc;
% this procedure should indicate, if no overlap, whether it over- or
% undershoots the reference CI.

%% importance sampling, plotting and looking for suitable g
close all
k_gam = 7;
theta_gam = 2;
% setting different parameters for Gamma dist (our g), 7 and 2 seemed to be
% the best fit. 

x = linspace(1,30);

phif = zeros(12,100);
fs = zeros(12,100);
for m = 1:12
    t = wblpdf(x,lambda(m),k(m));
    pwr = P(x);
    phif(m,:) = t'.*pwr;
    fs(m,:) = t;
end

% Plotting all relevant curves, to see how well g behaves. 
plot(x,1e8*fs,'b');
hold on
plot(x, P(x))
hold on 
plot(x,1e8*gampdf(x,k_gam,theta_gam))
plot(x,20*phif,'y')

%% Importance sampling 2.b)

tauN_imp = zeros(12,3);

for m = 1:12
wind = gamrnd(k_gam,theta_gam,[n,1]);
power = P(wind).*wblpdf(wind,lambda(m),k(m))./gampdf(wind,k_gam,theta_gam);
mean_power = mean(power);
sigma = std(power);
tauN_imp(m,2) = mean_power;
tauN_imp(m,3) = mean_power + z*sigma/sqrt(n);
tauN_imp(m,1) = mean_power - z*sigma/sqrt(n);
end
figure
plot(tauN_imp(:,2),'b')
hold on 
plot(tauN_imp(:,1),'-.r')
plot(tauN_imp(:,3), '-.r')
hold off
width_imp= tauN_imp(:,3) - tauN_imp(:,1);

inside_imp= zeros(12,1);
for m = 1:12
if tauN_imp(m,1) > tauN_stand(m,3)
inside_imp(m) = 1;
       
elseif tauN_imp(m,3) < tauN_stand(m,1)
inside_imp(m) = - 1;
    
end
   
end
% inside_imp

%% 2.c) Antithetic sampling
close all
% rng(1);
tauN_anti = zeros(12,3);

for m = 1:12
    u = rand(n,1);
    u_tilde = 1-u;
    ind = u~=1 & u~=0;
    wind = wblinv((wblcdf(b, lambda(m), k(m)) - wblcdf(a,lambda(m), k(m)))*u(ind) + wblcdf(a,lambda(m),k(m)),lambda(m),k(m));
    wind_tilde = wblinv((wblcdf(b, lambda(m), k(m)) - wblcdf(a,lambda(m), k(m)))*u_tilde(ind) + wblcdf(a,lambda(m),k(m)),lambda(m),k(m));
    power = P(wind);
    power_tilde = P(wind_tilde);
    W = (power + power_tilde)/2;
    mean_power = mean(W).*(wblcdf(b,lambda(m),k(m)) - wblcdf(a,lambda(m),k(m)));
    sigma = std(W);
    tauN_anti(m,2) = mean_power;
    tauN_anti(m,3) = mean_power + z*sigma/sqrt(n);
    tauN_anti(m,1) = mean_power - z*sigma/sqrt(n);% using n consistently,  as length(W) = n = 10000;
end

plot(tauN_anti(:,2),'b')
hold on 
plot(tauN_anti(:,1),'-.r')
plot(tauN_anti(:,3), '-.r')
width_anti= tauN_anti(:,3) - tauN_anti(:,1)

inside_anti= zeros(12,1);
for m = 1:12
if tauN_anti(m,1) > tauN_stand(m,3)
inside_anti(m) = 1;
       
elseif tauN_anti(m,3) < tauN_stand(m,1)
inside_anti(m) = - 1;
    
end
   
end
% inside_anti
%% 2.d)
prob_power = zeros(12,1);
for m = 1:12
   prob_power(m) = wblcdf(b,lambda(m),k(m)) - wblcdf(a,lambda(m),k(m)); 
end
% prob_power
% P(V in I) = F(b) - F(a)
%% 2.e) Power coefficient
close all
rho = 1.225;
d = 164;
B = (1/8)*rho*pi*d.^2; % collecting the constants
% rng(1);

tauN_ratio = zeros(12,3);
for m = 1:12 
   ptot = B*lambda(m).^3*gamma(1+3/k(m));
   tauN_ratio(m,2) = tauN_anti(m,2)/ptot;
   tauN_ratio(m,1) = tauN_anti(m,1)/ptot;
   tauN_ratio(m,3) = tauN_anti(m,3)/ptot;
   % scaling width of CI by the same division of deterministic constant.
end

plot(tauN_ratio(:,2),'b')
hold on 
plot(tauN_ratio(:,1),'-.r')
plot(tauN_ratio(:,3), '-.r')


%% 2.f) 
% rng(1);
capacity = zeros(12,1);

for m = 1:12
   wind = wblrnd(lambda(m), k(m),n,1);
   power = P(wind);
   capacity(m,1) = sum(power)/(9.5*1e6*n);
end
mean(capacity)
mean(prob_power) % availability over 12 months. 
% quite reasonable to build here, although a bit too low availability. 

%% Some input for table in report
tab_width = [width_stand width_trunc width_imp width_anti]
tab_inside = [inside_stand inside_trunc inside_imp inside_anti] 

% All methods CI:s overlap with the reference. At least with rng(1), but
% also when trying a couple of times without seed.


%% Problem 3 setup
close all
% rng(1);
n = 10000;
lambda = 9.13;
k = 1.96;
alfa = 0.638;
p = 3;
q = 1.5;

% setting up relevant functions. 
f = @(v) wblpdf(v,lambda,k);
F = @(v) wblcdf(v,lambda,k);
f_joint = @(v1,v2) f(v1).*f(v2).*(1 + alfa*(1-F(v1).^p).^(q-1).*(1-F(v2).^p).^(q-1).*(F(v1).^p.*(1+p.*q)-1).*(F(v2).^p.*(1+p.*q)-1));
g = @(v) wblpdf(v,lambda,k);
g_joint = @(v1,v2) wblpdf(v1,lambda,k).*wblpdf(v2,lambda,k);
P_sum = @(v1,v2) P(v1) + P(v2);

% Checking the behaviour of g(v1,v2) together with Phi_prod and Phi_sum,
% see report.

x = linspace(0,30,100);
[X,Y] = meshgrid(x,x);

PSUM = zeros(100,100);
PPROD = zeros(100,100);
for i = 1:100
    for j = 1:100
        PSUM(i,j) = P_sum(x(i),x(j));
        PPROD(i,j) = P(x(i))*P(x(j));
    end
end
mesh(X,Y, PSUM.*f_joint(X,Y)./g_joint(X,Y))
figure
mesh(X,Y, PPROD.*f_joint(X,Y)./g_joint(X,Y))
%% 3.a) 
% this reduces to a 1 dim problem, so we just simulate from 
% one weibull distribution, see report for reasoning. 
v = wblrnd(lambda,k,[n, 1]);

tau_sum = 2*mean(P(v))

%% 3.b) 

v1 = wblrnd(lambda, k,[n 1]);
v2 = wblrnd(lambda,k,[n 1]);
P1 = P(v1); 
P2 = P(v2);
% note, for one dim, f/g = 1 once again.

% see Phi_prod in report.
P_prod = P1.*P2.*f_joint(v1,v2)./g_joint(v1,v2);
cov_star = mean(P_prod) - mean(P1)*mean(P2)
var(P1) 
% for comparison, seems to be that cov and var are reasonably 
% scaled in relation to each other

%% 3.c) 
% Simply taking variance and stdev from simulated Phi_sum.
plus = P_sum(v1,v2).*f_joint(v1,v2)./g_joint(v1,v2);

variance = var(plus)
dev = std(plus)

%% 3.d)
% 'success' being over 9.5
Ns = sum(plus > 9.5e6);
p_over = Ns/n

% 'fail' being under 9.5
Nf = sum(plus < 9.5e6);
p_under = Nf/n

sgm = sqrt(p_over.*p_under./n);
% CI:s
upr_over = p_over + z*sgm
lwr_over = p_over - z*sgm

upr_under = p_under + z*sgm
lwr_under = p_under - z*sgm

