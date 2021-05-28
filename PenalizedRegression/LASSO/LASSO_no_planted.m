%% RS saddle point equations for LASSO
%
%  Caution: Macroscopic quantities are normalized by alpha.
%           Rescale them appropriately.
%
clear variables;
close all;

THETA = 1e-10;

alpha = input('alpha = ');
lambda = input('lambda = ');
% mean of data
mu_y = 0;
% variance of data
sigma_y2 = 1;

% set appropriate initial condition
Q = 1;
chi = 200;  
% damping factor
dam = 0.1;

shusoku = 0;
while shusoku == 0
    chi_old = chi;
    Q_old = Q;
    
    Qh = 1/(1+chi);
    chih = (Q+sigma_y2+mu_y^2)/(1+chi)^2;
    theta = lambda/sqrt(2*chih);
    rhoh = erfc(theta);
    
    chi = rhoh/(alpha-rhoh);
    Q_tmp = (1+2*theta^2)*rhoh-2*theta*exp(-theta^2)/sqrt(pi);
    Q = (sigma_y2+mu_y^2)*Q_tmp/(alpha-Q_tmp);
    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA);
    
    Q = dam*rand(1)*(Q-Q_old)+Q_old;
    chi = dam*rand(1)*(chi-chi_old)+chi_old;
end
% training error
err_train = chih;
% expectation value of the L1 penalty
l1 = chi*chih-Q*Qh;
% generalized degree of freedom
df = rhoh/alpha;
