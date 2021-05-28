%% Linear regression penalized by Elastic Net
%
%  Caution: Macroscopic quantities are normalized by alpha.
%           Rescale them appropriately.
%
clear variables;
close all;

% convergence threshold
THETA = 1e-10;

alpha = input('alpha = ');
% coefficient for L1 term
lambda1 = input('lambda1 = ');
% coefficient for L2 term
lambda2 = input('lambda2 = ');
mu_y = input('mean of data = ');
sigma_y2 = input('variance of data = ');

Q = 0.99;
chi = 0.99;

shusoku = 0;
while shusoku == 0
    chi_old = chi;
    Q_old = Q;
    
    Qh = 1/(1+chi);
    chih = (Q+sigma_y2+mu_y^2)/(1+chi)^2;
    theta = lambda1/sqrt(sqrt(2*chih));
    rhoh = erfc(theta);
    
    chi = rhoh/alpha/(Qh+lambda2);
    Q_tmp = (1+2*theta^2)*rhoh-2*theta*exp(-theta^2)/sqrt(pi);
    Q = chih*Q_tmp/alpha/((Qh+lambda2)^2);
    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA);
end

% training error
err_train = chih;
% degrees of freedom
df = rhoh/alpha-chi*lambda2;
% AT instability
AT = rhoh*Qh^2/alpha/((Qh+lambda2)^2);

