%% Linear regression penalized by MCP
%
%  Caution: Macroscopic quantities are normalized by alpha.
%           Rescale them appropriately.
%
clear variables;
close all;

THETA = 1e-10;
sqrt_pi = sqrt(pi);
TOL = 1e-10;
ITMAX = 100;

sgn = 0;

alpha = input('alpha = ');

lambda = input('lambda = ');
a = input('a = ');

mu_y = input('mu_y = ');
sigma_y2 = 1;

% initial condition
Q = 1;
chi = 100;

% damping factor
dam = 0.1;

shusoku = 0;
while shusoku == 0
    Q_old = Q;
    chi_old = chi;
    
    Qh = 1/(1+chi);
    chih = (Q+sigma_y2+mu_y^2)/(1+chi)^2;
    
    theta1 = lambda/sqrt(2*chih);
    theta2 = a*lambda*Qh/sqrt(2*chih);
    
    erfc1 = erfc(theta1);
    erfc2 = erfc(theta2);
    
    e1 = exp(-theta1^2)/sqrt_pi;
    e2 = exp(-theta2^2)/sqrt_pi;
    
    rhoh = erfc1;

    pi1 = (chih+lambda^2)/(Qh-1/a)*(erfc1-erfc2)-sqrt(2*chih)*lambda/(Qh-1/a)*(e1-e2+(a*Qh-1)*e2);
    pi2 = chih/Qh*(2*theta2*e2+erfc2)-a*lambda^2*erfc2;
    
    Q = pi1/(Qh-1/a)/alpha+chih/(Qh^2)/alpha*(2*theta2*e2+erfc2);
    chi = (rhoh+1/a/(Qh-1/a)*(erfc1-erfc2))/alpha/Qh;

    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA);

    Q = dam*(Q-Q_old)+Q_old;
    chi = dam*(chi-chi_old)+chi_old;
    [Q, chi]
end
% training error
err_train = chih;
pi_all = pi1 + pi2;

% expectation value of the penalty
r_ave = alpha*(chih*chi-0.5*Q*Qh)-0.5*pi_all;
% free energy
f = 0.5*alpha*(Q+sigma_y2+mu_y^2)/(1+chi)-0.5*alpha*(Q*Qh-chi*chih)-0.5*pi_all;

% generalized degrees of freedom
df = sigma_y2*chi/(1+chi);

% AT instability (if > 1, the RS ansatz is unstable)
AT = (rhoh+((Qh/(Qh-1/a))^2-1)*(erfc1-erfc2))/alpha;
