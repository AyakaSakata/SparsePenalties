%% Saddle point equations for linear regression penalized by SCAD
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

alpha = input('alpha = ');

% SCAD parameters
lambda = input('lambda = ');
a = input('a = ');

% mean of data
mu_y = 0;
% variance of data
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
    theta2 = lambda*(Qh+1)/sqrt(2*chih);
    theta3 = a*lambda*Qh/sqrt(2*chih);
    
    erfc1 = erfc(theta1);
    erfc2 = erfc(theta2);
    erfc3 = erfc(theta3);
    
    e1 = exp(-theta1^2)/sqrt_pi;
    e2 = exp(-theta2^2)/sqrt_pi;
    e3 = exp(-theta3^2)/sqrt_pi;
    
    xi1 = chih/Qh*(-2*theta1*(e1+(Qh-1)*e2)+(1+2*theta1^2)*(erfc1-erfc2));
    xi2 = chih/(Qh-1/(a-1))*(2*((theta2-2*theta3/Qh/(a-1))*e2-(1-2/Qh/(a-1))*theta3*e3)...
        +(1+2*(theta3/Qh/(a-1))^2)*(erfc2-erfc3));
    xi3 = chih/Qh*(2*theta3*e3+erfc3);
    xi4 = erfc2-erfc3;
    
    Q = (xi1/Qh+xi2/(Qh-1/(a-1))+xi3/Qh)/alpha;
    chi = (1/(a-1)/(Qh-1/(a-1))*xi4+erfc1)/alpha/Qh;

    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA);

    Q = dam*(Q-Q_old)+Q_old;
    chi = dam*(chi-chi_old)+chi_old;
    [Q, chi]
end

% effective density of non-zero components
rhoh = erfc1;
% traning error
err_train = chih;

xi_all = xi1 + xi2 + xi3 + lambda^2*xi4/(a-1) - (a+1)*lambda^2*erfc3;
% expectation of regularization term
r_ave = -0.5*(alpha*(Q*Qh-chi*chih)+xi_all-alpha*chi*chih);
% free energy
f = 0.5*alpha*(Q+sigma_y2+mu_y^2)/(1+chi)-0.5*alpha*(Q*Qh-chi*chih)-0.5*xi_all;

sta1 = (rhoh+(1/(Qh*(a-1)))/(1-1/(Qh*(a-1)))*xi4)/alpha;
sta2 = (rhoh-xi4)/alpha;

% degrees of freedom
GDF = chi/(1+chi);
% CV error
CVE = Q+sigma_y2+mu_y^2;
% AT instability (if AT > 1, then RS ansatz is unstable)
AT = (rhoh/Qh^2+((Qh-1/(a-1))^(-2)-Qh^(-2))*(erfc2-erfc3))/alpha/((1+chi)^2);
