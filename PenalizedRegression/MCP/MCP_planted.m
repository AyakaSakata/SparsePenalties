%% Linear regression penalized by MCP with planted solution
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
rho = input('rho = ');
% variance of the planted signal
sigma_x2 = 0.5;
%  variance of the noise
sigma_y2 = 0.8;

% MCP parameters
lambda = input('lambda = ');
a = input('a = ');

% initial condition
Q = 1;
chi = 1;
m = 0.1*rho;

% damping factor
dam = 0.1;

shusoku = 0;
while shusoku == 0
    Q_old = Q;
    m_old = m;
    chi_old = chi;
    
    Qh = 1/(1+chi);
    mh = 1/(1+chi);
    chih = (Q-2*m+rho*sigma_x2/alpha+sigma_y2)/(1+chi)^2;

    sigma_n = sqrt(chih);
    sigma_p = sqrt(chih+mh^2*sigma_x2);
    
    theta1_n = lambda/sqrt(2)/sigma_n;
    theta2_n = a*lambda*Qh/sqrt(2)/sigma_n;
    
    theta1_p = lambda/sqrt(2)/sigma_p;
    theta2_p = a*lambda*Qh/sqrt(2)/sigma_p;

    erfc1_n = erfc(theta1_n);
    erfc2_n = erfc(theta2_n);
    
    erfc1_p = erfc(theta1_p);
    erfc2_p = erfc(theta2_p);
    
    e1_n = exp(-theta1_n^2)/sqrt_pi;
    e2_n = exp(-theta2_n^2)/sqrt_pi;
    
    e1_p = exp(-theta1_p^2)/sqrt_pi;
    e2_p = exp(-theta2_p^2)/sqrt_pi;

    rhoh = (1-rho)*erfc1_n + rho*erfc1_p;
    
    xi1_n = (sigma_n^2+lambda^2)/(Qh-1/a)*(erfc1_n-erfc2_n)...
        -2*sigma_n^2*theta1_n/(Qh-1/a)*(e1_n-e2_n+(a*Qh-1)*e2_n);
    xi2_n = sigma_n^2/Qh*(2*theta2_n*e2_n+erfc2_n)-a*lambda^2*erfc2_n;
    xi3_n = erfc1_n-erfc2_n;
    
    xi1_p = (sigma_p^2+lambda^2)/(Qh-1/a)*(erfc1_p-erfc2_p)...
        -2*sigma_p^2*theta1_p/(Qh-1/a)*(e1_p-e2_p+(a*Qh-1)*e2_p);
    xi2_p = sigma_p^2/Qh*(2*theta2_p*e2_p+erfc2_p)-a*lambda^2*erfc2_p;
    xi3_p = erfc1_p-erfc2_p;
    
    Q_n = xi1_n/(Qh-1/a)/alpha+sigma_n^2/(Qh^2)/alpha*(2*theta2_n*e2_n+erfc2_n);
    Q_p = xi1_p/(Qh-1/a)/alpha+sigma_p^2/(Qh^2)/alpha*(2*theta2_p*e2_p+erfc2_p);
    Q = (1-rho)*Q_n+rho*Q_p;
    
    xi3_ave = (1-rho)*xi3_n+rho*xi3_p;
    chi = (rhoh + xi3_ave/a/(Qh-1/a))/alpha/Qh;
    
    m = rho*sigma_x2/alpha*(erfc1_p+xi3_p/a/(Qh-1/a));

    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA && abs(m-m_old)<THETA);

    Q = dam*(Q-Q_old)+Q_old;
    chi = dam*(chi-chi_old)+chi_old;
    m = dam*(m-m_old)+m_old;
    [Q, chi, m]
end
% training error
err_train = chih;
% generalized degrees of freedom
GDF = chi/(1+chi);
% cross validation error
CVE = Q-2*m+rho*sigma_x2/alpha+sigma_y2;
% AT instability (if > 1, then the RS ansatz is unstable.)
AT = (rhoh+((Qh/(Qh-1/a))^2-1)*xi3_ave)/alpha;
