%% Saddle point equations for linear regression penalized by SCAD
%
%  Caution: Macroscopic quantities are normalized by alpha.
%           Rescale them appropriately.
%

clear variables;
close all;

THETA = 1e-10;
sqrt_pi = sqrt(pi);
sq2 = sqrt(2);

%alpha = input('alpha = ');
%rho = input('rho = ');
alpha = 0.5;
rho = 0.1;

% SCAD parameters
lambda = input('lambda = ');
a = input('a = ');

sigma_x2 = 1;
sigma_y2 = 1;

% initial condition
Q = rho;
chi = 30;
m = rho;
% damping factor
dam = 0.1;

    
while shusoku == 0

    Q_old = Q;
    chi_old = chi;
    m_old = m;
    
    Qh = 1/(1+chi);
    chih = (Q-2*m+rho*sigma_x2/alpha+sigma_y2)/(1+chi)^2;
    mh = 1/(1+chi);
    
    std_p = sqrt(chih+mh^2*sigma_x2);
    std_n = sqrt(chih);
    
    var_p = std_p^2;
    var_n = std_n^2;
    
    theta1_p = lambda/sq2/std_p;
    theta2_p = lambda*(Qh+1)/sq2/std_p;
    theta3_p = a*lambda*Qh/sq2/std_p;
    
    theta1_n = lambda/sq2/std_n;
    theta2_n = lambda*(Qh+1)/sq2/std_n;
    theta3_n = a*lambda*Qh/sq2/std_n;
    
    erfc1_p = erfc(theta1_p);
    erfc2_p = erfc(theta2_p);
    erfc3_p = erfc(theta3_p);
    
    erfc1_n = erfc(theta1_n);
    erfc2_n = erfc(theta2_n);
    erfc3_n = erfc(theta3_n);

    e1_p = exp(-theta1_p^2)/sqrt_pi;
    e2_p = exp(-theta2_p^2)/sqrt_pi;
    e3_p = exp(-theta3_p^2)/sqrt_pi;
    
    e1_n = exp(-theta1_n^2)/sqrt_pi;
    e2_n = exp(-theta2_n^2)/sqrt_pi;
    e3_n = exp(-theta3_n^2)/sqrt_pi;

    xi1_p = var_p/Qh*(-2*theta1_p*(e1_p+(Qh-1)*e2_p)+(1+2*theta1_p^2)*(erfc1_p-erfc2_p));
    xi2_p = var_p/(Qh-1/(a-1))*(2*((theta2_p-2*theta3_p/Qh/(a-1))*e2_p-(1-2/Qh/(a-1))*theta3_p*e3_p)...
        +(1+2*(theta3_p/Qh/(a-1))^2)*(erfc2_p-erfc3_p));
    xi3_p = var_p/Qh*(2*theta3_p*e3_p+erfc3_p);
    xi4_p = erfc2_p-erfc3_p;
    
    xi1_n = var_n/Qh*(-2*theta1_n*(e1_n+(Qh-1)*e2_n)+(1+2*theta1_n^2)*(erfc1_n-erfc2_n));
    xi2_n = var_n/(Qh-1/(a-1))*(2*((theta2_n-2*theta3_n/Qh/(a-1))*e2_n-(1-2/Qh/(a-1))*theta3_n*e3_n)...
        +(1+2*(theta3_n/Qh/(a-1))^2)*(erfc2_n-erfc3_n));
    xi3_n = var_n/Qh*(2*theta3_n*e3_n+erfc3_n);
    xi4_n = erfc2_n-erfc3_n;

    xi1 = (1-rho)*xi1_n + rho*xi1_p;
    xi2 = (1-rho)*xi2_n + rho*xi2_p;
    xi3 = (1-rho)*xi3_n + rho*xi3_p;
    xi4 = (1-rho)*xi4_n + rho*xi4_p;

    rhoh = (1-rho)*erfc1_n+rho*erfc1_p;

    Q = (xi1/Qh+xi2/(Qh-1/(a-1))+xi3/Qh)/alpha;
    chi = (rhoh+1/(a-1)/(Qh-1/(a-1))*xi4)/alpha/Qh;
    m = rho*mh*sigma_x2*(erfc1_p+1/(a-1)/(Qh-1/(a-1))*xi4_p)/alpha/Qh;
    
    shusoku = (abs(chi-chi_old)<THETA && abs(Q-Q_old)<THETA && abs(m-m_old)<THETA);

    Q = dam*(Q-Q_old)+Q_old;
    chi = dam*(chi-chi_old)+chi_old;
    m = dam*(m-m_old)+m_old;
    [Q, chi, m]
end

erfc2 = (1-rho)*erfc2_n+rho*erfc2_p;
erfc3 = (1-rho)*erfc3_n+rho*erfc3_p;

xi_all = xi1 + xi2 + xi3 + lambda^2*xi4/(a-1) - (a+1)*lambda^2*erfc3;

% AT instability (if AT > 1, then RS ansatz is unstable)
AT = (rhoh/Qh^2+((Qh-1/(a-1))^(-2)-Qh^(-2))*(erfc2-erfc3))/alpha/((1+chi)^2);

% training error
err_train = chih;
% cross validation error
CVE = Q-2*m+rho/alpha*sigma_x2+sigma_y2;
% generalized degrees of freedom
GDF = chi/(1+chi);
% in sample error
ISE = err_train + 2*sigma_y2*GDF;

