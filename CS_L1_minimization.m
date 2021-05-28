%% L1 minimization for compressed sensing
%
% Caution: Quantities in the loop are normalized by "alpha".
% 
clear variables;
close all;

THETA = 1e-10;
count_max = 10000000;
count = 0;

alpha = input('alpha = ');
rho = input('rho = ');
delta = rho/alpha;

lambda = 1;
sigma_x2 = 1;
sigma2 = 0;

shusoku = 0;
dam = 0.01;

% initial condition
%load('mem.mat');
Q = rho/alpha;
m = Q*rand(1);
chi = 100;

while shusoku == 0
    count = count+1;
    
    Q_old = Q;
    m_old = m;
    chi_old = chi;
    
    Qh = 1/chi;
    mh = Qh;
    E = abs(Q-2*m+delta*sigma_x2+sigma2);
    chih = E/chi^2;
    theta_n = lambda/sqrt(2*chih);
    theta_p = lambda/sqrt(2*(chih+(mh^2)*sigma_x2));
        
    chi = ((1-rho)*erfc(theta_n)+rho*erfc(theta_p))/Qh/alpha;
    Q_tmp_p = (1+2*theta_p^2)*erfc(theta_p)-2*theta_p*exp(-theta_p^2)/sqrt(pi);
    Q_tmp_n = (1+2*theta_n^2)*erfc(theta_n)-2*theta_n*exp(-theta_n^2)/sqrt(pi);
    Q = ((1-rho)*chih*Q_tmp_n+rho*(chih+(mh^2)*sigma_x2)*Q_tmp_p)/(Qh^2)/alpha;
    m = rho*sigma_x2*erfc(theta_p)/alpha;
    
    shusoku = (abs(Q-Q_old)<THETA) && (abs(m-m_old)<THETA) && (abs(chi-chi_old)<THETA);
    Q = dam*(Q-Q_old)+Q_old;
    m = dam*(m-m_old)+m_old;
    chi = dam*(chi-chi_old)+chi_old;
    
    if(count>=count_max)
        break;
    end
end

rhoh = (1-rho)*erfc(theta_n)+rho*erfc(theta_p);

str = sprintf('LASSO_alpha%.2f_saddle.dat', alpha);
fp = fopen(str, 'a');
fprintf(fp, '%f %f %f %f %f %f %f\n',...
    rho, Q, m, chi, Qh, mh, chih);
fclose(fp);

save('mem.mat','Q','m','chi');
