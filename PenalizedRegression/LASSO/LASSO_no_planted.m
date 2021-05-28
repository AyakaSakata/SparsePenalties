clear all;
close all;

THETA = 1e-10;

alpha = input('alpha = ');
lambda = input('lambda = ');
mu_y = 0;
sigma_y2 = 1;

% “Á‚Élambda‚ª¬‚³‚¢‚É
% ‰ŠúğŒ‚Ì’²ß‚ª•K—v
Q = 1;
chi = 200;  

dam = 0.01; % lambda‚ª¬‚³‚¢‚ÍA¬‚³‚­‚·‚é

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

err_train = chih;
l1 = chi*chih-Q*Qh;

df = rhoh/alpha;
cvo = Q+sigma_y2+mu_y^2-err_train;

str = sprintf('l1_alpha%.2f_mu%.2f_var%.2f.dat', alpha, mu_y, sigma_y2);
fp = fopen(str, 'a');
fprintf(fp, '%f %f %f %f %f %f %f %f\n', lambda, 2*df*sigma_y2, cvo, rhoh, err_train, l1, Q, chi);
fclose(fp);
