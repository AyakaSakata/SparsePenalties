%% Calculation of the complexity 
% 
% See http://dx.doi.org/10.1088/1742-5468/2016/06/063302
% or proceedings of ISIT 2014 
% "Statistical Mechanical Analysis of Overcomplete Basis Compression"
% (Japanese)
%
clear variables;
close all;

%alpha = input('alpha = ');
alpha = 0.8;
R = input('rate = ');
mu = input('mu = ');

TOL = 1e-10;
ITMAX = 100;

THETA = 1e-8;

% range of the integral
int_max = 10;
% bin of the integral
d_int = 1e-4;
% normalization constant of the integral
norm_int = d_int/sqrt(2*pi);

sigma_y = 1;
sigma_y2 = sigma_y^2;

d2 = (-int_max:d_int:int_max).^2;
% Gaussian weight
ker = exp(-0.5*d2);

param_csolve.alpha = alpha;
param_csolve.int_max = int_max;
param_csolve.d_int = d_int;
param_csolve.R = R;

% initial condition
% load 'tmp.mat' when you want to start from the fixed point of the
% previous run
%load 'tmp.mat';
Q = rand(1);
chi = rand(1);
q0 = rand(1);
Rh = rand(1);

shusoku = 0;

while shusoku == 0
    Q_old = Q;
    q0_old = q0;
    chi_old = chi;
    
    chih = mu^2*(Q-q0)/(1+chi)/(1+chi+mu*(Q-q0));
    q0h = mu^2*(q0+sigma_y2)/(1+chi+mu*(Q-q0))^2;
    Qh = mu/(1+chi);

    param_csolve.chih = chih;
    param_csolve.q0h = q0h;
    param_csolve.Qh = Qh;

    [Rh, rc] = csolve(@Rh_solve, Rh, @Rh_jacobi, TOL, ITMAX, param_csolve);

    W = 1+sqrt(Qh/(Qh-chih))*exp(-Rh+0.5*q0h*d2/(Qh-chih));
    weight = (W-1)./W;
    Q = R*chih/Qh/(Qh-chih)+q0h/(Qh-chih)^2*(weight.*d2)*ker'*norm_int/alpha;
    q0 = R/(Qh-chih)+(q0h/(Qh-chih)^2-1/(Qh-chih))*(weight.*d2)*ker'*norm_int/alpha;
    chi = mu*R/Qh;
    
    shusoku = (abs(Q_old-Q)<THETA)*(abs(q0_old-q0)<THETA)*(abs(chi_old-chi)<THETA);
    
    [Q - Q_old, q0 - q0_old, chi - chi_old]
    
    Q = rand(1)*(Q-Q_old)+Q_old;
    q0 = rand(1)*(q0-q0_old)+q0_old;
    chi = rand(1)*(chi-chi_old)+chi_old;
end

den = 1+chi+mu*(Q-q0);
% Generating function
phi = 0.5*alpha*(log((1+chi)/den)-mu*(q0+sigma_y2)/den+Q*Qh...
    -(Q+chi/mu)*(chih+q0h)+q0*q0h)+alpha*R*Rh+log(W)*ker'*norm_int;
% Distortion
D = -(-(Q-q0)/den-(q0+sigma_y2)/den+mu*(q0+sigma_y2)*(Q-q0)/(den^2)+chi*(chih+q0h)/mu^2);
% Complexity
S = phi+0.5*alpha*mu*D;

str = sprintf('OC_alpha%.2f_R%.2f.dat', alpha, R);
fp = fopen(str, 'a');
fprintf(fp, '%f %f %f %f %f %f\n', mu, D, S, Q, chi, q0);
save('tmp.mat','Q','chi','q0','Rh');
fclose(fp);
