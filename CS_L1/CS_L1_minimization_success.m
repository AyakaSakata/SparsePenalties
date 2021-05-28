%% Check the AT instability of the success solution
%
%  AT > 1 means the instability of the success solution under RS ansatz
clear variables;
close all;

alpha = input('alpha = ');
rho = input('rho = ');
drho = 1e-6;

TOL = 1e-10;
ITMAX = 100;

param_csolve.alpha = alpha;
param_csolve.rho = rho;

% initial condition
chih = 1e-5;

[chih, rc_2] = csolve(@chih_solve, chih, @chih_jacobi, TOL, ITMAX, param_csolve);
AT = ((1-rho)*erfc(1/sqrt(2*chih))+rho)/alpha;

