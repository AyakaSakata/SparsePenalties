%% Bayes optimal compressed sensing
%% with Bernoulli-Gaussian signal
%%
clear variables;
close all;

% Parameters
alpha = 0.5;
rho_0 = 0.32;
rho = rho_0;
sigma_x = 1;
sigma_x2 = sigma_x^2;

% threshold for convergence
THETA = 1e-8;
% maximum step for iteration
count_max = 1000;
% damping factor
dmp = 0.1;

% variables for integral
z_max = 10;
dz = 0.01;
z = (0:dz:z_max)';
ker_z = exp(-0.5*z.^2);
norm_int = dz/sqrt(2*pi)*2;

% initial condition
shusoku = 0;

Q = rho;
m = Q*rand(1);
d_Q = m*rand(1);
q = Q-d_Q;

while shusoku == 0 
    
    Q_old = Q;
    d_Q_old = d_Q;
    m_old = m;
    
    q = Q-d_Q;
    
    mh = alpha/d_Q;
    qh = alpha*(q-2*m+rho_0*sigma_x2)/(d_Q)^2;
    Qh = mh-qh;

    var_x = Qh+qh+1/sigma_x2;

    sigma_p = sqrt(qh+mh^2*sigma_x2);
    sigma_n = sqrt(qh);
    
    weight_p = exp(-0.5*(sigma_p*z).^2/var_x);
    weight_n = exp(-0.5*(sigma_n*z).^2/var_x);
    
    den_p = (1-rho)*sigma_x*sqrt(var_x)*weight_p+rho;
    den_n = (1-rho)*sigma_x*sqrt(var_x)*weight_n+rho;
    
    m_tmp = ker_z'*(z.^2./den_p)*norm_int;
    m = rho_0*rho*mh*sigma_x2*m_tmp/var_x;
    
    d_Q_tmp_n = ker_z'*(z.^2./den_n)*norm_int;
    d_Q = rho*((1-rho_0)*d_Q_tmp_n+rho_0*m_tmp)/var_x;
    
    Q_tmp_p = ker_z'*(((sigma_p*z).^2+var_x)./den_p)*norm_int;
    Q_tmp_n = ker_z'*(((sigma_n*z).^2+var_x)./den_n)*norm_int;
    Q = rho*((1-rho_0)*Q_tmp_n + rho_0*Q_tmp_p)/(var_x^2);
    
    shusoku = (abs(Q_old-Q)<THETA && abs(m-m_old)<THETA && abs(d_Q-d_Q_old)<THETA);

    % damping
    Q = dmp*(Q-Q_old)+Q_old;
    m = dmp*(m-m_old)+m_old;
    d_Q = dmp*(d_Q-d_Q_old)+d_Q_old;

end
% MSE between true and reconstructed signals
D = Q-2*m+rho_0*sigma_x2;

