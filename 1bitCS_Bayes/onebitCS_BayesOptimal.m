%% RS saddle point for 1-bit CS
%
% This code outputs q(SG order parameter) and generalization error 
% as a function of alpha,
% which is the number of the measurements divided by the system size.
% If you want to get them for a particular value of alpha,
% remove the loop with respect to alpha.
%
clear variables;
close all;

% number of measurements divided by the system size
alpha_all = 0.01:0.01:5;
ind_alpha = size(alpha_all,2);

q_all = zeros(ind_alpha,1);
phi_all = q_all;
err_gen = phi_all;
qh_all = q_all;

rho = 0.125;

dx = 0.005;
x = (-10:dx:10)';
x_size = size(x,1);
z2 = ((0:dx:20).^2)';

ker_x = exp(-0.5*x.^2);
ker_z = exp(-0.5*z2);
norm = dx/sqrt(2*pi);

% standard deviation of the signal
sigma = 1;
sigma2 = sigma^2;

% damping factor
dmp = 0.1;
% convergence threshold
THETA = 1e-5;

% initial condition
q = rho*rand(1);

for i_a = 1:ind_alpha

    alpha = alpha_all(i_a);
    shusoku = 0;
    count_loop = 0;
   
    while shusoku == 0
        count_loop = count_loop+1;
        q_old = q;
    
        % for small alpha
        %zeta_0 = x*sqrt(0.5*q/abs(rho-q));
        %erfc_0 = erfc(zeta_0);
        %num = ker_x'*(exp(-2*zeta_0.^2)./erfc_0);
        %qh = 2*alpha/pi/(rho-q)*num*norm;
    
        % for large alpha
        zeta_0 = x*sqrt(0.5*q/rho);
        erfc_0 = erfc(zeta_0);
        num = ker_x'*(exp(-zeta_0.^2)./erfc_0);
        qh = 2*alpha/pi/sqrt(rho*(rho-q))*num*norm;
        
        % for small alpha
        %w_z = exp(0.5*qh*sigma2*z2/(1+qh*sigma2));
        %xi = (1-rho)+rho/sqrt(1+qh*sigma2)*w_z;
        %q_tmp = (rho^2*sigma2)^2*qh/(1+qh*sigma2)^3;
        %q = 2*q_tmp*norm*ker_z'*(z2./xi.*exp(qh*sigma2*z2/(1+qh*sigma2)));
        
        % for large alpha
        w_z = exp(-0.5*qh*sigma2*(z2));
        xi = sqrt(1+qh*sigma2)*(1-rho)*w_z+rho;

        q_tmp = ((rho*sigma2)^2)*qh/(1+qh*sigma2);
        q = 2*q_tmp*ker_z'*(z2./xi)*norm;
    
        shusoku = (abs(q-q_old)<THETA);
    
        q = dmp*(q-q_old)+q_old;
    end
    q_all(i_a) = q;
    qh_all(i_a) = qh;
    phi_all(i_a) = -log(2)*alpha-0.5*q*qh...
        +ker_z'*(xi.*log(xi))*norm...
        +ker_x'*(erfc_0.*log(erfc_0))*norm*alpha;
    W_t = erfc(zeta_0).*erfc(-zeta_0);
    % generalization error
    err_gen(i_a) = ker_x'*W_t*0.5*norm; 
end

all = [alpha_all',q_all,err_gen];
str = sprintf('Random_1bitCS_rho%.2e.dat',rho);
save(str,'all','-ascii');

