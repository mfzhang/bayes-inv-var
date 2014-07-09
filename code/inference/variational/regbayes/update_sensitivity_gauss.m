function A = update_sensitivity_gauss(y, sigmay, lambda_reg, fwdfun, fwdparam, ...
    mu, Sigma)
% Update sensitiviy matrix given parameters of the posterior
% mu:  approximate posterior mean
% Sigma: approximate posterior covariance
fwd_val = feval(fwdfun, fwdparam{:}, mu);


% modified on 28/10/2013
C       = Sigma + mu*mu' + 2*lambda_reg*sigmay*(mu*mu');
A = (y*mu' + 2*lambda_reg*sigmay*fwd_val*mu')/C; % (.)*C^-1
% Convex formulation
%C       = (1-lambda_reg)*(Sigma + mu*mu') + 2*lambda_reg*sigmay*(mu*mu');
%A       = ((1-lambda_reg)*y*mu' + 2*lambda_reg*sigmay*fwd_val*mu')/C; % (.)*C^-1
return;

