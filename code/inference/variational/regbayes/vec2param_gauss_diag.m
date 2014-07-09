function param = vec2param_gauss_diag(theta)
% reads vector of parameters into a cell 
% for the variances: theta_i = log sigma_i --> theta_i = 0.5*log(Sigma_ii)
% --> Sigma_ii = exp(2*theta_i)

L          = length(theta);
mu1        = theta(1:L/2);
diagSigma1 = exp(2*theta(L/2+1:L));

param = {mu1, diagSigma1};

return;

