function f = logNormalFwdModel(theta, rho, z)
% rho: precision of Gaussian to be passed through exponential funct.
sigma  = (1/rho);
x      = theta + sqrt(sigma)*z;
f      = exp(x);

return;