%%  generates samples from the true GP
function y = sampleGP(xstar, covfunc, loghyper, MIN_NOISE)
n = size(xstar,1);
Ktilde  = feval(covfunc, loghyper, xstar) + MIN_NOISE*eye(n);
%z =  randn(n,1);
%y = chol(Ktilde)'*z;
mu  = zeros(n,1);
y  = sampleGauss(mu, Ktilde, 1);

return;


