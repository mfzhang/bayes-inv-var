function nelbo = varObjectiveSigmaFixed(mu_q, y, fwdFunc, mu_p, invSigma_p, ...
                                        logdetSigma_p, logdetSigma_q, sigmay)
%VAROBJECTIVESIGMAFIXED  Variational objective for Fixed Sigma
%
N      = size(y,1);
fwdval = feval(fwdFunc, mu_q); % fwd model and Jacobian

% Quadratic term coming from likelihood
lZ      =  - (N/2)*log(2*pi*sigmay);
lquad1  =  - (1/(2*sigmay))*(y-fwdval)'*(y-fwdval);

% Quadratic term from KL
lquad2  =  - 0.5*( (mu_p - mu_q)'*invSigma_p*(mu_p - mu_q) );

% log determinant term from KL
ldet = -0.5*( logdetSigma_p -  logdetSigma_q);

elbo  = lZ + lquad1 + lquad2 + ldet;

nelbo = - elbo;

return;



