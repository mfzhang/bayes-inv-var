function prior = opperRecomputePrior(covFunc, loghyper, prior, x)
K                     = feval(covFunc, loghyper, x);   
cholK                 = getCholSafe(K);
Kinv                  = getInverseChol(cholK); 
prior.Lambda    = Kinv;       % precision
prior.nu        = Kinv*prior.mu; % precison-adjusted mean
prior.Sigma     = K;
prior.cholSigma = cholK;

return;

