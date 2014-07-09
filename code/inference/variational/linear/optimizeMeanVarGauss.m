function [ mu_q nelbo] = optimizeMeanVarGauss(y, fwdFunc, mu_p, L_p, mu_q, L_q, lambday)
%OPTIMIZEMEANVARGAUSS Optimizes the menan of the variational Gauss
%   Detailed explanation goes here

sigmay = 1/lambday;

optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = 100;  
optconf.MaxFunEvals = 100000; 
optconf.progTol     = 0;
optconf.numDiff     = 1; % 0: user-supplied; 1: fwd-difference, 2: central-difference


invSigma_p    = getInverseChol(L_p);
logdetSigma_p = getLogDetChol(L_p);
logdetSigma_q = getLogDetChol(L_q);

        
        
 fobj = @(mu) varObjectiveSigmaFixed(mu, y, fwdFunc, mu_p, invSigma_p, ...
                                    logdetSigma_p, logdetSigma_q, sigmay);
 [mu_q, nelbo] = minFunc(fobj, mu_q, optconf);

return;


