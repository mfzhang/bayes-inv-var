function param  = opperOptimizeHyper( param, x, y, conf )
%OPPEROPTIMIZEHYPER Optimization of Hyper-parametes for O&A
%   Detailed explanation goes here

loghyper0 = param.loghyper;

%% Optimization configuration
optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = conf.MaxIter;  
optconf.MaxFunEvals = 500; 
optconf.progTol     = 1e-6;
optconf.optTol      = 1e-10;

optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION

meanFunc = param.meanFunc;
covFunc  = param.covFunc;


%fobj = @(logtheta) opperKLGaussGPReparam(meanFunc, covFunc, logtheta, x ,param.post);
%loghyper = minFunc(fobj, loghyper0, optconf);  
N = size(y,1);
S          = conf.S;
fwdFunc    = @(xx) param.fwdFunc(xx, param.fwd{:});
elogMethod = conf.elogMethod;
lambday    = param.like.lambda;
sigmay     = 1/lambday;
diagSigmay = sigmay*ones(N,1);

fobj = @(logtheta) opperNelboGaussFullGP(logtheta, meanFunc, covFunc, fwdFunc, ...
                                            param.prior, param.post, diagSigmay, ...
                                                    x, y, S, elogMethod);
loghyper = minFunc(fobj, loghyper0, optconf);  
                                                
param.loghyper = loghyper;

return;





