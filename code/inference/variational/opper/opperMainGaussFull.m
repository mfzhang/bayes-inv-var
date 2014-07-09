function [ param, muq, sigmaPred, nelbo ] = opperMainGaussFull( y, param, conf )
%OPPERMAINGAUSSFULL Variational inference of Opper and Archembau
%   (2009, The Variational Gaussian Approx. revisited, Neural Computation)

%% Read firm configuration
S          = conf.S;
fwdFunc    = @(xx) param.fwdFunc(xx, param.fwd{:});
elogMethod = conf.elogMethod;

%% Reads prior and likelihood parameters
N       = size(y,1);
nup     = param.prior.nu; 
Lambdap = param.prior.Lambda;
mup     = param.prior.mu;
Sigmap  = param.prior.Sigma;
%
%
lambday    = param.like.lambda;
sigmay     = 1/lambday;
diagSigmay = sigmay*ones(N,1);
 
 %% Parameter initialization
[omega0 eta0]         = opperInitParameters(Lambdap, param.post);
theta0                = opperPack( omega0, eta0 );

%% Optimization configuration
optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = conf.MaxIter;  
optconf.MaxFunEvals = 1000; 
optconf.progTol     = 0;
optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION



%% Optimizes nelbo
fobj = @(theta) opperNelboGaussFull( theta,  param.prior, ...
                                    fwdFunc, diagSigmay, y, S, elogMethod);
[theta, nelbo] = minFunc(fobj, theta0, optconf);


%% Gets mean parameters
[ omegaq, eta ] = opperUnpack( theta ); %% this omegaq is Opper parameterization 
[ muq, Sigmaq Lambdaq ] = opperGetMeanParam( omegaq, eta, Sigmap, Lambdap );
sigmaPred       = diag(Sigmaq);


param.post.nu     = Lambdaq*muq; 
param.post.Lambda = Lambdaq;
param.post.mu     = muq;
param.post.Sigma  = Sigmaq;

return;





 