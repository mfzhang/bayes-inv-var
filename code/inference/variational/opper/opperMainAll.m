function [ param, muq, sigmaPred, nelbo  ] = opperMainAll( x, y, param, conf )
%OPPERMAINALL Learns everything
%   Detailed explanation goes here
N              = size(x,1);
covFunc        = param.covFunc;
meanFunc       = param.meanFunc;
param.prior.mu = meanFunc*ones(N,1); % prior mean



%% Read firm configuration
S          = conf.S;
fwdFunc    = @(xx) param.fwdFunc(xx, param.fwd{:});
elogMethod = conf.elogMethod;

%% reads likelihood
lambday    = param.like.lambda;
sigmay     = 1/lambday;
diagSigmay = sigmay*ones(N,1);


%% Optimization configuration
optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = conf.MaxIter;  
optconf.MaxFunEvals = 1000; 
optconf.progTol     = 1e-6;
optconf.optTol      = 1e-6;
optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION


[ omega0, eta0 ] = opperInitParametersReparam({}, N );
thetaVar         = opperPack( omega0, eta0 );
loghyper         = opperInitHyper(covFunc, x, y);
theta0           = [thetaVar; loghyper];

%% Indicees of parameters for optimization
idxOmega   = 1:N;
idxEta     = N+1:2*N;
idxVar     = [idxOmega, idxEta];
idxHyper   = 2*N+1:length(theta0);

%% Here we optimize wrt everything
% fobj = @(theta) opperNelboAll( theta, meanFunc, covFunc, fwdFunc, ...
%                                           diagSigmay, x, y, S, elogMethod);                                      
% [theta, nelbo]          = minFunc(fobj, theta0, optconf);
% [post.omega,post.eta ]  = opperUnpack( theta(idxVar) );
% param.loghyper          = theta(idxHyper);                                      

iter = 0;
diffHyper  = Inf;
loghyperOld = log(loghyper);
nelboOld    = -Inf;
diffNelbo    = Inf;
while ( (iter < 100) && (diffHyper > 1e-3) && (diffNelbo > 1e-3) )
    %% here we optimize wrt variational parameters only
    optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
    optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION
    fobj = @(xx) wrapNelboSubset([xx;loghyper], meanFunc, covFunc, fwdFunc, ...
                                           diagSigmay, x, y, S, elogMethod, idxVar);                                      
    [thetaVar, nelbo]          = minFunc(fobj, thetaVar, optconf);
    fprintf('\n iter=%d *** Optimization of variational parameters done *** \n',iter);

    %% here we optimize wrt hyper-parameters using numerical optimization
    optconf.numDiff     = 1; % 0: user-supplied; 1: fwd-difference, 2: central-difference
    optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION
    fobj = @(xx) wrapNelboSubset([thetaVar;xx], meanFunc, covFunc, fwdFunc, ...
                                           diagSigmay, x, y, S, elogMethod, idxHyper);                                      
    [loghyper, nelboHyper]          = minFunc(fobj, loghyper, optconf);
    fprintf('\n iter=%d *** Optimization of hyperparameter  done *** \n', iter);
    
    iter        = iter + 1;
    diffHyper   = norm(loghyper - loghyperOld);
    diffNelbo   = abs(nelbo - nelboOld);
    loghyperOld = loghyper;
end
%% One Last update of variational parametes
optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION
fobj = @(xx) wrapNelboSubset([xx;loghyper], meanFunc, covFunc, fwdFunc, ...
                                 diagSigmay, x, y, S, elogMethod, idxVar);                                      
[thetaVar, nelbo]          = minFunc(fobj, thetaVar, optconf);
fprintf('\n iter=%d *** Final optimization of variational parameters done *** \n',iter);




%  Assing results from optimization
param.loghyper          = loghyper;
[post.omega,post.eta ]  = opperUnpack( thetaVar );


%% Return stuff
prior = opperRecomputePrior(covFunc, param.loghyper, param.prior, x);
[ muq, Sigmaq Lambdaq ] = opperGetMeanParamStruct( prior, post);
sigmaPred         = diag(Sigmaq);
post.nu     = Lambdaq*muq; 
post.Lambda = Lambdaq;
post.mu     = muq;
post.Sigma  = Sigmaq;
param.post  = post;

return;


  
%%
function[ nelbo gradTheta ] = wrapNelboSubset(theta, meanFunc, covFunc, fwdFunc, ...
                                          diagSigmay, x, y, S, elogMethod,idxOpt)
if (nargout == 1)
  nelbo = opperNelboAll( theta, meanFunc, covFunc, fwdFunc, ...
                                          diagSigmay, x, y, S, elogMethod); 
else
 [ nelbo gradTheta ] = opperNelboAll( theta, meanFunc, covFunc, fwdFunc, ...
                                          diagSigmay, x, y, S, elogMethod);    
gradTheta = gradTheta(idxOpt);
end

return;              