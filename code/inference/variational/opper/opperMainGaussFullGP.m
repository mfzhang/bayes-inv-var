function [ param, muq, sigmaPred, nelbo ] = opperMainGaussFullGP(x, y, param, conf)
%OOPERMAINGAUSSFULLGP Now using a GP prior
%   Detailed explanation goes here
% For now we assume zero-mean GP prior
% Simple support for constant mean function
N         = size(y,1);
covFunc   = param.covFunc;
meanFunc  = param.meanFunc;
mup       = meanFunc*ones(N,1); % prior mean


% for hyperparameter learning
param.loghyper      = opperInitHyper(covFunc, x, y);



%% We alternate between optimizing hyper and the posterior parmaeters
for i = 1 : 10%conf.MaxIter
    
    %% Recompute prior depending on hyper-parameters
    K                     = feval(covFunc, param.loghyper, x);   
    cholK                 = getCholSafe(K);
    Kinv                  = getInverseChol(cholK); 
    param.prior.Lambda    = Kinv;       % precision
    param.prior.nu        = Kinv*mup; % precison-adjusted mean
    param.prior.Sigma     = K;
    param.prior.mu        = mup;
    param.prior.cholSigma = cholK;

    %% Optimizes posterior parameters given fixed hyper-parametes    
    [param, muq, sigmaPred, nelbo]  = opperMainGaussFull( y, param, conf );
    fprintf('\n\n *** Posterior Parameters Optimization Done --> ');
    fprintf('NELBO=%.6f ***\n\n', nelbo); 
    
%     %% Optimizers hyper-parameters given posterior parameters
     param = opperOptimizeHyper(param, x, conf);
     fprintf('\n\n*** Optimization of all Hyper-parameters Done ***\n\n');    
    
end

%% Final parameter values
[ param, muq, sigmaPred, nelbo ] = opperMainGaussFull( y, param, conf );



return;



