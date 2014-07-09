function [ param, muq, sigmaPred, nelbo ] = opperMainGaussFullGPReparam(x, y, param, conf)
%OOPERMAINGAUSSFULLGP Now using a GP prior
%   Detailed explanation goes here
% For now we assume zero-mean GP prior
% Simple support for constant mean function
N         = size(y,1);
covFunc   = param.covFunc;
meanFunc  = param.meanFunc;
param.prior.mu= meanFunc*ones(N,1); % prior mean


% for hyperparameter learning
param.loghyper      = opperInitHyper(covFunc, x, y);

%% Recompute prior depending on hyper-parameters
param.prior  = opperRecomputePrior(covFunc, param.loghyper, param.prior, x);

%% THIS DIDNT WORK
% param.post.mu     = param.prior.mu;
% param.post.Lambda = param.prior.Lambda; 
% param.post.Sigma  = param.prior.Sigma; 
% param.post.omega  = param.prior.Lambda*param.post.mu;
% param.post.eta    = diag(param.post.Lambda - param.prior.Lambda);
% param             = opperOptimizeHyper(param, x, y, conf);


%% We alternate between optimizing hyper and the posterior parmaeters
for i = 1 : 1 %conf.MaxIter

    %% Optimizes posterior parameters given fixed hyper-parametes    
    [param, muq, sigmaPred, nelbo]  = opperMainGaussFullReparam( y, param, conf );
    fprintf('\n\n *** Posterior Parameters Optimization Done --> ');
    fprintf('NELBO=%.6f ***\n\n', nelbo); 
    
%      %% Optimizers hyper-parameters given posterior parameters
     param = opperOptimizeHyper(param, x, y, conf);
     fprintf('\n\n*** Optimization of all Hyper-parameters Done ***\n\n');    

     %% Recompute prior depending on hyper-parameters
    param.prior  = opperRecomputePrior(covFunc, param.loghyper, param.prior, x);
end


%% Final parameter values
[ param, muq, sigmaPred, nelbo ] = opperMainGaussFullReparam( y, param, conf );


return;






return;



    