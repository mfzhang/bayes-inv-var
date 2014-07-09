function [ param, muq, sigmaPred ] = evbMainGaussFull( y, param, conf )
%EVBMAINGAUSSFULL Summary of this function goes here
%   Detailed explanation goes here

%% Assign parameters from configuratione
niter   = conf.MaxIter;
lrate   = conf.alpha; % lerarning rate
fwdFunc = conf.fwdfunc;
S       = conf.S;

%% Reads prior and likelihood parameters
nup     = param.prior{1};
Lambdap = param.prior{2};
lambday = param.like{1};
[mup, Sigmap, cholLambdap] = getMeanFromNaturalGauss(nup,Lambdap);
Sigmay = 1/lambday;
D  = length(nup);

%% Parmeter initialization
[nuq, Lambdaq cholLq] = initVariationalGaussNatural(D,param.post);
[muq, Sigmaq, cholLq] = getMeanFromNaturalGauss(nuq,Lambdaq);


%% Negative lower bound: -LVAR = - (ELOG - KL)
elogFunc = @(muPost, SigmaPost) evbElogGaussFull(y, fwdFunc, muPost, SigmaPost, Sigmay, S);
klFunc   = @(muPost, SigmaPost) evbKlGaussFull( mup, Sigmap, muPost, SigmaPost);
val = zeros(1,niter);                      
for i = 1 : niter 
    [val(i), gradMuq, gradSigmaq] = evbVarObjectiveGaussFull(elogFunc, klFunc, muq, Sigmaq);
    muq        = muq   - lrate*gradMuq;
    Sigmaq    = Sigmaq - lrate*gradSigmaq;
    plot(i,val(i), '.');
    hold on;
    drawnow;
    i
end

sigmaPred  = diag(Sigmaq);
param.post = {}; % TODO: FILL IN
                      
                      
return;                      



%% Computes variational objective (nelbo)
function [val, gradMuq, gradSigmaq] = evbVarObjectiveGaussFull(elogFunc, klFunc, muq, Sigmaq)
[valElog gradMuElog gradSigmaElog]   = feval(elogFunc, muq, Sigmaq);
[valKl gradMuKl gradSigmaKl]   = feval(klFunc, muq, Sigmaq);
val        = - (valElog - valKl );
gradMuq    = - (gradMuElog - gradMuKl);
gradSigmaq = - (gradSigmaElog - gradSigmaKl);
    

return;    
