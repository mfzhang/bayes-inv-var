function [ pred ] = opperPredictLatent( meanFunc, covFunc, loghyper,  ...
                                        prior, post, ...
                                        xtrain, xtest )
%OPPERPREDICTLATENT Computes the predictive distribution for the latent
% functions at xstar
% post: the posterior distribution
prior    = opperRecomputePrior(covFunc, loghyper, prior, xtrain);

nTest  = size(xtest,1);
% K      = prior.Sigma;
% IT WAS LIKE THIS BEFORE
%Kinv   = prior.Lambda;
%
Kinv   = post.Lambda - diag(post.eta);

muq     = post.mu;
Sigmaq  = post.Sigma;

muVal     = meanFunc*ones(nTest,1);
[kss Kfs] = feval(covFunc, loghyper, xtrain, xtest);
Ksf = Kfs';

% Predictive Dsitribution
% IT WAS LIKE THIS BEFORE
% muStar    = muVal  + Ksf*Kinv*muq;
%
muStar    = muVal  + Ksf*post.omega;

Cval      = Ksf*Kinv;
SigmaStar =  Cval*Sigmaq*Cval' + diag(kss) - Cval*Kfs;


pred.mu    = muStar;
pred.Sigma = SigmaStar;

return;





