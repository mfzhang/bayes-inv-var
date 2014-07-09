function  kl  = klGaussStructChol( prior, post )
%KLGAUSSSTRUCT KL Gauss when the parameters are all structures
% Assumes Cholesky factorization

D = length(prior.mu);


%% trace term
% ltr = traceProd(invSigmap,Sigmaq);
ltr = traceProd(prior.Lambda,post.Sigma);

%% Quadratic term
% lquad = (mup - muq)'*invSigmap*(mup - muq);
v = prior.mu - post.mu;
lquad = v'*(prior.Lambda)*v;

%% Log determinant term log | Sp^-1 Sq|
% ldet = -  logdetSigmap +  logdetSigmaq;
logdetSigmap =  getLogDetChol(prior.cholSigma);
logdetSigmaq = getLogDetChol(post.cholSigma);
ldet         = -  logdetSigmap +  logdetSigmaq;


%% final result
kl = 0.5*(ltr + lquad - ldet - D);



return;
