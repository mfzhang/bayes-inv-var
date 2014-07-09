function  [kl invSigmap Lq]  = klGauss( mup, Sigmap, muq, Sigmaq, strDecomp )
%KLGAUSS KL divergence between two Full Gaussians KL(q | p)
% klGauss is useful whe computing a variationa objective
% with p being the prior and q being the posterior
% Note that that the first two arguments correspond to the parameters 
% of the RHS of the KL
% kl  = klGauss( mup, Sigmap, muq, Sigmaq, strDecomp )
% mup: The mean of the RHS disitribution
% Sigmap: The Covariance of the RHS 
% muq: The mean of the the LHS
% Sigmaq: The Covariance of the LHS
% srtDecomp: What decomposition to use 'chol' or 'eigen'

if (nargin == 4) % Cholesky decomposition by default
    strDecomp = 'chol';
end

D = length(mup);

switch lower(strDecomp)
    case 'eigen'
        [Ep Vp]      = getEigenSymmetric(Sigmap);
        [Eq Vq]      = getEigenSymmetric(Sigmaq);
        invSigmap    = getInverseEigenSymmetric(Ep, Vp);
        logdetSigmap = getLogDetEigenSymmetric(Vp);     
        logdetSigmaq = getLogDetEigenSymmetric(Vq);
    case 'chol'
        Lp           = getChol(Sigmap);
        Lq           = getChol(Sigmaq);
        invSigmap    = getInverseChol(Lp);
        logdetSigmap = getLogDetChol(Lp);
        logdetSigmaq = getLogDetChol(Lq);
    otherwise
        dips('Unkown decomposition');
end
%% trace term
ltr = traceProd(invSigmap,Sigmaq);

%% Quadratic term
lquad = (mup - muq)'*invSigmap*(mup - muq);

%% Log determinant term log | Sp^-1 Sq|
ldet = -  logdetSigmap +  logdetSigmaq;


%% final result
kl = 0.5*(ltr + lquad - ldet - D);


return;

