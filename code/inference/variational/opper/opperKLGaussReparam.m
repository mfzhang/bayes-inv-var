function  kl  = opperKLGaussReparam( prior, post )
%OPPERKLGAUSSREPARAM KL Gauss with new reparametrization
%   Detailed explanation goes here

%% Reads prior covariance
K = prior.Sigma;

%% Reads posterior parameters in Oppers parametrization
omegaq  = post.omega;
etaq    = post.eta;
N       = length(omegaq);


%% Useful matrix for inverses
sqrtD   = diag(sqrt(etaq)); 
I       = eye(N);
E       = I + sqrtD*K*sqrtD;
cholE   = getChol(E);
Einv    = getInverseChol(cholE);


%% Trace Term
ltrace = trace(I - sqrtD*Einv*sqrtD*K);

%% Quadratic Term
lquad  = omegaq'*K*omegaq;

%% Determinant term
ldet   = getLogDetChol(cholE);

%% KL divergence
kl    = 0.5*(ltrace + lquad + ldet);


return;


