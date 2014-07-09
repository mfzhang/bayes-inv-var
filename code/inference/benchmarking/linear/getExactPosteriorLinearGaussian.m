function [mu, Sigma] = getExactPosteriorLinearGaussian(param,y)

N = size(y,1);

%% prior
nu_p            = param.prior.nu;
Lambda_p        = param.prior.Lambda;
[mu_p, Sigma_p] = getMeanFromNaturalGauss(nu_p, Lambda_p);

%% Likelihood
A       = param.fwd{1};
[cols rows]   = size(A);
if (cols> 1 && rows==1) % a vector has been specified
    A = diag(A);
end

lambda_y = param.like.lambda;
sigma_y  = 1/lambda_y;
Sigma_y  = sigma_y*eye(N);
Lambda_y = lambda_y*eye(N); 

%% Posterior mean
K     = A*Sigma_p*A' + Sigma_y;
cholK = chol(K, 'lower');
mu    = mu_p + (A*Sigma_p)'*solveLinChol(cholK,y); 

%% Posterior covariance
Kinv  = getInverseChol(cholK);
Sigma = Sigma_p - (A*Sigma_p)'*Kinv*(A*Sigma_p);
% Sigma  = inv(Sigma_p + A'*Lambda_y*A;)



return;


 