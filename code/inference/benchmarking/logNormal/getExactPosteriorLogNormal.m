function [ mu, sigma ] = getExactPosteriorLogNormal( param, y )
%GETEXACTPOSTERIORLOGNORMAL Exact posterior for Gaussian prior and
%logNormal transformation
%   This is as if the likelihood was a lognormal distribution

%% Gets prior parameters
N               = length(y);
nu_p            = param.prior{1};
Lambda_p        = param.prior{2};
lambda_p        = Lambda_p(1,1); % assumes all the same
[mu_p, Sigma_p] = getMeanFromNaturalGauss(nu_p, Lambda_p);

%% Parameters of logNormal
rho = param.fwd{1}; % rho


%% Posterior
fbar = mean(log(y)); % = mean(log(f)) as the noise is v small 
mu     = (mu_p*lambda_p + N*rho*fbar)/(lambda_p + N*rho);
lambda = lambda_p + N*rho;
sigma  = 1/lambda;

return;