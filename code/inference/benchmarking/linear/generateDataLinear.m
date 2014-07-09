function [y param] = generateDataLinear( M, N, diagFlag )
%GENERATEDATA Generate Data for Gaussian Prior, Linear Fwd Models
% M: Number of parameters
% N: Number of observations

%% Assigns number of parameters and observations
if (nargin == 0)
    M       = 10; 
    N       = 100; 
    diagFlag = 0;
end

%% Prior parametes 
mu_p    = zeros(M,1);
Sigma_p = eye(M);

%% likelihood Parameters
fwdFunc  = @linearFwdModel;
fwdParam =  rand(N,M);
if (diagFlag)
    fwdParam = diag(diag(fwdParam));
end
%% Parameters of conditional likelihood
sigma_y = 1e-5;
lambda_y = 1/sigma_y;
Sigma_y = sigma_y*eye(N);


%% Draws one sample from M-dimensional prior
theta   = sampleGauss(mu_p,Sigma_p,1); 


%% Draws from conditional likelihhood
f  = feval(fwdFunc,theta, fwdParam);
y  = sampleGauss(f, Sigma_y, 1);
plot(y);

%% Get Natural Parameters
[ nu_p, Lambda_p, cholSigma_p] = getMeanFromNaturalGauss( mu_p, Sigma_p );


param.prior{1} = nu_p;
param.prior{2} = Lambda_p;
param.like{1}  = lambda_y;
param.fwd{1}   = fwdParam;
 
return;

