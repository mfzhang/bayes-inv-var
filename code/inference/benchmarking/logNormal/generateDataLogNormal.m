function [ y, param ] = generateDataLogNormal( N )
%GENERATEDATALOGNORMAL Generates data for normal prior
% and log normal fwd model
% In this simple case N = M
%
% from M=1 parameter to N observations
M = 1;
if (nargin == 0)
    N = 1000;    
end

%% Prior parameters
mu_p    = zeros(M,1);
sigma_p = 1;
Sigma_p = sigma_p*eye(M);

%% Parameters of fwd model
fwdParam{1} = 1000; % rho: Precision of Gaussian to pass through exp
%  N Number of observations
z           = randn(N,1); % Needs to keep this fixed so that the fwd
                          % model does not change randomlya
fwdParam{2} = z; % sample Gaussian to pass through exp

%% Draws one sample from N-dimensional prior
theta   = sampleGauss(mu_p,Sigma_p,1); 

%% Parameters of conditional likelihood: "identity function"
% used only by inference algorithm
sigma_y  = 1e-10; % very small variance
lambda_y = 1/sigma_y;
Sigma_y  = sigma_y*eye(N);


%% Fwd Model Parameters
f = logNormalFwdModel(theta, fwdParam{:});

%% We don't sample from conditional likelihood as it is the identity
y = f;
%hist(log(y)); title('log(y)');

%% Get Natural Parameters
[ nu_p, Lambda_p, cholSigmap] = getMeanFromNaturalGauss( mu_p, Sigma_p );

param.prior{1} = nu_p;
param.prior{2} = Lambda_p;
param.like{1}  = lambda_y;
param.fwd      = fwdParam;

return;





    