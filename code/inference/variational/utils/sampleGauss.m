function X = sampleGauss(mu, Sigma, N)
% x = gauss_sample(mu, Sigma, N)
% Generate N values from a multivariate Gaussian distribution with specified
% mean vector mu and covariance matrix Sigma.
%
% eg.  N = 1000; %Number of samples to generate
% mu = [1;2]; % This is the mean vector that you saved for each problem
% Sigma = [1 0.5; 0.5 2]; %Covariance matrix you saved for each problem
% mu,x: column vectors
% Each column of X is a different sample
D = length(mu);

%if ( isDiag(Sigma) )
%    L = sqrt(Sigma);
%else
%    L = chol(Sigma, 'lower'); % 
%end
L = getChol(Sigma);

z =  randn(D,N); % D N-dimensional independent Gaussian vectors
X = repmat(mu,1,N) + L*z; % Correlated Gaussian vectors 


return;

