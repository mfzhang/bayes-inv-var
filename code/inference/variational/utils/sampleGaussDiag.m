function  X  = sampleGaussDiag( mu, diagSigma, N )
%SAMPLEGAUSSDIAG Sample form a diagonal Gaussian
%   mu: mean
% diagSigma: vector with variances
% N: Number of observations

D = length(mu);
z = randn(D,N); % D N-dimensional independent Gaussian vectors
S = repmat(diagSigma, 1, N);
X = repmat(mu,1,N) + sqrt(S).*z; % Correlated Gaussian vectors 



return;



